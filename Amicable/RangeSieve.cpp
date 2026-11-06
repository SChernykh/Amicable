#include "stdafx.h"
#include "PrimeTables.h"
#include "Engine.h"
#include "RangeGen.h"
#include "RangeSieve.h"
#include "sprp64.h"
#include <atomic>
#include <fstream>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <mutex>

// (x86-64 only: SSE2, and AVX2 / AVX-512 chosen at runtime; the portable code elsewhere)
#if AMICABLE_X64
#include <emmintrin.h>
#include <immintrin.h>
#define RANGE_SIEVE_SSE2 1
#else
#define RANGE_SIEVE_SSE2 0
#endif

// (MSVC compiles AVX intrinsics in any function; GCC needs the target attribute)
#if defined(_MSC_VER)
#define RANGE_SIEVE_TARGET_AVX2
#define RANGE_SIEVE_TARGET_AVX512
#else
#define RANGE_SIEVE_TARGET_AVX2 __attribute__((target("avx2")))
#define RANGE_SIEVE_TARGET_AVX512 __attribute__((target("avx512f")))
#endif

// See RangeSieve.h for the idea. Implementation notes:
//
// Layout: only p coprime to 30 are sieved, one row per residue class mod 30 (8 rows of SEG entries, 32-bit fixed-point sums of
// natural logarithms scaled by 2^26). Position j of class c is p = segStart + 30 * j + W30[c].
//
// Exactness: every increment is a difference of one rounded table entry T(q, e) = round(log(sigma(q^e) / q^e) * 2^26), so the sum
// over all powers of q telescopes to exactly T(q, v_q(N)). The total error is therefore at most 0.5 per distinct prime factor,
// far below MARGIN. If a prime power is longer than the whole range it can divide N(p) for at most one p of the range; that
// position gets the exact exponent of q, computed directly ("sparse" hits), so no power of q is ever missed.
//
// Prime powers are handled by size:
//  - 2^k: stride 2^(k-1) in every row; strides <= PRE2_PERIOD are summed into a periodic pattern used to initialize the rows
//  - 3^k, 5^k (k >= 2) and odd moduli <= M_SMALL: per-row strides, processed row by row while the row is in L1 cache
//  - other moduli coprime to 30, up to the range length: stepping through their hits that are coprime to 30 (wheel pattern)
//  - moduli longer than the range: at most one hit, exact exponent computed directly

bool g_UseRangeSieve = true;
num64 g_RangeSieveMinPrimes = 512; // measured crossover point for the 10^22 search limit (Zen 5)
int g_RangeSieveTierOverride = -1;
int g_RangeSieveScanLevel = -1;
int g_RangeSieveLayout = -1;
// (measured on Zen 5 with 32 threads: k9, t3, k11, k1)
int g_RangeSieveMaxTierByW[8] = { 5, 5, 5, 4, 3, 2, 2, 2 };
std::atomic<num64> g_RangeSieveStats[12];
int g_RangeSieveBand = 1;
int g_RangeSieveBandMinSpanLog2 = 14;
double g_RangeSieveTestBandL0 = 0.0;
double g_RangeSieveTestBandHi = 0.0;
std::atomic<num64> g_RangeSieveNumbersByW[8];
std::atomic<num64> g_RangeSieveTimeByW[8];

// time of the calls of SearchRangeSieve (all threads), and the slowest ranges
namespace
{
	std::atomic<num64> g_RangeSieveTimeNs(0);
	struct SlowRange
	{
		num64 ns;
		num128 M0;
		num64 p_lo, p_hi, survivors;
	};
	enum { NUM_SLOW_RANGES = 8 };
	SlowRange g_slowRanges[NUM_SLOW_RANGES] = {};
	std::atomic<num64> g_slowRangeMin(0);
	std::mutex g_slowRangeLock;
}
std::atomic<num64> g_RangeSieveVerifyCount(0);
std::atomic<num64> g_RangeSieveVerifyErrors(0);

namespace
{
	enum : unsigned
	{
		SEG = 8192,             // entries per row and segment: 8 rows * 8192 * 4 bytes = 256 KB
		SEG_STRIDE = SEG + 32,  // row stride (32 extra entries for the vectorized scan)
		PRE2_PERIOD = 1024,     // period (in entries) of the 2-adic pattern, must divide SEG
		M_SMALL = 2048,         // odd moduli up to M_SMALL use per-row strides
		MAX_Q = 65521,          // largest sieving prime
		T_MAX_E = 10,           // tables of T(q, e) for e <= T_MAX_E

		// Sieve depth by range size, measured for the 10^22 search limit (Zen 5): sieving primes up to 1021 for ranges with fewer
		// than TIER0_MAX_PRIMES primes (cheaper setup), up to 4093 below TIER1_MAX_PRIMES, up to 32749 for longer ranges (fewer
		// exact checks)
		NUM_TIERS = 6,
		TIER0_MAX_PRIMES = 4096,
		TIER1_MAX_PRIMES = 12288,
		NUM_W_CLASSES = 8,
	};

	const long double FX = 67108864.0L; // fixed point: 2^26 per unit of natural logarithm
	const int MARGIN = 64;              // fixed-point slack (actual rounding error is at most ~10)

	const unsigned int W30[8] = { 1, 7, 11, 13, 17, 19, 23, 29 };
	signed char CLS30[30];

	// Wheel stepping for a modulus m coprime to 30, u = index of (m mod 30) in W30:
	// from a hit in class c, the next hit coprime to 30 is DK[u][c] steps of m later, in class CN[u][c], and its position index
	// grows by DK[u][c] * (m / 30) + E30[u][c].
	unsigned char DK[8][8], CN[8][8], E30[8][8];
	// From residue s (mod 30), T0[u][s] steps of m reach the first residue coprime to 30, which is in class T0C[u][s]
	unsigned char T0[8][30], T0C[8][30];

	struct Tier
	{
		unsigned int Q;
		unsigned int numPrimes;
		double rho[100]; // rho[k] = max sigma(x)/x for x < 2^k with all prime factors > Q (rounded up)
	};
	Tier g_tiers[NUM_TIERS];

	std::vector<unsigned int> g_q;   // odd primes 3, 5, 7, ..., MAX_Q
	std::vector<num64> g_c64;        // 2^64 mod q
	std::vector<num64> g_c64sq;      // 2^64 mod q^2
	std::vector<num64> g_R;          // floor(2^64 / q), for fast reduction
	std::vector<num64> g_Rsq;        // floor(2^64 / q^2)
	std::vector<int> g_T[T_MAX_E + 1];
	int g_T2[128];
	num64 g_inv15;                   // 15^-1 mod 2^64
	bool g_initialized = false;

	// Inverses of the setup (a^-1 mod q for every sieving prime of a range): Fermat (a^(q-2)) with Montgomery multiplications in
	// 16-bit lanes (R = 2^16, q < 2^15), 8 primes per SSE2 vector, instead of an extended Euclid with a division per step for every
	// prime. The tables are padded to a multiple of MONT_GROUP with q = 3.
	enum : unsigned { MONT_GROUP = 8, MAX_TIER_PRIMES = 3584 };
	alignas(16) uint16_t g_mq[MAX_TIER_PRIMES];   // q
	alignas(16) uint16_t g_mqn[MAX_TIER_PRIMES];  // -q^-1 mod 2^16
	alignas(16) uint16_t g_mone[MAX_TIER_PRIMES]; // 2^16 mod q (1 in Montgomery form)
	alignas(16) uint16_t g_mr2[MAX_TIER_PRIMES];  // 2^32 mod q
	alignas(16) uint16_t g_mexp[MAX_TIER_PRIMES]; // q - 2
	unsigned char g_mtop[MAX_TIER_PRIMES / MONT_GROUP]; // the highest bit of the exponents of each group

	// Moduli m = q^k <= M_SMALL (q >= 7) of the rows: the first position in row c is jb - off[c] (mod m) with
	// jb = (r - P0) * 30^-1 mod m and off[c] = W30[c] * 30^-1 mod m (no divisions in AddModulus)
	struct RowMod
	{
		num64 R; // floor(2^64 / m)
		unsigned int m;
		unsigned int inv30;
		unsigned short off[8];
		// the same for the 210-wheel (q != 7): 210^-1 mod m and W210[s7][c] * 210^-1 mod m
		unsigned int inv210;
		unsigned short off210[6][8];
	};
	std::vector<RowMod> g_rowMod;
	unsigned short g_firstIdx[32768];
	std::vector<int> g_rowModIdx; // per prime index (q <= M_SMALL): index of q^1 in g_rowMod, -1 for 3 and 5

	FORCEINLINE long double LogSigmaRatio(long double q, int e)
	{
		if (e <= 0)
		{
			return 0.0L;
		}
		// sigma(q^e) / q^e = (q - q^-e) / (q - 1)
		return logl((q - powl(q, static_cast<long double>(-e))) / (q - 1.0L));
	}

	FORCEINLINE int FixedPoint(long double x)
	{
		return static_cast<int>(llroundl(x * FX));
	}

	FORCEINLINE int T(unsigned int qi, int e)
	{
		return (e <= static_cast<int>(T_MAX_E)) ? g_T[e][qi] : FixedPoint(LogSigmaRatio(g_q[qi], e));
	}

	// x mod d for d < 2^32, R = floor(2^64 / d)
	FORCEINLINE num64 FastMod(num64 x, num64 d, num64 R)
	{
		num64 q;
		_umul128(x, R, &q);
		num64 r = x - q * d;
		if (r >= d)
		{
			r -= d;
		}
		return r;
	}

	// x mod d for 128-bit x and d < 2^32, c64 = 2^64 mod d
	FORCEINLINE num64 Mod128(const num128& x, num64 d, num64 c64, num64 R)
	{
		return FastMod(FastMod(HighWord(x), d, R) * c64 + FastMod(LowWord(x), d, R), d, R);
	}

	// a^-1 mod m, gcd(a, m) = 1, 0 < a < m < 2^32
	FORCEINLINE num64 InvMod32(num64 a, num64 m)
	{
		int64_t t = 0, nt = 1;
		unsigned int r = static_cast<unsigned int>(m), nr = static_cast<unsigned int>(a);
		while (nr)
		{
			const unsigned int q = r / nr;
			const int64_t tmp = t - static_cast<int64_t>(q) * nt;
			t = nt;
			nt = tmp;
			const unsigned int tmp2 = r - q * nr;
			r = nr;
			nr = tmp2;
		}
		return static_cast<num64>((t < 0) ? (t + static_cast<int64_t>(m)) : t);
	}

	// a^-1 mod m for any odd m < 2^63 coprime to a < m
	FORCEINLINE num64 InvMod64(num64 a, num64 m)
	{
		int64_t t = 0, nt = 1;
		num64 r = m, nr = a;
		while (nr)
		{
			const num64 q = r / nr;
			const int64_t tmp = t - static_cast<int64_t>(q) * nt;
			t = nt;
			nt = tmp;
			const num64 tmp2 = r - q * nr;
			r = nr;
			nr = tmp2;
		}
		return static_cast<num64>((t < 0) ? (t + static_cast<int64_t>(m)) : t);
	}

	FORCEINLINE num64 MulMod64(num64 a, num64 b, num64 m)
	{
		num64 hi;
		const num64 lo = _umul128(a, b, &hi);
		num64 rem;
		udiv128(hi, lo, m, &rem);
		return rem;
	}

#if RANGE_SIEVE_SSE2
	// x * y * 2^-16 mod q for x, y < q < 2^15 (q odd, qn = -q^-1 mod 2^16), result in [0, q)
	FORCEINLINE __m128i MontMul16(const __m128i x, const __m128i y, const __m128i q, const __m128i qn)
	{
		const __m128i lo = _mm_mullo_epi16(x, y);
		const __m128i hi = _mm_mulhi_epu16(x, y);
		const __m128i m = _mm_mullo_epi16(lo, qn);
		const __m128i mh = _mm_mulhi_epu16(m, q);
		// x y + m q == 0 (mod 2^16): the low halves add up to 2^16 (a carry) unless lo == 0
		const __m128i carry = _mm_andnot_si128(_mm_cmpeq_epi16(lo, _mm_setzero_si128()), _mm_set1_epi16(1));
		const __m128i u = _mm_add_epi16(_mm_add_epi16(hi, mh), carry); // < 1.5 q + 1 < 2 q
		// u - q is negative as a signed 16-bit number iff u < q (q < 2^15)
		const __m128i w = _mm_sub_epi16(u, q);
		return _mm_add_epi16(w, _mm_and_si128(_mm_srai_epi16(w, 15), q));
	}
#endif

	// For the first n primes (n a multiple of MONT_GROUP): inv[i] = a[i]^-1 mod q and r[i] = b[i] * a[i]^-1 mod q (a[i] < q, b[i] < q;
	// a[i] == 0 gives inv[i] == r[i] == 0)
	void BatchInverse(const uint16_t* a, const uint16_t* b, uint16_t* inv, uint16_t* r, unsigned int n)
	{
#if RANGE_SIEVE_SSE2
		const __m128i one = _mm_set1_epi16(1);
		unsigned int i = 0;
		// 4 groups at a time: each exponentiation is a chain of dependent multiplications (one group alone was latency bound)
		for (; i + 4 * MONT_GROUP <= n; i += 4 * MONT_GROUP)
		{
			__m128i am[4];
			__m128i x[4];
			int top = 0;
			for (unsigned int k = 0; k < 4; ++k)
			{
				const unsigned int j = i + k * MONT_GROUP;
				const __m128i q = _mm_load_si128(reinterpret_cast<const __m128i*>(g_mq + j));
				const __m128i qn = _mm_load_si128(reinterpret_cast<const __m128i*>(g_mqn + j));
				am[k] = MontMul16(_mm_load_si128(reinterpret_cast<const __m128i*>(a + j)), _mm_load_si128(reinterpret_cast<const __m128i*>(g_mr2 + j)), q, qn);
				x[k] = _mm_load_si128(reinterpret_cast<const __m128i*>(g_mone + j));
				top = std::max(top, static_cast<int>(g_mtop[j / MONT_GROUP]));
			}
			for (int bit = top; bit >= 0; --bit)
			{
				const __m128i bm = _mm_set1_epi16(static_cast<short>(1 << bit));
				for (unsigned int k = 0; k < 4; ++k)
				{
					const unsigned int j = i + k * MONT_GROUP;
					const __m128i q = _mm_load_si128(reinterpret_cast<const __m128i*>(g_mq + j));
					const __m128i qn = _mm_load_si128(reinterpret_cast<const __m128i*>(g_mqn + j));
					const __m128i e = _mm_load_si128(reinterpret_cast<const __m128i*>(g_mexp + j));
					const __m128i sq = MontMul16(x[k], x[k], q, qn);
					const __m128i y = MontMul16(sq, am[k], q, qn);
					const __m128i mask = _mm_cmpeq_epi16(_mm_and_si128(e, bm), bm);
					x[k] = _mm_or_si128(_mm_and_si128(mask, y), _mm_andnot_si128(mask, sq));
				}
			}
			for (unsigned int k = 0; k < 4; ++k)
			{
				const unsigned int j = i + k * MONT_GROUP;
				const __m128i q = _mm_load_si128(reinterpret_cast<const __m128i*>(g_mq + j));
				const __m128i qn = _mm_load_si128(reinterpret_cast<const __m128i*>(g_mqn + j));
				_mm_store_si128(reinterpret_cast<__m128i*>(inv + j), MontMul16(x[k], one, q, qn));
				_mm_store_si128(reinterpret_cast<__m128i*>(r + j), MontMul16(x[k], _mm_load_si128(reinterpret_cast<const __m128i*>(b + j)), q, qn));
			}
		}
		for (; i < n; i += MONT_GROUP)
		{
			const __m128i q = _mm_load_si128(reinterpret_cast<const __m128i*>(g_mq + i));
			const __m128i qn = _mm_load_si128(reinterpret_cast<const __m128i*>(g_mqn + i));
			const __m128i e = _mm_load_si128(reinterpret_cast<const __m128i*>(g_mexp + i));
			const __m128i am = MontMul16(_mm_load_si128(reinterpret_cast<const __m128i*>(a + i)), _mm_load_si128(reinterpret_cast<const __m128i*>(g_mr2 + i)), q, qn);
			__m128i x = _mm_load_si128(reinterpret_cast<const __m128i*>(g_mone + i));
			for (int bit = g_mtop[i / MONT_GROUP]; bit >= 0; --bit)
			{
				x = MontMul16(x, x, q, qn);
				const __m128i y = MontMul16(x, am, q, qn);
				const __m128i bm = _mm_set1_epi16(static_cast<short>(1 << bit));
				const __m128i mask = _mm_cmpeq_epi16(_mm_and_si128(e, bm), bm);
				x = _mm_or_si128(_mm_and_si128(mask, y), _mm_andnot_si128(mask, x));
			}
			// x = a^-1 * 2^16 mod q
			_mm_store_si128(reinterpret_cast<__m128i*>(inv + i), MontMul16(x, one, q, qn));
			_mm_store_si128(reinterpret_cast<__m128i*>(r + i), MontMul16(x, _mm_load_si128(reinterpret_cast<const __m128i*>(b + i)), q, qn));
		}
#else
		for (unsigned int i = 0; i < n; ++i)
		{
			const num64 q = g_mq[i];
			const num64 ai = (a[i] == 0) ? 0 : InvMod32(a[i], q);
			inv[i] = static_cast<uint16_t>(ai);
			r[i] = static_cast<uint16_t>((ai * b[i]) % q);
		}
#endif
	}

	// The scan: the sums v = row[j] + pat[j % PRE2_PERIOD] (pat: the 2-adic pattern of the row with the constants) of the positions
	// j < len (len a multiple of 32), the ones with lo <= v <= hi written to out (positions) and outVal (sums), their number returned;
	// it sets row[0..len) back to 0 (the rows are 0 outside the scan: the next segment and range start from 0, instead of a copy of
	// the pattern in every row of every segment). SSE2, AVX2 and AVX-512 versions, chosen at RangeSieveInit by the CPU (the build
	// targets SSE2; g_RangeSieveScanLevel can choose a lower one). (x - lo) <= (hi - lo) as unsigned 32-bit numbers <=> lo <= x <= hi.
	typedef unsigned int (*ScanRowFn)(int* row, const int* pat, unsigned int len, int lo, int hi, unsigned short* out, int* outVal);

	FORCEINLINE unsigned int ScanBlock(const int* row, const int* pp, unsigned int j, unsigned int cnt, int lo, int hi, unsigned short* out, int* outVal, unsigned int n)
	{
		for (unsigned int t = 0; t < cnt; ++t)
		{
			const int v = row[j + t] + pp[t];
			if ((v >= lo) && (v <= hi))
			{
				out[n] = static_cast<unsigned short>(j + t);
				outVal[n] = v;
				++n;
			}
		}
		return n;
	}

	unsigned int ScanRowSSE2(int* row, const int* pat, unsigned int len, int lo, int hi, unsigned short* out, int* outVal)
	{
		unsigned int n = 0;
		const bool empty = (hi < lo);
#if RANGE_SIEVE_SSE2
		const __m128i vlo = _mm_set1_epi32(lo - 1);
		const __m128i vhi = _mm_set1_epi32(hi + 1);
		const __m128i zero = _mm_setzero_si128();
#endif
		for (unsigned int j = 0; j < len; j += 16)
		{
			const int* pp = pat + (j & (PRE2_PERIOD - 1));
#if RANGE_SIEVE_SSE2
			__m128i* r = reinterpret_cast<__m128i*>(row + j);
			if (!empty)
			{
				const __m128i a = _mm_add_epi32(_mm_loadu_si128(r), _mm_loadu_si128(reinterpret_cast<const __m128i*>(pp)));
				const __m128i b = _mm_add_epi32(_mm_loadu_si128(r + 1), _mm_loadu_si128(reinterpret_cast<const __m128i*>(pp + 4)));
				const __m128i d = _mm_add_epi32(_mm_loadu_si128(r + 2), _mm_loadu_si128(reinterpret_cast<const __m128i*>(pp + 8)));
				const __m128i e = _mm_add_epi32(_mm_loadu_si128(r + 3), _mm_loadu_si128(reinterpret_cast<const __m128i*>(pp + 12)));
				const __m128i ma = _mm_and_si128(_mm_cmpgt_epi32(a, vlo), _mm_cmplt_epi32(a, vhi));
				const __m128i mb = _mm_and_si128(_mm_cmpgt_epi32(b, vlo), _mm_cmplt_epi32(b, vhi));
				const __m128i md = _mm_and_si128(_mm_cmpgt_epi32(d, vlo), _mm_cmplt_epi32(d, vhi));
				const __m128i me = _mm_and_si128(_mm_cmpgt_epi32(e, vlo), _mm_cmplt_epi32(e, vhi));
				if (_mm_movemask_epi8(_mm_or_si128(_mm_or_si128(ma, mb), _mm_or_si128(md, me))) != 0)
				{
					n = ScanBlock(row, pp, j, 16, lo, hi, out, outVal, n);
				}
			}
			_mm_storeu_si128(r, zero);
			_mm_storeu_si128(r + 1, zero);
			_mm_storeu_si128(r + 2, zero);
			_mm_storeu_si128(r + 3, zero);
#else
			if (!empty)
			{
				n = ScanBlock(row, pp, j, 16, lo, hi, out, outVal, n);
			}
			for (unsigned int t = 0; t < 16; ++t)
			{
				row[j + t] = 0;
			}
#endif
		}
		return n;
	}

#if AMICABLE_X64
	RANGE_SIEVE_TARGET_AVX2 unsigned int ScanRowAVX2(int* row, const int* pat, unsigned int len, int lo, int hi, unsigned short* out, int* outVal)
	{
		unsigned int n = 0;
		const bool empty = (hi < lo);
		const __m256i vlo = _mm256_set1_epi32(lo);
		const __m256i vw = _mm256_set1_epi32(static_cast<int>(static_cast<unsigned int>(hi) - static_cast<unsigned int>(lo)));
		const __m256i zero = _mm256_setzero_si256();
		for (unsigned int j = 0; j < len; j += 32)
		{
			const int* pp = pat + (j & (PRE2_PERIOD - 1));
			__m256i* r = reinterpret_cast<__m256i*>(row + j);
			if (!empty)
			{
				const __m256i a = _mm256_sub_epi32(_mm256_add_epi32(_mm256_loadu_si256(r), _mm256_loadu_si256(reinterpret_cast<const __m256i*>(pp))), vlo);
				const __m256i b = _mm256_sub_epi32(_mm256_add_epi32(_mm256_loadu_si256(r + 1), _mm256_loadu_si256(reinterpret_cast<const __m256i*>(pp + 8))), vlo);
				const __m256i d = _mm256_sub_epi32(_mm256_add_epi32(_mm256_loadu_si256(r + 2), _mm256_loadu_si256(reinterpret_cast<const __m256i*>(pp + 16))), vlo);
				const __m256i e = _mm256_sub_epi32(_mm256_add_epi32(_mm256_loadu_si256(r + 3), _mm256_loadu_si256(reinterpret_cast<const __m256i*>(pp + 24))), vlo);
				const __m256i ma = _mm256_cmpeq_epi32(_mm256_min_epu32(a, vw), a);
				const __m256i mb = _mm256_cmpeq_epi32(_mm256_min_epu32(b, vw), b);
				const __m256i md = _mm256_cmpeq_epi32(_mm256_min_epu32(d, vw), d);
				const __m256i me = _mm256_cmpeq_epi32(_mm256_min_epu32(e, vw), e);
				const __m256i any = _mm256_or_si256(_mm256_or_si256(ma, mb), _mm256_or_si256(md, me));
				if (!_mm256_testz_si256(any, any))
				{
					n = ScanBlock(row, pp, j, 32, lo, hi, out, outVal, n);
				}
			}
			_mm256_storeu_si256(r, zero);
			_mm256_storeu_si256(r + 1, zero);
			_mm256_storeu_si256(r + 2, zero);
			_mm256_storeu_si256(r + 3, zero);
		}
		_mm256_zeroupper();
		return n;
	}

	RANGE_SIEVE_TARGET_AVX512 unsigned int ScanRowAVX512(int* row, const int* pat, unsigned int len, int lo, int hi, unsigned short* out, int* outVal)
	{
		unsigned int n = 0;
		const bool empty = (hi < lo);
		const __m512i vlo = _mm512_set1_epi32(lo);
		const __m512i vw = _mm512_set1_epi32(static_cast<int>(static_cast<unsigned int>(hi) - static_cast<unsigned int>(lo)));
		const __m512i zero = _mm512_setzero_si512();
		alignas(64) int vals[32];
		for (unsigned int j = 0; j < len; j += 32)
		{
			const int* pp = pat + (j & (PRE2_PERIOD - 1));
			const __m512i a = _mm512_add_epi32(_mm512_loadu_si512(row + j), _mm512_loadu_si512(pp));
			const __m512i b = _mm512_add_epi32(_mm512_loadu_si512(row + j + 16), _mm512_loadu_si512(pp + 16));
			_mm512_storeu_si512(row + j, zero);
			_mm512_storeu_si512(row + j + 16, zero);
			if (!empty)
			{
				const unsigned int ka = _mm512_cmple_epu32_mask(_mm512_sub_epi32(a, vlo), vw);
				const unsigned int kb = _mm512_cmple_epu32_mask(_mm512_sub_epi32(b, vlo), vw);
				unsigned int k = ka | (kb << 16);
				if (k)
				{
					_mm512_store_si512(vals, a);
					_mm512_store_si512(vals + 16, b);
					do
					{
						unsigned long t;
						_BitScanForward(&t, k);
						out[n] = static_cast<unsigned short>(j + t);
						outVal[n] = vals[t];
						++n;
						k &= k - 1;
					} while (k);
				}
			}
		}
		_mm256_zeroupper();
		return n;
	}

#endif

	ScanRowFn g_scanRow = ScanRowSSE2;

	// 0: SSE2 (or the portable code), 1: AVX2, 2: AVX-512 (F), each with OS support of the registers
	int DetectScanLevel()
	{
#if !AMICABLE_X64
		return 0;
#elif defined(_MSC_VER)
		int r[4];
		__cpuid(r, 0);
		const int maxLeaf = r[0];
		__cpuid(r, 1);
		const bool osxsave = ((r[2] >> 27) & 1) != 0;
		const bool avx = ((r[2] >> 28) & 1) != 0;
		if (!osxsave || !avx || (maxLeaf < 7))
		{
			return 0;
		}
		const unsigned long long xcr0 = _xgetbv(0);
		__cpuidex(r, 7, 0);
		const bool avx2 = ((r[1] >> 5) & 1) != 0;
		const bool avx512f = ((r[1] >> 16) & 1) != 0;
		if (avx512f && ((xcr0 & 0xE6) == 0xE6))
		{
			return 2;
		}
		if (avx2 && ((xcr0 & 6) == 6))
		{
			return 1;
		}
		return 0;
#else
		__builtin_cpu_init();
		if (__builtin_cpu_supports("avx512f"))
		{
			return 2;
		}
		if (__builtin_cpu_supports("avx2"))
		{
			return 1;
		}
		return 0;
#endif
	}

	FORCEINLINE num64 PopCount64(num64 x)
	{
		x = x - ((x >> 1) & 0x5555555555555555ULL);
		x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
		x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
		return (x * 0x0101010101010101ULL) >> 56;
	}

	struct NumberAndSum
	{
		num128 N;
		num128 S;
	};

	// max sigma(x)/x over x < maxN made of consecutive primes starting after p (non-increasing exponents), as in PrimeTables.cpp
	void MaxSumRatio(PrimeIterator p, num64 maxPower, const num128& maxN, NumberAndSum cur, double& best)
	{
		++p;
		const num128 startSum = cur.S;
		for (num64 k = 1; k <= maxPower; ++k)
		{
			const num128 nextN = cur.N * p.Get();
			if (nextN >= maxN)
			{
				const double r = Num128ToDouble(cur.S) / Num128ToDouble(cur.N);
				if (r > best)
				{
					best = r;
				}
				return;
			}
			cur.N = nextN;
			cur.S = cur.S * p.Get() + startSum;
			MaxSumRatio(p, k, maxN, cur, best);
		}
	}

	// Moduli of a row: fixed capacity (RangeSieveInit checks it: the moduli <= M_SMALL, at most 64 powers of 2 and 64 of 3 and 5)
	enum : unsigned { ROW_CAP = 512 };
	struct RowList
	{
		num64 next[ROW_CAP];
		num64 stride[ROW_CAP];
		int inc[ROW_CAP];
		unsigned int count;

		FORCEINLINE void Add(const num64 n, const num64 s, const int v)
		{
			next[count] = n;
			stride[count] = s;
			inc[count] = v;
			++count;
		}
	};

	// A modulus m > M_SMALL coprime to 30 (wheel-pattern stepping): its hits coprime to 30 repeat with period 30 m, 8 hits per
	// period, one in each row; the 8 hits are stored in order (row offset and position increment after the hit), and k is the
	// next one. (Before: the next row and step from tables indexed by the current row, a chain of dependent loads per hit.)
	struct WheelMod
	{
		unsigned short row[8];  // row offset (row * SEG_STRIDE) of hit k of the period
		unsigned short step[8]; // position increment after hit k
		unsigned int j;         // position of the next hit in its row
		int inc;
		unsigned int k;
	};
	static_assert(7 * SEG_STRIDE < 65536, "WheelMod::row");
	enum : unsigned { WHEEL_MOD_MAX = 1U << 18 }; // WheelMod for m < WHEEL_MOD_MAX (16-bit steps: 6 * (m / 30) + 6 < 2^16)

	// larger moduli (prime powers: fewer than 0.25 hits per segment): the same with 64-bit positions and steps
	struct BigWheelMod
	{
		num64 step[8];
		unsigned int row[8];
		num64 j;
		int inc;
		unsigned int k;
	};

	// ---- layouts ----
	// 30-wheel (s7 = 0): row c of a segment holds p = segStart + 30 j + W30[c]. 210-wheel (s7 = 1..6, long ranges): six sub-sieves,
	// one per p mod 7 = s7, each with the rows of the 30-wheel and the step 210: row c holds p = segStart + 210 j + W210[s7][c], where
	// W210[s7][c] = W30[c] + 30 tau == s7 (mod 7). A sub-sieve segment has as many positions as a 30-wheel segment and covers 7 times
	// as many numbers (6 sub-sieve segments instead of 7 segments): the positions with 7 | p are gone, and so are the hits of 7 (its
	// first power is a constant of the sub-sieve p == r (mod 7), its higher powers are row strides 7^(k-1) there). A modulus m
	// coprime to 210 has the hits p == r (mod m), p == s7 (mod 7) in a sub-sieve: the wheel stepping of 7 m.
	unsigned int W210[7][8]; // W210[0] = W30
	const unsigned int INV7[7] = { 0, 1, 4, 5, 2, 3, 6 };
	num64 g_inv35pow[2][2][48]; // [3 or 5][30- or 210-wheel][k]: (step / q)^-1 mod q^(k-1)
	num64 g_inv7pow[24];        // [k]: 30^-1 mod 7^(k-1)
	num64 g_inv105;             // 105^-1 mod 2^64
	enum : unsigned { LAYOUT210_MIN_SPAN = 4 * 210 * SEG };

	// a modulus q^k <= span (q odd, not q = 3, 5 with k = 1) with its root r (q^k | A1 p + B1 <=> p == r (mod q^k)) and P0 mod q^k
	struct ModRoot
	{
		num64 m;
		num64 r;
		num64 P0modm;
		const RowMod* rm; // the entry of g_rowMod if q >= 7 and m <= M_SMALL
		int inc;
		unsigned int q;
		int k;
	};

	// a power of 2 (2^k <= span, k >= 2) with its root
	struct TwoRoot
	{
		num64 m;
		num64 r;
		int inc;
	};

	struct SparseHit
	{
		num64 p;
		int inc;

		bool operator<(const SparseHit& other) const { return p < other.p; }
	};

	struct SieveState
	{
		RowList rows[8];
		std::vector<WheelMod> wheel;
		std::vector<BigWheelMod> big;
		std::vector<SparseHit> sparse;
		std::vector<ModRoot> mods;
		std::vector<TwoRoot> two;
		std::vector<unsigned char> bandBytes;
		CofactorBand band;
		int pre2[8 * PRE2_PERIOD];
		int S[8 * SEG_STRIDE];
		// setup: A1 mod q, -B1 mod q, A1^-1 mod q and the root -B1 / A1 mod q for every sieving prime
		alignas(16) uint16_t a1[MAX_TIER_PRIMES];
		alignas(16) uint16_t nb1[MAX_TIER_PRIMES];
		alignas(16) uint16_t inv[MAX_TIER_PRIMES];
		alignas(16) uint16_t root[MAX_TIER_PRIMES];
		// positions of a row inside the window and their sums (the scan)
		unsigned short cand[SEG_STRIDE];
		int candVal[SEG_STRIDE];
	};

	THREAD_LOCAL SieveState* t_state = nullptr;
	THREAD_LOCAL int t_lastTier = -1;
	THREAD_LOCAL double t_lastCheckShare = 0.0;

	int DirectSum(num128 N, unsigned int numQ, int* sums = nullptr);

	struct SieveRange
	{
		num128 M0, sumM0, A, B;
		unsigned int mask;
	};

	FORCEINLINE void CheckCandidate(const SieveRange& R, num64 p)
	{
		if ((CandidatesDataMask[Mod385(p + 1)] & R.mask) == 0)
		{
			CheckPair128NoInline(R.M0 * p, R.sumM0 * (p + 1));
		}
	}

	// A prime index (of g_q) at most the first one with T(q, 1) <= delta (numQ if there is none): T(q, 1) = round(2^26 log(1 + 1/q)) and
	// log(1 + 1/q) >= 1 / (q + 1/2), so T(q, 1) <= delta needs q >= 2^26 / (delta + 1/2) - 1/2. g_firstIdx[h]: the first index with
	// g_q[index] >= 2 h. (A binary search over the T(q, 1) was 5% of k9.)
	FORCEINLINE unsigned int FirstIndexT1AtMost(const int delta, const unsigned int numQ)
	{
		const double x = static_cast<double>(FX) / (static_cast<double>(delta) + 0.5) - 0.5 - 1e-3;
		if (x >= 65536.0)
		{
			return numQ;
		}
		const unsigned int idx = (x > 0.0) ? g_firstIdx[static_cast<unsigned int>(x) >> 1] : 0U;
		return (idx < numQ) ? idx : numQ;
	}

	// Exact check of the sieve survivor p (prime) with the sieve sum S: the same result as CheckCandidate. The sum is exactly
	// T2(v_2(N)) + the sum of T(q, v_q(N)) over the odd primes q <= Q (the first numQ primes of g_q) that divide N = A * p + B, so the
	// small factors are found without trial division by all of them: while delta (the part of S not found yet) is > 0, every prime
	// left has T(q, 1) <= T(q, e) <= delta, the next one can't be below the first prime with T(q, 1) <= delta. N without its factors
	// <= Q goes to CheckCofactor. (Sums that don't fit are checked by CheckCandidate: g_RangeSieveStats[4], never seen.)
	NOINLINE void CheckSurvivor(const SieveRange& R, const num64 p, const int S, const unsigned int numQ, const CofactorBand* band)
	{
		if ((CandidatesDataMask[Mod385(p + 1)] & R.mask) != 0)
		{
			return;
		}

		const num128 n1 = R.M0 * p;
		num128 t = R.sumM0 * (p + 1);
		num128 n = t - n1;

		// powers of 2 (as in InitialCheck128)
		unsigned long v;
		if (LowWord(n))
		{
			_BitScanForward64(&v, LowWord(n));
			n >>= v;
		}
		else
		{
			_BitScanForward64(&v, HighWord(n));
			n = HighWord(n) >> v;
			v += 64;
		}
		if (v > 0)
		{
			t *= PowersOf2_128DivisibilityData[v].first;
			if (t > PowersOf2_128DivisibilityData[v].second)
			{
				return;
			}
		}
		if (n >= t)
		{
			return;
		}

		int delta = S - g_T2[v];
		unsigned int qi = (delta > 0) ? FirstIndexT1AtMost(delta, numQ) : 0;
		while (delta > 0)
		{
			if (qi >= numQ)
			{
				break;
			}

			// q = g_q[qi] (PrimeInverses128 has the same index: 0 is 3); the tests are exact for n < 2^100 (see InitialCheck128)
			const num128 inv = PrimeInverses128[qi];
			num128 q = n * inv;
			if (q <= n)
			{
				const InverseData128* data = PowersOfP_128DivisibilityData[qi];
				int e = 1;
				for (;;)
				{
					n = q;
					q *= inv;
					if (HighWord(q) > HighWord(n))
					{
						break;
					}
					++data;
					++e;
				}

				// sigma(q^e) must divide t
				if (data->shift)
				{
					if (LowWord(t) & data->shift_bits)
					{
						return;
					}
					t >>= data->shift;
				}
				t *= data->inverse;
				if (t > data->max_value)
				{
					return;
				}
				if ((n >= t) && (n != 1))
				{
					return;
				}
				delta -= T(qi, e);
				++qi;
				if (delta > 0)
				{
					// (only after a factor: the bound changes with delta)
					const unsigned int first = FirstIndexT1AtMost(delta, numQ);
					if (first > qi)
					{
						qi = first;
					}
				}
				continue;
			}
			++qi;
		}

		if (delta != 0)
		{
			g_RangeSieveStats[4] += 1;
			CheckCandidate(R, p);
			return;
		}

		CheckCofactor(n1, t, n, numQ, band, p, R.sumM0 * (p + 1));
	}

	// Registers the modulus e for the layout s7 (cls: the constants of the rows)
	FORCEINLINE void AddModulusLayout(SieveState& st, int (&cls)[8], const ModRoot& e, const num64 P0, const unsigned int s7)
	{
		const unsigned int* w = W210[s7];
		const num64 m = e.m;
		const num64 r = e.r;
		const unsigned int q = e.q;
		const int inc = e.inc;

		if ((q == 3) || (q == 5) || ((q == 7) && s7))
		{
			if (q == 7)
			{
				// (210-wheel) only the sub-sieve p == r (mod 7) has hits of 7
				if ((r % 7) != s7)
				{
					return;
				}
				if (e.k == 1)
				{
					for (unsigned int c = 0; c < 8; ++c)
					{
						cls[c] += inc;
					}
					return;
				}
			}

			// q | step: p = P0 + step j + w[c] == r (mod q^k)  <=>  (step / q) j == (r - P0 - w[c]) / q (mod q^(k-1))
			const num64 md = m / q;
			const num64 inv = (q == 7) ? g_inv7pow[e.k] : g_inv35pow[(q == 5) ? 1 : 0][s7 ? 1 : 0][e.k];
			for (unsigned int c = 0; c < 8; ++c)
			{
				num64 t = e.P0modm + w[c];
				while (t >= m)
				{
					t -= m;
				}
				const num64 rhs = (r >= t) ? (r - t) : (r + m - t);
				if (rhs % q)
				{
					continue;
				}
				st.rows[c].Add(MulMod64(rhs / q, inv, md), md, inc);
			}
			return;
		}

		if (m <= M_SMALL)
		{
			// first position of row c: (r - P0 - w[c]) * step^-1 mod m
			const RowMod* rm = e.rm;
			const num64 d = (r >= e.P0modm) ? (r - e.P0modm) : (r + m - e.P0modm);
			const num64 jb = FastMod(d * (s7 ? rm->inv210 : rm->inv30), m, rm->R);
			const unsigned short* off = s7 ? rm->off210[s7 - 1] : rm->off;
			for (unsigned int c = 0; c < 8; ++c)
			{
				st.rows[c].Add((jb >= off[c]) ? (jb - off[c]) : (jb + m - off[c]), m, inc);
			}
			return;
		}

		// wheel-pattern stepping of M = m (30-wheel) or 7 m (210-wheel): the first p >= P0 with p == r (mod m) (and p == s7 (mod 7)),
		// then the first one coprime to 30
		num64 x = r + m - e.P0modm;
		if (x >= m)
		{
			x -= m;
		}
		x += P0;
		num64 M = m;
		const num64 step = s7 ? 210 : 30;
		if (s7)
		{
			x += ((s7 + 7 - static_cast<unsigned int>(x % 7)) * INV7[m % 7] % 7) * m;
			M = 7 * m;
		}
		const unsigned int u = static_cast<unsigned int>(CLS30[M % 30]);
		const unsigned int s = static_cast<unsigned int>(x % 30);
		x += T0[u][s] * M;
		unsigned int c = T0C[u][s];
		const num64 j = (x - P0 - w[c]) / step;
		if (m < WHEEL_MOD_MAX)
		{
			// (16-bit steps: (6 M + 209) / step < 2^16)
			WheelMod wm;
			wm.j = static_cast<unsigned int>(j);
			for (unsigned int i = 0; i < 8; ++i)
			{
				const unsigned int cn = CN[u][c];
				wm.row[i] = static_cast<unsigned short>(c * SEG_STRIDE);
				wm.step[i] = static_cast<unsigned short>((DK[u][c] * M + w[c] - w[cn]) / step);
				c = cn;
			}
			wm.inc = inc;
			wm.k = 0;
			st.wheel.push_back(wm);
		}
		else
		{
			BigWheelMod bm;
			bm.j = j;
			for (unsigned int i = 0; i < 8; ++i)
			{
				const unsigned int cn = CN[u][c];
				bm.row[i] = c * SEG_STRIDE;
				bm.step[i] = (DK[u][c] * M + w[c] - w[cn]) / step;
				c = cn;
			}
			bm.inc = inc;
			bm.k = 0;
			st.big.push_back(bm);
		}
	}

	// A range, set up for the layouts
	struct RangeSetup
	{
		const SieveRange* R;
		const Tier* tier;
		num64 p_lo, p_hi, P0;
		unsigned int numQ;
		int base;      // constant part of every sum
		int clsAdd[8]; // first powers of 3 and 5 (W210[s7][c] == W30[c] (mod 30): the same in all layouts)
		num64 numSurvivors, numChecked;
		num64 checkCycles; // (ReadCycleCounter) spent on the survivors
		const CofactorBand* band; // band sieve of the range, or nullptr
	};

	// Sieves the range in the layout s7 (0: the 30-wheel, 1..6: the sub-sieve p == s7 (mod 7) of the 210-wheel) and checks the
	// survivors
	void SieveLayout(SieveState& st, RangeSetup& rs, const unsigned int s7)
	{
		const SieveRange& R = *rs.R;
		const Tier& tier = *rs.tier;
		const unsigned int* w = W210[s7];
		const num64 step = s7 ? 210 : 30;
		const num64 P0 = rs.P0;
		const num64 p_lo = rs.p_lo;
		const num64 p_hi = rs.p_hi;
		const unsigned int numQ = rs.numQ;

		for (unsigned int c = 0; c < 8; ++c)
		{
			st.rows[c].count = 0;
		}
		st.wheel.clear();
		st.big.clear();
		int cls[8];
		for (unsigned int c = 0; c < 8; ++c)
		{
			cls[c] = rs.base + rs.clsAdd[c];
		}
		for (unsigned int i = 0; i < 8 * PRE2_PERIOD; ++i)
		{
			st.pre2[i] = 0;
		}

		// powers of 2: p = P0 + step j + w[c] == r (mod 2^k)  <=>  (step / 2) j == (r - P0 - w[c]) / 2 (mod 2^(k-1))
		const num64 invHalf = s7 ? g_inv105 : g_inv15;
		for (const TwoRoot& e : st.two)
		{
			const num64 m = e.m;
			const num64 stride = m >> 1;
			for (unsigned int c = 0; c < 8; ++c)
			{
				const num64 j0 = ((((e.r - P0 - w[c]) & (m - 1)) >> 1) * invHalf) & (stride - 1);
				if (stride <= PRE2_PERIOD)
				{
					int* pat = st.pre2 + c * PRE2_PERIOD;
					for (num64 j = j0; j < PRE2_PERIOD; j += stride)
					{
						pat[j] += e.inc;
					}
				}
				else
				{
					st.rows[c].Add(j0, stride, e.inc);
				}
			}
		}

		for (const ModRoot& e : st.mods)
		{
			AddModulusLayout(st, cls, e, P0, s7);
		}

		// the constants, folded into the 2-adic pattern
		for (unsigned int c = 0; c < 8; ++c)
		{
			const int add = cls[c];
			int* pat = st.pre2 + c * PRE2_PERIOD;
			for (unsigned int j = 0; j < PRE2_PERIOD; ++j)
			{
				pat[j] += add;
			}
		}

		// ---- segments ----
		size_t sp = 0;
		for (num64 segStart = P0; segStart <= p_hi; segStart += step * SEG)
		{
			const num64 remaining = (p_hi - segStart) / step + 1;
			const unsigned int segLen = (remaining < SEG) ? static_cast<unsigned int>(remaining) : static_cast<unsigned int>(SEG);
			const unsigned int scanLen = (segLen + 31) & ~31U;

			// (the rows are 0 here; the scan adds the 2-adic pattern, which also contains the constant terms)

			// per-row moduli, row by row (the row stays in L1 cache)
			for (unsigned int c = 0; c < 8; ++c)
			{
				int* row = st.S + c * SEG_STRIDE;
				RowList& L = st.rows[c];
				const size_t cnt = L.count;
				num64* next = L.next;
				const num64* stride = L.stride;
				const int* inc = L.inc;
				for (size_t i = 0; i < cnt; ++i)
				{
					num64 j = next[i];
					if (j < segLen)
					{
						const num64 s = stride[i];
						const int v = inc[i];
						const num64 s4 = s * 4;
						while (j + s * 3 < segLen)
						{
							row[j] += v;
							row[j + s] += v;
							row[j + s * 2] += v;
							row[j + s * 3] += v;
							j += s4;
						}
						while (j < segLen)
						{
							row[j] += v;
							j += s;
						}
					}
					next[i] = j - segLen;
				}
			}

			// wheel-pattern moduli
			{
				WheelMod* const wm = st.wheel.data();
				const size_t cnt = st.wheel.size();
				int* const Sbase = st.S;
				const unsigned int sl = segLen;
				for (size_t i = 0; i < cnt; ++i)
				{
					WheelMod& wh = wm[i];
					unsigned int j = wh.j;
					if (j < sl)
					{
						unsigned int k = wh.k;
						const int v = wh.inc;
						for (;;)
						{
							Sbase[wh.row[k] + j] += v;
							j += wh.step[k];
							if (j >= sl)
							{
								k = (k + 1) & 7;
								break;
							}
							const unsigned int k1 = (k + 1) & 7;
							Sbase[wh.row[k1] + j] += v;
							j += wh.step[k1];
							k = (k + 2) & 7;
							if (j >= sl)
							{
								break;
							}
						}
						wh.k = k;
					}
					wh.j = j - sl;
				}

				BigWheelMod* const bm = st.big.data();
				const size_t numBig = st.big.size();
				for (size_t i = 0; i < numBig; ++i)
				{
					BigWheelMod& b = bm[i];
					num64 j = b.j;
					if (j < segLen)
					{
						unsigned int k = b.k;
						const int v = b.inc;
						do
						{
							Sbase[b.row[k] + j] += v;
							j += b.step[k];
							k = (k + 1) & 7;
						} while (j < segLen);
						b.k = k;
					}
					b.j = j - segLen;
				}
			}

			// sparse hits (of this layout)
			const num64 segEnd = segStart + step * segLen;
			while ((sp < st.sparse.size()) && (st.sparse[sp].p < segEnd))
			{
				const num64 p = st.sparse[sp].p;
				if ((s7 == 0) || ((p % 7) == s7))
				{
					const unsigned int c = static_cast<unsigned int>(CLS30[p % 30]);
					st.S[c * SEG_STRIDE + (p - segStart - w[c]) / step] += st.sparse[sp].inc;
				}
				++sp;
			}

#ifdef RANGE_SIEVE_VERIFY
			// debug: compare every sieved sum with trial division
			for (PrimeIterator it(segStart); it.Get() < segEnd; ++it)
			{
				const num64 p = it.Get();
				if ((p < 7) || (p < p_lo) || (p > p_hi) || (s7 && ((p % 7) != s7)))
				{
					continue;
				}
				const unsigned int c = static_cast<unsigned int>(CLS30[p % 30]);
				const num64 j = (p - segStart - w[c]) / step;
				++g_RangeSieveVerifyCount;
				if (st.S[c * SEG_STRIDE + j] + st.pre2[c * PRE2_PERIOD + (j & (PRE2_PERIOD - 1))] != DirectSum(R.A * p + R.B, numQ))
				{
					++g_RangeSieveVerifyErrors;
				}
			}
#endif

			// window for this segment: S(p) must be in [log(T/N) - log(rho), log(T/N)]
			const num64 lo_p = (segStart > p_lo) ? segStart : p_lo;
			const num64 hi_p = ((segEnd - 1) < p_hi) ? (segEnd - 1) : p_hi;
			int lo_fx = 1;
			int hi_fx = 0;
			if (lo_p <= hi_p)
			{
				const num128 Nlo = R.A * lo_p + R.B;
				const num128 Nhi = R.A * hi_p + R.B;
				const long double tmin = Num128ToDouble(R.B * (lo_p + 1)) / Num128ToDouble(Nlo);
				const long double tmax = Num128ToDouble(R.B * (hi_p + 1)) / Num128ToDouble(Nhi);
				unsigned int L = 0;
				{
					num128 t = Nhi;
					while (t > 1)
					{
						t >>= 1;
						++L;
					}
				}
				const long double rho = tier.rho[(L + 1 < 100) ? (L + 1) : 99];
				lo_fx = static_cast<int>(floorl((logl(tmin) - logl(rho)) * FX)) - MARGIN;
				hi_fx = static_cast<int>(ceill(logl(tmax) * FX)) + MARGIN;
			}

			// scan for positions inside the window (rare)
			for (unsigned int c = 0; c < 8; ++c)
			{
				int* row = st.S + c * SEG_STRIDE;
				const unsigned int numInWindow = g_scanRow(row, st.pre2 + c * PRE2_PERIOD, scanLen, lo_fx, hi_fx, st.cand, st.candVal);
				const num64 c0 = numInWindow ? ReadCycleCounter() : 0;
				for (unsigned int i = 0; i < numInWindow; ++i)
				{
					const unsigned int t = st.cand[i];
					const int v = st.candVal[i];
					const num64 p = segStart + step * t + w[c];
					if ((lo_p <= p) && (p <= hi_p))
					{
						++rs.numSurvivors;
						if (IsPrime(p))
						{
							++rs.numChecked;
#ifdef RANGE_SIEVE_CHECK_COMPARE
							const unsigned int f0 = GetNumFoundPairsInThisThread();
							CheckCandidate(R, p);
							const unsigned int f1 = GetNumFoundPairsInThisThread();
							CheckSurvivor(R, p, v, numQ, rs.band);
							const unsigned int f2 = GetNumFoundPairsInThisThread();
							if ((f1 - f0) != (f2 - f1))
							{
								g_RangeSieveStats[5] += 1;
								std::cerr << "CheckSurvivor differs for p = " << p << ", M0 = " << R.M0 << " (old " << (f1 - f0) << ", new " << (f2 - f1) << ")" << std::endl;
							}
#else
							CheckSurvivor(R, p, v, numQ, rs.band);
#endif
						}
					}
				}
				if (numInWindow)
				{
					rs.checkCycles += ReadCycleCounter() - c0;
				}
			}
		}
	}
}

void RangeSieveInit()
{
	if (g_initialized)
	{
		return;
	}

	{
		const int detected = DetectScanLevel();
		const int level = ((g_RangeSieveScanLevel >= 0) && (g_RangeSieveScanLevel < detected)) ? g_RangeSieveScanLevel : detected;
#if AMICABLE_X64
		g_scanRow = (level == 2) ? ScanRowAVX512 : ((level == 1) ? ScanRowAVX2 : ScanRowSSE2);
#else
		g_scanRow = ScanRowSSE2;
#endif
		g_RangeSieveScanLevel = level;
		CheckCofactorInit(level);
	}

	for (int i = 0; i < 30; ++i)
	{
		CLS30[i] = -1;
	}
	for (unsigned int c = 0; c < 8; ++c)
	{
		CLS30[W30[c]] = static_cast<signed char>(c);
	}

	// 210-wheel: W210[s7][c] = W30[c] + 30 tau == s7 (mod 7) (30^-1 == 4 (mod 7))
	for (unsigned int s7 = 0; s7 < 7; ++s7)
	{
		for (unsigned int c = 0; c < 8; ++c)
		{
			W210[s7][c] = s7 ? (W30[c] + 30 * (((s7 + 7 - W30[c] % 7) * 4) % 7)) : W30[c];
		}
	}
	{
		const unsigned int qs[2] = { 3, 5 };
		const num64 steps[2] = { 30, 210 };
		for (unsigned int qi = 0; qi < 2; ++qi)
		{
			for (unsigned int li = 0; li < 2; ++li)
			{
				num64 md = 1;
				for (unsigned int k = 0; k < 48; ++k)
				{
					g_inv35pow[qi][li][k] = 0;
				}
				for (unsigned int k = 2; k < 48; ++k)
				{
					if (md > num64(-1) / 8 / qs[qi])
					{
						break;
					}
					md *= qs[qi];
					g_inv35pow[qi][li][k] = InvMod64((steps[li] / qs[qi]) % md, md);
				}
			}
		}
		num64 md = 1;
		for (unsigned int k = 0; k < 24; ++k)
		{
			g_inv7pow[k] = 0;
		}
		for (unsigned int k = 2; k < 24; ++k)
		{
			md *= 7;
			g_inv7pow[k] = InvMod64(30 % md, md);
		}
		g_inv105 = 105;
		for (int i = 0; i < 6; ++i)
		{
			g_inv105 *= 2 - 105 * g_inv105;
		}
	}

	for (unsigned int u = 0; u < 8; ++u)
	{
		const unsigned int w = W30[u];
		for (unsigned int c = 0; c < 8; ++c)
		{
			unsigned int dk = 1;
			while (CLS30[(W30[c] + dk * w) % 30] < 0)
			{
				++dk;
			}
			const unsigned int cn = static_cast<unsigned int>(CLS30[(W30[c] + dk * w) % 30]);
			DK[u][c] = static_cast<unsigned char>(dk);
			CN[u][c] = static_cast<unsigned char>(cn);
			E30[u][c] = static_cast<unsigned char>((dk * w + W30[c] - W30[cn]) / 30);
		}
		for (unsigned int s = 0; s < 30; ++s)
		{
			unsigned int t = 0;
			while (CLS30[(s + t * w) % 30] < 0)
			{
				++t;
			}
			T0[u][s] = static_cast<unsigned char>(t);
			T0C[u][s] = static_cast<unsigned char>(CLS30[(s + t * w) % 30]);
		}
	}

	g_q.clear();
	for (PrimeIterator it(3); it.Get() <= MAX_Q; ++it)
	{
		g_q.push_back(static_cast<unsigned int>(it.Get()));
	}

	const size_t n = g_q.size();
	g_c64.resize(n);
	g_c64sq.resize(n);
	g_R.resize(n);
	g_Rsq.resize(n);
	for (unsigned int e = 0; e <= T_MAX_E; ++e)
	{
		g_T[e].resize(n);
	}
	for (size_t i = 0; i < n; ++i)
	{
		const num64 q = g_q[i];
		const num64 q2 = q * q;
		g_c64[i] = (num64(-1) % q + 1) % q;
		g_c64sq[i] = (num64(-1) % q2 + 1) % q2;
		g_R[i] = num64(-1) / q;
		g_Rsq[i] = num64(-1) / q2;
		for (unsigned int e = 0; e <= T_MAX_E; ++e)
		{
			g_T[e][i] = FixedPoint(LogSigmaRatio(static_cast<long double>(q), static_cast<int>(e)));
		}
	}
	for (int v = 0; v < 128; ++v)
	{
		g_T2[v] = FixedPoint(LogSigmaRatio(2.0L, v));
	}

	g_inv15 = 15;
	for (int i = 0; i < 6; ++i)
	{
		g_inv15 *= 2 - 15 * g_inv15;
	}

	// FirstIndexT1AtMost
	{
		unsigned int qi = 0;
		for (unsigned int h = 0; h < 32768; ++h)
		{
			while ((qi < n) && (g_q[qi] < 2 * h))
			{
				++qi;
			}
			g_firstIdx[h] = static_cast<unsigned short>(qi);
		}
	}

	// Montgomery tables of BatchInverse (padding: q = 3)
	for (unsigned int i = 0; i < MAX_TIER_PRIMES; ++i)
	{
		const unsigned int q = ((i < n) && (g_q[i] < 32768)) ? g_q[i] : 3;
		unsigned int qinv = q;
		for (int j = 0; j < 4; ++j)
		{
			qinv *= 2 - q * qinv;
		}
		g_mq[i] = static_cast<uint16_t>(q);
		g_mqn[i] = static_cast<uint16_t>(0U - qinv);
		g_mone[i] = static_cast<uint16_t>(65536U % q);
		g_mr2[i] = static_cast<uint16_t>((num64(1) << 32) % q);
		g_mexp[i] = static_cast<uint16_t>(q - 2);
	}
	for (unsigned int g = 0; g < MAX_TIER_PRIMES / MONT_GROUP; ++g)
	{
		unsigned int e = 0;
		for (unsigned int i = g * MONT_GROUP; i < (g + 1) * MONT_GROUP; ++i)
		{
			e |= g_mexp[i];
		}
		unsigned int top = 0;
		while ((e >> top) > 1)
		{
			++top;
		}
		g_mtop[g] = static_cast<unsigned char>(top);
	}

	// moduli of the rows
	g_rowMod.clear();
	g_rowModIdx.assign(n, -1);
	for (size_t i = 0; (i < n) && (g_q[i] <= M_SMALL); ++i)
	{
		const unsigned int q = g_q[i];
		if ((q == 3) || (q == 5))
		{
			continue;
		}
		g_rowModIdx[i] = static_cast<int>(g_rowMod.size());
		for (num64 m = q; m <= M_SMALL; m *= q)
		{
			RowMod rm;
			rm.R = num64(-1) / m;
			rm.m = static_cast<unsigned int>(m);
			rm.inv30 = static_cast<unsigned int>(InvMod32(30 % m, m));
			for (unsigned int c = 0; c < 8; ++c)
			{
				rm.off[c] = static_cast<unsigned short>((W30[c] * rm.inv30) % m);
			}
			rm.inv210 = 0;
			for (unsigned int s7 = 1; s7 < 7; ++s7)
			{
				for (unsigned int c = 0; c < 8; ++c)
				{
					rm.off210[s7 - 1][c] = 0;
				}
			}
			if (q != 7)
			{
				rm.inv210 = static_cast<unsigned int>(InvMod32(210 % m, m));
				for (unsigned int s7 = 1; s7 < 7; ++s7)
				{
					for (unsigned int c = 0; c < 8; ++c)
					{
						rm.off210[s7 - 1][c] = static_cast<unsigned short>((W210[s7][c] * rm.inv210) % m);
					}
				}
			}
			g_rowMod.push_back(rm);
		}
	}
	if (g_rowMod.size() + 128 > ROW_CAP)
	{
		std::cerr << "RangeSieveInit: ROW_CAP is too small" << std::endl;
		abort();
	}

	const unsigned int tierQ[NUM_TIERS] = { 1021, 2039, 4093, 8191, 16381, 32749 };
	for (unsigned int t = 0; t < NUM_TIERS; ++t)
	{
		Tier& tier = g_tiers[t];
		tier.Q = tierQ[t];
		tier.numPrimes = 0;
		while ((tier.numPrimes < n) && (g_q[tier.numPrimes] <= tier.Q))
		{
			++tier.numPrimes;
		}
		for (unsigned int k = 0; k < 100; ++k)
		{
			double best = 1.0;
			if (k > 0)
			{
				NumberAndSum cur;
				cur.N = 1;
				cur.S = 1;
				MaxSumRatio(PrimeIterator(tier.Q), num64(-1), num128(1) << k, cur, best);
			}
			tier.rho[k] = best * (1.0 + 1e-12);
		}
		// (BatchInverse: 16-bit Montgomery lanes, q < 2^15)
		if ((((tier.numPrimes + MONT_GROUP - 1) & ~(MONT_GROUP - 1)) > MAX_TIER_PRIMES) || (g_q[tier.numPrimes - 1] >= 32768))
		{
			std::cerr << "RangeSieveInit: tier " << t << " has too many or too large primes for BatchInverse" << std::endl;
			abort();
		}
	}

	g_initialized = true;
}

num64 RangeSieveCountPrimes(num64 lo, num64 hi)
{
	if (hi < lo)
	{
		return 0;
	}

	num64 count = 0;
	// ("small" is a macro in the Windows headers)
	const num64 smallPrimes[4] = { 2, 3, 5, 7 };
	for (unsigned int i = 0; i < 4; ++i)
	{
		if ((lo <= smallPrimes[i]) && (smallPrimes[i] <= hi))
		{
			++count;
		}
	}
	if (hi < 11)
	{
		return count;
	}
	if (lo < 11)
	{
		lo = 11;
	}

	while (bitOffset[lo % PrimeTableParameters::Modulo] >= PrimeTableParameters::NumOffsets)
	{
		++lo;
	}
	while (bitOffset[hi % PrimeTableParameters::Modulo] >= PrimeTableParameters::NumOffsets)
	{
		--hi;
	}
	if (hi < lo)
	{
		return count;
	}

	const num64 b0 = (lo / PrimeTableParameters::Modulo) * PrimeTableParameters::NumOffsets + bitOffset[lo % PrimeTableParameters::Modulo];
	const num64 b1 = (hi / PrimeTableParameters::Modulo) * PrimeTableParameters::NumOffsets + bitOffset[hi % PrimeTableParameters::Modulo];
	const num64* words = reinterpret_cast<const num64*>(MainPrimeTable);
	const num64 w0 = b0 / 64;
	const num64 w1 = b1 / 64;
	if (w0 == w1)
	{
		return count + PopCount64((words[w0] >> (b0 % 64)) & (num64(-1) >> (63 - (b1 - b0))));
	}
	count += PopCount64(words[w0] >> (b0 % 64));
	for (num64 w = w0 + 1; w < w1; ++w)
	{
		count += PopCount64(words[w]);
	}
	count += PopCount64(words[w1] & (num64(-1) >> (63 - (b1 % 64))));
	return count;
}

void SearchRangeSieve(const RangeData& r, num64 p_lo, num64 p_hi, unsigned int is_over_abundant_mask, num64 num_primes, int tierHint)
{
	num64 primes[MaxPrimeFactors];
	unsigned int exponents[MaxPrimeFactors];
	const unsigned int n = static_cast<unsigned int>(r.last_factor_index + 1);
	for (unsigned int i = 0; i < n; ++i)
	{
		primes[i] = r.factors[i].p.Get();
		exponents[i] = r.factors[i].k;
	}
	SearchRangeSieve(r.value, r.sum, primes, exponents, n, p_lo, p_hi, is_over_abundant_mask, num_primes, tierHint);
}

void SearchRangeSieve(num128 M0, num128 sumM0, const num64* M0_primes, const unsigned int* M0_exponents, unsigned int M0_num_factors, num64 p_lo, num64 p_hi, unsigned int is_over_abundant_mask, num64 num_primes, int tierHint)
{
	(void)M0_exponents;
	t_lastTier = -1;

	if (p_hi < p_lo)
	{
		return;
	}

	if (!t_state)
	{
		t_state = new SieveState();
	}
	SieveState& st = *t_state;
	const auto startTime = std::chrono::steady_clock::now();

	SieveRange R;
	R.M0 = M0;
	R.sumM0 = sumM0;
	R.A = sumM0 - M0;
	R.B = sumM0;
	R.mask = is_over_abundant_mask;

	// Primes below 7 are not in the wheel
	const num64 wheelPrimes[3] = { 2, 3, 5 };
	for (const num64 p : wheelPrimes)
	{
		if ((p_lo <= p) && (p <= p_hi))
		{
			CheckCandidate(R, p);
		}
	}

	// N(p) = A * p + B = g * (A1 * p + B1), gcd(A1, B1) = 1
	num128 g;
	{
		num128 a = R.A;
		num128 b = R.B;
		while (b != 0)
		{
			const num128 t = a % b;
			a = b;
			b = t;
		}
		g = a;
	}
	const num128 A1 = R.A / g;
	const num128 B1 = R.B / g;

	// The window height W above the known part K of N (as in the GPU version): K = the prime powers of g, times 2 if A1 * p + B1 is
	// always even. sigma(N) / N >= sigma(K) / K, and T / N grows with p: if log(T / N) at p_hi is below log(sigma(K) / K), no p of
	// the range gives an amicable number. With a large W the survivors hardly depend on the sieve depth (a smaller Q is faster).
	// exponents of the odd prime factors of M0 in g
	int gExp[MaxPrimeFactors];
	for (unsigned int i = 0; i < M0_num_factors; ++i)
	{
		const num64 q = M0_primes[i];
		gExp[i] = 0;
		if (q & 1)
		{
			// x / q exactly if q | x (x * q^-1 mod 2^128 <= x), see InitialCheck128
			num128 qinv = q;
			for (int k = 0; k < 7; ++k)
			{
				qinv *= num128(2) - qinv * q;
			}
			num128 x = g;
			for (;;)
			{
				const num128 y = x * qinv;
				if ((y > x) || (y * q != x))
				{
					break;
				}
				x = y;
				++gExp[i];
			}
		}
	}

	unsigned int wClass;
	long double logHi;
	long double logK2;
	unsigned int v2K;
	{
		unsigned int v2g = 0;
		while (((LowWord(g) >> v2g) & 1) == 0)
		{
			++v2g;
		}
		v2K = v2g + (((LowWord(A1) & 1) && (LowWord(B1) & 1)) ? 1 : 0);
		logK2 = LogSigmaRatio(2.0L, static_cast<int>(v2K));
		long double logK = logK2;
		for (unsigned int i = 0; i < M0_num_factors; ++i)
		{
			logK += LogSigmaRatio(static_cast<long double>(M0_primes[i]), gExp[i]);
		}
		logHi = logl(static_cast<long double>(Num128ToDouble(R.B * (p_hi + 1)))) - logl(static_cast<long double>(Num128ToDouble(R.A * p_hi + R.B)));
		const double W = static_cast<double>(logHi - logK);
		if (W < -1e-9)
		{
			g_RangeSieveStats[6] += 1;
			return;
		}
		const double wb[NUM_W_CLASSES - 1] = { 1e-4, 1e-3, 1e-2, 0.03, 0.1, 0.3, 1.0 };
		wClass = 0;
		while ((wClass < NUM_W_CLASSES - 1) && (W >= wb[wClass]))
		{
			++wClass;
		}
	}

	// Sieve depth: fewer primes for short ranges (cheaper setup), more for long ranges (fewer exact checks), at most
	// g_RangeSieveMaxTierByW[wClass]
	unsigned int tierIndex = (num_primes < TIER0_MAX_PRIMES) ? 0U : ((num_primes < TIER1_MAX_PRIMES) ? 2U : 5U);
	if (tierIndex > static_cast<unsigned int>(g_RangeSieveMaxTierByW[wClass]))
	{
		tierIndex = static_cast<unsigned int>(g_RangeSieveMaxTierByW[wClass]);
	}
	if (tierHint >= 0)
	{
		tierIndex = (tierHint < static_cast<int>(NUM_TIERS)) ? static_cast<unsigned int>(tierHint) : (NUM_TIERS - 1);
	}
	if ((g_RangeSieveTierOverride >= 0) && (g_RangeSieveTierOverride < static_cast<int>(NUM_TIERS)))
	{
		tierIndex = static_cast<unsigned int>(g_RangeSieveTierOverride);
	}
	const Tier& tier = g_tiers[tierIndex];

	// (a multiple of 210 for both layouts)
	const num64 P0 = (p_lo / 210) * 210;
	const num64 span = p_hi - P0 + 1;
	const num128 Nmax1 = A1 * p_hi + B1;

	st.sparse.clear();
	st.mods.clear();
	st.two.clear();

	RangeSetup rs;
	rs.R = &R;
	rs.tier = &tier;
	rs.p_lo = p_lo;
	rs.p_hi = p_hi;
	rs.P0 = P0;
	rs.numQ = tier.numPrimes;
	rs.numSurvivors = 0;
	rs.numChecked = 0;
	rs.checkCycles = 0;
	const num64 startCycles = ReadCycleCounter();
	int (&clsAdd)[8] = rs.clsAdd;
	for (unsigned int c = 0; c < 8; ++c)
	{
		clsAdd[c] = 0;
	}
	int base = 0;

	// ---- 2-adic part ----
	int v2g = 0;
	{
		num128 x = g;
		while ((LowWord(x) & 1) == 0)
		{
			x >>= 1;
			++v2g;
		}
	}
	if ((LowWord(A1) & 1) && (LowWord(B1) & 1))
	{
		// every odd p makes A1 * p + B1 even
		base += g_T2[v2g + 1];

		num64 inv = LowWord(A1);
		for (int i = 0; i < 6; ++i)
		{
			inv *= 2 - LowWord(A1) * inv;
		}
		const num64 r64 = num64(0) - LowWord(B1) * inv; // -B1 / A1 mod 2^64

		for (int k = 2; k < 64; ++k)
		{
			const num64 m = num64(1) << k;
			const num64 r = r64 & (m - 1);
			const int inc = g_T2[v2g + k] - g_T2[v2g + k - 1];
			if (m > span)
			{
				const num64 x = P0 + ((r - P0) & (m - 1));
				if ((p_lo <= x) && (x <= p_hi) && (CLS30[x % 30] >= 0))
				{
					num128 Nx = A1 * x + B1;
					int v = 0;
					while ((LowWord(Nx) & 1) == 0)
					{
						Nx >>= 1;
						++v;
					}
					SparseHit h;
					h.p = x;
					h.inc = g_T2[v2g + v] - g_T2[v2g + k - 1];
					st.sparse.push_back(h);
				}
				break;
			}
			TwoRoot e;
			e.m = m;
			e.r = r;
			e.inc = inc;
			st.two.push_back(e);
			if (num128(m) * 2 > Nmax1)
			{
				break;
			}
		}
	}
	else
	{
		base += g_T2[v2g];
	}

	// ---- odd primes ----
	// odd prime factors of g (they are among the prime factors of M0)
	num64 gp[MaxPrimeFactors];
	int ge[MaxPrimeFactors];
	unsigned int ng = 0;
	for (unsigned int i = 0; i < M0_num_factors; ++i)
	{
		const num64 q = M0_primes[i];
		if ((q & 1) && (q <= tier.Q))
		{
			const int e = gExp[i];
			if (e > 0)
			{
				// insertion sort: the factorization of M0 doesn't have to be in increasing order
				unsigned int pos = ng;
				while ((pos > 0) && (gp[pos - 1] > q))
				{
					gp[pos] = gp[pos - 1];
					ge[pos] = ge[pos - 1];
					--pos;
				}
				gp[pos] = q;
				ge[pos] = e;
				++ng;
			}
		}
	}

	unsigned int gi = 0;
	const unsigned int numQ = tier.numPrimes;

	// A1^-1 and the roots -B1 / A1 mod q of all primes at once
	{
		const unsigned int numQv = (numQ + MONT_GROUP - 1) & ~(MONT_GROUP - 1);
		for (unsigned int qi = 0; qi < numQ; ++qi)
		{
			const num64 q = g_q[qi];
			const num64 Rq = g_R[qi];
			const num64 b1 = Mod128(B1, q, g_c64[qi], Rq);
			st.a1[qi] = static_cast<uint16_t>(Mod128(A1, q, g_c64[qi], Rq));
			st.nb1[qi] = static_cast<uint16_t>(b1 ? (q - b1) : 0);
		}
		for (unsigned int qi = numQ; qi < numQv; ++qi)
		{
			st.a1[qi] = 1;
			st.nb1[qi] = 0;
		}
		BatchInverse(st.a1, st.nb1, st.inv, st.root, numQv);
	}

	for (unsigned int qi = 0; qi < numQ; ++qi)
	{
		const unsigned int q = g_q[qi];
		int eg = 0;
		if ((gi < ng) && (gp[gi] == q))
		{
			eg = ge[gi];
			++gi;
			base += T(qi, eg);
		}

		const num64 Rq = g_R[qi];
		if (st.a1[qi] == 0)
		{
			continue; // q divides A1 but not B1, so it never divides A1 * p + B1
		}
		const num64 inv = st.inv[qi];
		num64 r = st.root[qi];
		num64 m = q;
		num64 P0modm = FastMod(P0, q, Rq);
		const RowMod* rm = ((q <= M_SMALL) && (g_rowModIdx[qi] >= 0)) ? &g_rowMod[static_cast<size_t>(g_rowModIdx[qi])] : nullptr;

		for (int k = 1;; ++k)
		{
			const int inc = T(qi, eg + k) - T(qi, eg + k - 1);
			if (m > span)
			{
				// at most one p in the range: find it and use the exact exponent of q
				num64 x = r + m - P0modm;
				if (x >= m)
				{
					x -= m;
				}
				x += P0;
				if ((p_lo <= x) && (x <= p_hi) && (CLS30[x % 30] >= 0))
				{
					num128 Nx = A1 * x + B1;
					int e = 0;
					while ((Nx % static_cast<num64>(q)) == 0)
					{
						Nx /= q;
						++e;
					}
					SparseHit h;
					h.p = x;
					h.inc = T(qi, eg + e) - T(qi, eg + k - 1);
					st.sparse.push_back(h);
				}
				break;
			}

			if (((q == 3) || (q == 5)) && (k == 1))
			{
				// constant on each residue class mod 30
				for (unsigned int c = 0; c < 8; ++c)
				{
					if (W30[c] % q == r)
					{
						clsAdd[c] += inc;
					}
				}
			}
			else
			{
				ModRoot e;
				e.m = m;
				e.r = r;
				e.P0modm = P0modm;
				e.rm = rm;
				e.inc = inc;
				e.q = q;
				e.k = k;
				st.mods.push_back(e);
			}
			if (rm)
			{
				// (the entry of the next power of q, if it is <= M_SMALL)
				rm = (m * q <= M_SMALL) ? (rm + 1) : nullptr;
			}

			if (num128(m) * q > Nmax1)
			{
				break;
			}

			// lift r to the next power of q: A1 * (r + t * m) + B1 == 0 (mod m * q)
			const num64 mq = m * q;
			num64 digit; // (A1 * r + B1) / m mod q
			if (k == 1)
			{
				const num64 Rsq = g_Rsq[qi];
				const num64 A1m = Mod128(A1, mq, g_c64sq[qi], Rsq);
				const num64 B1m = Mod128(B1, mq, g_c64sq[qi], Rsq);
				const num64 val = FastMod(A1m * r + B1m, mq, Rsq); // divisible by q
				_umul128(val, Rq, &digit);
				if (val - digit * q >= q)
				{
					++digit;
				}
				P0modm = FastMod(P0, mq, Rsq);
			}
			else
			{
				const num64 A1m = LowWord(A1 % num128(mq));
				const num64 B1m = LowWord(B1 % num128(mq));
				const num64 val = (MulMod64(A1m, r, mq) + B1m) % mq;
				digit = (val / m) % q;
				P0modm = P0 % mq;
			}
			const num64 t = FastMod((digit ? (q - digit) : 0) * inv, q, Rq);
			r += t * m;
			m = mq;
		}
	}

	std::sort(st.sparse.begin(), st.sparse.end());
	rs.base = base;

	// Band sieve (GPU version: SieveBand, DecideByBand) for the near-degenerate ranges: every number of the range has
	// T'/n2 <= e^W_Q for the part n2 of N without its prime factors <= Q (W_Q = log of the largest T/N minus the log of the abundancy of
	// the part of the known part K with primes <= Q: they divide the part of N that the survivor check removes), so D = T' - n2 - 1 gives
	// L = n2 / D > L0 = 0.999 / (e^W_Q - 1). With W_Q < 1e-4 the primes of (max(L0, Q), hi] are sieved over the range: a survivor whose
	// smallest prime factor would have to be in (L, K L + K + 2] within that band and whose bit is 0 is not amicable. (A prime of g in
	// the band divides every N: no band then.) hi: 4 L0 (the numbers of t1 have L ~ L0 and K = 3), at most 1.7 N^(1/3) + 8 (smaller D
	// is decided by DecideAtMost2Factors) and within the tables of CheckCofactor.
	rs.band = nullptr;
	if (g_RangeSieveBand)
	{
		double L0 = 0.0;
		double hi = 0.0;
		bool allDecidable = false;
		if (g_RangeSieveTestBandHi > 0.0)
		{
			L0 = g_RangeSieveTestBandL0;
			hi = g_RangeSieveTestBandHi;
		}
		else
		{
			// (KQ: the part of K with the primes <= Q, a divisor of every N of the range that CheckSurvivor removes)
			long double logKQ = logK2;
			double KQ = ldexp(1.0, static_cast<int>(v2K));
			for (unsigned int i = 0; i < M0_num_factors; ++i)
			{
				if ((M0_primes[i] & 1) && (M0_primes[i] <= tier.Q))
				{
					logKQ += LogSigmaRatio(static_cast<long double>(M0_primes[i]), gExp[i]);
					KQ *= pow(static_cast<double>(M0_primes[i]), gExp[i]);
				}
			}
			const double WQ = static_cast<double>(logHi - logKQ);
			const double Nmax = Num128ToDouble(R.A * p_hi + R.B);
			if ((WQ > 0.0) && (WQ < 1e-4))
			{
				L0 = 0.999 / expm1(WQ);
				// L > n2^(1/3) / 3 means D < 3 n2^(2/3): DecideAtMost2Factors decides those numbers without trial division (n2 <= Nmax / KQ)
				const double Lthr = cbrt(Nmax / KQ) * (1.01 / 3.0);
				allDecidable = (L0 > Lthr);
				// K for L >= L0 (K L + K + 2 bounds E for all numbers of the range)
				const num128 NmaxI = R.A * p_hi + R.B;
				unsigned long t;
				unsigned int bitsN;
				if (HighWord(NmaxI))
				{
					_BitScanReverse64(&t, HighWord(NmaxI));
					bitsN = static_cast<unsigned int>(t) + 65;
				}
				else
				{
					_BitScanReverse64(&t, LowWord(NmaxI));
					bitsN = static_cast<unsigned int>(t) + 1;
				}
				_BitScanReverse64(&t, static_cast<num64>(L0));
				double K = static_cast<double>(bitsN / static_cast<unsigned int>(t));
				if (K < 2.0)
				{
					K = 2.0;
				}
				// The numbers whose prime factors <= Q are those of K have L <= L1 = 1 / (e^W_Q(p_lo) - 1 - 1 / n2) (T / N grows with p;
				// n2 >= (A1 p_lo + B1) / 2). Any other prime power <= Q (q^2 for a prime q of K, a new prime) changes W by much more and
				// gives a larger L. Where W crosses 0 in the range (W_Q(p_lo) <= 0), L has no bound near that point.
				const long double logLo = logl(static_cast<long double>(Num128ToDouble(R.B * (p_lo + 1)))) - logl(static_cast<long double>(Num128ToDouble(R.A * p_lo + R.B)));
				const double WQlo = static_cast<double>(logLo - logKQ);
				const double den = expm1(WQlo) - 2.0 / Num128ToDouble(A1 * p_lo + B1);
				const double L1 = ((WQlo > 1e-9) && (den > 0.5 * expm1(WQlo))) ? (1.001 / den) : std::numeric_limits<double>::infinity();
				if (allDecidable)
				{
					// the band is only a faster test than DecideAtMost2Factors (see the cost check below): up to E of L1, at most 4 L0 (t1:
					// L ~ L0, K = 3)
					hi = std::min(4.0 * L0, (K * L1 + K + 2.0) * 1.0001 + 10.0);
				}
				else
				{
					// the band covers every number that DecideAtMost2Factors can't decide: L <= min(L1, Lthr). (Up to 4 L0 only, the ranges of
					// t1 where W crosses 0 left 43M numbers to a trial division up to E ~ 10^7 (t1b2: 224 s).)
					hi = (K * std::min(L1, Lthr) + K + 2.0) * 1.0001 + 10.0;
				}
			}
		}
		const double lo = std::max(L0, static_cast<double>(tier.Q));
		hi = std::min(hi, CofactorTableMaxPrime() * 0.999);
		// (shortest ranges with a band: 2^20 numbers where DecideAtMost2Factors decides all numbers (the cost check below decides there),
		// 2^14 where it can't (each number it can't decide would scan the primes up to E: 2^20 left t1b4, the end of the family
		// 2*5*7*11*23*41*2083*m, at 272.7 s, 2^18 24.1 s, 2^14 4.2 s; k4s9 2.58 -> 2.27 s; 2^0 .. 2^14 the same))
		// (at most 2^27 numbers: the band holds a byte per number prime to 30 of the range; shared ranges come in chunks of 34.4M, a range up to
		// twice that is checked whole, a longer one only when all slots for shared ranges are taken)
		const num64 minSpan = num64(1) << (allDecidable ? 20 : g_RangeSieveBandMinSpanLog2);
		bool useBand = (L0 > 0.0) && (hi > lo * 1.5) && ((span >= minSpan) || (g_RangeSieveTestBandHi > 0.0)) && (span <= (num64(1) << 27));
		for (unsigned int i = 0; useBand && (i < M0_num_factors); ++i)
		{
			const double q = static_cast<double>(M0_primes[i]);
			if ((gExp[i] > 0) && (q > lo) && (q <= hi))
			{
				useBand = false;
			}
		}
		if (useBand && allDecidable)
		{
			// Without the band, DecideAtMost2Factors decides each number (~250 cycles); the band costs ~32 cycles per prime (its root
			// and offset) and per hit, and rejects most numbers right away. With ~0.006 numbers to check per number of the range (t1):
			// the band only if (primes + hits) < 0.046 span. (Measured on t1: its chunks of M0 = 1512487130 have 0.019 (the band saves
			// ~150 cycles per number checked), the ranges of M0 ~ 7.5e14 0.081 (it costs ~680 more).)
			const double hits = static_cast<double>(span) * (8.0 / 30.0) * log(log(hi) / log(lo));
			useBand = (CofactorTablePrimeCount(lo, hi) + hits < 0.046 * static_cast<double>(span));
		}
		if (useBand)
		{
			st.bandBytes.assign(static_cast<size_t>(((p_hi - P0) / 30 + 1) * 8), 0);
			st.band.P0 = P0;
			st.band.L0 = L0;
			st.band.lo = lo;
			st.band.hi = hi;
			const unsigned int numBandPrimes = BandSieve(A1, B1, P0, p_hi, st.band, st.bandBytes.data());
			if (numBandPrimes > 0)
			{
				st.band.subBands = st.bandBytes.data();
				rs.band = &st.band;
				g_RangeSieveStats[7] += 1;
				g_RangeSieveStats[8] += numBandPrimes;
			}
		}
	}

	// the 210-wheel for long ranges (6 setups of the moduli instead of one); p = 7 is not in its layout
	const bool use210 = (g_RangeSieveLayout == 210) || ((g_RangeSieveLayout != 30) && (span >= LAYOUT210_MIN_SPAN));
	if (use210)
	{
		if ((p_lo <= 7) && (7 <= p_hi))
		{
			CheckCandidate(R, 7);
		}
		for (unsigned int s7 = 1; s7 < 7; ++s7)
		{
			SieveLayout(st, rs, s7);
		}
	}
	else
	{
		SieveLayout(st, rs, 0);
	}
	const num64 numSurvivors = rs.numSurvivors;
	const num64 numChecked = rs.numChecked;
	{
		num64 bs[3];
		CofactorBandStats(bs);
		g_RangeSieveStats[9] += bs[0];
		g_RangeSieveStats[10] += bs[1];
		g_RangeSieveStats[11] += bs[2];
	}
	{
		const num64 totalCycles = ReadCycleCounter() - startCycles;
		t_lastTier = static_cast<int>(tierIndex);
		t_lastCheckShare = (totalCycles > 0) ? (static_cast<double>(rs.checkCycles) / static_cast<double>(totalCycles)) : 0.0;
	}

	g_RangeSieveStats[0] += 1;
	g_RangeSieveStats[1] += p_hi - p_lo + 1;
	g_RangeSieveStats[2] += numSurvivors;
	g_RangeSieveStats[3] += numChecked;

	const num64 ns = static_cast<num64>(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - startTime).count());
	g_RangeSieveTimeNs += ns;
	g_RangeSieveNumbersByW[wClass] += p_hi - p_lo + 1;
	g_RangeSieveTimeByW[wClass] += ns;
	if (ns > g_slowRangeMin)
	{
		std::lock_guard<std::mutex> lock(g_slowRangeLock);
		unsigned int k = 0;
		for (unsigned int i = 1; i < NUM_SLOW_RANGES; ++i)
		{
			if (g_slowRanges[i].ns < g_slowRanges[k].ns)
			{
				k = i;
			}
		}
		if (ns > g_slowRanges[k].ns)
		{
			g_slowRanges[k].ns = ns;
			g_slowRanges[k].M0 = M0;
			g_slowRanges[k].p_lo = p_lo;
			g_slowRanges[k].p_hi = p_hi;
			g_slowRanges[k].survivors = numChecked;
			num64 m = g_slowRanges[0].ns;
			for (unsigned int i = 1; i < NUM_SLOW_RANGES; ++i)
			{
				m = std::min(m, g_slowRanges[i].ns);
			}
			g_slowRangeMin = m;
		}
	}
}

void RangeSieveLastCall(int& tier, double& checkShare)
{
	tier = t_lastTier;
	checkShare = t_lastCheckShare;
}

void RangeSievePrintStats()
{
	std::cerr << "Range sieve: " << g_RangeSieveStats[0] << " ranges, " << g_RangeSieveStats[1] << " numbers, " << g_RangeSieveStats[2] << " survivors, " << g_RangeSieveStats[3] << " with a prime p, " << g_RangeSieveStats[4] << " inconsistent sums";
#ifdef RANGE_SIEVE_CHECK_COMPARE
	std::cerr << ", " << g_RangeSieveStats[5] << " different results of the old and new checks";
#endif
	std::cerr << std::endl;
	std::cerr << "Range sieve: band sieve in " << g_RangeSieveStats[7] << " ranges (" << g_RangeSieveStats[8] << " band primes), survivors rejected by the band " << g_RangeSieveStats[9] << ", with a band prime " << g_RangeSieveStats[10] << ", band not covered " << g_RangeSieveStats[11] << std::endl;
#ifdef RANGE_SIEVE_VERIFY_BAND
	std::cerr << "Band verification: " << g_bandVerifyCount << " lookups checked (every 1024th), " << g_bandVerifyErrors << " mismatches" << std::endl;
#endif
	std::cerr << "Range sieve: " << g_RangeSieveStats[6] << " ranges skipped (no amicable numbers possible); by W (below 1e-4, 1e-3, 1e-2, 0.03, 0.1, 0.3, 1, more): numbers";
	for (unsigned int i = 0; i < NUM_W_CLASSES; ++i)
	{
		std::cerr << " " << static_cast<double>(g_RangeSieveNumbersByW[i]);
	}
	std::cerr << ", seconds";
	for (unsigned int i = 0; i < NUM_W_CLASSES; ++i)
	{
		std::cerr << " " << (static_cast<double>(g_RangeSieveTimeByW[i]) * 1e-9);
	}
	std::cerr << std::endl;
	std::cerr << "Range sieve: " << (static_cast<double>(g_RangeSieveTimeNs) * 1e-9) << " s in SearchRangeSieve (all threads), slowest ranges:";
	std::sort(g_slowRanges, g_slowRanges + NUM_SLOW_RANGES, [](const SlowRange& a, const SlowRange& b) { return a.ns > b.ns; });
	for (unsigned int i = 0; (i < NUM_SLOW_RANGES) && (g_slowRanges[i].ns > 0); ++i)
	{
		std::cerr << " " << (static_cast<double>(g_slowRanges[i].ns) * 1e-9) << " s (M0 " << g_slowRanges[i].M0 << ", p " << g_slowRanges[i].p_lo << " - " << g_slowRanges[i].p_hi << ", " << g_slowRanges[i].survivors << " checked)";
	}
	std::cerr << std::endl;
}

// ------------------------------------------------------------------------------------------------------------------------------
// Self-tests
// ------------------------------------------------------------------------------------------------------------------------------

namespace
{
	// Direct value of the sieve sum S(p) (natural log of sigma(F)/F for the part F of N(p) with prime factors <= Q), in fixed point
	// If sums != nullptr, sums[t] receives the value for tier t (numQ is ignored then)
	int DirectSum(num128 N, unsigned int numQ, int* sums)
	{
		if (sums)
		{
			numQ = g_tiers[NUM_TIERS - 1].numPrimes;
		}
		unsigned int nextTier = 0;
		int s = 0;
		int v = 0;
		while ((LowWord(N) & 1) == 0)
		{
			N >>= 1;
			++v;
		}
		s += g_T2[v];
		for (unsigned int qi = 0; qi < numQ; ++qi)
		{
			if (sums)
			{
				while ((nextTier < NUM_TIERS) && (g_tiers[nextTier].numPrimes == qi))
				{
					sums[nextTier++] = s;
				}
			}
			const num64 q = g_q[qi];
			if (Mod128(N, q, g_c64[qi], g_R[qi]) == 0)
			{
				int e = 0;
				do
				{
					N /= q;
					++e;
				} while (Mod128(N, q, g_c64[qi], g_R[qi]) == 0);
				s += T(qi, e);
			}
		}
		if (sums)
		{
			while (nextTier < NUM_TIERS)
			{
				sums[nextTier++] = s;
			}
		}
		return s;
	}
}

bool TestRangeSieve()
{
	RangeSieveInit();

	// For every known pair m = M0 * p below the search limit whose largest prime factor p has exponent 1 and is <= the main prime table
	// bound (and 6 doesn't divide m: RangeGen never generates such M0) - the pairs of the first stage:
	// 1) The window never rejects it: sigma(F)/F of the small-prime part F of n satisfies the window of every sieve depth.
	// 2) CheckSurvivor, the check of a sieve survivor, reports it at every sieve depth (with the exact sieve sum of n): without a band,
	//    with the over-abundance mask of the range (as SearchRange computes it), and with four bands around the bounds L and E of the
	//    cofactor of n (one that covers it, a tight one, one that doesn't reach E, one above L), the band byte of p set from the
	//    factorization of n (BandSieve itself is checked by the sieve runs of 3) and by RANGE_SIEVE_VERIFY_BAND).
	// 3) Full sieve runs on a short interval around p find it (sampled: each run has a setup cost), also with a band forced around L and
	//    E (tier 2). Every number they find that is not in the c2 files is printed: a new pair (to be verified) or a false report.
	const char* files[] = { "c2_3.txt", "c2_4.txt", "c2_5.txt", "c2_6.txt", "c2_7.txt", "c2_8.txt", "c2_9.txt", "c2_10.txt", "c2_11.txt", "c2_12.txt", "c2_13.txt", "c2_14.txt", "c2_15.txt", "c2_16.txt", "c2_17.txt", "c2_18.txt", "c2_19.txt", "c2_20.txt", "c2_20_1.txt", "c2_20_2.txt", "c2_21.txt", "c2_21_1.txt", "c2_21_2.txt", "c2_21_3.txt", "c2_21_4.txt", "c2_22.txt" };

	// factorization "x=p1^e1*p2*..." of a line of the c2 files, sorted (the files don't always sort them); false if malformed
	auto parseFactorization = [](const std::string& line, num64* primes, unsigned int* exps, unsigned int& nf) -> bool
	{
		nf = 0;
		const char* s = strchr(line.c_str(), '=');
		if (!s)
		{
			return false;
		}
		++s;
		while ((*s >= '0') && (*s <= '9'))
		{
			if (nf >= 64)
			{
				return false;
			}
			num64 q = 0;
			while ((*s >= '0') && (*s <= '9'))
			{
				q = q * 10 + static_cast<num64>(*s - '0');
				++s;
			}
			unsigned int e = 1;
			if (*s == '^')
			{
				++s;
				e = 0;
				while ((*s >= '0') && (*s <= '9'))
				{
					e = e * 10 + static_cast<unsigned int>(*s - '0');
					++s;
				}
			}
			primes[nf] = q;
			exps[nf] = e;
			++nf;
			if (*s == '*')
			{
				++s;
			}
			else
			{
				break;
			}
		}
		for (unsigned int i = 1; i < nf; ++i)
		{
			for (unsigned int j = i; (j > 0) && (primes[j - 1] > primes[j]); --j)
			{
				std::swap(primes[j - 1], primes[j]);
				std::swap(exps[j - 1], exps[j]);
			}
		}
		return nf > 0;
	};

	// all members of all known pairs (for the numbers found by the sieve runs)
	std::vector<num128> known;
	for (const char* name : files)
	{
		std::ifstream f(name);
		std::string line;
		while (std::getline(f, line))
		{
			const char* eq = strchr(line.c_str(), '=');
			if (eq && (line[0] >= '0') && (line[0] <= '9'))
			{
				known.push_back(atoi128(line.c_str()));
			}
		}
	}
	std::sort(known.begin(), known.end());

	num64 numFiles = 0;
	num64 numAll = 0;
	num64 numPairs = 0;
	num64 numWindowChecks = 0;
	num64 numSurvivorChecks = 0;
	num64 numMaskChecks = 0;
	num64 numMaskSkipped = 0;
	num64 numBandChecks = 0;
	num64 numSieveRuns = 0;
	num64 numBandRuns = 0;
	num64 numUnknown = 0;
	std::vector<num128> found;
	for (const char* name : files)
	{
		std::ifstream f(name);
		if (!f)
		{
			continue; // (c2_20.txt and c2_21.txt are split into parts)
		}
		++numFiles;
		std::string line1, line2, line3, line4;
		while (std::getline(f, line1))
		{
			if (line1.empty())
			{
				continue; // (extra empty lines between the records)
			}
			if (!std::getline(f, line2) || !std::getline(f, line3))
			{
				std::cerr << "TestRangeSieve: incomplete record in " << name << ": " << line1 << std::endl;
				return false;
			}
			std::getline(f, line4);

			num64 primes[64];
			unsigned int exps[64];
			unsigned int nf = 0;
			num64 nprimes[64];
			unsigned int nexps[64];
			unsigned int nnf = 0;
			if (!parseFactorization(line2, primes, exps, nf) || !parseFactorization(line3, nprimes, nexps, nnf))
			{
				std::cerr << "TestRangeSieve: bad record in " << name << ": " << line1 << " / " << line2 << " / " << line3 << std::endl;
				return false;
			}
			num128 m = atoi128(line2.c_str());
			num128 n = atoi128(line3.c_str());
			if (n < m)
			{
				std::swap(m, n);
				std::swap(primes, nprimes);
				std::swap(exps, nexps);
				std::swap(nf, nnf);
			}
			if (m >= SearchLimit::value)
			{
				continue;
			}
			++numAll;
			if ((nf < 2) || (exps[nf - 1] != 1) || (primes[nf - 1] > SearchLimit::MainPrimeTableBound) || ((m % 6) == 0))
			{
				continue; // not a first-stage pair
			}
			++numPairs;

			const num64 p = primes[nf - 1];
			num128 M0 = 1;
			num128 sumM0 = 1;
			for (unsigned int i = 0; i + 1 < nf; ++i)
			{
				num128 pk = 1;
				num128 sp = 1;
				for (unsigned int k = 0; k < exps[i]; ++k)
				{
					pk *= primes[i];
					sp += pk;
				}
				M0 *= pk;
				sumM0 *= sp;
			}
			if ((M0 * p != m) || (sumM0 * (p + 1) != m + n))
			{
				std::cerr << "TestRangeSieve: bad input for " << m << std::endl;
				return false;
			}

			// 1) window check at every tier
			++numWindowChecks;
			int sums[NUM_TIERS];
			DirectSum(n, 0, sums);
			for (unsigned int t = 0; t < NUM_TIERS; ++t)
			{
				const int S = sums[t];
				const num128 T = m + n;
				const long double tn = Num128ToDouble(T) / Num128ToDouble(n);
				unsigned int L = 0;
				{
					num128 x = n;
					while (x > 1)
					{
						x >>= 1;
						++L;
					}
				}
				const int lo_fx = static_cast<int>(floorl((logl(tn) - logl(static_cast<long double>(g_tiers[t].rho[L + 1]))) * FX)) - MARGIN;
				const int hi_fx = static_cast<int>(ceill(logl(tn) * FX)) + MARGIN;
				if ((S < lo_fx) || (S > hi_fx))
				{
					std::cerr << "TestRangeSieve: window rejects " << m << " (tier " << t << ")" << std::endl;
					return false;
				}
			}

			// 2) CheckSurvivor at every tier (p >= 7: 2, 3 and 5 are checked by CheckCandidate)
			if (p >= 7)
			{
				SieveRange R;
				R.M0 = M0;
				R.sumM0 = sumM0;
				R.A = sumM0 - M0;
				R.B = sumM0;
				R.mask = 0;

				// the mask of the range, from the factors of M0 as RangeGen sets them up (PrimeIterator needs the primes in the table)
				bool haveMask = (primes[nf - 2] < (num64(1) << 34));
				unsigned int rangeMask = 0;
				if (haveMask)
				{
					RangeData r;
					r.value = M0;
					r.sum = sumM0;
					r.last_factor_index = static_cast<int>(nf) - 2;
					r.start_prime = 0;
					r.index_start_prime = 0;
					const num64 smallPrimes[5] = { 2, 3, 5, 7, 11 };
					for (unsigned int i = 0; i + 1 < nf; ++i)
					{
						Factor& fa = r.factors[i];
						fa.p = PrimeIterator(primes[i]);
						fa.k = exps[i];
						fa.index = 1000000;
						for (int s = 0; s < 5; ++s)
						{
							if (primes[i] == smallPrimes[s])
							{
								fa.index = s;
							}
						}
						fa.p_inv = 0;
						fa.p_inv128 = -modular_inverse128(primes[i]);
					}
					rangeMask = RangeOverAbundantMask(r);
				}

				// cofactor of n and its target for each tier (the part with prime factors <= Q removed), for the bands
				num64 gPrimes[64];
				unsigned int gCount = 0;
				{
					num128 a = R.A;
					num128 b = R.B;
					while (b != 0)
					{
						const num128 tt = a % b;
						a = b;
						b = tt;
					}
					for (unsigned int i = 0; i < nnf; ++i)
					{
						if ((a % num128(nprimes[i])) == 0)
						{
							gPrimes[gCount++] = nprimes[i];
						}
					}
				}

				for (unsigned int t = 0; t < NUM_TIERS; ++t)
				{
					const unsigned int numQ = g_tiers[t].numPrimes;
					auto expectFound = [&](const char* what, const CofactorBand* band) -> bool
					{
						const unsigned int found0 = GetNumFoundPairsInThisThread();
						CheckSurvivor(R, p, sums[t], numQ, band);
						if (GetNumFoundPairsInThisThread() != found0 + 1)
						{
							std::cerr << "TestRangeSieve: CheckSurvivor missed " << m << " (tier " << t << ", " << what << ")" << std::endl;
							return false;
						}
						return true;
					};

					R.mask = 0;
					if (!expectFound("no band", nullptr))
					{
						return false;
					}
					++numSurvivorChecks;
					if (haveMask)
					{
						R.mask = rangeMask;
						if (!expectFound("mask of the range", nullptr))
						{
							return false;
						}
						R.mask = 0;
						++numMaskChecks;
					}
					else if (t == 0)
					{
						++numMaskSkipped;
					}

					// the cofactor n2 of n and its target t2 after the primes <= Q
					const num64 Q = g_tiers[t].Q;
					num128 n2 = 1;
					num128 sigmaSmall = 1;
					for (unsigned int i = 0; i < nnf; ++i)
					{
						num128 pk = 1;
						num128 sp = 1;
						for (unsigned int k = 0; k < nexps[i]; ++k)
						{
							pk *= nprimes[i];
							sp += pk;
						}
						if (nprimes[i] <= Q)
						{
							sigmaSmall *= sp;
						}
						else
						{
							n2 *= pk;
						}
					}
					const num128 t2 = (m + n) / sigmaSmall;
					double L = 0.0;
					double E = 0.0;
					if ((t2 <= n2 + 1) || !CofactorBandBounds(n2, t2, L, E))
					{
						continue; // (D = 0: no band bound)
					}
					const double variants[4][2] = { { 0.5, 1.5 }, { 1.0 - 1e-6, 1.0 + 1e-6 }, { 0.5, -1.0 }, { 1.5, 1.5 } };
					for (unsigned int v = 0; v < 4; ++v)
					{
						CofactorBand band;
						band.P0 = (p / 210) * 210;
						band.L0 = L * variants[v][0];
						band.lo = std::max(band.L0, static_cast<double>(Q));
						band.hi = (variants[v][1] > 0.0) ? (E * variants[v][1]) : ((L + E) * 0.5);
						band.hi = std::min(band.hi, CofactorTableMaxPrime() * 0.999);
						if (!CofactorBandSetup(band) || (band.sub[0] >= band.sub[BAND_SUBBANDS]))
						{
							continue;
						}
						bool gInBand = false;
						for (unsigned int i = 0; i < gCount; ++i)
						{
							if ((static_cast<double>(gPrimes[i]) > band.lo) && (static_cast<double>(gPrimes[i]) <= band.hi))
							{
								gInBand = true;
							}
						}
						if (gInBand)
						{
							continue; // (no band there: the primes of g divide every number of the range)
						}
						// the byte of p: the sub-band of the smallest prime factor of n in the band
						unsigned int sb = 0;
						for (unsigned int i = 0; i < nnf; ++i)
						{
							if ((static_cast<double>(nprimes[i]) > band.lo) && (static_cast<double>(nprimes[i]) <= band.hi))
							{
								sb = CofactorBandSubBandOf(band, nprimes[i]);
								break;
							}
						}
						static const int CLS[30] = { -1, 0, -1, -1, -1, -1, -1, 1, -1, -1, -1, 2, -1, 3, -1, -1, -1, 4, -1, 5, -1, -1, -1, 6, -1, -1, -1, -1, -1, 7 };
						const num64 d = p - band.P0;
						std::vector<unsigned char> bytes(static_cast<size_t>((d / 30 + 1) * 8), 0);
						bytes[static_cast<size_t>((d / 30) * 8 + static_cast<num64>(CLS[d % 30]))] = static_cast<unsigned char>(sb);
						band.subBands = bytes.data();
						if (!expectFound((v == 0) ? "band (L/2, 1.5 E]" : ((v == 1) ? "tight band" : ((v == 2) ? "band below E" : "band above L")), &band))
						{
							std::cerr << "  (L " << L << ", E " << E << ", band " << band.lo << " .. " << band.hi << ", sub-band " << sb << ")" << std::endl;
							return false;
						}
						++numBandChecks;
					}
				}
			}

			// 3) full sieve runs on a short interval around p (every 256th pair; "/test" initializes the prime table up to 2.5e10)
			if (((numPairs % 256) == 1) && (p < (num64(1) << 32) - (num64(1) << 16)))
			{
				const num64 W = 1 << 14;
				num128 limit = (SearchLimit::value - 1) / M0;
				if (limit > SearchLimit::MainPrimeTableBound)
				{
					limit = SearchLimit::MainPrimeTableBound;
				}
				num64 p_lo = (p > W) ? (p - W) : 2;
				if (p_lo <= primes[nf - 2])
				{
					p_lo = primes[nf - 2] + 1;
				}
				num64 p_hi = p + W;
				if (p_hi > LowWord(limit))
				{
					p_hi = LowWord(limit);
				}
				// the numbers a run finds: m, and any other must be a known pair
				auto checkFound = [&](const char* what) -> bool
				{
					bool haveM = false;
					for (const num128& x : found)
					{
						if (x == m)
						{
							haveM = true;
						}
						else if (!std::binary_search(known.begin(), known.end(), x))
						{
							char buf[64];
							std::cerr << "Range sieve test found " << itoa128(x, buf, sizeof(buf)) << " which is not in the known pairs" << std::endl;
							++numUnknown;
						}
					}
					if (!haveM)
					{
						std::cerr << "TestRangeSieve: " << what << " missed " << m << std::endl;
					}
					return haveM;
				};
				found.clear();
				g_FoundNumbersLog = &found;
				SearchRangeSieve(M0, sumM0, primes, exps, nf - 1, p_lo, p_hi, 0, RangeSieveCountPrimes(p_lo, p_hi));
				g_FoundNumbersLog = nullptr;
				if (!checkFound("sieve"))
				{
					return false;
				}
				++numSieveRuns;

				// the same run with tier 2 and a band forced around the bounds L, E of the cofactor of n (n without its prime factors
				// <= Q, as in CheckSurvivor): band primes of (max(L / 2, Q), 1.5 E], so the pair goes through the band check
				{
					const unsigned int tb = 2;
					num128 n2 = n;
					num128 t = m + n;
					while ((LowWord(n2) & 1) == 0)
					{
						n2 >>= 1;
					}
					{
						num128 x = n / n2;
						num128 s2 = 1;
						while (x > 1)
						{
							x >>= 1;
							s2 = s2 * 2 + 1;
						}
						t = t / s2;
					}
					for (unsigned int qi = 0; qi < g_tiers[tb].numPrimes; ++qi)
					{
						const num128 q = static_cast<num64>(g_q[qi]);
						if ((n2 % q) == 0)
						{
							num128 sq = 1;
							num128 pk = 1;
							while ((n2 % q) == 0)
							{
								n2 = n2 / q;
								pk = pk * q;
								sq = sq + pk;
							}
							t = t / sq;
						}
					}
					double L = 0.0;
					double E = 0.0;
					if ((t > n2 + 1) && CofactorBandBounds(n2, t, L, E) && (E * 1.5 < CofactorTableMaxPrime() * 0.99))
					{
						g_RangeSieveTierOverride = static_cast<int>(tb);
						g_RangeSieveTestBandL0 = L * 0.5;
						g_RangeSieveTestBandHi = E * 1.5;
						g_testBandP = p;
						g_testBandPasses = 0;
						const num64 bandRanges0 = g_RangeSieveStats[7];
						found.clear();
						g_FoundNumbersLog = &found;
						SearchRangeSieve(M0, sumM0, primes, exps, nf - 1, p_lo, p_hi, 0, RangeSieveCountPrimes(p_lo, p_hi));
						g_FoundNumbersLog = nullptr;
						const bool used = (g_RangeSieveStats[7] != bandRanges0);
						const num64 passes = g_testBandPasses;
						g_RangeSieveTierOverride = -1;
						g_RangeSieveTestBandL0 = 0.0;
						g_RangeSieveTestBandHi = 0.0;
						g_testBandP = 0;
						if (!checkFound("sieve with the band"))
						{
							std::cerr << "  (L " << L << ", E " << E << ")" << std::endl;
							return false;
						}
						if (used)
						{
							if (passes == 0)
							{
								std::cerr << "TestRangeSieve: " << m << " did not pass the band check (L " << L << ", E " << E << ")" << std::endl;
								return false;
							}
							++numBandRuns;
						}
					}
				}
			}
		}
	}

	std::cout << numFiles << " files, " << numAll << " pairs below the search limit, " << numPairs << " first-stage pairs: " << numWindowChecks << " window checks, " << numSurvivorChecks << " survivor checks, " << numMaskChecks << " with the mask of the range (" << numMaskSkipped << " pairs without), " << numBandChecks << " with a band, " << numSieveRuns << " sieve runs (" << numBandRuns << " through the band), " << numUnknown << " unknown numbers found...";
	return (numPairs > 0);
}
