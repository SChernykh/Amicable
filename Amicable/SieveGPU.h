#pragma once

#include <string>
#include <vector>

// GPU arithmetic-progression sieve for the first search stage (ranges of m = M0 * p, p prime, p > largest prime factor of M0).
//
// For a fixed M0, N(p) = sigma(M0) * (p + 1) - M0 * p = A * p + B with A = sigma(M0) - M0 and B = sigma(M0) is linear in p,
// so every prime power q^k divides N(p) for exactly one residue class of p modulo q^k (or for none). Instead of trial-dividing
// every N(p), the sieve adds log(sigma(q^e) / q^e) for all prime powers q^e || N(p) with q <= Q, for a whole segment of p at once.
// sigma(N) = T = sigma(M0) * (p + 1) is only possible if the accumulated value S(p) satisfies
//
//     log(T / N) - log(rho_Q) - slack <= S(p) <= log(T / N)
//
// where rho_Q is the largest possible sigma(x) / x for the part of N without prime factors <= Q, and slack covers prime powers
// that are too large to be sieved (q^k > 65535). Only primes p that pass this test (a fraction of a percent) are checked exactly,
// with the same code (CheckPairPhase1 and phases 2-3) as all numbers in the original kernel.
//
// Implementation:
// - N = g * (A1 * p + B1), g = gcd(A, B) = gcd(M0, sigma(M0)) (computed on the host from the factorization of M0, so g is known
//   exactly); the GPU sieves A1 * p + B1 and adds the constant contribution of g
// - only p coprime to 210 are sieved, in 6 sub-sieves (p mod 7 = 1 .. 6) of one row per residue class modulo 30 each: in the sub-sieve
//   s7, p = Pseg + 210 * j + w(c) with w(c) == W30[c] (mod 30), w(c) == s7 (mod 7) (8 rows); a chunk runs all its segments once per
//   sub-sieve (the segments keep the layout of 8 rows, but cover 7 times as many numbers as with the wheel of 30)
// - sums are 32-bit fixed-point natural logarithms (2^26 per unit); every increment is a difference of rounded table entries
//   T(q, e), so the sum over all powers of q telescopes to exactly T(q, v_q(N)) (rounding error <= 0.5 per prime factor)
// - q = 2: the 2-adic valuation of A1 * p + B1 is computed directly for every position (from its 32 low bits)
// - q = 3, 5: the first power is constant in a row; higher powers are "row" moduli with stride q^(k-1)
// - q = 7: the first power is constant in a sub-sieve; higher powers are "row" moduli with stride 7^(k-1) in the sub-sieve p == r (mod 7)
// - other q^k <= 65535 (q <= Q): "row" moduli (an arithmetic progression with stride q^k in each row) if q^k is small,
//   "wheel" moduli (a single state that steps through the hits coprime to 30) if q^k is large
// - one work-group processes one chunk (several consecutive segments of 8 * SegJ positions) of one range; the roots of all
//   moduli are computed by the work-group itself at the beginning of the chunk (a few hundred instructions per work-item)
// - the hits of a segment are split into tasks (one row of a modulus or a slice of it, or one wheel modulus), statically
//   assigned to work-items to balance the work
// - survivors (positions inside the window and inside [p_lo, p_hi], without prime factors below 100) are stored as
//   (p << 27) | range index and checked by CheckSieveSurvivors: Miller-Rabin test for p, then CheckPairPhase1

struct RangeData;

enum SieveGPUParams
{
	SieveMaxTiers = 4,
	SieveRangeStride = 24,       // ulongs per range record
	SieveRangeIndexBits = 27,    // survivors store (p << 27) | range index, so p < 2^37 and at most 2^27 ranges per batch
	SieveMaxModulus = 65535,     // all moduli are prime powers <= 65535 (32-bit arithmetic on the GPU)
	SieveMaxOverrides = 64,      // wheel moduli whose prime divides g (see kernel_sieve.cl)
	SieveTStride = 64,           // T(q, e) is stored for e < 64
	SieveT2Size = 128,           // T(2, e) for e < 128
	SieveMargin = 1024,          // safety margin of the window in fixed-point units (about 1.5e-5)
	SieveMaxGFactors = 16,
	SieveMaxL = 65535,           // largest progression modulus (kernel: 16-bit row offsets)
	SieveKillQMax = 16381,       // wheel primes up to this bound kill positions where q + 1 can't divide T (kernel_sieve.cl)
	SieveUndoCap = 512,          // capacity of the list of q^2 hits that undo the kill, per chunk
	SieveKill = 0x10000000,
	SieveSurvLocal = 1024,       // capacity of the local survivor list of one segment (kernel_sieve.cl), 2 lists per work-group
	SieveWheel = 210,            // positions: p coprime to 210 (p = P0 + 210 * j + w), SievePositionsPerJ of them for every j
	SievePositionsPerJ = 48,
	SieveSubSieves = 6,          // p mod 7 = 1 .. 6, 8 rows each
	SieveMaxChunkSegments = 128, // (kernel_sieve.cl: the undo list stores sub-sieve * 128 + segment)
	SieveWarpWidth = 32,         // SIMD width of the GPU (tasks are balanced per warp)
	SieveStatsHeader = 4,        // SieveLaunchDone statistics: total survivors (64-bit), 2 unused, then one overflow flag per launch
	SieveMaxLaunchesPerRound = 4096,
};

// Fixed-point scale of the sums: 2^26 per unit of natural logarithm
#define SIEVE_FX 67108864.0

struct SieveGPUConfig
{
	SieveGPUConfig();

	unsigned int wgSize;          // work-group size of SieveRanges
	unsigned int segJ;            // positions per row in one segment (power of 2)
	unsigned int chunkSegments;   // segments per chunk (one work-group)
	unsigned int hitsPerTask;     // target number of hits per row task and segment
	unsigned int wheelThreshold;  // moduli >= this value (and coprime to 30) use wheel stepping, 0 = automatic (segJ * 3 / 4)
	unsigned int numTiers;
	unsigned int tierQ[SieveMaxTiers];            // largest sieving prime of each tier (increasing)
	num64 tierMinPrimes[SieveMaxTiers];           // tier t is used for ranges with at least this many (estimated) primes
	double highTierW;             // see ChooseTier
	unsigned int killQMax;        // divisibility kill of the wheel primes up to this bound (0 = off, the kernel code is left out)
	unsigned int survLocal;       // capacity of the local survivor list of one segment
};

#include "RangeReduction.h"
typedef RangeReduction SieveReduction;

class SieveGPUTables
{
public:
	// Computes the reduction for a range (see SieveReduction)
	void ComputeReduction(num128 M0, num128 sumM0, const num64* primes, const unsigned int* exponents, unsigned int numFactors, num64 p_lo, num64 p_hi, num64 maxL, SieveReduction& r) const;

	// Must be called after PrimeTablesInit (needs primes up to the largest tierQ). Returns false if the configuration is invalid.
	bool Init(const SieveGPUConfig& config);

	// Chooses the tier for a range with the given (estimated) number of primes. With W >= 0 (height of the window above the known
	// part, see RangeReduction), the last tier (largest Q) is only used for ranges with W < cfg.highTierW: they have many
	// survivors (the window contains every N whose cofactor has no small prime factors), which a larger Q reduces.
	unsigned int ChooseTier(num64 numPrimes, double W = -1.0) const;

	// Fills a range record (SieveRangeStride ulongs) for M0 given by its factorization (primes don't have to be sorted), for
	// p in [p_lo, p_hi]. Returns false if the range can't be sieved on the GPU (p_hi >= 2^37).
	// If adaptiveTier is true, a higher tier than "tier" is used when the constant part of the sum (from g and the 2-adic part)
	// is already inside the window of this tier: then every p where A1 * p + B1 has no small prime factors would survive.
	bool PrepareRange(num128 M0, num128 sumM0, const num64* primes, const unsigned int* exponents, unsigned int numFactors, num64 p_lo, num64 p_hi, unsigned int tier, bool adaptiveTier, num64* out) const;
	bool PrepareRange(const RangeData& r, num64 p_lo, num64 p_hi, unsigned int tier, bool adaptiveTier, num64* out) const;

	// Exact sieve sum S(p) for the self-test (same fixed-point values as the GPU), and the window of a range record for p
	int DirectSum(num128 N, unsigned int tier) const;
	// The sum the GPU computes for position p of a range record (prime powers above SieveMaxModulus are truncated)
	int ExpectedSum(const num64* rangeRecord, num64 p) const;
	static bool InsideWindow(const num64* rangeRecord, num64 p, int S);

	// Analysis: does every prime power q^v || N (q = 2 or odd q <= Q) satisfy sigma(q^v) | T = sigma(M0) * (p + 1)?
	bool DivisibilityOK(num128 M0, num128 sumM0, num64 p, unsigned int Q) const;
	// Restricted variants that a sieve can apply at its hits: mode 1 = odd q <= maxQ not dividing g, exponent exactly 1 in N;
	// mode 2 = mode 1 and the power of 2 (exact)
	bool DivisibilityOK2(num128 M0, num128 sumM0, num128 g, num64 p, unsigned int maxQ, int mode) const;

	// ln(rho) + slack + margin in fixed point for every tier and bit length of N (128 entries per tier), for the large primes
	// stage where the GPU prepares the ranges itself
	void LoExtraTable(std::vector<int>& out) const;

	// Estimated number of primes in [lo, hi]
	static double EstimatePrimes(num64 lo, num64 hi);

	// Number of positions (p coprime to 30) of a range
	num64 NumPositions(num64 p_lo, num64 p_hi) const;

	// Number of chunks of a range, and the chunk descriptor (range index, first j, number of segments)
	num64 NumChunks(num64 p_lo, num64 p_hi) const;

	SieveGPUConfig cfg;

	// Data for the GPU buffers (see kernel_sieve.cl)
	std::vector<unsigned int> primeRecs;   // uint4: q, R(q), first modulus | number of moduli << 16, offset of T(q, 0) in Ttab
	std::vector<unsigned int> modRecs;     // uint4: m, R(m), stride, inverse
	std::vector<unsigned int> modRecs2;    // uint4: (m / q)^-1 mod 2^32, INC index (row moduli), type | class << 8, R(stride)
	std::vector<int> Ttab;                 // T(2, e) for e < 128, then T(q, e) for e < 64 for each sieving prime
	std::vector<unsigned int> tasks;       // uint2 task descriptors, per tier [slot * wgSize + lid]
	std::vector<unsigned int> tierInfo;    // 2 x uint4 per tier: (number of primes, number of moduli, number of slots, offset of
	                                       // tasks), (offset of modTaskIdx, offset of modTaskList, number of row moduli, 0)
	std::vector<unsigned int> wheelTabs;   // 6 x 64 entries (DK | CN << 8 | E << 16, one table per sub-sieve), then 240 entries (T0 | T0C << 8)
	std::vector<unsigned int> modTaskIdx;  // uint2 per tier and modulus: first entry and number of entries in modTaskList
	std::vector<unsigned int> modTaskList; // per tier: slot << 16 | lane of every task of every modulus

	unsigned int maxSlots;
	unsigned int maxSlotsMain;              // over all tiers but the last (the program without the kill, kernel_sieve.cl)
	unsigned int maxModuli;
	unsigned int maxRowModuli;
	unsigned int localMemSize;
	unsigned int killQMax = 0;              // = cfg.killQMax

private:
	std::vector<unsigned int> myPrimes;     // odd sieving primes 3, 5, 7, ..., largest tierQ
	double myLnRho[SieveMaxTiers][128];     // ln(rho) for numbers < 2^L
	double myMaxTrunc;                      // largest ln(sigma(q^inf) / q^inf) - ln(sigma(q^K) / q^K) over sieved q
	std::vector<unsigned int> myMaxPower;   // K_q: largest sieved power of each prime
public:
	double myMaxWarpLoad = 0.0;             // estimated sieve work per segment of the busiest warp and of an average warp (last tier)
	double myAvgWarpLoad = 0.0;
};

// Known amicable pairs of the first stage (m = M0 * p, p = largest prime factor of m with exponent 1, p <= MainPrimeTableBound,
// m not divisible by 6) below SearchLimit::value, read from the c2_*.txt files
struct KnownPair
{
	num128 m;
	num128 n;
	num128 M0;
	num128 sumM0;
	num64 p;
	unsigned int numFactors;     // number of prime factors of M0
	num64 primes[32];            // prime factors of M0 (increasing)
	unsigned int exponents[32];
};

// Keeps every stride-th first-stage pair; allAmicable (if not null) receives both members of all known pairs (sorted)
bool LoadKnownPairs(unsigned int stride, std::vector<KnownPair>& pairs, std::vector<num128>* allAmicable);
// The first-stage pair (largePrimes: the pair of the large primes stage, largest prime factor > LinearLimit) of the two lines
// "m=factorization", "n=factorization" of a record of the c2 files; false if it is none (bad: the factorization doesn't match)
bool ParseKnownPair(const std::string& line2, const std::string& line3, KnownPair& pair, bool& bad, bool largePrimes = false);

// Self-test (used by "/test"): every given pair must be inside the window of every tier
bool TestSieveGPUWindow(const SieveGPUTables& t, const std::vector<KnownPair>& pairs);

// Self-test: every given pair satisfies the necessary conditions of SieveReduction
bool TestSieveReduction(const SieveGPUTables& t, const std::vector<KnownPair>& pairs);
