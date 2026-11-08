#pragma once

// OpenCL kernels for the arithmetic-progression sieve (see SieveGPU.h for the algorithm and the data layouts).
// This source is compiled in the same program as kernel.cl, after it (CheckSieveSurvivors calls CheckPairPhase1).
// MSVC limits a single string literal to 16380 characters, so the source is split into several literals.

static const char* kernel_sieve_cl = R"===(
#ifdef SIEVE_WG_SIZE
// Build options set by the host (the sieve kernels are only compiled when the sieve is used):
//   SIEVE_WG_SIZE      work-group size of SieveRanges
//   SIEVE_SEG_J        positions per row in one segment (power of 2), SIEVE_SEG_SHIFT = log2(SIEVE_SEG_J)
//   SIEVE_MAX_SLOTS    maximum number of sieve tasks per work-item (over all tiers)
//   SIEVE_MAX_ROW_MODULI maximum number of row moduli (over all tiers)
//   SIEVE_RANGE_STRIDE number of ulongs per range record
//   SIEVE_MARGIN       safety margin of the window in fixed-point units
//   SIEVE_SURV_LOCAL   capacity of the local survivor list of one segment

// Positions are the p coprime to 210, in 6 sub-sieves (p mod 7 = 1 .. 6) of 8 rows (p mod 30) each: row c of the sub-sieve s7 holds
// p = Pseg + 210 * j + w(c) for 0 <= j < SIEVE_SEG_J, w(c) = W30[c] + 30 * tau(c) == s7 (mod 7). A chunk runs its segments once per
// sub-sieve, so a segment has the layout of the 30-wheel (8 rows of SIEVE_SEG_J positions), but it covers 7 times as many numbers.
#define SIEVE_ROWS 8
#define SIEVE_SUBS 6
#define SIEVE_SEG_POS (SIEVE_ROWS * SIEVE_SEG_J)
#define SIEVE_POS_PER_ITEM (SIEVE_SEG_POS / SIEVE_WG_SIZE)
#define SIEVE_NEVER 0xFFFFFFFFu
#define SIEVE_ROW_NEVER 0xFFFFu
#define SIEVE_MAX_OVERRIDES 64
#define SIEVE_TASK_EMPTY 0xFFFFFFFFu
#define SIEVE_TASK_WHEEL 0x8000u
#define SIEVE_GEO_WHEEL 0x80000000u
// Divisibility kill (see SieveSetupModulus): positions where a wheel prime q || N has sigma(q) = q + 1 not dividing T get +SIEVE_KILL
#define SIEVE_KILL 0x10000000
// sums and survivor lists; during the setup the per-modulus tables (the host checks that they fit)
#define SIEVE_SHARED_WORDS (SIEVE_SEG_POS + 2 * SIEVE_SURV_LOCAL)
#define SIEVE_J0_STRIDE (SIEVE_ROWS + 1)
// steps of 4 positions per work-item in the scan of a segment (at most 8: the candidate bits are one uint)
#define SIEVE_SCAN_ITERS (((SIEVE_SEG_POS / (SIEVE_WG_SIZE * 4)) > 0) ? (SIEVE_SEG_POS / (SIEVE_WG_SIZE * 4)) : 1)
#define SIEVE_FX 67108864.0f
#define SIEVE_T_STRIDE 64
#define SIEVE_TBC_STRIDE 33
#define SIEVE_RANGE_INDEX_BITS 27

__constant uint SieveW30[8] = { 1, 7, 11, 13, 17, 19, 23, 29 };
__constant uchar SieveCls30[30] = { 255, 0, 255, 255, 255, 255, 255, 1, 255, 255, 255, 2, 255, 3, 255, 255, 255, 4, 255, 5, 255, 255, 255, 6, 255, 255, 255, 255, 255, 7 };

// (hi:lo) mod m for m < 2^32 (slow, only used in the chunk setup)
static inline ulong SieveMod128_64(const ulong lo, const ulong hi, const ulong m)
{
	const ulong c64 = (0UL - m) % m; // 2^64 mod m
	const ulong x = ((hi % m) * c64) % m;
	return (x + lo % m) % m;
}


#ifdef SIEVE_WARP_TIMING
// Development aid (NVIDIA only, profiling harness): cycles of the sieve and scan phases of every warp, summed over all work-groups
// into the last 2 * (SIEVE_WG_SIZE / 32) ulongs of the survivor buffer
static inline ulong SieveClock()
{
	ulong t;
	asm volatile("mov.u64 %0, %%clock64;" : "=l"(t));
	return t;
}
#endif

// W30[c] without a memory access
#define SIEVE_W30(c) ((uint)(0x1D1713110D0B0701UL >> ((c) << 3)) & 0xFFu)
// x^-1 mod 7 for x = 1 .. 6
#define SIEVE_INV7(x) ((0x6325410u >> ((x) << 2)) & 0xFu)
// class of 7 * W30[u] mod 30 (the wheel tasks step by 7 m)
#define SIEVE_U7(u) ((0xC9F129u >> ((u) * 3)) & 7u)
// w(c) of a sub-sieve from its offsets packed into a ulong (8 bits per row)
#define SIEVE_WSUB(wpack, c) ((uint)((wpack) >> ((c) << 3)) & 0xFFu)
// uints per entry of the band table (SieveBand)
#define SIEVE_BAND_ENTRY 8

// x mod m for 32-bit x, R = floor((2^32 - 1) / m), m < 2^16
static inline uint SieveMod32(const uint x, const uint m, const uint R)
{
	const uint r = x - mul_hi(x, R) * m;
	return (r >= m) ? (r - m) : r;
}

// a^e mod m, a < m < 2^16
static inline uint SievePowMod(uint a, uint e, const uint m, const uint R)
{
	uint result = 1;
	while (e)
	{
		if (e & 1)
		{
			result = SieveMod32(result * a, m, R);
		}
		a = SieveMod32(a * a, m, R);
		e >>= 1;
	}
	return result;
}

// 2^32 mod m
static inline uint SieveC32(const uint m, const uint R)
{
	return SieveMod32(0u - m, m, R);
}

// (hi:lo) mod m for a 128-bit number, m < 2^16 (all intermediate values stay below m * m < 2^32)
static inline uint SieveMod128(const ulong lo, const ulong hi, const uint m, const uint R)
{
	const uint c32 = SieveC32(m, R);
	uint x = SieveMod32((uint)(hi >> 32), m, R);
	x = SieveMod32(x * c32 + SieveMod32((uint)hi, m, R), m, R);
	x = SieveMod32(x * c32 + SieveMod32((uint)(lo >> 32), m, R), m, R);
	x = SieveMod32(x * c32 + SieveMod32((uint)lo, m, R), m, R);
	return x;
}

// D(q) = (q + 1) / gcd(q + 1, sigma(M0)) for a prime q: if q || N, then q + 1 | T = sigma(M0) * (p + 1) needs D(q) | p + 1 (0: D(q) < 7,
// no kill). The kill of a position p (SieveRanges, the undo list in SieveSetupModulus, ExpectedSum on the host): D(q) does not divide
// p + 1, and if 7 | D(q), the sub-sieve p == 6 (mod 7) kills only if D(q) >= 49.
static inline uint SieveKillD(const uint q, const ulong sumM0lo, const ulong sumM0hi)
{
	const uint m = q + 1;
	uint a = m;
	uint b = SieveMod128(sumM0lo, sumM0hi, m, 0xFFFFFFFFu / m);
	while (b)
	{
		const uint t = a % b;
		a = b;
		b = t;
	}
	const uint D = m / a;
	return (D >= 7) ? D : 0;
}

// Writes the n survivors of a segment (positions c << SIEVE_SEG_SHIFT | j in sb, row offsets w) to the global list at index base
static inline void SieveFlushSurvivors(__local const uint* sb, const uint n, const uint base, const ulong Pseg, __local const uint* w, const uint rangeIndex, __global ulong* survivors, const uint survivorCapacity, const uint lid)
{
	for (uint i = lid; i < n; i += SIEVE_WG_SIZE)
	{
		const uint pos = sb[i];
		const ulong p = Pseg + 210UL * (pos & (SIEVE_SEG_J - 1)) + w[pos >> SIEVE_SEG_SHIFT];
		const uint k = base + i;
		if (k < survivorCapacity)
		{
#ifdef SIEVE_LARGE_PRIMES
			survivors[k * 2] = p;
			survivors[k * 2 + 1] = rangeIndex;
#else
			survivors[k] = (p << SIEVE_RANGE_INDEX_BITS) | rangeIndex;
#endif
		}
	}
}

// Exponent of the sieving prime q in g = gcd(A, B): from the factorization of g (ranges), or by division (large primes, g < 2^40)
static uint SieveExponentInG(const uint q, const uint Rq, const uint ngf, __local const uint* GF, const ulong gOdd)
{
	uint eg = 0;
	for (uint f = 0; f < ngf; ++f)
	{
		const uint gf = GF[f];
		if ((gf & 0xFFFFFFu) == q)
		{
			eg = gf >> 24;
		}
	}
	if (gOdd > 1)
	{
		ulong g = gOdd;
		const uint c32q = SieveC32(q, Rq);
		while (SieveMod32(SieveMod32((uint)(g >> 32), q, Rq) * c32q + SieveMod32((uint)g, q, Rq), q, Rq) == 0)
		{
			g /= q;
			++eg;
		}
	}
	return eg;
}

// The parts of the setup that belong to a sieving prime: the constant T(q, e_g) of q in g, for q = 3, 5 the first power, which is the
// same in a whole row, and for q = 7 the first power, which is the same in a whole sub-sieve (their higher powers are moduli of
// SieveSetupModulus).
static void SieveSetupPrimeConst(
	const uint pi,
	__global const uint4* primeRecs,
	__global const int* Ttab,
	const ulong A1lo, const ulong A1hi, const ulong B1lo, const ulong B1hi,
	const uint ngf,
	__local const uint* GF,
	const ulong gOdd,
	__local int* BASE,
	__local uint* SUB7)
{
	const uint4 pr = primeRecs[pi]; // q, R(q), first modulus | number of moduli << 16, offset of T(q, 0) in Ttab
	const uint q = pr.x;
	const uint Rq = pr.y;
	__global const int* Tq = Ttab + pr.w;
	const uint eg = SieveExponentInG(q, Rq, ngf, GF, gOdd);
	if (eg)
	{
		atomic_add(BASE + SIEVE_ROWS, Tq[eg]);
	}
	if (q <= 5)
	{
		// q = 3, 5: p mod q is the same for the whole row, so q | A1 * p + B1 holds for complete rows
		const uint a = SieveMod128(A1lo, A1hi, q, Rq);
		const uint b = SieveMod128(B1lo, B1hi, q, Rq);
		const uint inv = a ? SievePowMod(a, q - 2, q, Rq) : 0;
		const uint r = SieveMod32((b ? (q - b) : 0) * inv, q, Rq);
		if (a && r)
		{
			const int add = Tq[eg + 1] - Tq[eg];
			for (uint c = 0; c < SIEVE_ROWS; ++c)
			{
				if (SIEVE_W30(c) % q == r)
				{
					atomic_add(BASE + c, add);
				}
			}
		}
	}
	else if (q == 7)
	{
		// 7 | A1 * p + B1 for all positions of the sub-sieve p == r (mod 7) (r = 0: for none, p is coprime to 7)
		const uint a = SieveMod128(A1lo, A1hi, 7, Rq);
		const uint b = SieveMod128(B1lo, B1hi, 7, Rq);
		const uint r = SieveMod32((b ? (7 - b) : 0) * SievePowMod(a, 5, 7, Rq), 7, Rq);
		if (a && r)
		{
			SUB7[0] = r;
			SUB7[1] = (uint)(Tq[eg + 1] - Tq[eg]);
		}
	}
}

// Computes the root of one modulus m = q^k for the current chunk: for a wheel modulus mi, STW[mi] = x, the root is p == Pchunk + x
// (mod m) (SIEVE_NEVER: no hits); for a row modulus with row index ri, J0TAB[ri * 9 + c] = first hit j in row c of the 30-wheel,
// p = Pchunk + 30 j + W30[c] (SIEVE_ROW_NEVER: none) and J0TAB[ri * 9 + 8] = its stride in that row. The work-items then compute the
// roots of their own tasks from these tables (a row modulus has up to 8 * 255 tasks), and every sub-sieve the states from the roots.
// One work-item per modulus: the root of A1 * p + B1 == 0 (mod m) is -B1 * A1^(phi(m) - 1) (Euler). (With one work-item per prime that
// lifted the root through all powers of q and computed their rows, the work-items of the small primes - 3^2 .. 3^10, 5^2 .. 5^6, 7 .. 7^5
// with 8 rows each, dependent loads of the records in every step - made all others wait: k9, small chunks: 36% of the kernel time.)
static void SieveSetupModulus(
	const uint mi,
	__global const uint4* primeRecs,
	__global const uint4* modRecs,
	__global const uint4* modRecs2,
	__global const int* Ttab,
	const ulong A1lo, const ulong A1hi, const ulong B1lo, const ulong B1hi,
	const ulong Pchunk,
	const uint ngf,
	__local const uint* GF,
	const ulong gOdd,
	__local uint* STW,
	__local uint* J0TAB,
	__local int* INC,
	__local uint* OVR,
	const ulong sumM0lo, const ulong sumM0hi,
	const uint nseg,
	const uint jLast,
	const bool killTier,
	__local uint* UNDO,
	__local uint* UNDON)
{
	const uint4 mr = modRecs[mi];   // m, R(m), stride, inverse
	// q, row moduli: INC index | increment (q not dividing g) << 8 (wheel moduli: all ones), type | class << 8 | prime index << 16 | k << 28,
	// q > 5: R(q), q = 3, 5: R(stride). (The prime record and T(q, e) are loaded only if q divides g: dependent loads were the critical
	// path of the setup.)
	const uint4 mr2 = modRecs2[mi];
	const uint m = mr.x;
	const uint R = mr.y;
	const uint q = mr2.x;
	const uint Rq = (q > 5) ? mr2.w : ((q == 3) ? 0x55555555u : 0x33333333u);
	const uint k = mr2.z >> 28;     // m = q^k
	const uint eg = SieveExponentInG(q, Rq, ngf, GF, gOdd);

	const uint A1m = SieveMod128(A1lo, A1hi, m, R);
	const uint B1m = SieveMod128(B1lo, B1hi, m, R);
	const uint a = SieveMod32(A1m, q, Rq);

	const bool wheel = (mr2.z & 1) != 0;
	if (eg)
	{
		const uint pi = (mr2.z >> 16) & 0xFFFu;
		__global const int* Tq = Ttab + primeRecs[pi].w;
		const int inc = Tq[eg + k] - Tq[eg + k - 1];
		if (!wheel)
		{
			INC[mr2.y & 0xFFu] = inc;
		}
		else if (a)
		{
			// wheel tasks have the increment for q not dividing g in their descriptor, override it
			const uint o = atomic_inc(OVR);
			if (o < SIEVE_MAX_OVERRIDES)
			{
				OVR[1 + o * 2] = mi;
				OVR[2 + o * 2] = (uint)inc;
			}
		}
	}
	else if (!wheel)
	{
		INC[mr2.y & 0xFFu] = (int)(mr2.y >> 8);
	}

	if (!a)
	{
		// q divides A1 (but not B1): q never divides A1 * p + B1
		if (wheel)
		{
			STW[mi] = SIEVE_NEVER;
		}
		else
		{
			for (uint c = 0; c < SIEVE_ROWS; ++c)
			{
				J0TAB[(mr2.y & 0xFFu) * SIEVE_J0_STRIDE + c] = SIEVE_ROW_NEVER;
			}
		}
		return;
	}

	// root: p == r (mod m)
	const uint inv = SievePowMod(A1m, (m / q) * (q - 1) - 1, m, R);
	const uint r = SieveMod32((B1m ? (m - B1m) : 0) * inv, m, R);

	// p == r (mod m)  <=>  p - Pchunk == x (mod m)
	const uint c32 = SieveC32(m, R);
	const uint Pm = SieveMod32(SieveMod32((uint)(Pchunk >> 32), m, R) * c32 + SieveMod32((uint)Pchunk, m, R), m, R);
	uint x = r + m - Pm;
	if (x >= m)
	{
		x -= m;
	}

	if (killTier && wheel && (k == 1) && !eg && (q <= SIEVE_KILL_QMAX))
	{
		// Divisibility kill: the wheel task adds SIEVE_KILL where D(q) does not divide p + 1 (the exponent of q in N would be 1).
		// Positions with q^2 | N (q^2 is not sieved) get the kill back from the list UNDO, made here: sub-sieve * 128 + segment << 16 |
		// c << SIEVE_SEG_SHIFT | j. (m = q: inv and r are mod q; the positions up to p_hi: the last segment ends at jLast)
		const uint D = SieveKillD(q, sumM0lo, sumM0hi);
		if (D)
		{
			const ulong q2 = (ulong)q * q;
			const ulong A1m2 = SieveMod128_64(A1lo, A1hi, q2);
			const ulong B1m2 = SieveMod128_64(B1lo, B1hi, q2);
			const ulong val = (A1m2 * r + B1m2) % q2; // divisible by q
			const uint digit = (uint)((val / q) % q);
			const uint t2 = SieveMod32((digit ? (q - digit) : 0) * inv, q, Rq);
			const ulong r2 = r + (ulong)t2 * q;
			const ulong span = 210UL * (((ulong)(nseg - 1) << SIEVE_SEG_SHIFT) + jLast);
			const ulong P1 = (Pchunk + 1) % D;
			for (ulong off = (r2 + q2 - Pchunk % q2) % q2; off < span; off += q2)
			{
				const uint cls = SieveCls30[off % 30];
				const uint s7 = (uint)(off % 7); // p mod 7 (Pchunk == 0 (mod 210))
				if ((cls != 255) && (s7 != 0) && ((s7 != 6) || (D % 7) || (D >= 49)) && ((P1 + off) % D))
				{
					const uint tau = ((s7 + 7 - SIEVE_W30(cls) % 7) * 4) % 7;
					const ulong jg = (off - SIEVE_W30(cls) - 30 * tau) / 210;
					const uint k2 = atomic_inc(UNDON);
					if (k2 < SIEVE_UNDO_CAP)
					{
						UNDO[k2] = (((s7 - 1) * 128 + (uint)(jg >> SIEVE_SEG_SHIFT)) << 16) | (cls << SIEVE_SEG_SHIFT) | (uint)(jg & (SIEVE_SEG_J - 1));
					}
				}
			}
		}
	}

	if (wheel)
	{
		// (the first hit depends on the sub-sieve)
		STW[mi] = x;
		return;
	}

	// row modulus: first position in each row, 30 j == x - W30[c] (mod m)
	const uint ri = mr2.y & 0xFFu;
	J0TAB[ri * SIEVE_J0_STRIDE + SIEVE_ROWS] = mr.z;
	for (uint c = 0; c < SIEVE_ROWS; ++c)
	{
		const uint d = SieveMod32(x + m - SieveMod32(SIEVE_W30(c), m, R), m, R);
		uint j0;
		if (q <= 5)
		{
			// (30 / q) j == (x - W30[c]) / q (mod q^(k-1)), only rows where q | x - W30[c]
			const uint dq = (q == 3) ? (d / 3) : (d / 5);
			j0 = (dq * q == d) ? SieveMod32(dq * mr.w, mr.z, mr2.w) : SIEVE_ROW_NEVER;
		}
		else
		{
			j0 = SieveMod32(d * mr.w, m, R);
		}
		J0TAB[ri * SIEVE_J0_STRIDE + c] = j0;
	}
}
)==="
R"===(
// Scan of a segment (SieveRanges) in two passes: (1) every work-item reads SIEVE_SCAN_ITERS groups of 4 consecutive positions of one row
// (every row is a multiple of 4 * 32 positions, so the row is the same for a whole warp), clears them and marks the ones inside the
// window [lo, lo + width] (bit it * 4 + k of the result); (2) the caller records the marked positions (rare). Without the survivor code
// inside the first pass, its loads are not serialized. partial: the last segment of a chunk, the positions from jEnd on were not sieved
// (their sums are 0) and are not scanned (the full segments, with partial = false, don't test it).
static inline uint SieveScan(__local int* S, __local const int* TBC, __local const uint* subW, const uint lid, const uint twoAdic, const uint A1lo32, const uint B1lo32, const uint Pseg32, const int lo, const uint width, const uint jEnd, const bool partial
#ifdef SIEVE_DUMP
	, __global int* dump, const int c7
#endif
	)
{
	const uint dn = 210u * A1lo32;
	uint cand = 0;
	#pragma unroll
	for (uint it = 0; it < SIEVE_SCAN_ITERS; ++it)
	{
		const uint i = (it * SIEVE_WG_SIZE + lid) * 4;
		if (i >= SIEVE_SEG_POS)
		{
			// (segments with fewer than 4 positions per work-item)
			break;
		}
		const uint c = i >> SIEVE_SEG_SHIFT;
		const uint j = i & (SIEVE_SEG_J - 1);
		if (partial && (j >= jEnd))
		{
			continue;
		}
		__local int4* p4 = (__local int4*)(S + i);
		int4 sum = *p4;
		*p4 = (int4)(0);
		__local const int* tbc = TBC + c * SIEVE_TBC_STRIDE;
		if (twoAdic)
		{
			uint n = A1lo32 * (Pseg32 + 210u * j + subW[c]) + B1lo32;
			sum.x += tbc[popcount((n & (0u - n)) - 1u)];
			n += dn;
			sum.y += tbc[popcount((n & (0u - n)) - 1u)];
			n += dn;
			sum.z += tbc[popcount((n & (0u - n)) - 1u)];
			n += dn;
			sum.w += tbc[popcount((n & (0u - n)) - 1u)];
		}
		else
		{
			sum += (int4)(tbc[0]);
		}
#ifdef SIEVE_DUMP
		// self-test: store all sums (with the power of 7 of the sub-sieve), don't record survivors
		*(__global int4*)(dump + i) = sum + (int4)(c7);
#else
		// vector comparisons give -1 for true
		const int4 inside = as_uint4(sum - (int4)(lo)) <= (uint4)(width);
		cand |= (((uint)inside.x & 1u) | (((uint)inside.y & 1u) << 1) | (((uint)inside.z & 1u) << 2) | (((uint)inside.w & 1u) << 3)) << (it * 4);
#endif
	}
	return cand;
}
)==="
R"===(
// p with a prime factor 11 <= q <= 127 (q | p <=> p * q^-1 mod 2^64 <= (2^64 - 1) / q; the positions are prime to 210): p is composite,
// the position is not a number of the range. CheckSieveSurvivors rejects these positions before it reads their band word (the first pass,
// p not known to be prime), so SieveBand skips their hits (v73: half of its atomic_min; a wrong word of a composite p could only reject
// that position, never report it). (The same test in the scan of SieveRanges, to record half of the survivors of the ranges with
// W ~ 0: the survivor checks 1.53 -> 0.84 s on x1, but SieveRanges 0.69 -> 1.21 s, k1 3.05 -> 3.45 s: not used.)
__constant ulong SieveBandSmallInv[27] = { 0x2E8BA2E8BA2E8BA3UL, 0x4EC4EC4EC4EC4EC5UL, 0xF0F0F0F0F0F0F0F1UL, 0x86BCA1AF286BCA1BUL, 0xD37A6F4DE9BD37A7UL, 0x34F72C234F72C235UL, 0xEF7BDEF7BDEF7BDFUL, 0x14C1BACF914C1BADUL, 0x8F9C18F9C18F9C19UL, 0x82FA0BE82FA0BE83UL, 0x51B3BEA3677D46CFUL, 0x21CFB2B78C13521DUL, 0xCBEEA4E1A08AD8F3UL, 0x4FBCDA3AC10C9715UL, 0xF0B7672A07A44C6BUL, 0x193D4BB7E327A977UL, 0x7E3F1F8FC7E3F1F9UL, 0x9B8B577E613716AFUL, 0xA3784A062B2E43DBUL, 0xF47E8FD1FA3F47E9UL, 0xA3A0FD5C5F02A3A1UL, 0x3A4C0A237C32B16DUL, 0xDAB7EC1DD3431B57UL, 0x77A04C8F8D28AC43UL, 0xA6C0964FDA6C0965UL, 0x90FDBC090FDBC091UL, 0x7EFDFBF7EFDFBF7FUL };
__constant ulong SieveBandSmallLim[27] = { 0x1745D1745D1745D1UL, 0x13B13B13B13B13B1UL, 0x0F0F0F0F0F0F0F0FUL, 0x0D79435E50D79435UL, 0x0B21642C8590B216UL, 0x08D3DCB08D3DCB08UL, 0x0842108421084210UL, 0x06EB3E45306EB3E4UL, 0x063E7063E7063E70UL, 0x05F417D05F417D05UL, 0x0572620AE4C415C9UL, 0x04D4873ECADE304DUL, 0x0456C797DD49C341UL, 0x04325C53EF368EB0UL, 0x03D226357E16ECE5UL, 0x039B0AD12073615AUL, 0x0381C0E070381C0EUL, 0x033D91D2A2067B23UL, 0x03159721ED7E7534UL, 0x02E05C0B81702E05UL, 0x02A3A0FD5C5F02A3UL, 0x0288DF0CAC5B3F5DUL, 0x027C45979C95204FUL, 0x02647C69456217ECUL, 0x02593F69B02593F6UL, 0x0243F6F0243F6F02UL, 0x0204081020408102UL };

static inline bool BandSmallFactor(const ulong p)
{
	bool f = false;
	#pragma unroll
	for (int k = 0; k < 27; ++k)
	{
		f |= (p * SieveBandSmallInv[k] <= SieveBandSmallLim[k]);
	}
	return f && (p > 127);
}

// One work-group sieves one chunk (up to SIEVE_CHUNK_SEGMENTS segments) of one range: all p = Pchunk + 210 * j + w(c) coprime to 210
// with 0 <= j < nseg * SIEVE_SEG_J, in the 6 sub-sieves one after the other, and records the positions whose sum is inside the window
// (and p_lo <= p <= p_hi) as survivors: (p << 27) | range index
__kernel
__attribute__((reqd_work_group_size(SIEVE_WG_SIZE, 1, 1)))
void SieveRanges(
#ifdef SIEVE_LARGE_PRIMES
	__global const uchar* candidates,  // AmicableCandidate records (10 bytes)
	const ulong4 lpParams,             // P1, P2, P0 = P1 rounded down to a multiple of 210, chunks per candidate | tier << 32 | segments per chunk << 40
	const float limitF,                // SearchLimit::value as float
	__global const int* loExtraTab,    // per tier and bit length of N: ln(rho) + slack + margin in fixed point
	__global const uint* primeBits,    // prime p = P0 + 30 * j + W30[c] in [P1, P2]: bit j of row c (wordsPerRow words per row)
	const uint wordsPerRow,
#else
	__global const ulong* ranges,
	__global const uint4* chunks,
#endif
	const uint chunkOffset,
	__global const uint4* primeRecs,
	__global const uint4* modRecs,
	__global const uint4* modRecs2,
	__global const int* Ttab,
	__global const uint2* tasks,
	__global const uint4* tierInfo,
	__global const uint* wheelTabs,
	__global const uint2* modTaskIdx,
	__global const uint* modTaskList,
	__global uint* survivorCount,
	__global ulong* survivors,
	const uint survivorCapacity)
{
	// ST (initial task states, written by the setup and then moved to registers) shares its memory with S (sums of the current segment)
	// and SB (survivors of the current and of the previous segment: positions c << SIEVE_SEG_SHIFT | j)
	__local uint4 MEM4[(SIEVE_SHARED_WORDS + 3) / 4];
	__local int* S = (__local int*)MEM4;
	__local uint* SB = (__local uint*)MEM4 + SIEVE_SEG_POS;
	// during the setup: J0TAB (SIEVE_J0_STRIDE words per row modulus), then STW (one word per modulus)
	__local uint* J0TAB = (__local uint*)MEM4;
	__local uint* STW = (__local uint*)MEM4 + SIEVE_MAX_ROW_MODULI * SIEVE_J0_STRIDE;
	__local int INC[SIEVE_MAX_ROW_MODULI];
	__local int TBC[SIEVE_ROWS * SIEVE_TBC_STRIDE];
	__local int BASE[SIEVE_ROWS + 1];
	// steps of the wheel tasks, one table per sub-sieve (class of m mod 30, row: DK | CN << 8 | E << 16, see SieveGPUTables::Init), and the
	// first hits coprime to 30 (class, residue mod 30: T | class << 8)
	__local uint WHEEL[SIEVE_SUBS * 64];
	__local uint T0TAB[8 * 30];
	__local uint GF[16];
	__local uint OVR[1 + SIEVE_MAX_OVERRIDES * 2];
	// numbers of survivors in the two lists of SB, and the global index of the previous segment's survivors (reserved by work-item 0
	// with one atomic per segment)
	__local uint SBN[3];
	// the sub-sieve (p mod 7) where 7 divides every A1 * p + B1 (0: none), and the increment T(7, e_g + 1) - T(7, e_g)
	__local uint SUB7[2];
	// row offsets w(c) of the sub-sieves, alternating (the survivors of the last segment of a sub-sieve are written in the next one)
	__local uint SUBW[2 * SIEVE_ROWS];
	__local int WINS[2 * SIEVE_CHUNK_SEGMENTS]; // window of every segment of the chunk
#if SIEVE_KILL_QMAX > 0
	__local uint UNDO[SIEVE_UNDO_CAP]; // q^2 hits of killing wheel primes: sub-sieve * 128 + segment << 16 | c << SIEVE_SEG_SHIFT | j
#else
	__local uint UNDO[1];
#endif
	__local uint UNDON[1];

	const uint lid = get_local_id(0);
#ifdef SIEVE_LARGE_PRIMES
	// One work-group = one chunk of one candidate c: m = c * p, p in [P1, min(P2, SearchLimit::value / c)]
	const uint chunksPerCandidate = (uint)lpParams.w;
	const uint tier = (uint)(lpParams.w >> 32) & 0xFFu;
	const uint chunkSegments = (uint)(lpParams.w >> 40);
	const uint gid = get_group_id(0) + chunkOffset;
	const uint rangeIndex = gid / chunksPerCandidate;
	const uint chunkIndex = gid - rangeIndex * chunksPerCandidate;
	__global const uchar* cd = candidates + (ulong)(rangeIndex) * 10;
	const ulong cv = (ulong)cd[0] | ((ulong)cd[1] << 8) | ((ulong)cd[2] << 16) | ((ulong)cd[3] << 24) | ((ulong)cd[8] << 32);
	const ulong cs = ((ulong)cd[4] | ((ulong)cd[5] << 8) | ((ulong)cd[6] << 16) | ((ulong)cd[7] << 24) | ((ulong)cd[9] << 32)) + cv * 2;

	// N = (sigma(c) - c) * p + sigma(c) = g * (A1 * p + B1), g = gcd(sigma(c) - c, sigma(c)) (binary GCD)
	const ulong A = cs - cv;
	const ulong B = cs;
	const uint v2g = (uint)popcount(((A | B) & (0 - (A | B))) - 1);
	ulong ga = A >> popcount((A & (0 - A)) - 1);
	ulong gb = B;
	do
	{
		gb >>= popcount((gb & (0 - gb)) - 1);
		const ulong t = min(ga, gb);
		gb = max(ga, gb) - t;
		ga = t;
	} while (gb);
	const ulong gOdd = ga;
	ulong ginv = gOdd;
	for (int k = 0; k < 5; ++k)
	{
		ginv *= 2 - gOdd * ginv;
	}
	const ulong A1lo = (A >> v2g) * ginv;
	const ulong B1lo = (B >> v2g) * ginv;
	const ulong A1hi = 0;
	const ulong B1hi = 0;
	const uint twoAdic = (uint)(A1lo & B1lo & 1);
	const uint ngf = 0;
	const ulong sumM0lo = cs;
	const ulong sumM0hi = 0;

	// window: ln(B1 / A1) in float is accurate to about 30 fixed-point units, far below SIEVE_MARGIN
	const int lnr_fx = (int)(log((float)B1lo / (float)A1lo) * SIEVE_FX + 0.5f);
	const float rm1 = (float)cv / (float)A;

	// p range of this candidate (the upper bound only needs to be >= (SearchLimit::value - 1) / c)
	const ulong P0 = lpParams.z;
	const ulong p_lo = lpParams.x;
	// (for c < 543 the quotient does not fit in 64 bits, then the only limit is P2)
	const float qff = limitF / (float)cv;
	const ulong qf = (qff < 1.8e19f) ? (ulong)qff : 0xFFFFFFFFFFFFFFFFUL;
	const ulong p_hi = (qf >= 0xFFF0000000000000UL) ? lpParams.y : min(lpParams.y, qf + (qf >> 20) + 1);
	const ulong Pchunk = P0 + 210UL * (((ulong)(chunkIndex) * chunkSegments) << SIEVE_SEG_SHIFT);
	if ((Pchunk > p_hi) || (p_lo > p_hi))
	{
		return;
	}
	const uint nseg = (uint)min((((p_hi - Pchunk) / 210) >> SIEVE_SEG_SHIFT) + 1, (ulong)chunkSegments);

	// rho and slack for the largest N = A * p_hi + B, divided by the known power of 2 in N (see PrepareRange on the host)
	const ulong Nlo = A * p_hi + B;
	const ulong Nhi = mul_hi(A, p_hi) + ((Nlo < B) ? 1 : 0);
	const uint L = (Nhi ? (128 - (uint)clz(Nhi)) : (64 - (uint)clz(Nlo))) - v2g - twoAdic;
	const int loExtra = loExtraTab[tier * 128 + min(L, 127u)];
#else
	const uint4 ch = chunks[get_group_id(0) + chunkOffset];
	const uint rangeIndex = ch.x;
	__global const ulong* rg = ranges + (ulong)(rangeIndex) * SIEVE_RANGE_STRIDE;

	const ulong A1lo = rg[0];
	const ulong A1hi = rg[1];
	const ulong B1lo = rg[2];
	const ulong B1hi = rg[3];
	const ulong P0 = rg[8];
	const ulong p_lo = rg[9];
	const ulong p_hi = rg[10];
	const int lnr_fx = (int)((uint)(rg[11]));
	const int loExtra = (int)((uint)(rg[11] >> 32));
	const float rm1 = as_float((uint)(rg[12]));
	const uint tier = (uint)(rg[12] >> 32) & 0xFFu;
	const uint twoAdic = (uint)(rg[12] >> 40) & 1u;
	const uint v2g = (uint)(rg[12] >> 48) & 0xFFu;
	const uint ngf = (uint)(rg[12] >> 56);
	const ulong Pchunk = P0 + 210UL * ch.y;
	const uint nseg = ch.z;
	const ulong gOdd = 0;
	const ulong sumM0lo = rg[6];
	const ulong sumM0hi = rg[7];
#endif
	// positions of the last segment: j < jLast (the ones up to p_hi; the others are neither sieved nor scanned)
	const uint jLast = (uint)min((ulong)SIEVE_SEG_J, (p_hi - Pchunk - 210UL * ((ulong)(nseg - 1) << SIEVE_SEG_SHIFT)) / 210 + 1);
	const uint4 ti = tierInfo[tier * 2];      // number of primes, number of moduli, number of task slots, offset of the task table
	// the divisibility kill is only used in the last tier (ranges with a narrow window and many survivors), where it pays off
	const bool killTier = (SIEVE_KILL_QMAX > 0) && (tier == SIEVE_KILL_TIER);
	__global const uint2* tierTasks = tasks + ti.w;

#ifndef SIEVE_LARGE_PRIMES
	if (lid < 16)
	{
		GF[lid] = (lid < ngf) ? ((__global const uint*)(rg + 13))[lid] : 0;
	}
#endif
	if (lid <= SIEVE_ROWS)
	{
		BASE[lid] = 0;
	}
	if (lid == 0)
	{
		OVR[0] = 0;
		UNDON[0] = 0;
		SBN[0] = 0;
		SBN[1] = 0;
		SBN[2] = 0;
		SUB7[0] = 0;
		SUB7[1] = 0;
	}
	for (uint i = lid; i < SIEVE_SUBS * 64 + 8 * 30; i += SIEVE_WG_SIZE)
	{
		const uint w = wheelTabs[i];
		if (i < SIEVE_SUBS * 64)
		{
			WHEEL[i] = w;
		}
		else
		{
			T0TAB[i - SIEVE_SUBS * 64] = w;
		}
	}
	const uint numSlots = ti.z;
	barrier(CLK_LOCAL_MEM_FENCE);

	// roots of all moduli (one work-item per modulus, consecutive moduli in a warp: spreading the expensive row moduli over the warps left
	// only a third of the lanes of each warp busy) and the constants of the primes
	for (uint k = lid; k < ti.x; k += SIEVE_WG_SIZE)
	{
		SieveSetupPrimeConst(k, primeRecs, Ttab, A1lo, A1hi, B1lo, B1hi, ngf, GF, gOdd, BASE, SUB7);
	}
	for (uint mi = lid; mi < ti.y; mi += SIEVE_WG_SIZE)
	{
		SieveSetupModulus(mi, primeRecs, modRecs, modRecs2, Ttab, A1lo, A1hi, B1lo, B1hi, Pchunk, ngf, GF, gOdd, STW, J0TAB, INC, OVR, sumM0lo, sumM0hi, nseg, jLast, killTier, UNDO, UNDON);
	}

	// window of every segment: ln(T/N) = ln(B1/A1) - ln(1 + (B1/A1 - 1) / (p + 1)) grows with p
	for (uint s = lid; s < nseg; s += SIEVE_WG_SIZE)
	{
		const ulong Ps = Pchunk + 210UL * ((ulong)(s) << SIEVE_SEG_SHIFT);
		const ulong pa = max(Ps, p_lo);
		const ulong pb = min(Ps + 210UL * SIEVE_SEG_J, p_hi);
		const float xa = rm1 / ((float)(pa) + 1.0f);
		const float xb = rm1 / ((float)(pb) + 1.0f);
		const float ca = (xa < 1e-3f) ? (xa * (1.0f - xa * (0.5f - xa * (1.0f / 3.0f)))) : log(1.0f + xa);
		const float cb = (xb < 1e-3f) ? (xb * (1.0f - xb * (0.5f - xb * (1.0f / 3.0f)))) : log(1.0f + xb);
		WINS[s * 2] = lnr_fx - (int)(ca * SIEVE_FX) - 1 - loExtra;
		WINS[s * 2 + 1] = lnr_fx - (int)(cb * SIEVE_FX) + SIEVE_MARGIN;
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	// constant part of the sum for each row and 2-adic valuation of A1 * p + B1
	__global const int* T2 = Ttab + v2g;
	for (uint i = lid; i < SIEVE_ROWS * SIEVE_TBC_STRIDE; i += SIEVE_WG_SIZE)
	{
		const uint c = i / SIEVE_TBC_STRIDE;
		const uint v = i - c * SIEVE_TBC_STRIDE;
		TBC[i] = BASE[SIEVE_ROWS] + BASE[c] + T2[v];
	}
	const uint numOverrides = min(OVR[0], (uint)SIEVE_MAX_OVERRIDES);
	// the kill needs the complete undo list
	const uint numUndo = UNDON[0];
	const bool killEnabled = killTier && (numUndo <= SIEVE_UNDO_CAP);

)==="
R"===(
	// The tasks of this work-item live in registers: state (j << 3 | c; it decreases by SIEVE_SEG_J << 3 per segment, so SIEVE_NEVER
	// stays out of reach), geometry and increment:
	//   row task:   geometry = (c << SIEVE_SEG_SHIFT) << 16 | step, increment
	//   wheel task: geometry = SIEVE_GEO_WHEEL | class of m mod 30 << 28 | D << 14 | kill mode << 12 | m / 30, increment | x << 18 (kill,
	//               see SieveKillD: x = (p + 1) / 7 mod D for the next hit p, the kill applies if x != 0; D = 0: no kill)
	// The states at the start of the 6 sub-sieves are computed here from the root of the task (rt) and kept in tsub (with | x << 18 for
	// the wheel tasks with the kill): an index the compiler can't resolve keeps tsub in private memory instead of registers. (Computed
	// at the start of every sub-sieve, they cost 9% of the instructions of a sub-segment on k9.)
	uint tst[SIEVE_MAX_SLOTS];
	uint tgeo[SIEVE_MAX_SLOTS];
	int tinc[SIEVE_MAX_SLOTS];
	uint tsub[SIEVE_MAX_SLOTS * SIEVE_SUBS];
	const uint tsubZero = survivorCapacity >> 31;
	#pragma unroll
	for (uint s = 0; s < SIEVE_MAX_SLOTS; ++s)
	{
		tst[s] = SIEVE_NEVER;
		tgeo[s] = 0;
		tinc[s] = 0;
		// root: row task: first hit j of the task in row c of the 30-wheel (p = Pchunk + 30 j + W30[c]; | 0x80000000 for q = 7); wheel task:
		// x | (Pchunk + 1) mod D << 16, where p == Pchunk + x (mod m); SIEVE_NEVER: no hits
		uint rt = SIEVE_NEVER;
		if (s < numSlots)
		{
			const uint2 td = tierTasks[s * SIEVE_WG_SIZE + lid];
			if (td.x != SIEVE_TASK_EMPTY)
			{
				if (td.x & SIEVE_TASK_WHEEL)
				{
					const uint mi = td.x & 0xFFFu;
					const uint xr = STW[mi];
					int v = (int)(td.y >> 14);
					bool overridden = false;
					for (uint o = 0; o < numOverrides; ++o)
					{
						if (OVR[1 + o * 2] == mi)
						{
							v = (int)OVR[2 + o * 2];
							overridden = true;
						}
					}
					// kill (SieveKillD): mode 0, 7 does not divide D(q): D = D(q); mode 1, 7 | D(q) >= 49: D = D(q) / 7; mode 2, 7 | D(q) < 49:
					// D = 7, and 0 in the sub-sieve p == 6 (mod 7)
					uint D = 0;
					uint mode = 0;
					uint pc1 = 0;
#if SIEVE_KILL_QMAX > 0
					const uint m = modRecs[mi].x;
					if (killEnabled && !overridden && (m <= SIEVE_KILL_QMAX) && ((modRecs2[mi].z >> 28) == 1) && (xr != SIEVE_NEVER))
					{
						// m is prime (its first power)
						const uint Dq = SieveKillD(m, sumM0lo, sumM0hi);
						if (Dq)
						{
							pc1 = (uint)((Pchunk + 1) % Dq);
							mode = (Dq % 7) ? 0 : ((Dq >= 49) ? 1 : 2);
							D = (mode == 0) ? Dq : ((mode == 1) ? (Dq / 7) : 7);
						}
					}
#endif
					tinc[s] = v;
					rt = (xr == SIEVE_NEVER) ? SIEVE_NEVER : (xr | (pc1 << 16));
					tgeo[s] = SIEVE_GEO_WHEEL | (((td.x >> 12) & 7) << 28) | (D << 14) | (mode << 12) | (td.y & 0xFFFu);
				}
				else
				{
					// row task: row c, slice u of row modulus ri
					const uint ri = td.x & 0xFFFu;
					uint c = (td.x >> 12) & 7;
					if (td.x & (1u << 24))
					{
						// q = 3, 5: c is the index among the rows that have hits in this chunk
						uint k = c;
						c = SIEVE_ROWS;
						for (uint r = 0; r < SIEVE_ROWS; ++r)
						{
							if ((c == SIEVE_ROWS) && (J0TAB[ri * SIEVE_J0_STRIDE + r] != SIEVE_ROW_NEVER))
							{
								if (k == 0)
								{
									c = r;
								}
								--k;
							}
						}
					}
					const uint jc = (c < SIEVE_ROWS) ? J0TAB[ri * SIEVE_J0_STRIDE + c] : SIEVE_ROW_NEVER;
					rt = (jc == SIEVE_ROW_NEVER) ? SIEVE_NEVER : ((jc + ((td.x >> 16) & 0xFFu) * J0TAB[ri * SIEVE_J0_STRIDE + SIEVE_ROWS]) | (((td.x >> 25) & 1u) << 31));
					tinc[s] = INC[ri];
					tgeo[s] = (((c & 7) << SIEVE_SEG_SHIFT) << 16) | td.y;
				}
			}
		}

		// states at the start of the sub-sieves s7 = sub + 1
		const uint geo = tgeo[s];
		if (geo & SIEVE_GEO_WHEEL)
		{
			// hits p = Pchunk + y: y == x (mod m), y == s7 (mod 7) (Pchunk == 0 (mod 210)): y0 = x + k m with k = (s7 - x) m^-1 mod 7, which
			// grows by m^-1 mod 7 from one sub-sieve to the next; then the first y = y0 + 7 m t coprime to 30 (the hits step by 7 m)
			const uint x = rt & 0xFFFFu;
			const uint u = (geo >> 28) & 7u;
			const uint m = 30u * (geo & 0xFFFu) + SIEVE_W30(u);
			const uint m7 = 7u * m;
			const uint mu = SIEVE_INV7(m % 7u);
			const uint dy = mu * m;
			uint y0 = x + (((8u - x % 7u) * mu) % 7u) * m;
			__local const uint* t0 = T0TAB + SIEVE_U7(u) * 30;
#if SIEVE_KILL_QMAX > 0
			const uint D = (geo >> 14) & 0x3FFFu;
			const uint mode = (geo >> 12) & 3u;
			// the modulus of the kill state and its reciprocal for SieveMod32 (mode 0: D, mode 1: D(q) = 7 D in the sub-sieve p == 6 (mod 7))
			const uint Dk = (mode == 0) ? D : (7u * D);
			const uint Rk = D ? (0xFFFFFFFFu / Dk) : 0;
			const uint pc1 = rt >> 16;
			// Mode 0: the kill state of the first hit p = Pchunk + y, y = x + (k + 7 t) m, is (p + 1) 7^-1 = Z0 - k 7^-1 - t (mod D) (m == -1
			// (mod D)) with Z0 = (Pchunk + 1 + x) 7^-1; KI = k 7^-1 grows by m^-1 7^-1 from one sub-sieve to the next (k by m^-1 mod 7),
			// and by -1 where k wraps around (7 7^-1 == 1).
			const uint inv7 = (D && (mode == 0)) ? ((1u + ((7u - SIEVE_INV7(D % 7u)) % 7u) * D) / 7u) : 0;
			const uint Z0 = (mode == 0) ? SieveMod32(SieveMod32(pc1 + x, Dk, Rk) * inv7, Dk, Rk) : 0;
			uint KI = (mode == 0) ? SieveMod32((((8u - x % 7u) * mu) % 7u) * inv7, Dk, Rk) : 0;
			const uint dKI = (mode == 0) ? SieveMod32(mu * inv7, Dk, Rk) : 0;
#endif
			for (uint sub = 0; sub < SIEVE_SUBS; ++sub)
			{
				uint st = SIEVE_NEVER;
				if (rt != SIEVE_NEVER)
				{
					const uint tt = t0[y0 % 30u];
					const uint y = y0 + (tt & 0xFFu) * m7;
					// (y == w(c0) (mod 210), the row offset of the sub-sieve: j = y / 210)
					st = ((y / 210u) << 3) | (tt >> 8);
#if SIEVE_KILL_QMAX > 0
					if (D)
					{
						// kill state of the first hit p = Pchunk + y (st < 2^18): mode 0: (p + 1) 7^-1 mod D; sub-sieve p == 6 (mod 7), mode 1:
						// ((p + 1) mod D(q)) / 7 = ((p + 1) / 7) mod D (mode 2: no kill there); the other sub-sieves of modes 1 and 2: 7 does
						// not divide p + 1, every hit is killed: x = 16383 (it decreases by at most 6 per hit and doesn't reach 0 in a chunk)
						uint z;
						if (mode == 0)
						{
							int zz = (int)Z0 - (int)KI;
							zz += (zz < 0) ? (int)D : 0;
							zz -= (int)(tt & 0xFFu);
							zz += (zz < 0) ? (int)D : 0;
							z = (uint)zz;
						}
						else
						{
							z = (sub != 5) ? 16383u : ((mode == 1) ? (SieveMod32(pc1 + SieveMod32(y, Dk, Rk), Dk, Rk) / 7u) : 0u);
						}
						st |= z << 18;
					}
#endif
				}
				tsub[s * SIEVE_SUBS + (sub ^ tsubZero)] = st;
				y0 += dy;
				const bool wrap = (y0 >= m7);
				y0 -= wrap ? m7 : 0u;
#if SIEVE_KILL_QMAX > 0
				KI += dKI;
				KI -= (KI >= D) ? D : 0u;
				KI = wrap ? ((KI == 0) ? (D - 1) : (KI - 1)) : KI;
#endif
			}
		}
		else
		{
			// row task (and empty slots): the hits j30 = rt + i step in row c of the 30-wheel; those of the sub-sieve s7 have j30 == tau (mod 7)
			// (the step is prime to 7), from the i0 with rt + i0 step == tau (mod 7): i0 grows by 4 step^-1 mod 7 as tau grows by 4. q = 7:
			// all hits of 7^k are in the sub-sieve with p == r (mod 7).
			const uint c = (geo >> (16 + SIEVE_SEG_SHIFT)) & 7u;
			const uint step = geo & 0xFFFFu;
			const uint inv = SIEVE_INV7(step % 7u);
			uint tau = ((8u - SIEVE_W30(c) % 7u) * 4u) % 7u;
			uint i0 = ((tau + 7u - rt % 7u) * inv) % 7u;
			const uint di = (4u * inv) % 7u;
			for (uint sub = 0; sub < SIEVE_SUBS; ++sub)
			{
				uint st = SIEVE_NEVER;
				if (rt != SIEVE_NEVER)
				{
					if (rt & 0x80000000u)
					{
						const uint j30 = rt & 0x7FFFFFFFu;
						st = (j30 % 7u == tau) ? (((j30 - tau) * 0xB6DB6DB7u) << 3) : SIEVE_NEVER;
					}
					else
					{
						// (exact division by 7)
						st = ((rt + i0 * step - tau) * 0xB6DB6DB7u) << 3;
					}
				}
				tsub[s * SIEVE_SUBS + (sub ^ tsubZero)] = st;
				tau += 4u;
				tau -= (tau >= 7u) ? 7u : 0u;
				i0 += di;
				i0 -= (i0 >= 7u) ? 7u : 0u;
			}
		}
	}
	// the task states of sub-sieve 0
	#pragma unroll
	for (uint s = 0; s < SIEVE_MAX_SLOTS; ++s)
	{
		if (s < numSlots)
		{
			const uint v = tsub[s * SIEVE_SUBS + (0 ^ tsubZero)];
			tst[s] = v;
#if SIEVE_KILL_QMAX > 0
			if ((tgeo[s] & SIEVE_GEO_WHEEL) && (tgeo[s] & (0x3FFFu << 14)) && (v != SIEVE_NEVER))
			{
				tst[s] = v & 0x3FFFFu;
				tinc[s] = (int)(((uint)tinc[s] & 0x3FFFFu) | (v & 0xFFFC0000u));
				if ((0 == 5) && (((tgeo[s] >> 12) & 3u) == 2u))
				{
					// kill mode 2: no kill in the sub-sieve p == 6 (mod 7) (the last one)
					tgeo[s] &= ~(0x3FFFu << 14);
				}
			}
#endif
		}
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	// J0TAB and STW are no longer needed: clear the sums
	for (uint i = lid * 4; i < SIEVE_SEG_POS; i += SIEVE_WG_SIZE * 4)
	{
		*(__local int4*)(S + i) = (int4)(0);
	}
	barrier(CLK_LOCAL_MEM_FENCE);

)==="
R"===(
	const uint A1lo32 = (uint)A1lo;
	const uint B1lo32 = (uint)B1lo;

	uint pendingBase = 0; // work-item 0: global index of the previous segment's survivors
	ulong PsegPrev = 0;
	uint subPrev = 0;
#ifdef SIEVE_WARP_TIMING
	ulong tSieve = 0;
	ulong tScan = 0;
#endif
	for (uint sub = 0; sub < SIEVE_SUBS; ++sub)
	{
		// sub-sieve of the p == s7 (mod 7): row c holds p = Pseg + 210 j + w(c), w(c) = W30[c] + 30 tau(c) == s7 (mod 7)
		const uint s7 = sub + 1;
		__local const uint* wheelSub = WHEEL + sub * 64;
		__local const uint* subW = SUBW + (sub & 1) * SIEVE_ROWS;
		if (lid < SIEVE_ROWS)
		{
			SUBW[(sub & 1) * SIEVE_ROWS + lid] = SIEVE_W30(lid) + 30u * (((s7 + 7u - SIEVE_W30(lid) % 7u) * 4u) % 7u);
		}
#ifdef SIEVE_LARGE_PRIMES
		uint taupack = 0;
		#pragma unroll
		for (uint c = 0; c < SIEVE_ROWS; ++c)
		{
			taupack |= (((s7 + 7u - SIEVE_W30(c) % 7u) * 4u) % 7u) << (c * 3);
		}
#endif
		// the first power of 7 is the same for the whole sub-sieve: it moves the window instead of being added to the sums
		const int c7 = (SUB7[0] == s7) ? (int)SUB7[1] : 0;

		// (the task states of the sub-sieve: loaded in the last segment of the previous one, see below)

		for (uint seg = 0; seg < nseg; ++seg)
		{
			// (g: segment counter over all sub-sieves, for the survivor lists)
			const uint gseg = sub * nseg + seg;
			const uint b = gseg & 1;
			const uint jEnd = (seg + 1 < nseg) ? SIEVE_SEG_J : jLast;
#ifdef SIEVE_WARP_TIMING_SERIAL
			// measurement only: the warps sieve one after the other, so that the time of a warp is its own work
			for (uint tw = 0; tw < SIEVE_WG_SIZE / 32; ++tw)
			{
			barrier(CLK_LOCAL_MEM_FENCE);
			if ((lid >> 5) == tw)
			{
#endif
#ifdef SIEVE_WARP_TIMING
			const ulong tw0 = SieveClock();
#endif

			// sieve
			#pragma unroll
			for (uint s = 0; s < SIEVE_MAX_SLOTS; ++s)
			{
				if (s < numSlots)
				{
					uint st = tst[s];
					if (st < (jEnd << 3))
					{
						const uint geo = tgeo[s];
						uint j = st >> 3;
						if (geo & SIEVE_GEO_WHEEL)
						{
							uint c = st & 7;
							const uint mq = geo & 0xFFFu;
							const int D = (int)((geo >> 14) & 0x3FFFu);
							__local const uint* wt = wheelSub + ((geo >> 25) & 0x38u);
							const int v = (int)((uint)tinc[s] & 0x3FFFFu);
#if SIEVE_KILL_QMAX > 0
		// (the ranges stage runs this variant only for chunks of the last tier, the large primes stage for all chunks: D = 0 elsewhere)
		const int killv = D ? SIEVE_KILL : 0;
		int x = (int)((uint)tinc[s] >> 18);
		do
		{
			atomic_add(S + (c << SIEVE_SEG_SHIFT) + j, (x != 0) ? (v + killv) : v);
			const uint w = wt[c];
			const uint dk = w & 0xFFu;
			j += dk * mq + (w >> 16);
			c = (w >> 8) & 7;
			// p += 7 dk q and q == -1 (mod D): (p + 1) / 7 mod D decreases by dk
			x -= (int)dk;
			x += (x < 0) ? D : 0;
		} while (j < jEnd);
		st = (j << 3) | c;
		tinc[s] = (int)((uint)v | ((uint)x << 18));
#else
		(void)D;
		do
		{
			atomic_add(S + (c << SIEVE_SEG_SHIFT) + j, v);
			const uint w = wt[c];
			j += (w & 0xFFu) * mq + (w >> 16);
			c = (w >> 8) & 7;
		} while (j < jEnd);
		st = (j << 3) | c;
#endif
						}
						else
						{
							const int v = tinc[s];
							__local int* row = S + (geo >> 16);
							const uint step = geo & 0xFFFFu;
							do
							{
								atomic_add(row + j, v);
								j += step;
							} while (j < jEnd);
							st = j << 3;
						}
					}
					tst[s] = st - (SIEVE_SEG_J << 3);
				}
			}
#ifdef SIEVE_WARP_TIMING
			tSieve += SieveClock() - tw0;
#endif
#ifdef SIEVE_WARP_TIMING_SERIAL
			}
			}
#endif
			if ((seg + 1 == nseg) && (sub + 1 < SIEVE_SUBS))
			{
				// the task states of the next sub-sieve (the task registers are not used before its first sieve phase: the scan and the
				// barriers of this segment hide the latency of the loads)
				#pragma unroll
				for (uint s = 0; s < SIEVE_MAX_SLOTS; ++s)
				{
					if (s < numSlots)
					{
						const uint v = tsub[s * SIEVE_SUBS + ((sub + 1) ^ tsubZero)];
						tst[s] = v;
#if SIEVE_KILL_QMAX > 0
						if ((tgeo[s] & SIEVE_GEO_WHEEL) && (tgeo[s] & (0x3FFFu << 14)) && (v != SIEVE_NEVER))
						{
							tst[s] = v & 0x3FFFFu;
							tinc[s] = (int)(((uint)tinc[s] & 0x3FFFFu) | (v & 0xFFFC0000u));
							if (((sub + 1) == 5) && (((tgeo[s] >> 12) & 3u) == 2u))
							{
								// kill mode 2: no kill in the sub-sieve p == 6 (mod 7) (the last one)
								tgeo[s] &= ~(0x3FFFu << 14);
							}
						}
#endif
					}
				}
			}
#if SIEVE_KILL_QMAX > 0
			if (killEnabled)
			{
				for (uint i = lid; i < numUndo; i += SIEVE_WG_SIZE)
				{
					const uint e = UNDO[i];
					if ((e >> 16) == sub * 128 + seg)
					{
						atomic_add(S + (e & 0xFFFFu), -SIEVE_KILL);
					}
				}
			}
#endif

			const ulong Pseg = Pchunk + 210UL * ((ulong)(seg) << SIEVE_SEG_SHIFT);
			if ((lid == 0) && (gseg > 0))
			{
				SBN[2] = pendingBase;
			}
			barrier(CLK_LOCAL_MEM_FENCE);

			// survivors of the previous segment (it can belong to the previous sub-sieve)
			if (gseg > 0)
			{
				SieveFlushSurvivors(SB + (b ^ 1) * SIEVE_SURV_LOCAL, min(SBN[b ^ 1], (uint)SIEVE_SURV_LOCAL), SBN[2], PsegPrev, SUBW + (subPrev & 1) * SIEVE_ROWS, rangeIndex, survivors, survivorCapacity, lid);
			}

			const int lo = WINS[seg * 2] - c7;
			const int hi = WINS[seg * 2 + 1] - c7;
			const uint width = (uint)(hi - lo);
			const uint Pseg32 = (uint)Pseg;
#ifdef SIEVE_LARGE_PRIMES
			// j of the segment start in the 30-wheel layout of the bitmap of primes
			const uint segJ0 = (uint)((Pseg - P0) / 30);
#endif

#ifdef SIEVE_WARP_TIMING
			const ulong tw1 = SieveClock();
#endif
			// scan (SieveScan; the last segment of a chunk without its unsieved positions)
#ifdef SIEVE_DUMP
			__global int* dump = (__global int*)survivors + ((ulong)(get_group_id(0)) * SIEVE_SUBS + sub) * SIEVE_SEG_POS;
			if ((lid == 0) && (sub == 0))
			{
				// self-test: the window (every chunk has one segment)
				__global int* dw = (__global int*)survivors + (ulong)(get_num_groups(0)) * SIEVE_SUBS * SIEVE_SEG_POS + get_group_id(0) * 2;
				dw[0] = lo + c7;
				dw[1] = hi + c7;
			}
			uint cand = (jEnd == SIEVE_SEG_J) ? SieveScan(S, TBC, subW, lid, twoAdic, A1lo32, B1lo32, Pseg32, lo, width, jEnd, false, dump, c7) : SieveScan(S, TBC, subW, lid, twoAdic, A1lo32, B1lo32, Pseg32, lo, width, jEnd, true, dump, c7);
#else
			uint cand = (jEnd == SIEVE_SEG_J) ? SieveScan(S, TBC, subW, lid, twoAdic, A1lo32, B1lo32, Pseg32, lo, width, jEnd, false) : SieveScan(S, TBC, subW, lid, twoAdic, A1lo32, B1lo32, Pseg32, lo, width, jEnd, true);
#endif
)==="
R"===(
#ifdef EXP_NO_SURVIVORS
			cand = 0;
#endif
			while (cand)
			{
				const uint bit = 31 - clz(cand & (0u - cand));
				cand &= cand - 1;
				const uint i = ((bit >> 2) * SIEVE_WG_SIZE + lid) * 4;
				const uint c = i >> SIEVE_SEG_SHIFT;
				const uint jk = (i & (SIEVE_SEG_J - 1)) + (bit & 3);
				const ulong p = Pseg + 210UL * jk + subW[c];
#ifdef SIEVE_LARGE_PRIMES
				// only prime p are checked (the host sieved all primes of [P1, P2]; the bitmap has the 30-wheel layout)
				const uint jg = segJ0 + 7u * jk + ((taupack >> (c * 3)) & 7u);
				if ((p >= p_lo) && (p <= p_hi) && ((primeBits[c * wordsPerRow + (jg >> 5)] >> (jg & 31)) & 1))
#else
				if ((p >= p_lo) && (p <= p_hi))
#endif
				{
					const uint q = atomic_inc(SBN + b);
					if (q < SIEVE_SURV_LOCAL)
					{
						SB[b * SIEVE_SURV_LOCAL + q] = (c << SIEVE_SEG_SHIFT) | jk;
					}
					else
					{
						// the local list is full (very rare): write it directly
						const uint g = atomic_inc(survivorCount);
						if (g < survivorCapacity)
						{
#ifdef SIEVE_LARGE_PRIMES
							survivors[g * 2] = p;
							survivors[g * 2 + 1] = rangeIndex;
#else
							survivors[g] = (p << SIEVE_RANGE_INDEX_BITS) | rangeIndex;
#endif
						}
					}
				}
			}
#ifdef SIEVE_WARP_TIMING
			tScan += SieveClock() - tw1;
#endif
			barrier(CLK_LOCAL_MEM_FENCE);

			// reserve space for this segment's survivors (the result is only needed after the next sieve phase), and reuse the list of the
			// previous segment (all work-items copied it before the barrier)
			if (lid == 0)
			{
				SBN[b ^ 1] = 0;
				pendingBase = atomic_add(survivorCount, min(SBN[b], (uint)SIEVE_SURV_LOCAL));
			}
			PsegPrev = Pseg;
			subPrev = sub;
		}
	}

	// survivors of the last segment
	if (lid == 0)
	{
		SBN[2] = pendingBase;
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	const uint bLast = (SIEVE_SUBS * nseg - 1) & 1;
	SieveFlushSurvivors(SB + bLast * SIEVE_SURV_LOCAL, min(SBN[bLast], (uint)SIEVE_SURV_LOCAL), SBN[2], PsegPrev, SUBW + (subPrev & 1) * SIEVE_ROWS, rangeIndex, survivors, survivorCapacity, lid);
#ifdef SIEVE_WARP_TIMING
	if ((lid & 31) == 0)
	{
		__global ulong* timing = (__global ulong*)survivors + (ulong)survivorCapacity - 2 * (SIEVE_WG_SIZE / 32);
		atom_add(timing + (lid >> 5) * 2, tSieve);
		atom_add(timing + (lid >> 5) * 2 + 1, tScan);
	}
#endif
}
)==="
R"===(
#ifdef SIEVE_LARGE_PRIMES
// Exact check of the sieve survivors of the large primes stage: (p, candidate index) pairs, p is prime
__kernel
__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
void CheckSieveSurvivors(
	__global const ulong* survivors,
	__global const uint* survivorCount,
	const uint survivorCapacity,
	__global const uchar* candidates,
	__global const uint* smallPrimes,
	__global const ulong2* primeInverses,
	__global const ulong2* PQ,
	__constant ulong2* PowerOf2SumInverses,
	__global const uint* PowersOfP_128SumInverses_offsets,
	__global const ulong4* PowersOfP_128SumInverses,
	__global const ulong2* PrimeInverses_128,
	__global const ulong* SumEstimates_128,
	__global ulong* phase1_offset_to_resume_after_overflow,
	__global uint* phase2_numbers_count,
	__global ulong4* phase2_numbers_data,
	__global uint* amicable_numbers_count,
	__global ulong4* amicable_numbers_data)
{
	// grid-stride loop: this kernel is launched with a fixed size right after the sieve, without reading the count on the host
	const uint count = min(*survivorCount, survivorCapacity);
	for (uint i = get_global_id(0); i < count; i += get_global_size(0))
	{
		if (*phase2_numbers_count > PHASE2_MAX_COUNT)
		{
			return;
		}

		const ulong p = survivors[i * 2]; // prime (from the bitmap of primes)

		__global const uchar* cd = candidates + survivors[i * 2 + 1] * 10;
		const ulong cv = (ulong)cd[0] | ((ulong)cd[1] << 8) | ((ulong)cd[2] << 16) | ((ulong)cd[3] << 24) | ((ulong)cd[8] << 32);
		const ulong cs = ((ulong)cd[4] | ((ulong)cd[5] << 8) | ((ulong)cd[6] << 16) | ((ulong)cd[7] << 24) | ((ulong)cd[9] << 32)) + cv * 2;

		CheckPairPhase1(
			smallPrimes,
			primeInverses,
			PQ,
			PowerOf2SumInverses,
			PowersOfP_128SumInverses_offsets,
			PowersOfP_128SumInverses,
			PrimeInverses_128,
			SumEstimates_128,
			cv * p, mul_hi(cv, p), cs * (p + 1), mul_hi(cs, p + 1),
			0,
			phase1_offset_to_resume_after_overflow,
			phase2_numbers_count,
			phase2_numbers_data,
			amicable_numbers_count,
			amicable_numbers_data,
			0, phase2_numbers_count, (__global ulong*)phase2_numbers_data, 0
		);
	}
}
#else
// Exact check of the sieve survivors: primality of p, then the same checks as for every number in SearchMultipleRanges
__kernel
__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
void CheckSieveSurvivors(
	__global const ulong* survivors,
	__global const uint* survivorCount,
	const uint survivorCapacity,
	__global const ulong* ranges,
	__global const uint* smallPrimes,
	__global const ulong2* primeInverses,
	__global const ulong2* PQ,
	__constant ulong2* PowerOf2SumInverses,
	__global const uint* PowersOfP_128SumInverses_offsets,
	__global const ulong4* PowersOfP_128SumInverses,
	__global const ulong2* PrimeInverses_128,
	__global const ulong* SumEstimates_128,
	__global ulong* phase1_offset_to_resume_after_overflow,
	__global uint* phase2_numbers_count,
	__global ulong4* phase2_numbers_data,
		__global uint* amicable_numbers_count,
		__global ulong4* amicable_numbers_data,
		__global uint* deferred_count,
		__global ulong* deferred_list,
		const uint defer,
		__global const uint* bandTable,
		const uint bandFirst,
		const uint bandCount,
		__global const uint* bandBitmap)
	{
		// grid-stride loop: this kernel is launched with a fixed size right after the sieve, without reading the count on the host.
		// defer = 1: the survivors of the sieve (p not known to be prime), defer = 0: the survivors with a prime p (SieveTestDeferred).
		// The work-items of a work-group take the survivors of an iteration together, and a barrier at its end makes the warps reconverge
		// there (with a loop over i = get_global_id(0) + k * get_global_size(0), the work-items of a warp went through the survivors at
		// their own pace and the warps never reconverged: 5.9 active threads per warp instruction on average in k4, and the loop of
		// this form without the barrier was not reconverged by the compiler either; with the barrier the kernel took 2.2x less time).
		// All work-items of a work-group have the same number of iterations (base is the same for them), so they all reach the barrier.
		const uint count = min(*survivorCount, survivorCapacity);
		for (uint base = get_group_id(0) * get_local_size(0); base < count; base += get_global_size(0))
		{
		const uint i = base + get_local_id(0);
		do
		{
		if (i >= count)
		{
			break;
		}

		// p is not known to be prime: CheckPairPhase1 tests it only if the number passes all of its checks
		const ulong s = survivors[i];
		const ulong p = s >> SIEVE_RANGE_INDEX_BITS;

		// Most composite p (the positions are all numbers prime to 210) have a factor 11 <= q <= 127: 27 tests with the 64-bit inverses
		// (primeInverses[k] belongs to smallPrimes[k], k = 4..30) are much cheaper than phase 1 and the deferred test of p. In the
		// ranges with W ~ 0 almost every survivor passes phase 1, and about 2/3 of them have a composite p.
		if (defer && (p > 127))
		{
			bool smallFactor = false;
			#pragma unroll
			for (int k = 4; k <= 30; ++k)
			{
				smallFactor |= (p * primeInverses[k].x <= primeInverses[k].y);
			}
			if (smallFactor)
			{
				break;
			}
		}

		__global const ulong* rg = ranges + (s & ((1UL << SIEVE_RANGE_INDEX_BITS) - 1)) * SIEVE_RANGE_STRIDE;
		const ulong M0lo = rg[4];
		const ulong M0hi = rg[5];
		const ulong sumM0lo = rg[6];
		const ulong sumM0hi = rg[7];

		// band sieve (SieveBand, kernel.cl: DecideByBand): the smallest band prime of this position if its chunk is in the band table of the
		// launch (and the smallest prime of g = gcd(A, B) in the band: it divides N, but not A1 p + B1)
		ulong band = 0;
		if (bandCount != 0)
		{
			const uint rangeIndex = (uint)(s & ((1UL << SIEVE_RANGE_INDEX_BITS) - 1));
			const ulong J = (p - rg[8]) / 30;
			const uint w = (uint)(p - rg[8] - J * 30);
			// the entries of a launch are sorted by (range index, first j): the last one <= (rangeIndex, J)
			uint lo = 0;
			uint hi = bandCount;
			while (hi - lo > 1)
			{
				const uint mid = (lo + hi) >> 1;
				__global const uint* bm = bandTable + (bandFirst + mid) * SIEVE_BAND_ENTRY;
				if ((bm[0] < rangeIndex) || ((bm[0] == rangeIndex) && (bm[1] <= J)))
				{
					lo = mid;
				}
				else
				{
					hi = mid;
				}
			}
			__global const uint* be = bandTable + (bandFirst + lo) * SIEVE_BAND_ENTRY;
			const uint c = (w > 1) + (w > 7) + (w > 11) + (w > 13) + (w > 17) + (w > 19) + (w > 23);
			// (p = 2, 3, 5 are not positions of the sieve)
			if ((be[0] == rangeIndex) && (J >= be[1]) && (J - be[1] < be[2]) && (SIEVE_W30(c) == w))
			{
				const uint rmin = min(min(bandBitmap[be[3] + (uint)(J - be[1]) * 8 + c], be[6]), be[5]);
				band = (ulong)rmin | ((ulong)be[4] << 21) | ((ulong)be[5] << 42);
			}
		}

		CheckPairPhase1(
			smallPrimes,
			primeInverses,
			PQ,
			PowerOf2SumInverses,
			PowersOfP_128SumInverses_offsets,
			PowersOfP_128SumInverses,
			PrimeInverses_128,
			SumEstimates_128,
			M0lo * p, mul_hi(M0lo, p) + M0hi * p, sumM0lo * (p + 1), mul_hi(sumM0lo, p + 1) + sumM0hi * (p + 1),
			0,
			phase1_offset_to_resume_after_overflow,
			phase2_numbers_count,
			phase2_numbers_data,
						amicable_numbers_count,
						amicable_numbers_data,
						defer ? s : 0,
						deferred_count,
						deferred_list,
						band
					);
		} while (0);
		barrier(CLK_LOCAL_MEM_FENCE);
				}
			}

// x mod r for x < 2^64, r < 2^32, m = floor((2^64 - 1) / r) (Barrett: the quotient estimate is at most 2 too small)
static inline uint BandMod64(const ulong x, const uint r, const ulong m)
{
	ulong rem = x - mul_hi(x, m) * r;
	rem = (rem >= r) ? (rem - r) : rem;
	rem = (rem >= r) ? (rem - r) : rem;
	return (uint)rem;
}

static inline uint BandMulMod(const uint a, const uint b, const uint r, const ulong m)
{
	return BandMod64((ulong)a * b, r, m);
}

static uint BandPowMod(uint a, uint e, const uint r, const ulong m)
{
	uint x = 1;
	while (e)
	{
		const uint y = BandMulMod(x, a, r, m);
		x = (e & 1) ? y : x;
		a = BandMulMod(a, a, r, m);
		e >>= 1;
	}
	return x;
}

// Band sieve for the ranges with W ~ 0 (the host: RunSieveLaunches, kernel.cl: DecideByBand). Entry e of the launch (bandTable +
// SIEVE_BAND_ENTRY * (bandFirst + e)): range index, first j, number of j, offset in the bitmap, prime indexes [begin, end) of the band,
// the index of the smallest prime of g in the band (or more), 0. One work-item per prime r = smallPrimes[i]: for every entry whose primes
// contain r, bitmap[offset + (j - first j) * 8 + c] of the positions p = P0 + 30 j + W30(c) with r | A1 p + B1 gets min(itself, i) (the
// host fills it with 0xFFFFFFFF). A1 p + B1 = 0 (mod r) <=> j = -(A1 P0 + B1) / (30 A1) - W30(c) / 30 (mod r); the inverses of
// A1 mod r of 16 entries are computed with one modular inverse (Montgomery's batch inversion). gcd(A1, B1) = 1, so r | A1 means no
// position.
#define SIEVE_BAND_GROUP 16

__kernel
__attribute__((reqd_work_group_size(256, 1, 1)))
void SieveBand(
	__global const ulong* ranges,
	__global const uint* bandTable,
	const uint bandFirst,
	const uint bandCount,
	const uint primeBegin,
	const uint primeEnd,
	__global const uint* smallPrimes,
	__global uint* bitmap)
{
	const uint pi = primeBegin + get_global_id(0);
	if (pi >= primeEnd)
	{
		return;
	}
	const uint r = smallPrimes[pi];
	const ulong m = 0xFFFFFFFFFFFFFFFFUL / r;
	// 2^64 mod r
	uint c64 = (uint)(0xFFFFFFFFFFFFFFFFUL - m * r) + 1;
	c64 = (c64 == r) ? 0 : c64;
	// The inverse of 30 comes with the first batch inversion (inv(30 P) P = inv(30), inv(30 P) 30 = inv(P)): no exponentiation of its
	// own (the host runs this kernel over slices of the entries of a launch, every slice computed it again)
	uint inv30 = 0;
	uint wv[8];
	bool have30 = false;

	for (uint e0 = 0; e0 < bandCount; e0 += SIEVE_BAND_GROUP)
	{
		const uint ne = min((uint)SIEVE_BAND_GROUP, bandCount - e0);
		uint a[SIEVE_BAND_GROUP];
		uint bm[SIEVE_BAND_GROUP];
		uint pre[SIEVE_BAND_GROUP];
		uint prod = 1;
		uint used = 0;
		for (uint k = 0; k < ne; ++k)
		{
			__global const uint* be = bandTable + (bandFirst + e0 + k) * SIEVE_BAND_ENTRY;
			uint av = 1;
			bm[k] = 0;
			if ((pi >= be[4]) && (pi < be[5]))
			{
				__global const ulong* rg = ranges + (ulong)be[0] * SIEVE_RANGE_STRIDE;
				// A1 mod r, B1 mod r (128-bit: hi * 2^64 + lo)
				const uint ah = BandMod64(rg[1], r, m);
				uint am = BandMod64((ulong)ah * c64, r, m) + BandMod64(rg[0], r, m);
				am = (am >= r) ? (am - r) : am;
				if (am != 0)
				{
					const uint bh = BandMod64(rg[3], r, m);
					uint b = BandMod64((ulong)bh * c64, r, m) + BandMod64(rg[2], r, m);
					bm[k] = (b >= r) ? (b - r) : b;
					av = am;
					used |= 1U << k;
				}
			}
			a[k] = av;
			pre[k] = prod;
			prod = BandMulMod(prod, av, r, m);
		}
		if (used == 0)
		{
			continue;
		}
		uint inv;
		if (!have30)
		{
			const uint inv30prod = BandPowMod(BandMulMod(prod, 30, r, m), r - 2, r, m);
			inv30 = BandMulMod(inv30prod, prod, r, m);
			inv = BandMulMod(inv30prod, 30, r, m);
			for (uint c = 0; c < 8; ++c)
			{
				wv[c] = BandMulMod(SIEVE_W30(c), inv30, r, m);
			}
			have30 = true;
		}
		else
		{
			inv = BandPowMod(prod, r - 2, r, m);
		}
		for (int k = (int)ne - 1; k >= 0; --k)
		{
			const uint ainv = BandMulMod(inv, pre[k], r, m);
			inv = BandMulMod(inv, a[k], r, m);
			if ((used & (1U << k)) == 0)
			{
				continue;
			}
			__global const uint* be = bandTable + (bandFirst + e0 + (uint)k) * SIEVE_BAND_ENTRY;
			__global const ulong* rg = ranges + (ulong)be[0] * SIEVE_RANGE_STRIDE;
			// j0 = -(A1 P0 + B1) / (30 A1) mod r
			uint t = BandMulMod(a[k], BandMod64(rg[8], r, m), r, m) + bm[k];
			t = (t >= r) ? (t - r) : t;
			const uint jb = BandMulMod((t == 0) ? 0 : (r - t), BandMulMod(ainv, inv30, r, m), r, m);
			const uint js = be[1] % r;
			const uint count = be[2];
			const uint offset = be[3];
			for (uint c = 0; c < 8; ++c)
			{
				const uint jr = (jb >= wv[c]) ? (jb - wv[c]) : (jb + r - wv[c]);
				// (the test of p only in the entries with at least 2^20 j, 32 MB of bitmap: the host slices the launches at 32 MB, an entry
				// above it is a slice of its own whose hits miss the L2 cache. In the small entries ("t1", 2.6 chunks per range) the hits are
				// cheap, and the test everywhere made SieveBand 17.6 -> 18.5 s)
				const ulong pc = rg[8] + 30 * (ulong)be[1] + SIEVE_W30(c);
				const bool filter = (count >= (1U << 20));
				for (uint d = (jr >= js) ? (jr - js) : (jr + r - js); d < count; d += r)
				{
					if (!filter || !BandSmallFactor(pc + 30 * (ulong)d))
					{
						atomic_min(&bitmap[offset + d * 8 + c], pi);
					}
				}
			}
		}
	}
}

// Primality test of p for the survivors that passed CheckPairPhase1 (see Phase1Defer), all in parallel: the ones with a prime p are
// written to the start of the survivor buffer (it has been checked already) for the second CheckSieveSurvivors pass
__kernel
__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
// (spsp2: the base-2 strong pseudoprimes below 2^32, sorted: the test of p < 2^32 is one Miller-Rabin round and a binary search in it,
// kernel.cl: SieveIsPrimeList)
void SieveTestDeferred(
	__global const uint* deferred_count,
	__global const ulong* deferred_list,
	__global ulong* survivors,
	__global uint* primeCount,
	__global const uint* spsp2,
	const uint spsp2Count)
{
	const uint count = min(*deferred_count, (uint)SIEVE_DEFERRED_MAX);
	for (uint i = get_global_id(0); i < count; i += get_global_size(0))
	{
		const ulong s = deferred_list[i];
		if (SieveIsPrimeList(s >> SIEVE_RANGE_INDEX_BITS, spsp2, spsp2Count))
		{
			survivors[atomic_inc(primeCount)] = s;
		}
	}
}
#endif

// Runs after every sieve launch and the survivor checks: counts the survivors and the deferred numbers, flags the launch if the
// survivor buffer (1) or the phase 2 buffer (2) overflowed (then the host runs it again), resets the counters
__kernel
__attribute__((reqd_work_group_size(1, 1, 1)))
void SieveLaunchDone(__global uint* survivorCount, __global uint* stats, const uint survivorCapacity, const uint launchIndex, __global const uint* phase2_numbers_count, __global uint* deferredCount, __global uint* primeCount)
{
	const uint count = *survivorCount;
	const uint deferred = *deferredCount;
	*primeCount = 0;
	// survivors of the launch (for the estimates of the host) and the flags
	const uint countBits = min(count, 0x3FFFFFFFu) << 2;
	if (count > survivorCapacity)
	{
		stats[4 + launchIndex] = countBits | 1;
	}
	else if (*phase2_numbers_count > PHASE2_MAX_COUNT)
	{
		stats[4 + launchIndex] = countBits | 2;
	}
	else
	{
		stats[4 + launchIndex] = countBits;
		const ulong total = ((ulong)stats[1] << 32) + stats[0] + count;
		stats[0] = (uint)total;
		stats[1] = (uint)(total >> 32);
		const ulong totalDeferred = ((ulong)stats[3] << 32) + stats[2] + deferred;
		stats[2] = (uint)totalDeferred;
		stats[3] = (uint)(totalDeferred >> 32);
	}
	*survivorCount = 0;
	*deferredCount = 0;
}
#endif
)===";
