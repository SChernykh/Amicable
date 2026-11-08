#pragma once

static const char* kernel_cl = R"===(
#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable

#ifdef cl_clang_storage_class_specifiers
#pragma OPENCL EXTENSION cl_clang_storage_class_specifiers : enable
#endif

#define PHASE1_DEPTH 256

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Check amicability for a range of numbers of the form "M0*p" where p is prime
// p runs through all values in primes[p_offset]...primes[p_offset + get_global_size(0) - 1] (inclusive)
// primes stores prime numbers in compact form: 4 primes per each 8 bytes
// primeInverses is an array of multiplicative inverses of first 192725 prime numbers, this is used for checking divisibility by doing just one multiplication
// primeReciprocals is an array of reciprocals of first 192725 prime numbers, this is used for calculating N / p (not just checking divisibility)
// PQ corresponds to tables P and Q in lemma 2.1 from "Computation of All the Amicable Pairs Below 10^10 By H.J.J.te Riele": http://www.ams.org/journals/mcom/1986-47-175/S0025-5718-1986-0842142-3/S0025-5718-1986-0842142-3.pdf
// The only difference is that we calculate exact (hence better) upper bounds for S(m)/m instead of inexact estimates
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// A: smaller member of an amicable pair
// B: if it's not zero, it must be prime in order for A to be amicable number
// amicable_numbers_count[1] is the capacity of the buffer (the host reports an error if the counter is larger)
static void NumberFound(__global uint* amicable_numbers_count, __global ulong4* amicable_numbers_data, ulong A, ulong A_high, ulong B)
{
	const uint index = atom_inc(amicable_numbers_count);
	if (index < amicable_numbers_count[1])
	{
		amicable_numbers_data[index].x = A;
		amicable_numbers_data[index].y = A_high;
		amicable_numbers_data[index].z = B;
	}
}

// 32-bit floats don't have enough precision
// 64-bit floats are not supported by all OpenCL devices
// So we have to use integer arithmetic to calculate square roots

// Returns number x such that x^2 <= n < (x+1)^2
static uint IntegerSquareRoot(const ulong n)
{
	uint result = 1;
	result <<= ((63 - clz(n)) >> 1);
	for (uint cur_bit = result >> 1; cur_bit > 0; cur_bit >>= 1)
	{
		const uint result1 = result | cur_bit;
		const ulong k = result1;
		if (k * k <= n)
		{
			result = result1;
		}
	}
	return result;
}

// Returns number x such that x^2 <= n < (x+1)^2
static ulong IntegerSquareRoot128(const ulong2 n)
{
	ulong result = 1;
	const ulong highest_bit_index = n.y ? (127 - clz(n.y)) : (63 - clz(n.x));
	result <<= (highest_bit_index >> 1);
	for (ulong cur_bit = result >> 1; cur_bit > 0; cur_bit >>= 1)
	{
		const ulong k = result | cur_bit;
		ulong2 k2;
		k2.x = k * k;
		k2.y = mul_hi(k, k);
		if ((k2.y < n.y) || ((k2.y == n.y) && (k2.x <= n.x)))
		{
			result = k;
		}
	}
	return result;
}

// p | N for N = N_high * 2^64 + N with N_high < p (so N / p < 2^64), an odd prime p < 2^32 and inv = p^-1 mod 2^64: with q = N * inv mod 2^64,
// q * p = N mod 2^64 always, so p | N <=> q * p = N <=> the high 64 bits of q * p are N_high. Two 32x32-bit products and one 64-bit product
// instead of the three 64-bit products (one of them a high part) of the test with the 128-bit inverse.
static inline bool Divides128Fast(const ulong N, const uint N_high, const ulong inv, const uint p)
{
	const ulong q = N * inv;
	const ulong t = (ulong)((uint)q) * p;
	const ulong u = (ulong)((uint)(q >> 32)) * p + (t >> 32);
	return (uint)(u >> 32) == N_high;
}

// The same test with 32-bit products (exact division from the low end), pinv = p^-1 mod 2^32: q0 = n0 / p mod 2^32 makes N - q0 p a
// multiple of 2^32, N1 = (N - q0 p) / 2^32 = (N >> 32) - hi(q0 p), and p | N <=> N1 = q1 p with q1 = N1 / p mod 2^32 (N / p < 2^64
// because N_high < p) <=> hi(N1) = hi(q1 p) (lo(q1 p) = lo(N1) by construction; N1 < 0 gives hi(N1) = 2^32 - 1 > hi(q1 p)).
// Two products and two high halves of 32x32 bits instead of two 32x32->64-bit products and a 64-bit product: phase 2 was limited by the
// pipe of the integer multiplications (FMA-heavy, 80% busy).
static inline bool Divides128Fast32(const ulong N, const uint N_high, const uint pinv, const uint p)
{
	const uint q0 = (uint)N * pinv;
	const ulong N1 = ((N >> 32) | ((ulong)N_high << 32)) - (ulong)mul_hi(q0, p);
	const uint q1 = (uint)N1 * pinv;
	return (uint)(N1 >> 32) == mul_hi(q1, p);
}

// Returns number x such that x^3 <= n < (x+1)^3
static ulong IntegerCubeRoot128(const ulong2 n)
{
	ulong result = 1;
	const uint highest_bit_index = n.y ? (127 - clz(n.y)) : (63 - clz(n.x));
	result <<= (highest_bit_index * ((65536 / 3) + 1)) >> 16;
	for (ulong cur_bit = result >> 1; cur_bit > 0; cur_bit >>= 1)
	{
		const ulong k = result | cur_bit;
		ulong2 k2;
		k2.x = (k * k) * k;
		k2.y = mul_hi(k * k, k);
		if ((k2.y < n.y) || ((k2.y == n.y) && (k2.x <= n.x)))
		{
			result = k;
		}
	}
	return result;
}
)==="
R"===(
// Montgomery multiplication: a * b / 2^64 mod n for a, b < n < 2^63, ninv = -n^-1 mod 2^64
static inline ulong SieveMontMul(const ulong a, const ulong b, const ulong n, const ulong ninv)
{
	const ulong lo = a * b;
	const ulong m = lo * ninv;
	const ulong t = mul_hi(a, b) + mul_hi(m, n) + ((lo != 0) ? 1 : 0);
	return (t >= n) ? (t - n) : t;
}

// Montgomery multiplication: a * b / 2^32 mod n for a, b < n < 2^32 (n odd), ninv = -n^-1 mod 2^32. The low 32 bits of a * b + m * n are 0,
// so the sum is (hi(a * b) + hi(m * n) + 1) * 2^32 if the low 32 bits of a * b are not 0; the result before the subtraction is < 2n.
static inline uint SieveMontMul32(const uint a, const uint b, const uint n, const uint ninv)
{
	const ulong t = (ulong)a * b;
	const uint m = (uint)t * ninv;
	const ulong u = (t >> 32) + (ulong)mul_hi(m, n) + (((uint)t != 0) ? 1 : 0);
	return (u >= n) ? (uint)(u - n) : (uint)u;
}

// Miller-Rabin round to base a for odd n < 2^32 in 32-bit Montgomery arithmetic: ninv = -n^-1 mod 2^32, one = 2^32 mod n, r2 = 2^64 mod n,
// n - 1 = d 2^s with d odd (SIEVE_MONT32_SETUP)
static bool SieveSprp32(const uint n, const uint ninv, const uint one, const uint r2, const uint s, const uint d, const uint a)
{
	// x = a^d in Montgomery form
	uint base = SieveMontMul32(a, r2, n, ninv);
	uint x = one;
	uint e = d;
	while (e)
	{
		if (e & 1)
		{
			x = SieveMontMul32(x, base, n, ninv);
		}
		base = SieveMontMul32(base, base, n, ninv);
		e >>= 1;
	}
	const uint minus_one = n - one;
	if ((x == one) || (x == minus_one))
	{
		return true;
	}
	for (uint i = 1; i < s; ++i)
	{
		x = SieveMontMul32(x, x, n, ninv);
		if (x == minus_one)
		{
			return true;
		}
	}
	return false;
}

#define SIEVE_MONT32_SETUP(n) \
	uint ninv = n; /* n * n = 1 mod 8 for odd n: 3 bits, then 6, 12, 24, 48 */ \
	ninv *= 2 - n * ninv; \
	ninv *= 2 - n * ninv; \
	ninv *= 2 - n * ninv; \
	ninv *= 2 - n * ninv; \
	ninv = 0 - ninv; \
	const uint s = 31 - clz((n - 1) & (1 - n)); \
	const uint d = (n - 1) >> s; \
	const uint one = (0u - n) % n;                  /* 2^32 mod n */ \
	const uint r2 = (uint)(((ulong)one * one) % n); /* 2^64 mod n */

// Deterministic Miller-Rabin test for odd n < 2^32, n > 7, not divisible by 3, 5, 7: bases 2, 7, 61 (deterministic below 4759123141,
// Jaeschke), in 32-bit Montgomery arithmetic (the 64-bit test with 4 bases took 5.9 s of k4 in SieveTestDeferred: its p are < 2^32)
static bool SieveIsPrime32(const uint n)
{
	SIEVE_MONT32_SETUP(n)
	// (a loop over the bases: one copy of the round where the test is inlined, three copies made phase 1 slower)
	for (int b = 0; b < 3; ++b)
	{
		const uint a = (b == 0) ? 2 : ((b == 1) ? 7 : 61);
		if ((a != n) && !SieveSprp32(n, ninv, one, r2, s, d, a))
		{
			return false;
		}
	}
	return true;
}

// The same for the survivors of the sieve (SieveTestDeferred): base 2, then the list of the base-2 strong pseudoprimes below 2^32 (spsp2,
// sorted; tools/cuharness --spsp2 finds them with SieveSprp32 and SieveIsPrime32, which agrees with a sieve for all n < 2^32). With three
// bases almost every warp ran all three (about half of the numbers are prime), with the list one base and a binary search.
static bool SieveIsPrime32List(const uint n, __global const uint* spsp2, const uint spsp2Count)
{
	SIEVE_MONT32_SETUP(n)
	if (!SieveSprp32(n, ninv, one, r2, s, d, 2))
	{
		return false;
	}
	uint lo = 0;
	uint hi = spsp2Count;
	while (lo < hi)
	{
		const uint mid = (lo + hi) >> 1;
		if (spsp2[mid] < n)
		{
			lo = mid + 1;
		}
		else
		{
			hi = mid;
		}
	}
	return (lo == spsp2Count) || (spsp2[lo] != n);
}

)==="
R"===(
// Deterministic Miller-Rabin test for n < 1122004669633 (bases 2, 13, 23, 1662803; n < 2^32: SieveIsPrime32)
static bool SieveIsPrime(const ulong n)
{
	if (n < 2)
	{
		return false;
	}
	if (n < (1UL << 32))
	{
		const uint n32 = (uint)n;
		if ((n32 == 2) || (n32 == 3) || (n32 == 5) || (n32 == 7))
		{
			return true;
		}
		if (((n32 & 1) == 0) || (n32 % 3 == 0) || (n32 % 5 == 0) || (n32 % 7 == 0))
		{
			return false;
		}
		return SieveIsPrime32(n32);
	}
	if ((n == 2) || (n == 3) || (n == 5) || (n == 7))
	{
		return true;
	}
	if (((n & 1) == 0) || (n % 3 == 0) || (n % 5 == 0) || (n % 7 == 0))
	{
		return false;
	}

	ulong ninv = n;
	for (int i = 0; i < 5; ++i)
	{
		ninv *= 2 - n * ninv;
	}
	ninv = 0 - ninv;

	ulong d = n - 1;
	uint s = 0;
	while ((d & 1) == 0)
	{
		d >>= 1;
		++s;
	}

	const ulong one = (0 - n) % n; // 2^64 mod n
	const ulong minus_one = n - one;
	const ulong bases[4] = { 2, 13, 23, 1662803 };
	for (int b = 0; b < 4; ++b)
	{
		const ulong a = bases[b] % n;
		if (a == 0)
		{
			continue;
		}

		// x = a^d in Montgomery form
		const ulong am = (a * one) % n; // a < 2^21, one < n < 2^41
		ulong x = one;
		ulong base = am;
		ulong e = d;
		while (e)
		{
			if (e & 1)
			{
				x = SieveMontMul(x, base, n, ninv);
			}
			base = SieveMontMul(base, base, n, ninv);
			e >>= 1;
		}

		if ((x == one) || (x == minus_one))
		{
			continue;
		}

		bool composite = true;
		for (uint i = 1; i < s; ++i)
		{
			x = SieveMontMul(x, x, n, ninv);
			if (x == minus_one)
			{
				composite = false;
				break;
			}
		}
		if (composite)
		{
			return false;
		}
	}
	return true;
}

)==="
R"===(
// SieveIsPrime with the list of base-2 strong pseudoprimes for n < 2^32 (SieveIsPrime32List; SieveTestDeferred)
static bool SieveIsPrimeList(const ulong n, __global const uint* spsp2, const uint spsp2Count)
{
	if (n < (1UL << 32))
	{
		const uint n32 = (uint)n;
		if ((n32 == 2) || (n32 == 3) || (n32 == 5) || (n32 == 7))
		{
			return true;
		}
		if ((n32 < 2) || ((n32 & 1) == 0) || (n32 % 3 == 0) || (n32 % 5 == 0) || (n32 % 7 == 0))
		{
			return false;
		}
		return SieveIsPrime32List(n32, spsp2, spsp2Count);
	}
	return SieveIsPrime(n);
}

// Miller-Rabin test without false negatives (false only for a composite n): deterministic for n < 2^63 (the 7 bases of Jim Sinclair
// are deterministic below 2^64); n >= 2^63 returns true (SieveMontMul needs n < 2^63, the host tests these)
static bool SieveMayBePrime64(const ulong n)
{
	if (n < 1122004669633UL)
	{
		return SieveIsPrime(n);
	}
	if (n >= (1UL << 63))
	{
		return true;
	}
	if ((n & 1) == 0)
	{
		return false;
	}

	ulong ninv = n;
	for (int i = 0; i < 5; ++i)
	{
		ninv *= 2 - n * ninv;
	}
	ninv = 0 - ninv;

	ulong d = n - 1;
	uint s = 0;
	while ((d & 1) == 0)
	{
		d >>= 1;
		++s;
	}

	const ulong one = (0 - n) % n; // 2^64 mod n
	// 2^128 mod n (n < 2^63, so the doubling doesn't overflow)
	ulong r2 = one;
	for (int i = 0; i < 64; ++i)
	{
		r2 <<= 1;
		r2 = (r2 >= n) ? (r2 - n) : r2;
	}
	const ulong minus_one = n - one;
	const ulong bases[7] = { 2, 325, 9375, 28178, 450775, 9780504, 1795265022 };
	for (int b = 0; b < 7; ++b)
	{
		const ulong a = bases[b] % n;
		if (a == 0)
		{
			continue;
		}

		ulong x = one;
		ulong base = SieveMontMul(a, r2, n, ninv); // a * 2^64 mod n
		ulong e = d;
		while (e)
		{
			if (e & 1)
			{
				x = SieveMontMul(x, base, n, ninv);
			}
			base = SieveMontMul(base, base, n, ninv);
			e >>= 1;
		}

		if ((x == one) || (x == minus_one))
		{
			continue;
		}

		bool composite = true;
		for (uint i = 1; i < s; ++i)
		{
			x = SieveMontMul(x, x, n, ninv);
			if (x == minus_one)
			{
				composite = false;
				break;
			}
		}
		if (composite)
		{
			return false;
		}
	}
	return true;
}

// The cofactor B must be prime for A to be amicable: it is only reported if it may be prime (the host tests it again). In some families
// of ranges (sigma(K) * (R + 1) == T for all p, with N = K * R) every number with a rough R gets here, so without this test the host
// had to test millions of numbers per work unit (and the buffer of found numbers overflowed).
static void NumberFoundIfPrime(__global uint* amicable_numbers_count, __global ulong4* amicable_numbers_data, ulong A, ulong A_high, ulong B)
{
	if (SieveMayBePrime64(B))
	{
		NumberFound(amicable_numbers_count, amicable_numbers_data, A, A_high, B);
	}
}

// N = p^2 (D = p) or N = p * q (D = p + q, p < q) with primes p and q (D = targetSum - N - 1 < 2^64): 2 if N is one of them (*B: a factor
// for the host to test again, all factors passed SieveMayBePrime64), 1 if N has such a divisor pair with a composite factor (it is not the
// other member of a pair then), 0 if neither form fits
static uint DecideTwoFactors(const ulong N, const ulong N_high, const ulong D, ulong* B)
{
	if ((D * D == N) && (mul_hi(D, D) == N_high))
	{
		*B = D;
		return SieveMayBePrime64(D) ? 2 : 1;
	}
	// (q - p)^2 = D^2 - 4N
	ulong2 E;
	E.x = D * D;
	E.y = mul_hi(D, D);
	const ulong N4 = N << 2;
	const ulong N4_high = (N_high << 2) | (N >> 62);
	if ((E.y > N4_high) || ((E.y == N4_high) && (E.x > N4)))
	{
		E.y -= N4_high + ((E.x < N4) ? 1 : 0);
		E.x -= N4;
		const ulong sq = IntegerSquareRoot128(E);
		if ((sq * sq == E.x) && (mul_hi(sq, sq) == E.y))
		{
			const ulong p = (D - sq) >> 1;
			const ulong q = (D + sq) >> 1;
			*B = q;
			return (SieveMayBePrime64(p) && SieveMayBePrime64(q)) ? 2 : 1;
		}
	}
	return 0;
}

// The band sieve (kernel_sieve.cl: SieveBand, for the ranges with W ~ 0): band = rmin | begin << 21 | end << 42 (indexes in smallPrimes, 0:
// no band), the primes of the band are smallPrimes[begin, end), rmin is the index of the smallest of them that divides N (end or more:
// none). If N is the other member of a pair with D = targetSum - N - 1 > 0, every prime factor of N is > L = N / D (D > N / r), N has at
// most K = floor(ln N / ln L) of them, and its smallest one r is in (L, K L + K + 2]: with k >= 3 factors r - 1 <= k L / (1 - 1/(2L))
// (D < N ((1 + 1/(r - 1))^k - 1)), with 2 (p q, p < q, or p^2) L = p q / (p + q) is in [p/2, p). So: a band prime of N <= L means no
// match; if the band covers (L, ...] (L >= the prime below smallPrimes[begin]), the smallest prime factor of N above L is rmin, and if
// it is above K L + K + 2 (or there is none and the band reaches beyond K L + K + 2), N is not a match. Otherwise the trial division can
// start at rmin (*startIndex). Returns 0 (phase 1 continues as before) or 1 (N is not the other member of a pair).
static uint DecideByBand(const ulong N, const ulong N_high, const ulong targetSum, const ulong targetSumHigh, const ulong band, __global const uint* smallPrimes, uint* startIndex)
{
	const ulong D = targetSum - N - 1;
	const ulong D_high = targetSumHigh - N_high - (((targetSum < N + 1) || (N + 1 == 0)) ? 1 : 0);
	if ((band == 0) || (D_high != 0) || (D == 0))
	{
		return 0;
	}
	const uint rmin = (uint)(band & 0x1FFFFF);
	const uint begin = (uint)((band >> 21) & 0x1FFFFF);
	const uint end = (uint)(band >> 42);
	// x <= L <=> x D <= N (128-bit)
	#define BAND_AT_MOST_L(x) ((mul_hi((ulong)(x), D) < N_high) || ((mul_hi((ulong)(x), D) == N_high) && ((ulong)(x) * D <= N)))
	if (rmin < end)
	{
		const uint r = smallPrimes[rmin];
		if (BAND_AT_MOST_L(r))
		{
			return 1;
		}
	}
	const uint below = smallPrimes[begin - 1];
	if (!BAND_AT_MOST_L(below))
	{
		return 0;
	}
	#undef BAND_AT_MOST_L
	// K L + K + 2, rounded up (the errors of the float operations)
	const float Nf = convert_float(N_high) * 18446744073709551616.0f + convert_float(N);
	const float L = Nf / convert_float(D);
	const uint K = max((uint)(log(Nf) / log(L) * 1.0001f + 0.0001f), 2U);
	const float bandHi = (convert_float(K) * L + convert_float(K) + 2.0f) * 1.0001f + 1.0f;
	if (rmin < end)
	{
		if (convert_float(smallPrimes[rmin]) > bandHi)
		{
			return 1;
		}
		*startIndex = rmin;
		return 0;
	}
	if (convert_float(smallPrimes[end]) > bandHi)
	{
		return 1;
	}
	*startIndex = end;
	return 0;
}

// N = r^2 s (DecideAtMost2Factors): the sign of f(r) - D for f(r) = (r + 1)(r + N / r^2), which is the sign of (r + 1)(r^3 + N) - D r^2
// (r < 2^27, N < 2^78, D < 2^64: the products are below 2^120)
static int R2S_Sign(const ulong r, const ulong N, const ulong N_high, const ulong D)
{
	const ulong r2 = r * r;
	ulong a = r2 * r;
	ulong a_high = mul_hi(r2, r);
	a += N;
	a_high += N_high + ((a < N) ? 1 : 0);
	const ulong r1 = r + 1;
	const ulong b = a * r1;
	const ulong b_high = mul_hi(a, r1) + a_high * r1;
	const ulong c = D * r2;
	const ulong c_high = mul_hi(D, r2);
	if (b_high != c_high)
	{
		return (b_high > c_high) ? 1 : -1;
	}
	return (b > c) ? 1 : ((b < c) ? -1 : 0);
}

// f'(r) > 0 for f(r) = (r + 1)(r + N / r^2): r^3 f'(r) = 2 r^4 + r^3 - N r - 2 N = r^3 (2 r + 1) - N (r + 2)
static bool R2S_Rising(const ulong r, const ulong N, const ulong N_high)
{
	const ulong r2 = r * r;
	const ulong r3 = r2 * r;
	const ulong r3_high = mul_hi(r2, r);
	const ulong m = 2 * r + 1;
	const ulong x = r3 * m;
	const ulong x_high = mul_hi(r3, m) + r3_high * m;
	const ulong n = r + 2;
	const ulong y = N * n;
	const ulong y_high = mul_hi(N, n) + N_high * n;
	return (x_high > y_high) || ((x_high == y_high) && (x > y));
}

// Numbers with few prime factors. D = targetSum - N - 1 must be sigma(N) - N - 1, the sum of the divisors of N between 1 and N. If N has
// 4 or more prime factors, D >= N / p >= N^(3/4) (p: the smallest prime factor, N / p is such a divisor). With 3 different prime factors,
// D > r1 r2 + r1 r3 + r2 r3 >= 3 N^(2/3) (AM-GM), with N = r^2 s (r != s), D = (r + 1)(r + s) > r^2 + N / r >= 1.8899 N^(2/3), and
// N = r^3 has D = r + r^2. So if D < 3 floor(N^(1/3))^2 (and N >= 2^20: N^(3/4) > 3 N^(2/3)), N can only be the other member of a pair if
// it is prime (D = 0), p^2 (D = p), p * q (D = p + q, p < q), r^3 (D = r + r^2) or r^2 s (D = (r + 1)(r + s)) with primes p, q, r, s, and
// that is decided here without trial divisions. In the ranges where sigma(N) is close to N + 1 for all candidates (W ~ 0: sigma(K) (R + 1)
// ~ T for N = K * R), most numbers that pass phase 1 are decided by this: the trial divisions up to N^(1/4) in phase 2 and N^(1/3) in
// phase 3 were most of the time of these work units (t1).
// withR2S: the form r^2 s and D up to 3 N^(2/3) (only in Phase2Prefilter: its code made the kernels of phases 1-3 use more registers);
// without it, D < 1.88 N^(2/3) (the other forms).
// Returns 0 if it doesn't apply (D too large, or D = 0 and N >= 2^63: no deterministic test), 1 if N is not the other member of a pair,
// 2 if it is (*B: a prime factor of N for the host to test again; all prime factors of N passed the deterministic test here).
static uint DecideAtMost2Factors(const ulong N, const ulong N_high, const ulong targetSum, const ulong targetSumHigh, ulong* B, const bool withR2S)
{
	const ulong D = targetSum - N - 1;
	const ulong D_high = targetSumHigh - N_high - (((targetSum < N + 1) || (N + 1 == 0)) ? 1 : 0);
	if (D_high != 0)
	{
		return 0;
	}
	if (D == 0)
	{
		if ((N_high != 0) || (N >= (1UL << 63)))
		{
			return 0;
		}
		*B = N;
		return SieveMayBePrime64(N) ? 2 : 1;
	}

	// A float test first: most numbers have D >> N^(2/3) (a wrong result of the float test can only make this return 0)
	const float rf = cbrt(convert_float(N_high) * 18446744073709551616.0f + convert_float(N));
	if (convert_float(D) > rf * rf * (withR2S ? 3.03f : 1.9f))
	{
		return 0;
	}
	// r = floor(N^(1/3)) < 2^26 (N < 2^78), then D < 3 r^2 exactly (D < 2^54 here)
	ulong r = (ulong)rf;
	for (;;)
	{
		const ulong r2 = r * r;
		const ulong c_high = mul_hi(r2, r);
		if ((c_high > N_high) || ((c_high == N_high) && (r2 * r > N)))
		{
			--r;
			continue;
		}
		const ulong r1 = r + 1;
		const ulong r12 = r1 * r1;
		const ulong c1_high = mul_hi(r12, r1);
		if ((c1_high < N_high) || ((c1_high == N_high) && (r12 * r1 <= N)))
		{
			r = r1;
			continue;
		}
		break;
	}
	const ulong r2 = r * r;
	if ((withR2S ? (D >= r2 * 3) : (D * 25 >= r2 * 47)) || ((N_high == 0) && (N < (1UL << 20))))
	{
		return 0;
	}

	// N = r^3: D = r + r^2 (a cube has no other of these forms)
	if ((r2 * r == N) && (mul_hi(r2, r) == N_high))
	{
		*B = r;
		return ((D == r + r2) && SieveMayBePrime64(r)) ? 2 : 1;
	}

	// N = p^2 or p * q (if N has such a pair of divisors, it is not r^2 s with D = (r + 1)(r + s))
	{
		const uint d2 = DecideTwoFactors(N, N_high, D, B);
		if (d2 != 0)
		{
			return d2;
		}
	}

	// N = r^2 s (D >= 1.8899 N^(2/3)): D = f(r) = (r + 1)(r + N / r^2). f is convex for r > 0, with its minimum near (N / 2)^(1/3), and
	// D < 3 N^(2/3) puts both solutions in (0.34, 1.54) N^(1/3). Exact bisections (R2S_Sign: the sign of f(r) - D, R2S_Rising: f'(r) > 0)
	// find the integers next to them.
	if (withR2S && (D * 25 >= r2 * 47))
	{
		ulong lo = max((r * 34) / 100, 2UL);
		ulong hi = (r * 154) / 100 + 2;
		// the minimum: the first r in [lo, hi] with f'(r) > 0
		ulong a = lo;
		ulong b = hi;
		while (a < b)
		{
			const ulong m = (a + b) >> 1;
			if (R2S_Rising(m, N, N_high))
			{
				b = m;
			}
			else
			{
				a = m + 1;
			}
		}
		const ulong r0 = a;
		ulong candidates[2];
		// f decreases in [lo, r0 - 1] (r0 - 1 <= the real minimum): the last r with f(r) >= D
		a = lo;
		b = (r0 > lo) ? (r0 - 1) : lo;
		while (a < b)
		{
			const ulong m = (a + b + 1) >> 1;
			if (R2S_Sign(m, N, N_high, D) >= 0)
			{
				a = m;
			}
			else
			{
				b = m - 1;
			}
		}
		candidates[0] = a;
		// f increases in [r0, hi]: the first r with f(r) >= D
		a = r0;
		b = hi;
		while (a < b)
		{
			const ulong m = (a + b) >> 1;
			if (R2S_Sign(m, N, N_high, D) >= 0)
			{
				b = m;
			}
			else
			{
				a = m + 1;
			}
		}
		candidates[1] = a;
		for (int j = 0; j < 2; ++j)
		{
			const ulong x = candidates[j];
			if (R2S_Sign(x, N, N_high, D) != 0)
			{
				continue;
			}
			// D = (x + 1)(x + s), N = x^2 s
			const ulong s = D / (x + 1) - x;
			const ulong x2 = x * x;
			if (((x + 1) * (x + s) == D) && (x2 * s == N) && (mul_hi(x2, s) == N_high) && (s != x))
			{
				*B = s;
				return (SieveMayBePrime64(x) && SieveMayBePrime64(s)) ? 2 : 1;
			}
		}
	}
	return 1;
}

// The primes that can divide N if it is the other member of a pair. If D = targetSum - N - 1 > 0, N is composite, and for every prime
// factor r of N, N / r is a divisor of N between 1 and N: D = sigma(N) - N - 1 > N / r, so r > N / D. Returns the index of a prime
// <= N / D in smallPrimes, or 0: the trial divisions of phases 2 and 3 start there (the primes below N / D are skipped, so a number that
// passes all tests after that may still have a small factor: every report of phases 2 and 3 tests all prime factors that it assumes).
// In the ranges where sigma(N) is close to N + 1 for all candidates (W ~ 0), N / D is 3e4 - 3e6, and the trial divisions below it were
// most of the time of these work units.
// bucketEdge (phase 2): N / D is rounded down to the lower edge 5393 * 2^(c/3) of its bucket of SortPhase2Numbers (Phase2Bucket), so that
// the work-items of a warp (the numbers are sorted by these buckets) start at the same prime and load the same inverses: with a different
// start in every work-item, phase 2 was limited by the L1 cache (each load used half of its sector).
static int FirstUsefulPrimeIndex(const ulong N, const ulong N_high, const ulong targetSum, const ulong targetSumHigh, const bool bucketEdge)
{
	const ulong D = targetSum - N - 1;
	const ulong D_high = targetSumHigh - N_high - (((targetSum < N + 1) || (N + 1 == 0)) ? 1 : 0);
	if ((D | D_high) == 0)
	{
		return 0;
	}
	float L = (convert_float(N_high) * 18446744073709551616.0f + convert_float(N)) / (convert_float(D_high) * 18446744073709551616.0f + convert_float(D));
	if (!(L >= 5393.0f))
	{
		return 0;
	}
	if (bucketEdge)
	{
		L = 5393.0f * exp2(floor(3.0f * (log2(L) - 12.397f)) * (1.0f / 3.0f));
	}
	// pi(x) >= x / (ln x - 1) for x >= 5393 (Dusart), 0.99 for the errors of the float operations; smallPrimes[i] is the (i+1)-th prime
	const float index = 0.99f * L / (log(L) - 1.0f) - 1.0f;
	return (int)min(index, 1048576.0f);
}

)==="
R"===(
// A band prime r = smallPrimes[index] of N (DecideByBand; index >= 1: PrimeInverses_128[index - 1] and PowersOfP_128SumInverses belong to it,
// as in the loops of phases 1 and 2): if N is the other member of a pair, sigma(r^e) divides targetSum and N / r^e < targetSum / sigma(r^e)
// (or both are 1).
// Most candidates of the band fail this (it is what phase 2 finds after the trial division up to r). Returns false if N is not a match.
static bool BandPrimeDivides(ulong N, ulong N_high, ulong targetSum, ulong targetSumHigh, const uint index, __global const ulong2* PrimeInverses_128,
	__global const uint* PowersOfP_128SumInverses_offsets, __global const ulong4* PowersOfP_128SumInverses)
{
	const ulong2 inv = PrimeInverses_128[index - 1];
	ulong q_high = mul_hi(N, inv.x) + N * inv.y + N_high * inv.x;
	if (q_high > N_high)
	{
		// (r doesn't divide N: no test)
		return true;
	}
	ulong q = N * inv.x;
	int powerOfP = 0;
	do
	{
		++powerOfP;
		N = q;
		N_high = q_high;
		q = N * inv.x;
		q_high = mul_hi(N, inv.x) + N * inv.y + N_high * inv.x;
	} while (q_high <= N_high);
	const ulong4 inv2 = PowersOfP_128SumInverses[PowersOfP_128SumInverses_offsets[index - 1] + powerOfP - 1];
	if (inv2.z != 0)
	{
		if ((targetSum & inv2.w) != 0)
		{
			return false;
		}
		targetSum = (targetSum >> inv2.z) | (targetSumHigh << (64 - inv2.z));
		targetSumHigh >>= inv2.z;
	}
	q_high = mul_hi(targetSum, inv2.x) + targetSum * inv2.y + targetSumHigh * inv2.x;
	if (q_high > targetSumHigh)
	{
		return false;
	}
	targetSum *= inv2.x;
	targetSumHigh = q_high;
	// (N = r^e: a match if targetSum is 1 now; the comparison below rejected it, 2026-11-08)
	if ((N_high == 0) && (N == 1))
	{
		return (targetSumHigh == 0) && (targetSum == 1);
	}
	return (N_high < targetSumHigh) || ((N_high == targetSumHigh) && (N < targetSum));
}

// deferValue: 0, or (sieve survivors, p not known to be prime) the survivor record (p << 27 | range index). Then a number that
// passes all checks of this phase is not reported or saved for phase 2 here: its survivor record goes to the deferred list, and
// SieveTestDeferred tests p for all of them in parallel (a primality test in one work-item would stall its whole warp); the
// survivors with a prime p are checked again without deferral. A survivor is deferred at most once, so the list (as big as the
// survivor buffer) can't overflow.
#ifndef SIEVE_DEFERRED_MAX
#define SIEVE_DEFERRED_MAX 1
#endif

static void Phase1Defer(__global uint* deferred_count, __global ulong* deferred_list, const ulong value)
{
	const uint k = atom_inc(deferred_count);
	if (k < SIEVE_DEFERRED_MAX)
	{
		deferred_list[k] = value;
	}
}

#define PHASE1_FOUND() do { \
		if (deferValue != 0) \
		{ \
			Phase1Defer(deferred_count, deferred_list, deferValue); \
		} \
		else \
		{ \
			NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, 0); \
		} \
	} while (0)

static void CheckPairPhase1(
	__global const uint* smallPrimes,
	__global const ulong2* primeInverses,
	__global const ulong2* PQ,
	__constant ulong2* PowerOf2SumInverses,
	__global const uint* PowersOfP_128SumInverses_offsets,
	__global const ulong4* PowersOfP_128SumInverses,
	__global const ulong2* PrimeInverses_128,
	__global const ulong* SumEstimates_128,
	const ulong M,
	const ulong M_high,
	ulong targetSum,
	ulong targetSumHigh,
	const ulong global_offset,
	__global ulong* phase1_offset_to_resume_after_overflow,
	__global uint* phase2_numbers_count,
	__global ulong4* phase2_numbers_data,
	__global uint* amicable_numbers_count,
	__global ulong4* amicable_numbers_data,
		const ulong deferValue,
		__global uint* deferred_count,
		__global ulong* deferred_list,
		const ulong band
	)
{
	ulong N = targetSum - M;
	ulong N_high = targetSumHigh - M_high - ((targetSum < M) ? 1 : 0);

	// Divide out power of 2 using fast bitwise operations
	// All numbers in the same range have the same parity, so threads won't diverge here
	if ((N & 1) == 0)
	{
		if (N == 0)
		{
			// N = 2^64*k, k >= 0. It's not amicable for k = 0, 1, 2, ..., 5 and too large for k > 5.
			return;
		}

		const uint powerOf2 = 63 - clz(N ^ (N - 1));

		// M and targetSum (128-bit) come in ranges where they both grow monotonically,
		// so there will be exactly one divergence in each range and one divergence on a border between ranges
		// Threads don't diverge 99.9% of the time here
		if (targetSumHigh)
		{
			const ulong2 inv = PowerOf2SumInverses[powerOf2 + 64];
			const ulong q_high = mul_hi(targetSum, inv.x) + targetSum * inv.y + targetSumHigh * inv.x;

			// The proper check would be "if (128-bit value of q > (2^128 - 1) / divisor)"
			// "if (q_high > targetSumHigh)" doesn't give false positives because the result of division (by divisor > 1) can't be larger than the original value
			// It can give very rare false negatives when the divisor is very large (larger than 2^64 / targetSumHigh), but false negatives don't jeopardize the correctness of the search here
			// Checking only highest 64 bits and not having to store or calculate "(2^128 - 1) / divisor" is much faster on average, so it's worth it
			if (q_high > targetSumHigh)
			{
				return;
			}

			targetSum *= inv.x;
			targetSumHigh = q_high;

			N = (N >> powerOf2) | (N_high << (64 - powerOf2));
			N_high >>= powerOf2;

			if ((N_high > targetSumHigh) || ((N_high == targetSumHigh) && (N >= targetSum)))
			{
				return;
			}
		}
		else
		{
			targetSum *= PowerOf2SumInverses[powerOf2].x;
			if (targetSum > PowerOf2SumInverses[powerOf2].y)
			{
				return;
			}
			N >>= powerOf2;

			if (N >= targetSum)
			{
				return;
			}
		}
	}

	int i = 0;

	// M and targetSum (128-bit) come in ranges where they both grow monotonically,
	// so there will be exactly one divergence in each range and one divergence on a border between ranges
	// Threads don't diverge 99.9% of the time here
	if (targetSumHigh)
	{
		for (int cur_index = 0; cur_index < PHASE1_DEPTH; cur_index += 16)
		{
			const int e = cur_index + 16;
			for (i = cur_index + 1; i <= e; ++i)
			{
				const ulong2 inv = PrimeInverses_128[i - 1];
				ulong q_high = mul_hi(N, inv.x) + N * inv.y + N_high * inv.x;
				if (q_high > N_high)
				{
					continue;
				}

				ulong q = N * inv.x;
				int powerOfP = 0;
				do
				{
					++powerOfP;

					N = q;
					N_high = q_high;

					q = N * inv.x;
					q_high = mul_hi(N, inv.x) + N * inv.y + N_high * inv.x;
				} while (q_high <= N_high);

				const ulong4 inv2 = PowersOfP_128SumInverses[PowersOfP_128SumInverses_offsets[i - 1] + powerOfP - 1];
				if (inv2.z != 0)
				{
					if ((targetSum & inv2.w) != 0)
					{
						return;
					}
					targetSum = (targetSum >> inv2.z) | (targetSumHigh << (64 - inv2.z));
					targetSumHigh >>= inv2.z;
				}

				q_high = mul_hi(targetSum, inv2.x) + targetSum * inv2.y + targetSumHigh * inv2.x;
				if (q_high > targetSumHigh)
				{
					return;
				}

				targetSum *= inv2.x;
				targetSumHigh = q_high;

				const ulong p = smallPrimes[i];
				if ((N_high == 0) && (p * p > N))
				{
					// (N_high: N can still be >= 2^64 here, the test of the low 64 bits alone made rare numbers return)
					const ulong sumN = ((N > 1) ? N : 0) + 1;
					if (sumN == targetSum)
					{
						PHASE1_FOUND();
					}
					return;
				}

				if ((N_high > targetSumHigh) || ((N_high == targetSumHigh) && (N >= targetSum)))
				{
					return;
				}

				if (!targetSumHigh)
				{
					const int i1 = i + 1;
					i = i1 - (i1 & 15);
#ifdef DISABLE_GOTO
					break;
#else
					goto small_numbers;
#endif
				}
			}

#ifdef DISABLE_GOTO
			if (!targetSumHigh)
			{
				break;
			}
#endif

			ulong2 max_sum;
			const ulong max_sum_ratio = SumEstimates_128[cur_index >> 4];

#if SUM_ESTIMATES_128_SHIFT == 64
			max_sum.x = mul_hi(N, max_sum_ratio) + N_high * max_sum_ratio;
			max_sum.y = 0;
#else
			max_sum.x = N * max_sum_ratio;
			max_sum.y = mul_hi(N, max_sum_ratio) + N_high * max_sum_ratio;
			max_sum.x = (max_sum.x >> SUM_ESTIMATES_128_SHIFT) | (max_sum.y << (64 - SUM_ESTIMATES_128_SHIFT));
			max_sum.y >>= SUM_ESTIMATES_128_SHIFT;
#endif

			const ulong carry = (max_sum.x + N < max_sum.x) ? 1 : 0;
			max_sum.x = max_sum.x + N;
			max_sum.y = max_sum.y + N_high + carry;

						if ((max_sum.y < targetSumHigh) || ((max_sum.y == targetSumHigh) && (max_sum.x < targetSum)))
						{
							return;
						}

						// If no other prime of phase 1 can divide N (FirstUsefulPrimeIndex), phase 2 starts at the first useful prime
						if ((cur_index == 0) && (FirstUsefulPrimeIndex(N, N_high, targetSum, targetSumHigh, false) >= PHASE1_DEPTH))
						{
							i = PHASE1_DEPTH;
							break;
						}
					}

#ifdef DISABLE_GOTO
					if (targetSumHigh)
#endif
		{
			// (a found number is deferred like the others if p is not known to be prime yet)
			ulong B = 0;
			uint decided = DecideAtMost2Factors(N, N_high, targetSum, targetSumHigh, &B, false);
			uint startIndex = 0;
			if (decided == 0)
			{
				decided = DecideByBand(N, N_high, targetSum, targetSumHigh, band, smallPrimes, &startIndex);
				// (startIndex is the smallest band prime of N if it is the rmin field of band)
				if ((decided == 0) && (startIndex != 0) && (startIndex == (uint)(band & 0x1FFFFF)) &&
					!BandPrimeDivides(N, N_high, targetSum, targetSumHigh, startIndex, PrimeInverses_128, PowersOfP_128SumInverses_offsets, PowersOfP_128SumInverses))
				{
					decided = 1;
				}
			}
			if (decided == 1)
			{
				return;
			}
			if (deferValue != 0)
			{
				Phase1Defer(deferred_count, deferred_list, deferValue);
				return;
			}
			if (decided == 2)
			{
				NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, B);
				return;
			}
			const uint phase2_number_index = atom_inc(phase2_numbers_count);
			if (phase2_number_index < PHASE2_MAX_COUNT)
			{
				// (the primes below startIndex can't divide N: DecideByBand; one less because the 128-bit loop of phase 2 starts after the
				// prime of its first index)
				const ulong i_to_save = max((uint)(i - (i & 15)), ((startIndex > 0) ? min(startIndex - 1, (1U << 22) - 64) : 0) & ~15U);
				phase2_numbers_data[phase2_number_index].x = M;
				phase2_numbers_data[phase2_number_index].y = N;
				phase2_numbers_data[phase2_number_index].z = targetSum;
				phase2_numbers_data[phase2_number_index].w = i_to_save | (M_high << 22) | (N_high << 36) | (targetSumHigh << 50);
			}
			else if (phase2_number_index == PHASE2_MAX_COUNT)
			{
				*phase1_offset_to_resume_after_overflow = global_offset;
			}
			return;
		}
	}

)==="
R"===(
#ifndef DISABLE_GOTO
	small_numbers:
#endif

	if (i < 16)
	{
		// Collect divisibility data for first 15 odd primes
		// No divergence here, compiler will use predicates for "A ? B : C" operator
		uint pattern =
			((N * primeInverses[ 1].x <= primeInverses[ 1].y) ? 16384 : 0) |
			((N * primeInverses[ 2].x <= primeInverses[ 2].y) ? 8192  : 0) |
			((N * primeInverses[ 3].x <= primeInverses[ 3].y) ? 4096  : 0) |
			((N * primeInverses[ 4].x <= primeInverses[ 4].y) ? 2048  : 0) |
			((N * primeInverses[ 5].x <= primeInverses[ 5].y) ? 1024  : 0) |
			((N * primeInverses[ 6].x <= primeInverses[ 6].y) ? 512   : 0) |
			((N * primeInverses[ 7].x <= primeInverses[ 7].y) ? 256   : 0) |
			((N * primeInverses[ 8].x <= primeInverses[ 8].y) ? 128   : 0) |
			((N * primeInverses[ 9].x <= primeInverses[ 9].y) ? 64    : 0) |
			((N * primeInverses[10].x <= primeInverses[10].y) ? 32    : 0) |
			((N * primeInverses[11].x <= primeInverses[11].y) ? 16    : 0) |
			((N * primeInverses[12].x <= primeInverses[12].y) ? 8     : 0) |
			((N * primeInverses[13].x <= primeInverses[13].y) ? 4     : 0) |
			((N * primeInverses[14].x <= primeInverses[14].y) ? 2     : 0) |
			((N * primeInverses[15].x <= primeInverses[15].y) ? 1     : 0);

		i = 16;

		// Do trial divisions by primes found
		if (pattern)
		{
			ulong sumN = 1;
			do
			{
				// Bit 0: divisibility for i == 15
				// Bit 1: divisibility for i == 14
				// Bit 2: divisibility for i == 13
				// and so on ...
				const uint num_leading_zeroes = clz(pattern);
				pattern -= (1U << (31 - num_leading_zeroes));
				const uint index = num_leading_zeroes - 16;

				const ulong prevSumN = sumN;
				const uint p1 = smallPrimes[index];
				ulong q = N * primeInverses[index].x;
				do
				{
					sumN = sumN * p1 + prevSumN;
					N = q;
					q *= primeInverses[index].x;
				} while (q <= primeInverses[index].y);
			} while (pattern);

			if (N < 59 * 59)
			{
				if (N > 1)
				{
					sumN *= N + 1;
				}
				if (sumN == targetSum)
				{
					PHASE1_FOUND();
				}
				return;
			}

			const ulong oldTargetSum = targetSum;
			targetSum /= sumN;
			if (targetSum * sumN != oldTargetSum)
			{
				return;
			}

			if (N >= targetSum)
			{
				return;
			}
		}
	}

	// Check maximal possible sum for N. If it's too small, it can't be an amicable number
	int k = 8;
	while (N <= PQ[k].x)
	{
		--k;
	}
	++k;
	if (mul_hi(N, PQ[k].y) < targetSum - N)
	{
		return;
	}

	// If no other prime of phase 1 can divide N (FirstUsefulPrimeIndex), phase 2 starts at the first useful prime
	if ((i < PHASE1_DEPTH) && (FirstUsefulPrimeIndex(N, N_high, targetSum, targetSumHigh, false) >= PHASE1_DEPTH))
	{
		i = PHASE1_DEPTH;
	}

	if (i < PHASE1_DEPTH)
	{
		do
		{
			// Collect divisibility data for the next 16 primes
			// No divergence here, compiler will use predicates for "A ? B : C" operator
			uint pattern =
				((N * primeInverses[i +  0].x <= primeInverses[i +  0].y) ? 32768 : 0) |
				((N * primeInverses[i +  1].x <= primeInverses[i +  1].y) ? 16384 : 0) |
				((N * primeInverses[i +  2].x <= primeInverses[i +  2].y) ? 8192  : 0) |
				((N * primeInverses[i +  3].x <= primeInverses[i +  3].y) ? 4096  : 0) |
				((N * primeInverses[i +  4].x <= primeInverses[i +  4].y) ? 2048  : 0) |
				((N * primeInverses[i +  5].x <= primeInverses[i +  5].y) ? 1024  : 0) |
				((N * primeInverses[i +  6].x <= primeInverses[i +  6].y) ? 512   : 0) |
				((N * primeInverses[i +  7].x <= primeInverses[i +  7].y) ? 256   : 0) |
				((N * primeInverses[i +  8].x <= primeInverses[i +  8].y) ? 128   : 0) |
				((N * primeInverses[i +  9].x <= primeInverses[i +  9].y) ? 64    : 0) |
				((N * primeInverses[i + 10].x <= primeInverses[i + 10].y) ? 32    : 0) |
				((N * primeInverses[i + 11].x <= primeInverses[i + 11].y) ? 16    : 0) |
				((N * primeInverses[i + 12].x <= primeInverses[i + 12].y) ? 8     : 0) |
				((N * primeInverses[i + 13].x <= primeInverses[i + 13].y) ? 4     : 0) |
				((N * primeInverses[i + 14].x <= primeInverses[i + 14].y) ? 2     : 0) |
				((N * primeInverses[i + 15].x <= primeInverses[i + 15].y) ? 1     : 0);

			i += 16;

			// Do trial divisions by prime factors found
			if (pattern)
			{
				ulong sumN = 1;
				do
				{
					// Bit 0: divisibility for i - 1
					// Bit 1: divisibility for i - 2
					// Bit 2: divisibility for i - 3
					// and so on ...
					const uint num_leading_zeroes = clz(pattern);
					pattern -= (1U << (31 - num_leading_zeroes));
					const uint index = i + num_leading_zeroes - 32;

					const ulong prevSumN = sumN;
					const uint p1 = smallPrimes[index];
					ulong q = N * primeInverses[index].x;
					do
					{
						sumN = sumN * p1 + prevSumN;
						N = q;
						q *= primeInverses[index].x;
					} while (q <= primeInverses[index].y);
				} while (pattern);

				const ulong p = smallPrimes[i];
				if (p * p > N)
				{
					if (N > 1)
					{
						sumN *= N + 1;
					}
					if (sumN == targetSum)
					{
						PHASE1_FOUND();
					}
					return;
				}

				const ulong oldTargetSum = targetSum;
				targetSum /= sumN;
				if (targetSum * sumN != oldTargetSum)
				{
					return;
				}

				if (N >= targetSum)
				{
					return;
				}
			}

			// Check maximal possible sum for N. If it's too small, it can't be an amicable number
			// Same i for all threads, no divergence
			uint index = k * PQ_STRIDE_SIZE + (i >> 4) + 16;
			while (N <= PQ[index].x)
			{
				--k;
				index -= PQ_STRIDE_SIZE;
			}
			if (mul_hi(N, PQ[index].y) < targetSum - N)
			{
				return;
			}
		} while (i != PHASE1_DEPTH);
	}
	else
	{
		uint index = k * PQ_STRIDE_SIZE + (i >> 4) + 16;
		while (N <= PQ[index].x)
		{
			--k;
			index -= PQ_STRIDE_SIZE;
		}
		if (mul_hi(N, PQ[index].y) < targetSum - N)
		{
			return;
		}
	}

	// Data needed to resume from arbitrary loop iteration in phase 2:
	// 1) M (8 bytes + 14 bits)
	// 2) N (8 bytes + 14 bits)
	// 3) targetSum (8 bytes + 14 bits)
	// 4) i (<= 192724, can fit in 18 bits)
	// 5) k (< 16, can fit in 4 bits)
	// It can all fit in 32 bytes
	// 8 bytes + 14 bits is enough to store numbers up to ~3.022 * 10^23, so 10^23 is the maximum search limit for this version

	// (a found number is deferred like the others if p is not known to be prime yet)
	ulong B = 0;
	uint decided = DecideAtMost2Factors(N, N_high, targetSum, targetSumHigh, &B, false);
	uint startIndex = 0;
	if (decided == 0)
	{
		decided = DecideByBand(N, N_high, targetSum, targetSumHigh, band, smallPrimes, &startIndex);
		// (startIndex is the smallest band prime of N if it is the rmin field of band)
		if ((decided == 0) && (startIndex != 0) && (startIndex == (uint)(band & 0x1FFFFF)) &&
			!BandPrimeDivides(N, N_high, targetSum, targetSumHigh, startIndex, PrimeInverses_128, PowersOfP_128SumInverses_offsets, PowersOfP_128SumInverses))
		{
			decided = 1;
		}
	}
	if (decided == 1)
	{
		return;
	}
	if (deferValue != 0)
	{
		Phase1Defer(deferred_count, deferred_list, deferValue);
		return;
	}
	if (decided == 2)
	{
		NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, B);
		return;
	}
	const uint phase2_number_index = atom_inc(phase2_numbers_count);
	if (phase2_number_index < PHASE2_MAX_COUNT)
	{
		// (the primes below startIndex can't divide N: DecideByBand; 18 bits; one less because the 128-bit loop of phase 2 starts after
		// the prime of its first index)
		const ulong i_to_save = max((uint)i, ((startIndex > 0) ? min(startIndex - 1, (1U << 18) - 64) : 0) & ~15U);
		const ulong k_to_save = k;
		phase2_numbers_data[phase2_number_index].x = M;
		phase2_numbers_data[phase2_number_index].y = N;
		phase2_numbers_data[phase2_number_index].z = targetSum;
		phase2_numbers_data[phase2_number_index].w = i_to_save | (k_to_save << 18) | (M_high << 22) | (N_high << 36) | (targetSumHigh << 50);
	}
	else if (phase2_number_index == PHASE2_MAX_COUNT)
	{
		*phase1_offset_to_resume_after_overflow = global_offset;
	}
}
)==="
R"===(
__kernel
__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
void CheckPairPhase2(
	__global const uint* smallPrimes,
	__global const ulong4* numbersToCheck,
	const uint numbersToCheckOffset,
	__global const ulong2* primeInverses,
	__global const uint* PowersOfP_128SumInverses_offsets,
	__global const ulong4* PowersOfP_128SumInverses,
	__global const ulong2* PrimeInverses_128,
	__global const ulong* SumEstimates_128,
	__global const ulong2* PQ,
	__global uint* filtered_numbers_count,
	__global ulong4* filtered_numbers_data,
	__global uint* amicable_numbers_count,
	__global ulong4* amicable_numbers_data
)
{
	const ulong4 cur_number = numbersToCheck[get_global_id(0) + numbersToCheckOffset];
	if (cur_number.x == 0)
	{
		return;
	}

	const ulong M = cur_number.x;
	ulong N = cur_number.y;
	ulong targetSum = cur_number.z;
	int i = cur_number.w & ((1 << 18) - 1);
	int k = (cur_number.w >> 18) & ((1 << 4) - 1);

	const ulong M_high = (cur_number.w >> 22) & ((1U << 14) - 1);
	ulong N_high = (cur_number.w >> 36) & ((1U << 14) - 1);
	ulong targetSumHigh = cur_number.w >> 50;

	if (targetSumHigh)
	{
		ulong2 N_128;
		N_128.x = N;
		N_128.y = N_high;
		uint root4_N_rounded_up = convert_uint_rtp(sqrt(convert_float_rtp(IntegerSquareRoot128(N_128))));

		// (the primes below N / D can't divide N if it is the other member of a pair: FirstUsefulPrimeIndex)
		i = max(i, FirstUsefulPrimeIndex(N, N_high, targetSum, targetSumHigh, true) & ~15);

		// N and D = targetSum - N - 1 as floats for the bound at the end of each group (Phase2Bound4)
		float Nf = convert_float(N_high) * 18446744073709551616.0f + convert_float(N);
		float Df = convert_float(targetSumHigh - N_high - (((targetSum < N + 1) || (N + 1 == 0)) ? 1 : 0)) * 18446744073709551616.0f + convert_float(targetSum - N - 1);

		for (int cur_index = i; (smallPrimes[i] <= root4_N_rounded_up) && targetSumHigh; cur_index += 16)
		{
			// The divisibility by the 16 primes smallPrimes[cur_index + 1 .. cur_index + 16] (PrimeInverses_128[cur_index .. cur_index + 15]) without
			// branches, then the divisions for the primes that divide N (rare): with the division inside the loop over the primes, the loop
			// had ~35 instructions per prime (reconvergence and control flow around the rare path), the tests alone have ~16.
			uint pattern = 0;
			if (N_high < smallPrimes[cur_index + 1])
			{
				// (all 16 primes are > N_high: Divides128Fast32 is exact; the low 32 bits of PrimeInverses_128[k].x are the inverse mod 2^32)
				__global const uint* inv32 = (__global const uint*)PrimeInverses_128;
				#pragma unroll
				for (int j = 0; j < 16; ++j)
				{
					pattern |= Divides128Fast32(N, (uint)N_high, inv32[(cur_index + j) * 4], smallPrimes[cur_index + 1 + j]) ? (1U << j) : 0U;
				}
			}
			else
			{
				#pragma unroll
				for (int j = 0; j < 16; ++j)
				{
					const ulong2 invj = PrimeInverses_128[cur_index + j];
					const ulong qh = mul_hi(N, invj.x) + N * invj.y + N_high * invj.x;
					pattern |= (qh <= N_high) ? (1U << j) : 0U;
				}
			}
			i = cur_index + 17;
			while (pattern != 0)
			{
				// the smallest prime first
				const int ii = cur_index + 1 + (31 - (int)clz(pattern & (0U - pattern)));
				pattern &= pattern - 1;
				const ulong2 inv = PrimeInverses_128[ii - 1];
				ulong q_high = mul_hi(N, inv.x) + N * inv.y + N_high * inv.x;
				if (q_high > N_high)
				{
					continue;
				}

				ulong q = N * inv.x;
				int powerOfP = 0;
				do
				{
					++powerOfP;

					N = q;
					N_high = q_high;

					q = N * inv.x;
					q_high = mul_hi(N, inv.x) + N * inv.y + N_high * inv.x;
				} while (q_high <= N_high);

				const ulong4 inv2 = PowersOfP_128SumInverses[PowersOfP_128SumInverses_offsets[ii - 1] + powerOfP - 1];
				if (inv2.z != 0)
				{
					if ((targetSum & inv2.w) != 0)
					{
						return;
					}
					targetSum = (targetSum >> inv2.z) | (targetSumHigh << (64 - inv2.z));
					targetSumHigh >>= inv2.z;
				}

				q_high = mul_hi(targetSum, inv2.x) + targetSum * inv2.y + targetSumHigh * inv2.x;
				if (q_high > targetSumHigh)
				{
					return;
				}

				targetSum *= inv2.x;
				targetSumHigh = q_high;

				const ulong p = smallPrimes[ii];
				if ((N_high == 0) && (p * p > N))
				{
					// N is 1 or prime (unless a factor below N / D was skipped: FirstUsefulPrimeIndex). (N_high: N can still be >= 2^64 here,
					// the test of the low 64 bits alone made rare numbers return)
					const ulong sumN = ((N > 1) ? N : 0) + 1;
					if (sumN == targetSum)
					{
						if (N > 1)
						{
							NumberFoundIfPrime(amicable_numbers_count, amicable_numbers_data, M, M_high, N);
						}
						else
						{
							NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, 0);
						}
					}
					return;
				}

				if ((N_high > targetSumHigh) || ((N_high == targetSumHigh) && (N >= targetSum)))
				{
					return;
				}

				N_128.x = N;
				N_128.y = N_high;
				root4_N_rounded_up = convert_uint_rtp(sqrt(convert_float_rtp(IntegerSquareRoot128(N_128))));
				Nf = convert_float(N_high) * 18446744073709551616.0f + convert_float(N);
				Df = convert_float(targetSumHigh - N_high - (((targetSum < N + 1) || (N + 1 == 0)) ? 1 : 0)) * 18446744073709551616.0f + convert_float(targetSum - N - 1);
			}

			// All prime factors of N are >= pn = smallPrimes[i] (the next prime to test; the primes before it were tested or are below N / D).
			// If pn^5 > N, N has at most 4 of them (with multiplicity), and sigma(N) / N <= the product of (1 + 1/q) over them, which is
			// largest for three factors pn and one N / pn^3 (log(1 + e^-t) is convex): sigma(N) - N <= N ((1 + 1/pn)^3 - 1) + (pn + 1)^3.
			// The table bound below (SumEstimates_128, for any number below the search limit) allows more factors: with it the loop went
			// to ~4 N / D instead of ~3.2 N / D (k4: 29% more primes). Float arithmetic with a margin of 1e-4 (the rounding errors are < 1e-6).
			const float pn = convert_float(smallPrimes[i]);
			const float pn2 = pn * pn;
			if (pn2 * pn2 * pn > Nf * 1.0001f)
			{
				const float x = 1.0f / pn;
				const float pn1 = pn + 1.0f;
				if (Df > (Nf * x * (3.0f + x * (3.0f + x)) + pn1 * pn1 * pn1) * 1.0001f)
				{
					return;
				}
				continue;
			}

			ulong2 max_sum;
			const ulong max_sum_ratio = SumEstimates_128[cur_index >> 4];

#if SUM_ESTIMATES_128_SHIFT == 64
			max_sum.x = mul_hi(N, max_sum_ratio) + N_high * max_sum_ratio;
			max_sum.y = 0;
#else
			max_sum.x = N * max_sum_ratio;
			max_sum.y = mul_hi(N, max_sum_ratio) + N_high * max_sum_ratio;
			max_sum.x = (max_sum.x >> SUM_ESTIMATES_128_SHIFT) | (max_sum.y << (64 - SUM_ESTIMATES_128_SHIFT));
			max_sum.y >>= SUM_ESTIMATES_128_SHIFT;
#endif

			const ulong carry = (max_sum.x + N < max_sum.x) ? 1 : 0;
			max_sum.x = max_sum.x + N;
			max_sum.y = max_sum.y + N_high + carry;

			if ((max_sum.y < targetSumHigh) || ((max_sum.y == targetSumHigh) && (max_sum.x < targetSum)))
			{
				return;
			}
		}

		if (targetSumHigh)
		{
			// All prime factors of N are >= p (the primes before it were tested, or are below N / D). If p^4 > N, N has at most 3 of them and
			// sigma(N) - N - 1 <= 2 * (N / p + p * (p + 1)): the bound of phase 3 (CheckPairPhase3), which rejected most of these numbers at its
			// first test (k4: 93% of the numbers of the expensive group of phase 2 went to phase 3 only for that). p^4 is checked exactly
			// (the loop above ends at a float estimate of N^(1/4)).
			{
				const ulong p = smallPrimes[i];
				const ulong p2 = p * p;
				const ulong p4_high = mul_hi(p2, p2);
				if ((p4_high > N_high) || ((p4_high == N_high) && (p2 * p2 > N)))
				{
					const ulong pmax = primeInverses[i].y; // floor((2^64 - 1) / p)
					const ulong q_upper = N_high * (pmax + 1) + mul_hi(N, pmax) + 2;
					const ulong bound = (q_upper + p * (p + 1)) << 1;
					const ulong D = targetSum - N - 1;
					const ulong D_high = targetSumHigh - N_high - ((targetSum < N + 1) || (N + 1 == 0) ? 1 : 0);
					if (D_high || (bound <= D))
					{
						return;
					}
				}
			}
			ulong B = 0;
			const uint decided = DecideAtMost2Factors(N, N_high, targetSum, targetSumHigh, &B, false);
			if (decided != 0)
			{
				if (decided == 2)
				{
					NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, B);
				}
				return;
			}
			const uint filtered_numbers_data_id = atom_inc(filtered_numbers_count);
			filtered_numbers_data[filtered_numbers_data_id].x = M;
			filtered_numbers_data[filtered_numbers_data_id].y = N;
			filtered_numbers_data[filtered_numbers_data_id].z = targetSum;
			const ulong i_to_save = i;
			filtered_numbers_data[filtered_numbers_data_id].w = i_to_save | (M_high << 22) | (N_high << 36) | (targetSumHigh << 50);
			return;
		}

			i -= (i & 15);

			k = 8;
			while (N <= PQ[k].x)
			{
				--k;
			}
			++k;
			// (i can be beyond the columns of PQ after FirstUsefulPrimeIndex: then the 64-bit loop below doesn't run, smallPrimes[i] > 2^16)
			if ((i >> 4) + 16 < PQ_STRIDE_SIZE)
			{
				uint index = k * PQ_STRIDE_SIZE + (i >> 4) + 16;
				while (k && (N <= PQ[index].x))
				{
					--k;
					index -= PQ_STRIDE_SIZE;
				}
			}
		}

	uint root4_N_rounded_up = convert_uint_rtp(sqrt(sqrt(convert_float_rtp(N))));

	// All threads start with the same value of "i" here, which is set to PHASE1_DEPTH at the end of phase 1, so there will be no divergence
	// because of different values of "i" (except for the numbers where the primes below N / D are skipped: FirstUsefulPrimeIndex)
	i = max(i, FirstUsefulPrimeIndex(N, N_high, targetSum, targetSumHigh, true) & ~15);
	if (smallPrimes[i] <= root4_N_rounded_up)
	{
		for (;;)
		{
			// Collect divisibility data for the next 16 primes
			// No divergence here, compiler will use predicates for "A ? B : C" operator
			uint pattern =
				((N * primeInverses[i +  0].x <= primeInverses[i +  0].y) ? 32768 : 0) |
				((N * primeInverses[i +  1].x <= primeInverses[i +  1].y) ? 16384 : 0) |
				((N * primeInverses[i +  2].x <= primeInverses[i +  2].y) ? 8192  : 0) |
				((N * primeInverses[i +  3].x <= primeInverses[i +  3].y) ? 4096  : 0) |
				((N * primeInverses[i +  4].x <= primeInverses[i +  4].y) ? 2048  : 0) |
				((N * primeInverses[i +  5].x <= primeInverses[i +  5].y) ? 1024  : 0) |
				((N * primeInverses[i +  6].x <= primeInverses[i +  6].y) ? 512   : 0) |
				((N * primeInverses[i +  7].x <= primeInverses[i +  7].y) ? 256   : 0) |
				((N * primeInverses[i +  8].x <= primeInverses[i +  8].y) ? 128   : 0) |
				((N * primeInverses[i +  9].x <= primeInverses[i +  9].y) ? 64    : 0) |
				((N * primeInverses[i + 10].x <= primeInverses[i + 10].y) ? 32    : 0) |
				((N * primeInverses[i + 11].x <= primeInverses[i + 11].y) ? 16    : 0) |
				((N * primeInverses[i + 12].x <= primeInverses[i + 12].y) ? 8     : 0) |
				((N * primeInverses[i + 13].x <= primeInverses[i + 13].y) ? 4     : 0) |
				((N * primeInverses[i + 14].x <= primeInverses[i + 14].y) ? 2     : 0) |
				((N * primeInverses[i + 15].x <= primeInverses[i + 15].y) ? 1     : 0);

			i += 16;

			// Do trial divisions by prime factors found
			if (pattern)
			{
				ulong sumN = 1;
				do
				{
					// Bit 0: divisibility for i - 1
					// Bit 1: divisibility for i - 2
					// Bit 2: divisibility for i - 3
					// and so on ...
					const uint num_leading_zeroes = clz(pattern);
					pattern -= (1U << (31 - num_leading_zeroes));
					const uint index = i + num_leading_zeroes - 32;

					const ulong prevSumN = sumN;
					const uint p1 = smallPrimes[index];
					ulong q = N * primeInverses[index].x;
					do
					{
						sumN = sumN * p1 + prevSumN;
						N = q;
						q *= primeInverses[index].x;
					} while (q <= primeInverses[index].y);
				} while (pattern);

				const ulong p = smallPrimes[i];
				if (p * p > N)
				{
					// N is 1 or prime (unless a factor below N / D was skipped: FirstUsefulPrimeIndex)
					if (N > 1)
					{
						sumN *= N + 1;
					}
					if (sumN == targetSum)
					{
						if (N > 1)
						{
							NumberFoundIfPrime(amicable_numbers_count, amicable_numbers_data, M, M_high, N);
						}
						else
						{
							NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, 0);
						}
					}
					return;
				}

				const ulong oldTargetSum = targetSum;
				targetSum /= sumN;
				if (targetSum * sumN != oldTargetSum)
				{
					return;
				}

				if (N >= targetSum)
				{
					return;
				}

				root4_N_rounded_up = convert_uint_rtp(sqrt(sqrt(convert_float_rtp(N))));
			}

			if (smallPrimes[i] > root4_N_rounded_up)
			{
				break;
			}

			uint index = k * PQ_STRIDE_SIZE + (i >> 4) + 16;
			while (N <= PQ[index].x)
			{
				--k;
				index -= PQ_STRIDE_SIZE;
			}
			if (mul_hi(N, PQ[index].y) < targetSum - N)
			{
				return;
			}
		}
	}

	// No more than 3 factors remain here
	
	// The case when N is prime (it has 1 factor)
	if (N + 1 == targetSum)
	{
		NumberFoundIfPrime(amicable_numbers_count, amicable_numbers_data, M, M_high, N);
		return;
	}

	const ulong root3_N_rounded_down = convert_ulong_rtn(cbrt(convert_float_rtn(N)));

	// The case when N can only have 2 factors
	if (N >= targetSum - root3_N_rounded_down * (root3_N_rounded_down + 1))
	{
		// The case when N is squared prime
		const ulong sqrt_N = IntegerSquareRoot(N);
		if ((sqrt_N * sqrt_N == N) && (N + sqrt_N + 1 == targetSum))
		{
			NumberFoundIfPrime(amicable_numbers_count, amicable_numbers_data, M, M_high, sqrt_N);
			return;
		}

		// The case when N = p * q - a product of two different primes
		// p + q == targetSum - N - 1
		// p * q == N
		const ulong B = targetSum - N - 1;
		ulong2 B2;
		B2.x = B * B;
		B2.y = mul_hi(B, B);
		const ulong carry = (B2.x < (N << 2)) ? 1 : 0;
		B2.x -= (N << 2);
		if (B2.y >= (N >> 62) + carry)
		{
			B2.y -= (N >> 62) + carry;
			const ulong sqrt_D = IntegerSquareRoot128(B2);
			const ulong p = (B - sqrt_D) >> 1;
			const ulong q = (B + sqrt_D) >> 1;
			if ((p < q) && (p + q == B) && (p * q == N) && SieveMayBePrime64(p))
			{
				NumberFoundIfPrime(amicable_numbers_count, amicable_numbers_data, M, M_high, q);
				return;
			}
		}
		return;
	}

	// N can have 3 factors
	// Check maximal possible sum of divisors for numbers with 3 factors
	uint p1 = smallPrimes[i];
	ulong q = N / p1;

	ulong p2 = p1;
	p2 *= p1 + 1;

	// If it's too small, N can't be amicable
	if (N < targetSum - ((q + p2) << 1))
	{
		return;
	}

	// Data needed to resume from arbitrary loop iteration in phase 3:
	// 1) M (8 bytes + 14 bits)
	// 2) N (8 bytes + 14 bits)
	// 3) targetSum (8 bytes + 14 bits)
	// 4) i (<= 2798576, can fit in 22 bits)
	// It can all fit in 32 bytes

	const uint filtered_numbers_data_id = atom_inc(filtered_numbers_count);
	filtered_numbers_data[filtered_numbers_data_id].x = M;
	filtered_numbers_data[filtered_numbers_data_id].y = N;
	filtered_numbers_data[filtered_numbers_data_id].z = targetSum;
	const ulong i_to_save = i;
	filtered_numbers_data[filtered_numbers_data_id].w = i_to_save | (M_high << 22) | (N_high << 36) | (targetSumHigh << 50);
}
)==="
R"===(
// Phase 3 for the late iterations (RunPhase3: few numbers, each with up to ~1e5 primes; with one work-item per number most of the GPU was
// idle and the loads were serial). One warp per number: round k covers the primes [i + 1024 k, i + 1024 (k + 1)) (i rounded down to a
// multiple of 32), lane l tests the primes i + 1024 k + l + 32 t (t = 0..31: the loads of a warp are consecutive). The checks of
// CheckPairPhase3 at the first prime p of each round are the same for all lanes: the bound (D >= 2 (N / p + p (p + 1)) for at most 3
// factors >= p) and p > N^(1/3) (from a float estimate rounded up: a too large estimate only tests more primes). wideStatus[n] (all ones
// before the launch) gets atomic_min(k << 2 | type) of the first round with an event: 0 a prime of the round divides N, 1 the bound
// rejects N at p, 2 p > N^(1/3) (or N_high >= p, or the end of the table). CheckPairPhase3 (wideMode) then continues every number at the
// first prime of that round (the factors, the closed forms and the decisions stay there), the numbers without an event at max_i.
#define PHASE3_WIDE_ROUND_SHIFT 10

__kernel
__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
void CheckPairPhase3Wide(
	__global const ulong4* numbersToCheck,
	const uint count,
	const int max_i,
	const int tableLimit,
	__global const uint* smallPrimes,
	__global const ulong2* primeInverses,
	__global const ulong2* PrimeInverses_128,
	__global uint* wideStatus)
{
	const uint n = get_global_id(0) >> 5;
	const uint lane = get_global_id(0) & 31;
	if (n >= count)
	{
		return;
	}
	const ulong4 cur_number = numbersToCheck[n];
	if (cur_number.x == 0)
	{
		return;
	}
	const ulong N = cur_number.y;
	const ulong targetSum = cur_number.z;
	const ulong N_high = (cur_number.w >> 36) & ((1U << 14) - 1);
	const ulong targetSumHigh = cur_number.w >> 50;
	const int base = (int)(cur_number.w & ((1 << 22) - 1)) & ~31;
	const ulong D = targetSum - N - 1;
	const ulong D_high = targetSumHigh - N_high - (((targetSum < N + 1) || (N + 1 == 0)) ? 1 : 0);
	const float Nf = convert_float(N_high) * 18446744073709551616.0f + convert_float(N);
	const float root3 = cbrt(Nf) * 1.0001f + 1.0f;
	__global const uint* inv32 = (__global const uint*)PrimeInverses_128;
	for (uint k = 0;; ++k)
	{
		const int s0 = base + (int)(k << PHASE3_WIDE_ROUND_SHIFT);
		if (s0 >= max_i)
		{
			return;
		}
		const uint p = smallPrimes[s0];
		if ((s0 + (1 << PHASE3_WIDE_ROUND_SHIFT) > tableLimit) || (N_high >= p) || (convert_float(p) > root3))
		{
			if (lane == 0)
			{
				atomic_min(&wideStatus[n], (k << 2) | 2);
			}
			return;
		}
		// (the bound of CheckPairPhase3; N < 2^78 and p > 2^16, so N / p and the bound fit in 64 bits)
		const ulong pmax = primeInverses[s0].y; // floor((2^64 - 1) / p)
		const ulong q_upper = N_high * (pmax + 1) + mul_hi(N, pmax) + 2;
		const ulong bound = (q_upper + (ulong)p * (p + 1)) << 1;
		if (D_high || (bound <= D))
		{
			if (lane == 0)
			{
				atomic_min(&wideStatus[n], (k << 2) | 1);
			}
			return;
		}
		uint found = 0;
		#pragma unroll
		for (int t = 0; t < 32; ++t)
		{
			const int idx = s0 + (int)lane + 32 * t;
			found |= Divides128Fast32(N, (uint)N_high, inv32[(idx - 1) * 4], smallPrimes[idx]) ? 1U : 0U;
		}
		if (found)
		{
			atomic_min(&wideStatus[n], k << 2);
			return;
		}
	}
}

)==="
R"===(
__kernel
__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
void CheckPairPhase3(
	__global const ulong4* numbersToCheck,
	const uint numbersToCheckOffset,
	__global const ulong2* primeInverses,
	__global const ulong2* primeReciprocals,
	__global uint* filtered_numbers_count,
	__global ulong4* filtered_numbers_data,
	__global uint* amicable_numbers_count,
	__global ulong4* amicable_numbers_data,
	const int max_i,
	__global const uint* smallPrimes,
	__global const uint* PowersOfP_128SumInverses_offsets,
	__global const ulong4* PowersOfP_128SumInverses,
	__global const ulong2* PrimeInverses_128,
	__global const uint* wideStatus,
	const uint wideMode
)
{
	const ulong4 cur_number = numbersToCheck[get_global_id(0) + numbersToCheckOffset];
	if (cur_number.x == 0)
	{
		return;
	}

	const ulong M = cur_number.x;
	ulong N = cur_number.y;
	ulong targetSum = cur_number.z;

	int i = cur_number.w & ((1 << 22) - 1);
	if (wideMode)
	{
		// after CheckPairPhase3Wide: no prime of [i, max_i) divides N (and no bound check returns there) if wideStatus is all ones, then
		// the number goes to the next iteration from max_i; otherwise it continues at the first prime of the round of its first event
		const uint code = wideStatus[get_global_id(0) + numbersToCheckOffset];
		i = (code == 0xFFFFFFFFU) ? max(i, max_i) : ((i & ~31) + (int)((code >> 2) << PHASE3_WIDE_ROUND_SHIFT));
	}

	const ulong M_high = (cur_number.w >> 22) & ((1U << 14) - 1);
	ulong N_high = (cur_number.w >> 36) & ((1U << 14) - 1);
	ulong targetSumHigh = cur_number.w >> 50;

	if (targetSumHigh)
	{
		ulong2 N_128;
		N_128.x = N;
		N_128.y = N_high;
		uint root3_N_rounded_up = IntegerCubeRoot128(N_128);

		// (phase 2 decided the numbers with at most 2 prime factors, DecideAtMost2Factors, it is checked again when N changes)
		for (;;)
		{
			if ((i & 31) == 0)
			{
				if (!((i < max_i) && (smallPrimes[i] <= root3_N_rounded_up) && targetSumHigh))
				{
					break;
				}
				// All prime factors of N are >= p > N^(1/4), so N has at most 3 prime factors and sigma(N) - N - 1 <= 2 * (N / p + p * (p + 1))
				// (the same bound as in the 64-bit loop below). If it can't reach targetSum - N - 1, N is not the other member of a pair.
				// N < 2^78 and p > 2^16, so N / p and the bound fit in 64 bits.
				const ulong p = smallPrimes[i];
				const ulong pmax = primeInverses[i].y; // floor((2^64 - 1) / p)
				const ulong q_upper = N_high * (pmax + 1) + mul_hi(N, pmax) + 2;
				const ulong bound = (q_upper + p * (p + 1)) << 1;
				const ulong D = targetSum - N - 1;
				const ulong D_high = targetSumHigh - N_high - ((targetSum < N + 1) || (N + 1 == 0) ? 1 : 0);
				if (D_high || (bound <= D))
				{
					return;
				}
			}

			// The primes i .. (the next multiple of 32) - 1 (32 of them except at the start): the divisibility without branches, then the
			// divisions for the primes that divide N (rare). With the division inside the loop over the primes, the loop had ~35
			// instructions per prime (reconvergence and control flow around the rare path).
			const int base = i;
			const int count = ((i | 31) + 1) - i;
			uint pattern = 0;
			if (N_high < smallPrimes[i])
			{
				// (all primes of the block are > N_high: Divides128Fast32 is exact; PrimeInverses_128[k - 1] is the inverse of smallPrimes[k], its
				// low 32 bits the inverse mod 2^32)
				__global const uint* inv32 = (__global const uint*)PrimeInverses_128;
				#pragma unroll
				for (int j = 0; j < 32; ++j)
				{
					if (j < count)
					{
						pattern |= Divides128Fast32(N, (uint)N_high, inv32[(base + j - 1) * 4], smallPrimes[base + j]) ? (1U << j) : 0U;
					}
				}
			}
			else
			{
				for (int j = 0; j < count; ++j)
				{
					const ulong2 invj = PrimeInverses_128[base + j - 1];
					const ulong qh = mul_hi(N, invj.x) + N * invj.y + N_high * invj.x;
					pattern |= (qh <= N_high) ? (1U << j) : 0U;
				}
			}
			i = base + count;
			while (pattern != 0)
			{
				// the smallest prime first
				const int ii = base + (31 - (int)clz(pattern & (0U - pattern)));
				pattern &= pattern - 1;
				const ulong2 inv = PrimeInverses_128[ii - 1];
				ulong q_high = mul_hi(N, inv.x) + N * inv.y + N_high * inv.x;
				if (q_high > N_high)
				{
					continue;
				}

				ulong q = N * inv.x;
				int powerOfP = 0;
				do
				{
					++powerOfP;

					N = q;
					N_high = q_high;

					q = N * inv.x;
					q_high = mul_hi(N, inv.x) + N * inv.y + N_high * inv.x;
				} while (q_high <= N_high);

				const ulong4 inv2 = PowersOfP_128SumInverses[PowersOfP_128SumInverses_offsets[ii - 1] + powerOfP - 1];
				if (inv2.z != 0)
				{
					if ((targetSum & inv2.w) != 0)
					{
						return;
					}
					targetSum = (targetSum >> inv2.z) | (targetSumHigh << (64 - inv2.z));
					targetSumHigh >>= inv2.z;
				}

				q_high = mul_hi(targetSum, inv2.x) + targetSum * inv2.y + targetSumHigh * inv2.x;
				if (q_high > targetSumHigh)
				{
					return;
				}

				targetSum *= inv2.x;
				targetSumHigh = q_high;

				const ulong p = smallPrimes[ii];
				if ((N_high == 0) && (p * p > N))
				{
					// N is 1 or prime (unless a factor below N / D was skipped: FirstUsefulPrimeIndex). (N_high: N can still be >= 2^64 here,
					// the test of the low 64 bits alone made rare numbers return)
					const ulong sumN = ((N > 1) ? N : 0) + 1;
					if (sumN == targetSum)
					{
						if (N > 1)
						{
							NumberFoundIfPrime(amicable_numbers_count, amicable_numbers_data, M, M_high, N);
						}
						else
						{
							NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, 0);
						}
					}
					return;
				}

				if ((N_high > targetSumHigh) || ((N_high == targetSumHigh) && (N >= targetSum)))
				{
					return;
				}

					N_128.x = N;
					N_128.y = N_high;
					root3_N_rounded_up = IntegerCubeRoot128(N_128);

					ulong B = 0;
					const uint decided = DecideAtMost2Factors(N, N_high, targetSum, targetSumHigh, &B, false);
					if (decided != 0)
					{
						if (decided == 2)
						{
							NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, B);
						}
						return;
					}
			}
		}

		if (targetSumHigh)
		{
			if (smallPrimes[i] <= root3_N_rounded_up)
				{
					const uint filtered_numbers_data_id = atom_inc(filtered_numbers_count);
					filtered_numbers_data[filtered_numbers_data_id].x = M;
					filtered_numbers_data[filtered_numbers_data_id].y = N;
					filtered_numbers_data[filtered_numbers_data_id].z = targetSum;
					// (i >= max_i: it can be larger in the first iterations, the numbers continue there)
					const ulong i_to_save = max(i, max_i);
				filtered_numbers_data[filtered_numbers_data_id].w = i_to_save | (M_high << 22) | (N_high << 36) | (targetSumHigh << 50);
			}
			else
			{
				const ulong N2 = N + 1;
				const ulong N2_high = N_high + ((N2 == 0) ? 1 : 0);
				if (targetSum & 1)
				{
					const ulong carry = (targetSum < N2) ? 1 : 0;
					const ulong p = targetSum - N2;
					const ulong p_high = targetSumHigh - N2_high - carry;
					if ((p_high == 0) && (p * p == N) && (mul_hi(p, p) == N_high))
					{
						NumberFoundIfPrime(amicable_numbers_count, amicable_numbers_data, M, M_high, p);
					}
				}
				else
				{
					if ((targetSum == N2) && (targetSumHigh == N2_high))
					{
						NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, 0);
					}
					else
					{
						const ulong carry = (targetSum < N2) ? 1 : 0;
						ulong2 B;
						B.x = targetSum - N2;
						B.y = targetSumHigh - N2_high - carry;
						B.x = (B.x >> 1) | (B.y << 63);
						B.y >>= 1;
						if (B.y == 0)
						{
							ulong2 D;
							D.x = B.x * B.x;
							D.y = mul_hi(B.x, B.x);
							if ((D.y > N_high) || ((D.y == N_high) && (D.x > N)))
							{
								const ulong carry2 = (D.x < N) ? 1 : 0;
								D.x -= N;
								D.y -= N_high + carry2;
								ulong sqrtD = IntegerSquareRoot128(D);
								if ((sqrtD * sqrtD == D.x) && (mul_hi(sqrtD, sqrtD) == D.y))
								{
									const ulong p = B.x - sqrtD;
									const ulong q = B.x + sqrtD;
									if ((p * q == N) && (mul_hi(p, q) == N_high) && SieveMayBePrime64(q))
									{
										NumberFoundIfPrime(amicable_numbers_count, amicable_numbers_data, M, M_high, p);
									}
								}
							}
						}
					}
				}
			}
			return;
		}
	}

	// No more than 3 factors remain here
	const ulong root3_N_rounded_up = convert_uint_rtp(cbrt(convert_float_rtp(N)));
	const ulong MaximumSumOfDivisors3_target = targetSum - N - 1;

	// The numbers that come from the 128-bit part above can have at most 2 prime factors (DecideAtMost2Factors; phase 2 checked the others)
	{
		ulong B = 0;
		const uint decided = DecideAtMost2Factors(N, N_high, targetSum, targetSumHigh, &B, false);
		if (decided != 0)
		{
			if (decided == 2)
			{
				NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, B);
			}
			return;
		}
	}

	// All threads start with the same value of "i" here except for "Phase 3, iteration 1"
	// so there will be no divergence because of different values of "i" on iteration 2, 3, etc.
	// "Phase 3, iteration 1" is deliberately made shorter because of this
	for (;;)
	{
		uint pattern =
			((N * primeInverses[i +  0].x <= primeInverses[i +  0].y) ? (1U << 31) : 0) |
			((N * primeInverses[i +  1].x <= primeInverses[i +  1].y) ? (1U << 30) : 0) |
			((N * primeInverses[i +  2].x <= primeInverses[i +  2].y) ? (1U << 29) : 0) |
			((N * primeInverses[i +  3].x <= primeInverses[i +  3].y) ? (1U << 28) : 0) |
			((N * primeInverses[i +  4].x <= primeInverses[i +  4].y) ? (1U << 27) : 0) |
			((N * primeInverses[i +  5].x <= primeInverses[i +  5].y) ? (1U << 26) : 0) |
			((N * primeInverses[i +  6].x <= primeInverses[i +  6].y) ? (1U << 25) : 0) |
			((N * primeInverses[i +  7].x <= primeInverses[i +  7].y) ? (1U << 24) : 0) |
			((N * primeInverses[i +  8].x <= primeInverses[i +  8].y) ? (1U << 23) : 0) |
			((N * primeInverses[i +  9].x <= primeInverses[i +  9].y) ? (1U << 22) : 0) |
			((N * primeInverses[i + 10].x <= primeInverses[i + 10].y) ? (1U << 21) : 0) |
			((N * primeInverses[i + 11].x <= primeInverses[i + 11].y) ? (1U << 20) : 0) |
			((N * primeInverses[i + 12].x <= primeInverses[i + 12].y) ? (1U << 19) : 0) |
			((N * primeInverses[i + 13].x <= primeInverses[i + 13].y) ? (1U << 18) : 0) |
			((N * primeInverses[i + 14].x <= primeInverses[i + 14].y) ? (1U << 17) : 0) |
			((N * primeInverses[i + 15].x <= primeInverses[i + 15].y) ? (1U << 16) : 0) |
			((N * primeInverses[i + 16].x <= primeInverses[i + 16].y) ? (1U << 15) : 0) |
			((N * primeInverses[i + 17].x <= primeInverses[i + 17].y) ? (1U << 14) : 0) |
			((N * primeInverses[i + 18].x <= primeInverses[i + 18].y) ? (1U << 13) : 0) |
			((N * primeInverses[i + 19].x <= primeInverses[i + 19].y) ? (1U << 12) : 0) |
			((N * primeInverses[i + 20].x <= primeInverses[i + 20].y) ? (1U << 11) : 0) |
			((N * primeInverses[i + 21].x <= primeInverses[i + 21].y) ? (1U << 10) : 0) |
			((N * primeInverses[i + 22].x <= primeInverses[i + 22].y) ? (1U <<  9) : 0) |
			((N * primeInverses[i + 23].x <= primeInverses[i + 23].y) ? (1U <<  8) : 0) |
			((N * primeInverses[i + 24].x <= primeInverses[i + 24].y) ? (1U <<  7) : 0) |
			((N * primeInverses[i + 25].x <= primeInverses[i + 25].y) ? (1U <<  6) : 0) |
			((N * primeInverses[i + 26].x <= primeInverses[i + 26].y) ? (1U <<  5) : 0) |
			((N * primeInverses[i + 27].x <= primeInverses[i + 27].y) ? (1U <<  4) : 0) |
			((N * primeInverses[i + 28].x <= primeInverses[i + 28].y) ? (1U <<  3) : 0) |
			((N * primeInverses[i + 29].x <= primeInverses[i + 29].y) ? (1U <<  2) : 0) |
			((N * primeInverses[i + 30].x <= primeInverses[i + 30].y) ? (1U <<  1) : 0) |
			((N * primeInverses[i + 31].x <= primeInverses[i + 31].y) ? (1U <<  0) : 0);

		i += 32;

		// Do trial divisions by prime factors found
		if (pattern)
		{
			ulong sumN = 1;
			do
			{
				// Bit 0: divisibility for i - 1
				// Bit 1: divisibility for i - 2
				// Bit 2: divisibility for i - 3
				// and so on ...
				const uint num_leading_zeroes = clz(pattern);
				pattern -= (1U << (31 - num_leading_zeroes));
				const uint index = i + num_leading_zeroes - 32;

				const ulong prevSumN = sumN;
				const uint p1 = primeReciprocals[index].y >> 32;
				ulong q = N * primeInverses[index].x;
				do
				{
					sumN = sumN * p1 + prevSumN;
					N = q;
					q *= primeInverses[index].x;
				} while (q <= primeInverses[index].y);
			} while (pattern);

			const ulong p2 = primeReciprocals[i].y >> 32;
			if (p2 * p2 > N)
			{
				// N is 1 or prime (unless a factor below N / D was skipped: FirstUsefulPrimeIndex)
				if (N > 1)
				{
					sumN *= N + 1;
				}
				if (sumN == targetSum)
				{
					if (N > 1)
					{
						NumberFoundIfPrime(amicable_numbers_count, amicable_numbers_data, M, M_high, N);
					}
					else
					{
						NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, 0);
					}
				}
				return;
			}

			const ulong oldTargetSum = targetSum;
			targetSum /= sumN;
			if (targetSum * sumN != oldTargetSum)
			{
				return;
			}

			if (N >= targetSum)
			{
				return;
			}

			// No more than 2 factors remain, we need to exit the loop now, otherwise MaximumSumOfDivisors3 will fail
			break;
		}

		const ulong2 curReciprocal = primeReciprocals[i];
		const ulong p1 = curReciprocal.y >> 32;
		if (p1 > root3_N_rounded_up)
		{
			break;
		}
		const ulong q = mul_hi(N + ((curReciprocal.y >> 8) & 1), curReciprocal.x) >> (curReciprocal.y & 63);

		ulong p2 = p1;
		p2 *= p1 + 1;
		if (((q + p2) << 1) <= MaximumSumOfDivisors3_target)
		{
			return;
		}

		if (i >= max_i)
		{
			const uint filtered_numbers_data_id = atom_inc(filtered_numbers_count);
			filtered_numbers_data[filtered_numbers_data_id].x = M;
			filtered_numbers_data[filtered_numbers_data_id].y = N;
			filtered_numbers_data[filtered_numbers_data_id].z = targetSum;
			// (i can be larger than max_i in the first iterations, the numbers continue there: multiples of 32 keep the same i in a warp)
			const ulong i_to_save = max(i & ~31, max_i);
			filtered_numbers_data[filtered_numbers_data_id].w = i_to_save | (M_high << 22) | (N_high << 36) | (targetSumHigh << 50);
			return;
		}
	}

	// No more than 2 factors remain here

	// 1) The case when N is prime
	const ulong minSumN = N + 1;
	if (minSumN >= targetSum)
	{
		if (minSumN == targetSum)
		{
			// M is an amicable number if N is prime
			// We don't check for primality here, let CPU do it
			NumberFoundIfPrime(amicable_numbers_count, amicable_numbers_data, M, M_high, N);
		}
		return;
	}

	// 2) The case when N is squared prime
	const ulong sqrt_N = IntegerSquareRoot(N);
	if ((sqrt_N * sqrt_N == N) && (N + sqrt_N + 1 == targetSum))
	{
		NumberFoundIfPrime(amicable_numbers_count, amicable_numbers_data, M, M_high, sqrt_N);
		return;
	}

	// 3) The case when N = p * q - a product of two different primes
	// p + q == targetSum - N - 1
	// p * q == N
	const ulong B = targetSum - N - 1;
	ulong2 B2;
	B2.x = B * B;
	B2.y = mul_hi(B, B);
	const ulong carry = (B2.x < (N << 2)) ? 1 : 0;
	B2.x -= (N << 2);
	if (B2.y >= (N >> 62) + carry)
	{
		B2.y -= (N >> 62) + carry;
		const ulong sqrt_D = IntegerSquareRoot128(B2);
		const ulong p = (B - sqrt_D) >> 1;
		const ulong q = (B + sqrt_D) >> 1;
		if ((p < q) && (p + q == B) && (p * q == N) && SieveMayBePrime64(p))
		{
			NumberFoundIfPrime(amicable_numbers_count, amicable_numbers_data, M, M_high, q);
			return;
		}
	}
}
)==="
R"===(
__kernel
__attribute__((reqd_work_group_size(1, 1, 1)))
void SaveCounter(__global uint* phase2_numbers_count)
{
	if (phase2_numbers_count[0] <= PHASE2_MAX_COUNT)
	{
		phase2_numbers_count[1] = phase2_numbers_count[0];
	}
}

// The numbers of phase 2 sorted by the primes they need (the host: SortPhase2Numbers). The primes of phases 2 and 3 start below N / D
// (FirstUsefulPrimeIndex), and the work of a number grows with N / D: in the ranges with W ~ 0, the N / D of the numbers in a warp
// differed by 4x - 30x (the warps of phase 2 took ~2.5x the average). Buckets: 64-bit and 128-bit numbers (different loops) x 32 classes
// of N / D (class 0: below 5393, no primes are skipped, then 3 classes per octave).
#define PHASE2_BUCKETS 64

static uint Phase2Bucket(const ulong4 x)
{
	const ulong N = x.y;
	const ulong N_high = (x.w >> 36) & ((1U << 14) - 1);
	const ulong targetSum = x.z;
	const ulong targetSumHigh = x.w >> 50;
	const ulong D = targetSum - N - 1;
	const ulong D_high = targetSumHigh - N_high - (((targetSum < N + 1) || (N + 1 == 0)) ? 1 : 0);
	uint c = 0;
	if ((D | D_high) != 0)
	{
		const float L = (convert_float(N_high) * 18446744073709551616.0f + convert_float(N)) / (convert_float(D_high) * 18446744073709551616.0f + convert_float(D));
		if (L >= 5393.0f)
		{
			c = 1 + (uint)min(30.0f, 3.0f * (log2(L) - 12.397f));
		}
	}
	return (N_high ? 32 : 0) + c;
}

__kernel
__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
void Phase2BucketCount(__global const ulong4* numbers, const uint count, __global uint* bucketCounts)
{
	__local uint localCounts[PHASE2_BUCKETS];
	for (uint b = get_local_id(0); b < PHASE2_BUCKETS; b += WORK_GROUP_SIZE)
	{
		localCounts[b] = 0;
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	const uint id = get_global_id(0);
	if (id < count)
	{
		atomic_inc(&localCounts[Phase2Bucket(numbers[id])]);
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	for (uint b = get_local_id(0); b < PHASE2_BUCKETS; b += WORK_GROUP_SIZE)
	{
		if (localCounts[b] != 0)
		{
			atomic_add(&bucketCounts[b], localCounts[b]);
		}
	}
}

// bucketOffsets: the first index of each bucket in "sorted" (the kernel advances them)
__kernel
__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
void Phase2BucketScatter(__global const ulong4* numbers, const uint count, __global uint* bucketOffsets, __global ulong4* sorted)
{
	__local uint localCounts[PHASE2_BUCKETS];
	__local uint localBase[PHASE2_BUCKETS];
	for (uint b = get_local_id(0); b < PHASE2_BUCKETS; b += WORK_GROUP_SIZE)
	{
		localCounts[b] = 0;
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	const uint id = get_global_id(0);
	ulong4 x = (ulong4)(0, 0, 0, 0);
	uint bucket = 0;
	uint slot = 0;
	if (id < count)
	{
		x = numbers[id];
		bucket = Phase2Bucket(x);
		slot = atomic_inc(&localCounts[bucket]);
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	for (uint b = get_local_id(0); b < PHASE2_BUCKETS; b += WORK_GROUP_SIZE)
	{
		if (localCounts[b] != 0)
		{
			localBase[b] = atomic_add(&bucketOffsets[b], localCounts[b]);
		}
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	if (id < count)
	{
		sorted[localBase[bucket] + slot] = x;
	}
}

// Before phase 2 (the host: ProcessNumbersPhases2_3): the numbers with D < 3 N^(2/3) that phase 1 could not decide (the form r^2 s,
// DecideAtMost2Factors) are decided here and removed (M = 0: the kernels of phase 2 skip them). In the ranges with W ~ 0 they are most
// numbers of phase 2 for some N (t1: D / N^(2/3) ~ 1.88 - 3), and they needed the trial division up to 2 N / D.
__kernel
__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
void Phase2Prefilter(__global ulong4* numbers, const uint count, __global uint* amicable_numbers_count, __global ulong4* amicable_numbers_data)
{
	const uint id = get_global_id(0);
	if (id >= count)
	{
		return;
	}
	const ulong4 x = numbers[id];
	if (x.x == 0)
	{
		return;
	}
	const ulong M_high = (x.w >> 22) & ((1U << 14) - 1);
	const ulong N_high = (x.w >> 36) & ((1U << 14) - 1);
	const ulong targetSumHigh = x.w >> 50;
	ulong B = 0;
	const uint decided = DecideAtMost2Factors(x.y, N_high, x.z, targetSumHigh, &B, true);
	if (decided != 0)
	{
		if (decided == 2)
		{
			NumberFound(amicable_numbers_count, amicable_numbers_data, x.x, M_high, B);
		}
		numbers[id].x = 0;
	}
}

)==="
R"===(
// primes is an array of pointers to chunks of data
// each chunk is 128 MB in size (2^24 uint2 elements) or bigger
static ulong GetNthPrime(uint n,
	__global uint2* primes0
#if NUM_DATA_CHUNKS > 1
	, __global uint2* primes1
#endif
#if NUM_DATA_CHUNKS > 2
	, __global uint2* primes2
#endif
#if NUM_DATA_CHUNKS > 3
	, __global uint2* primes3
#endif
#if NUM_DATA_CHUNKS > 4
	, __global uint2* primes4
#endif
)
{
	if (n < 4)
		return (0x07050302 >> (n << 3)) & 0xFF;

	n &= MAIN_BUFFER_MASK;

#if NUM_DATA_CHUNKS == 1
	const uint2 data = primes0[n >> 2];
#else
	const uint global_offset = n >> 2;
	__global const uint2* chunk;
	switch (global_offset >> (CHUNK_SIZE_SHIFT - 3))
	{
	default: chunk = primes0; break;
#if NUM_DATA_CHUNKS > 1
	case 1: chunk = primes1; break;
#endif
#if NUM_DATA_CHUNKS > 2
	case 2: chunk = primes2; break;
#endif
#if NUM_DATA_CHUNKS > 3
	case 3: chunk = primes3; break;
#endif
#if NUM_DATA_CHUNKS > 4
	case 4: chunk = primes4; break;
#endif
	}

	const uint chunk_offset = global_offset & ((1U << (CHUNK_SIZE_SHIFT - 3)) - 1);
	const uint2 data = chunk[chunk_offset];
#endif

	ulong base = data.y & 31;
	base = (base << 32) + data.x;

	const ulong offsets = data.y >> 5;

	return base + (((offsets >> ((~n & 3) * 9)) & 511) << 1);
}

__kernel
__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
void SearchMultipleRanges(
	__global const uint* smallPrimes,
	__global const ulong2* primeInverses,
	__global const ulong2* PQ,
	__constant ulong2* PowerOf2SumInverses,
	__global const uint* PowersOfP_128SumInverses_offsets,
	__global const ulong4* PowersOfP_128SumInverses,
	__global const ulong2* PrimeInverses_128,
	__global ulong* SumEstimates_128,
	__global const ulong* RangesTable,
	__global const uint2* RangeLookupTable,
	const uint lookup_shift,
	const ulong global_offset,
	const uint global_size,
	__global ulong* phase1_offset_to_resume_after_overflow,
	__global uint* phase2_numbers_count,
	__global ulong4* phase2_numbers_data,
	__global uint* amicable_numbers_count,
	__global ulong4* amicable_numbers_data,

	__global uint2* primes0
#if NUM_DATA_CHUNKS > 1
	, __global uint2* primes1
#endif
#if NUM_DATA_CHUNKS > 2
	, __global uint2* primes2
#endif
#if NUM_DATA_CHUNKS > 3
	, __global uint2* primes3
#endif
#if NUM_DATA_CHUNKS > 4
	, __global uint2* primes4
#endif
)
{
	if (*phase2_numbers_count > PHASE2_MAX_COUNT)
	{
		return;
	}

	const uint globalIndex = get_global_id(0) + global_offset;
	if (globalIndex >= global_size)
	{
		return;
	}

	ulong2 M0, sumM0;
	ulong p;

	uint2 rangeIndexAndOffset = RangeLookupTable[globalIndex >> lookup_shift];
	uint curIndex = globalIndex & (0xFFFFFFFFU << lookup_shift);
	for (;;)
	{
		__global const ulong* curRange = RangesTable + rangeIndexAndOffset.x * 5;
		const uint numbersRemainingInCurRange = (curRange[4] >> 32) - rangeIndexAndOffset.y;

		// curIndex + numbersRemainingInCurRange now points exactly at the beginning of the next range
		// If the beginning of the next range is > globalIndex, then globalIndex is within the current range
		if (curIndex + numbersRemainingInCurRange > globalIndex)
		{
			M0.x = curRange[0];
			M0.y = curRange[1];
			sumM0.x = curRange[2];
			sumM0.y = curRange[3];

			//      primes[ start_prime_index            + range_offset          + globalIndex - curIndex]
			p = GetNthPrime((curRange[4] & 0xFFFFFFFFUL) + rangeIndexAndOffset.y + globalIndex - curIndex,
				primes0
#if NUM_DATA_CHUNKS > 1
				, primes1
#endif
#if NUM_DATA_CHUNKS > 2
				, primes2
#endif
#if NUM_DATA_CHUNKS > 3
				, primes3
#endif
#if NUM_DATA_CHUNKS > 4
				, primes4
#endif
			);
			break;
		}

		// Move curIndex to the beginning of the next range
		curIndex += numbersRemainingInCurRange;
		++rangeIndexAndOffset.x;
		rangeIndexAndOffset.y = 0;
	}

#if LPP == 1
	ulong value = p;
	ulong value_sum = p + 1;
#elif LPP == 2
	const ulong value = p * p;
	const ulong value_sum = p * (p + 1) + 1;
#elif LPP == 3
	const ulong value = p * p * p;
	const ulong value_sum = p * (p * (p + 1) + 1) + 1;
#endif

	CheckPairPhase1(
		smallPrimes,
		primeInverses,
		PQ,
		PowerOf2SumInverses,
		PowersOfP_128SumInverses_offsets,
		PowersOfP_128SumInverses,
		PrimeInverses_128,
		SumEstimates_128,
		M0.x * value, mul_hi(M0.x, value) + M0.y * value, sumM0.x * value_sum, mul_hi(sumM0.x, value_sum) + sumM0.y * value_sum,
		global_offset,
		phase1_offset_to_resume_after_overflow,
		phase2_numbers_count,
		phase2_numbers_data,
		amicable_numbers_count,
		amicable_numbers_data,
		0, phase2_numbers_count, (__global ulong*)phase2_numbers_data, 0
	);
}

__kernel
__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
void SearchLargePrimes(
	__global const uint* smallPrimes,
	__global const ulong2* primeInverses,
	__global const ulong2* PQ,
	__constant ulong2* PowerOf2SumInverses,
	__global const uint* PowersOfP_128SumInverses_offsets,
	__global const ulong4* PowersOfP_128SumInverses,
	__global const ulong2* PrimeInverses_128,
	__global const ulong* SumEstimates_128,
	__global const ulong* largePrimes,
	const uint largePrimesCount,
	const ulong largePrimesCountReciprocal,
	const uint largePrimesCountIncrementAndShift,
	const ulong global_offset,
	const ulong global_size,
	__global ulong* phase1_offset_to_resume_after_overflow,
	__global uint* phase2_numbers_count,
	__global ulong4* phase2_numbers_data,
	__global uint* amicable_numbers_count,
	__global ulong4* amicable_numbers_data,

	__global ulong* amicableCandidates0
#if NUM_DATA_CHUNKS > 1
	, __global ulong* amicableCandidates1
#endif
#if NUM_DATA_CHUNKS > 2
	, __global ulong* amicableCandidates2
#endif
#if NUM_DATA_CHUNKS > 3
	, __global ulong* amicableCandidates3
#endif
#if NUM_DATA_CHUNKS > 4
	, __global ulong* amicableCandidates4
#endif
#if NUM_DATA_CHUNKS > 5
	, __global ulong* amicableCandidates5
#endif
#if NUM_DATA_CHUNKS > 6
	, __global ulong* amicableCandidates6
#endif
#if NUM_DATA_CHUNKS > 7
	, __global ulong* amicableCandidates7
#endif
#if NUM_DATA_CHUNKS > 8
	, __global ulong* amicableCandidates8
#endif
#if NUM_DATA_CHUNKS > 9
	, __global ulong* amicableCandidates9
#endif
#if NUM_DATA_CHUNKS > 10
	, __global ulong* amicableCandidates10
#endif
)
{
	if (*phase2_numbers_count > PHASE2_MAX_COUNT)
	{
		return;
	}

	const ulong globalIndex = global_offset + get_global_id(0);
	if (globalIndex >= global_size)
	{
		return;
	}

	ulong amicableCandidateIndex;
	uint largePrimeIndex;

	if (largePrimesCountReciprocal != 0)
	{
		amicableCandidateIndex = mul_hi(globalIndex + (largePrimesCountIncrementAndShift & 1), largePrimesCountReciprocal);
		amicableCandidateIndex >>= (largePrimesCountIncrementAndShift >> 1);
		largePrimeIndex = globalIndex - amicableCandidateIndex * largePrimesCount;
	}
	else
	{
		amicableCandidateIndex = globalIndex >> largePrimesCountIncrementAndShift;
		largePrimeIndex = globalIndex - (amicableCandidateIndex << largePrimesCountIncrementAndShift);
	}

	// AmicableCandidate is 10 bytes: 32-bit value, 32-bit sum, 8 high bits of value, 8 high bits of sum
	const ulong amicableCandidateByteIndex = amicableCandidateIndex * 10;

#if NUM_DATA_CHUNKS == 1
	ulong data0 = amicableCandidates0[amicableCandidateByteIndex / 8];
	ulong data1 = amicableCandidates0[amicableCandidateByteIndex / 8 + 1];
#else
	__global const ulong* chunk;
	__global const ulong* next_chunk;

#if   NUM_DATA_CHUNKS == 2
	__global const ulong* amicableCandidates2 = amicableCandidates0;
#elif NUM_DATA_CHUNKS == 3
	__global const ulong* amicableCandidates3 = amicableCandidates0;
#elif NUM_DATA_CHUNKS == 4
	__global const ulong* amicableCandidates4 = amicableCandidates0;
#elif NUM_DATA_CHUNKS == 5
	__global const ulong* amicableCandidates5 = amicableCandidates0;
#elif NUM_DATA_CHUNKS == 6
	__global const ulong* amicableCandidates6 = amicableCandidates0;
#elif NUM_DATA_CHUNKS == 7
	__global const ulong* amicableCandidates7 = amicableCandidates0;
#elif NUM_DATA_CHUNKS == 8
	__global const ulong* amicableCandidates8 = amicableCandidates0;
#elif NUM_DATA_CHUNKS == 9
	__global const ulong* amicableCandidates9 = amicableCandidates0;
#elif NUM_DATA_CHUNKS == 10
	__global const ulong* amicableCandidates10 = amicableCandidates0;
#elif NUM_DATA_CHUNKS == 11
	__global const ulong* amicableCandidates11 = amicableCandidates0;
#endif

	switch (amicableCandidateByteIndex >> CHUNK_SIZE_SHIFT)
	{
	default: chunk = amicableCandidates0; next_chunk = amicableCandidates1; break;
#if NUM_DATA_CHUNKS > 1
	case 1: chunk = amicableCandidates1; next_chunk = amicableCandidates2; break;
#endif
#if NUM_DATA_CHUNKS > 2
	case 2: chunk = amicableCandidates2; next_chunk = amicableCandidates3; break;
#endif
#if NUM_DATA_CHUNKS > 3
	case 3: chunk = amicableCandidates3; next_chunk = amicableCandidates4; break;
#endif
#if NUM_DATA_CHUNKS > 4
	case 4: chunk = amicableCandidates4; next_chunk = amicableCandidates5; break;
#endif
#if NUM_DATA_CHUNKS > 5
	case 5: chunk = amicableCandidates5; next_chunk = amicableCandidates6; break;
#endif
#if NUM_DATA_CHUNKS > 6
	case 6: chunk = amicableCandidates6; next_chunk = amicableCandidates7; break;
#endif
#if NUM_DATA_CHUNKS > 7
	case 7: chunk = amicableCandidates7; next_chunk = amicableCandidates8; break;
#endif
#if NUM_DATA_CHUNKS > 8
	case 8: chunk = amicableCandidates8; next_chunk = amicableCandidates9; break;
#endif
#if NUM_DATA_CHUNKS > 9
	case 9: chunk = amicableCandidates9; next_chunk = amicableCandidates10; break;
#endif
#if NUM_DATA_CHUNKS > 10
	case 10: chunk = amicableCandidates10; next_chunk = amicableCandidates11; break;
#endif
	}

	const uint offset = (uint)((amicableCandidateByteIndex & ((1UL << CHUNK_SIZE_SHIFT) - 1)) >> 3);

	ulong data0 = chunk[offset];
	ulong data1;

	if (offset + 1 < (1U << (CHUNK_SIZE_SHIFT - 3))) {
		data1 = chunk[offset + 1];
	} else {
		data1 = next_chunk[0];
	}
#endif

	const uint shift = (amicableCandidateByteIndex & 7) * 8;

	if (shift > 0) {
		data0 = (data0 >> shift) | (data1 << (64u - shift));
		data1 >>= shift;
	}

	ulong value_ulong = (data0 & 4294967295u) | ((data1 & 255u) << 32u);
	ulong sum_ulong = (data0 >> 32u) | (((data1 >> 8u) & 255u) << 32u);

	sum_ulong += (value_ulong << 1);

	const ulong p = largePrimes[largePrimeIndex];

	CheckPairPhase1(
		smallPrimes,
		primeInverses,
		PQ,
		PowerOf2SumInverses,
		PowersOfP_128SumInverses_offsets,
		PowersOfP_128SumInverses,
		PrimeInverses_128,
		SumEstimates_128,
		value_ulong * p, mul_hi(value_ulong, p), sum_ulong * (p + 1), mul_hi(sum_ulong, p + 1),
		global_offset,
		phase1_offset_to_resume_after_overflow,
		phase2_numbers_count,
		phase2_numbers_data,
		amicable_numbers_count,
		amicable_numbers_data,
		0, phase2_numbers_count, (__global ulong*)phase2_numbers_data, 0
	);
}


__kernel
__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
void CheckPairs(
	__global const uint* smallPrimes,
	__global const ulong2* primeInverses,
	__global const ulong2* PQ,
	__constant ulong2* PowerOf2SumInverses,
	__global const uint* PowersOfP_128SumInverses_offsets,
	__global const ulong4* PowersOfP_128SumInverses,
	__global const ulong2* PrimeInverses_128,
	__global const ulong* SumEstimates_128,
	__global const ulong4* pairsToCheck,
	__global const ulong* bands,
	__global uint* filtered_numbers_count,
	__global ulong4* filtered_numbers_data,
	__global uint* amicable_numbers_count,
	__global ulong4* amicable_numbers_data
)
{
	const ulong4 cur_pair = pairsToCheck[get_global_id(0)];
	if (cur_pair.x == 0)
	{
		return;
	}
	// bands: the band parameter of each number (DecideByBand; 0: no band, the test of BandPrimeDivides sets it)
	CheckPairPhase1(smallPrimes, primeInverses, PQ, PowerOf2SumInverses, PowersOfP_128SumInverses_offsets, PowersOfP_128SumInverses, PrimeInverses_128, SumEstimates_128, cur_pair.x, cur_pair.y, cur_pair.z, cur_pair.w, 0, 0, filtered_numbers_count, filtered_numbers_data, amicable_numbers_count, amicable_numbers_data, 0, filtered_numbers_count, (__global ulong*)filtered_numbers_data, bands[get_global_id(0)]);
}
)===";
