#pragma once

#define RUN_PERFORMANCE_TEST 1

enum SearchLimit : number
{
	// Search up to 10^3, 1 pair
	//value = 1000,

	// Search up to 10^4, 5 pairs
	//value = 10000,

	// Search up to 10^5, 13 pairs
	//value = 100000,

	// Search up to 10^6, 42 pairs
	//value = 1000000,

	// Search up to 10^7, 108 pairs
	//value = 10000000,

	// Search up to 10^8, 236 pairs
	//value = 100000000,

	// Search up to 10^9, 586 pairs
	//value = 1000000000,

	// Search up to 10^10, 1427 pairs
	//value = 10000000000,

	// Search up to 10^11, 3340 pairs
	value = 100000000000,

	// Search up to 10^12, 7642 pairs
	//value = 1000000000000,

	// Search up to 10^13, 17519 pairs
	//value = 10000000000000,

	// Search up to 10^14, 39374 pairs
	//value = 100000000000000,

	// Search up to 10^15, 87102 pairs
	//value = 1000000000000000,

	// Search up to 10^16, 190775 pairs
	//value = 10000000000000000,

	// Search up to 10^17, 415523 pairs
	//value = 100000000000000000,

	// Search up to 10^18, 901312 pairs
	//value = 1000000000000000000,

	// Search up to 10^19, ~1955600 pairs
	//value = 10000000000000000000,

	// Search up to 2^64, ~2401900 pairs
	//value = 18446744073709551615,
};

template<number Pred, number A, number B> struct Predicate;

template<number A, number B>
struct Predicate<1, A, B>
{
	enum Value : number
	{
		value = A,
	};
};

template<number A, number B>
struct Predicate<0, A, B>
{
	enum Value : number
	{
		value = B,
	};
};

template<number N, number A, number B>
struct CompileTimeSQRTImpl
{
	enum Value : number
	{
		C = (A / 2) + (B / 2) + ((A & B) & 1) + 1,
		Pred = (C * C > N) ? 1 : 0,
		value = CompileTimeSQRTImpl<N, Predicate<Pred, A, C>::value, Predicate<Pred, C - 1, B>::value>::value,
	};
};

template<number N, number C>
struct CompileTimeSQRTImpl<N, C, C>
{
	enum Value : number
	{
		value = C,
	};
};

template<number N>
struct CompileTimeSQRT
{
	enum Value : number
	{
		value = CompileTimeSQRTImpl<N, 1, 4294967295>::value,
	};
};

template<number A, number B>
struct Min
{
	enum X : number
	{
		value = (A < B) ? A : B,
	};
};

template<number A, number B>
struct Max
{
	enum X : number
	{
		value = (A > B) ? A : B,
	};
};

enum CompileTimeParams : number
{
	// Linear search starts when p >= SearchLimit::LinearLimit
	// It iterates over 1 * p, 2 * p, 3 * p, ... until it reaches SearchLimit::value
	// It assumes that:
	// 1) If k is abundant then k * p is also abundant and vice versa
	// 2) If k if deficient then k * p if also deficient and vice versa
	// We can guarantee these assumptions only if p / k > 2 which means linear limit must be > sqrt(limit * 2)
	// But we don't want to store too many values for linear search, so we only store values up to 10000000
	LinearLimit = ((CompileTimeSQRT<SearchLimit::value>::value + 1) * 1414213563) / 1000000000 + 1,
	MainPrimeTableBound = Max<LinearLimit, 1000>::value,
	PrimeInversesBound = CompileTimeSQRT<SearchLimit::value / 4>::value,
	ShiftMultiplier = 2,

	// Safe upper bound for the largest prime factor
	//
	// Lemma: m = k * p (p is the largest prime factor) can't be a smaller member of an amicable pair when k < 20
	//
	// Proof.
	//
	// n = S(m) - m must be > m, so let's check that for all k < 20 it's false:
	//
	// 1) For all p > 31:
	//
	// S(1 * p) - p = (p + 1) - p = 1 < p
	// S(2 * p) - 2 * p = 3 * (p + 1) - 2 * p = p + 3 < 2 * p when p > 3
	// S(3 * p) - 3 * p = 4 * (p + 1) - 3 * p = p + 4 < 3 * p when p > 2
	// S(4 * p) - 4 * p = 7 * (p + 1) - 4 * p = 3 * p + 7 < 4 * p when p > 7
	// S(5 * p) - 5 * p = 6 * (p + 1) - 5 * p = p + 6 < 5 * p when p > 1
	// S(6 * p) - 6 * p = 12 * (p + 1) - 6 * p = 6 * p + 12 > 6 * p
	// But! 6 * p = 0 (mod 6) and S(6 * p) = 12 * (p + 1) - even number, so it can't be an amicable number
	//
	// S(7 * p) - 7 * p = 8 * (p + 1) - 7 * p = p + 8 < 7 * p when p > 1
	// S(8 * p) - 8 * p = 15 * (p + 1) - 8 * p = 7 * p + 15 < 8 * p when p > 15
	// S(9 * p) - 9 * p = 13 * (p + 1) - 9 * p = 4 * p + 13 < 9 * p when p > 2
	// S(10 * p) - 10 * p = 18 * (p + 1) - 10 * p = 8 * p + 18 < 10 * p when p > 9
	// S(11 * p) - 11 * p = 12 * (p + 1) - 11 * p = p + 12 < 11 * p when p > 1
	// S(12 * p) - 12 * p = 28 * (p + 1) - 12 * p = 16 * p + 28 > 12 * p
	// But! 12 * p = 0 (mod 6) and S(12 * p) = 28 * (p + 1) - even number, so it can't be an amicable number 
	//
	// S(13 * p) - 13 * p = 14 * (p + 1) - 13 * p = p + 14 < 13 * p when p > 1
	// S(14 * p) - 14 * p = 24 * (p + 1) - 14 * p = 10 * p + 24 < 14 * p when p > 6
	// S(15 * p) - 15 * p = 24 * (p + 1) - 15 * p = 9 * p + 24 < 15 * p when p > 4
	// S(16 * p) - 16 * p = 31 * (p + 1) - 16 * p = 15 * p + 31 < 16 * p when p > 31
	// S(17 * p) - 17 * p = 18 * (p + 1) - 17 * p = p + 18 < 17 * p when p > 1
	// S(18 * p) - 18 * p = 39 * (p + 1) - 18 * p = 21 * p + 39 > 18 * p
	// But! 18 * p = 0 (mod 6) and S(18 * p) = 39 * (p + 1) - even number because (p + 1) is even, so it can't be an amicable number 
	//
	// S(19 * p) - 19 * p = 20 * (p + 1) - 19 * p = p + 20 < 19 * p when p > 2
	// S(20 * p) - 20 * p = 42 * (p + 1) - 20 * p = 22 * p + 42 > 20 * p, so we must stop here
	//
	// 2) Amicable pairs where smaller member has the form k * p, k < 20, p <= 31
	// They must be <= 19 * 31 = 589, so it can only be 220 = 20 * 11, so k = 20 in this case (m = k * p)
	//
	// Lemma is proved.
	//
	// We can use SearchLimit::value / 20 as a safe upper bound
	// because if largest prime factor p > SearchLimit::value / 20 then we'll only have numbers of the form k * p where k < 20
	// SearchLimit::value / 20 is also the best possible bound because 220 = 20 * 11
	SafeLimit = value / 20,
};

// Let S(n) be the sum of all divisors of n
// Factorization in general form has k distinct primes: n = P1^e1*P2^e2*...*Pk^ek, P1 < P2 < ... < Pk
//
// Smallest number with 16 distinct primes in factorization is:
// 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29 * 31 * 37 * 41 * 43 * 47 * 53 = 32589158477190044730 ~ 1.766662 * 2^64 > 2^64
//
// Therefore, if n <= 2^64, it'll have at most 15 distinct primes in factorization
//
// Let's estimate the maximum value of S(n) / n, n <= 2^64:
// S(n) / n = multiply(i = 1...k, (Pi^(ei + 1) - 1) / (Pi^ei * (Pi - 1)))
//
// Divide numerator and denominator by Pi^ei for each i:
// S(n) / n = multiply(i = 1...k, (Pi - Pi^-ei) / (Pi - 1))
//
// Since Pi^-ei > 0, we can estimate each multiplier to be < Pi / (Pi - 1):
// S(n) / n < F(n) = multiply(i = 1...k, Pi / (Pi - 1))
//
// f(x) = x / (x - 1) is a monotonic function:
// f'(x + 1) = ((x + 1) / x)' = (1 + 1 / x)' = -1 / x^2,
// f'(x) = -1 / (x - 1)^2 < 0, so it's monotonically decreasing for all x > 1: 1 < x1 <= x2 implies that f(x1) >= f(x2)
//
// Let's get back to F(n). If we replace P1, ..., Pk distinct primes with _first_ k primes, each multiplier will increase or remain the same, because:
// 2 <= P1, 3 <= P2, 5 <= P3, ..., kth prime <= Pk
//
// F(n) <= multiply(i = 1...k, (ith prime) / (ith prime - 1))
//
// We've come to this: S(n) / n < F(n) <= multiply(i = 1...k, (ith prime) / (ith prime - 1))
// And we know that if n <= 2^64 then k <= 15, so S(n) / n < 2 / 1 * 3 / 2 * 5 / 4 * ... * 47 / 46 = 7.2095926010299534105882364948723...
// S(n) < n * 7.2095926010299534105882364948723...
//
// S(n) must be less than 2^64 which means that n must be < 2^64 / 7.2095926010299534105882364948723... = 2.558638898829633279... * 10^18
// So 2.558638898829633279... * 10^18 is a theoretical limit
//
// But the practical limit is 7 * 10^18 - it's approximately where we start getting pairs with pair sum >= 2^64
// The first actually missed known pair is (7364096816271422535, 11119300665402788793) with pair sum 1.0019869852272140775539810420014 * 2^64
//static_assert(SearchLimit::value <= 7000000000000000000ULL, "Search limit is too large, some pairs can be missed");

template<number N> struct CompileTimePrimes;

template<> struct CompileTimePrimes<0>  { enum { value =  2 }; };
template<> struct CompileTimePrimes<1>  { enum { value =  3 }; };
template<> struct CompileTimePrimes<2>  { enum { value =  5 }; };
template<> struct CompileTimePrimes<3>  { enum { value =  7 }; };
template<> struct CompileTimePrimes<4>  { enum { value = 11 }; };
template<> struct CompileTimePrimes<5>  { enum { value = 13 }; };
template<> struct CompileTimePrimes<6>  { enum { value = 17 }; };
template<> struct CompileTimePrimes<7>  { enum { value = 19 }; };
template<> struct CompileTimePrimes<8>  { enum { value = 23 }; };
template<> struct CompileTimePrimes<9>  { enum { value = 29 }; };
template<> struct CompileTimePrimes<10> { enum { value = 31 }; };
template<> struct CompileTimePrimes<11> { enum { value = 37 }; };
template<> struct CompileTimePrimes<12> { enum { value = 41 }; };
template<> struct CompileTimePrimes<13> { enum { value = 43 }; };
template<> struct CompileTimePrimes<14> { enum { value = 47 }; };
template<> struct CompileTimePrimes<15> { enum { value = 53 }; };
template<> struct CompileTimePrimes<16> { enum { value = 59 }; };
template<> struct CompileTimePrimes<17> { enum { value = 61 }; };
template<> struct CompileTimePrimes<18> { enum { value = 67 }; };
template<> struct CompileTimePrimes<19> { enum { value = 71 }; };
template<> struct CompileTimePrimes<20> { enum { value = 73 }; };
template<> struct CompileTimePrimes<21> { enum { value = 79 }; };
template<> struct CompileTimePrimes<22> { enum { value = 83 }; };
template<> struct CompileTimePrimes<23> { enum { value = 89 }; };
template<> struct CompileTimePrimes<24> { enum { value = 97 }; };
template<> struct CompileTimePrimes<25> { enum { value = 101 }; };
template<> struct CompileTimePrimes<26> { enum { value = 103 }; };
template<> struct CompileTimePrimes<27> { enum { value = 107 }; };
template<> struct CompileTimePrimes<28> { enum { value = 109 }; };
template<> struct CompileTimePrimes<29> { enum { value = 113 }; };
template<> struct CompileTimePrimes<30> { enum { value = 127 }; };
template<> struct CompileTimePrimes<31> { enum { value = 131 }; };
template<> struct CompileTimePrimes<32> { enum { value = 137 }; };
template<> struct CompileTimePrimes<33> { enum { value = 139 }; };
template<> struct CompileTimePrimes<34> { enum { value = 149 }; };
template<> struct CompileTimePrimes<35> { enum { value = 151 }; };
template<> struct CompileTimePrimes<36> { enum { value = 157 }; };
template<> struct CompileTimePrimes<37> { enum { value = 163 }; };
template<> struct CompileTimePrimes<38> { enum { value = 167 }; };
template<> struct CompileTimePrimes<39> { enum { value = 173 }; };
template<> struct CompileTimePrimes<40> { enum { value = 179 }; };
template<> struct CompileTimePrimes<41> { enum { value = 181 }; };
template<> struct CompileTimePrimes<42> { enum { value = 191 }; };
template<> struct CompileTimePrimes<43> { enum { value = 193 }; };
template<> struct CompileTimePrimes<44> { enum { value = 197 }; };
template<> struct CompileTimePrimes<45> { enum { value = 199 }; };
template<> struct CompileTimePrimes<46> { enum { value = 211 }; };
template<> struct CompileTimePrimes<47> { enum { value = 223 }; };
template<> struct CompileTimePrimes<48> { enum { value = 227 }; };
template<> struct CompileTimePrimes<49> { enum { value = 229 }; };
template<> struct CompileTimePrimes<50> { enum { value = 233 }; };
template<> struct CompileTimePrimes<51> { enum { value = 239 }; };
template<> struct CompileTimePrimes<52> { enum { value = 241 }; };
template<> struct CompileTimePrimes<53> { enum { value = 251 }; };
template<> struct CompileTimePrimes<54> { enum { value = 257 }; };
template<> struct CompileTimePrimes<55> { enum { value = 263 }; };
template<> struct CompileTimePrimes<56> { enum { value = 269 }; };
template<> struct CompileTimePrimes<57> { enum { value = 271 }; };
template<> struct CompileTimePrimes<58> { enum { value = 277 }; };
template<> struct CompileTimePrimes<59> { enum { value = 281 }; };
template<> struct CompileTimePrimes<60> { enum { value = 283 }; };
template<> struct CompileTimePrimes<61> { enum { value = 293 }; };
template<> struct CompileTimePrimes<62> { enum { value = 307 }; };
template<> struct CompileTimePrimes<63> { enum { value = 311 }; };
template<> struct CompileTimePrimes<64> { enum { value = 313 }; };
template<> struct CompileTimePrimes<65> { enum { value = 317 }; };
template<> struct CompileTimePrimes<66> { enum { value = 331 }; };
template<> struct CompileTimePrimes<67> { enum { value = 337 }; };
template<> struct CompileTimePrimes<68> { enum { value = 347 }; };
template<> struct CompileTimePrimes<69> { enum { value = 349 }; };
template<> struct CompileTimePrimes<70> { enum { value = 353 }; };
template<> struct CompileTimePrimes<71> { enum { value = 359 }; };
template<> struct CompileTimePrimes<72> { enum { value = 367 }; };
template<> struct CompileTimePrimes<73> { enum { value = 373 }; };
template<> struct CompileTimePrimes<74> { enum { value = 379 }; };
template<> struct CompileTimePrimes<75> { enum { value = 383 }; };
template<> struct CompileTimePrimes<76> { enum { value = 389 }; };
template<> struct CompileTimePrimes<77> { enum { value = 397 }; };
template<> struct CompileTimePrimes<78> { enum { value = 401 }; };
template<> struct CompileTimePrimes<79> { enum { value = 409 }; };
template<> struct CompileTimePrimes<80> { enum { value = 419 }; };
template<> struct CompileTimePrimes<81> { enum { value = 421 }; };
template<> struct CompileTimePrimes<82> { enum { value = 431 }; };
template<> struct CompileTimePrimes<83> { enum { value = 433 }; };
template<> struct CompileTimePrimes<84> { enum { value = 439 }; };
template<> struct CompileTimePrimes<85> { enum { value = 443 }; };
template<> struct CompileTimePrimes<86> { enum { value = 449 }; };
template<> struct CompileTimePrimes<87> { enum { value = 457 }; };
template<> struct CompileTimePrimes<88> { enum { value = 461 }; };
template<> struct CompileTimePrimes<89> { enum { value = 463 }; };
template<> struct CompileTimePrimes<90> { enum { value = 467 }; };
template<> struct CompileTimePrimes<91> { enum { value = 479 }; };
template<> struct CompileTimePrimes<92> { enum { value = 487 }; };
template<> struct CompileTimePrimes<93> { enum { value = 491 }; };
template<> struct CompileTimePrimes<94> { enum { value = 499 }; };
template<> struct CompileTimePrimes<95> { enum { value = 503 }; };
template<> struct CompileTimePrimes<96> { enum { value = 509 }; };
template<> struct CompileTimePrimes<97> { enum { value = 521 }; };

enum { CompileTimePrimesCount = 96 };

enum PrimeTableParameters
{
	Modulo = 210,
	NumOffsets = 48,
};

static_assert(PrimeTableParameters::Modulo < 256, "Modulo is too large");

const unsigned int NumbersCoprimeToModulo[NumOffsets * 2] = {
	  1,  11,  13,  17,  19,  23,  29,  31,
	 37,  41,  43,  47,  53,  59,  61,  67,
	 71,  73,  79,  83,  89,  97, 101, 103,
	107, 109, 113, 121, 127, 131, 137, 139,
	143, 149, 151, 157, 163, 167, 169, 173,
	179, 181, 187, 191, 193, 197, 199, 209,

	  1 + 210,  11 + 210,  13 + 210,  17 + 210,  19 + 210,  23 + 210,  29 + 210,  31 + 210,
	 37 + 210,  41 + 210,  43 + 210,  47 + 210,  53 + 210,  59 + 210,  61 + 210,  67 + 210,
	 71 + 210,  73 + 210,  79 + 210,  83 + 210,  89 + 210,  97 + 210, 101 + 210, 103 + 210,
	107 + 210, 109 + 210, 113 + 210, 121 + 210, 127 + 210, 131 + 210, 137 + 210, 139 + 210,
	143 + 210, 149 + 210, 151 + 210, 157 + 210, 163 + 210, 167 + 210, 169 + 210, 173 + 210,
	179 + 210, 181 + 210, 187 + 210, 191 + 210, 193 + 210, 197 + 210, 199 + 210, 209 + 210,
};

enum Byte
{
	Bits = 8,
};

FORCEINLINE bool IsPerfectSquareCandidate(const number n)
{
	// If N is a perfect square, then N mod 64 must be one of the following: 0, 1, 4, 9, 16, 17, 25, 33, 36, 41, 49, 57
	enum Modulo64SquareCheck : number
	{
		value =
		(number(1) << 0) |
		(number(1) << 1) |
		(number(1) << 4) |
		(number(1) << 9) |
		(number(1) << 16) |
		(number(1) << 17) |
		(number(1) << 25) |
		(number(1) << 33) |
		(number(1) << 36) |
		(number(1) << 41) |
		(number(1) << 49) |
		(number(1) << 57)
	};
	return (((number(1) << n) & Modulo64SquareCheck::value) != 0);
}

#pragma warning(push)

// Disable this warning here because number overflows are the essence of this algorithm: it computes everything modulo 2^64.
#pragma warning(disable : 4307) // warning C4307: '*' : integral constant overflow

template<number x1, number x2, number v1, number v2>
struct ExtendedEuclidCompileTime
{
	enum Variables : number
	{
		q = v1 / v2,
		x3 = x1 - q * x2,
		v3 = v1 - q * v2,
		value = ExtendedEuclidCompileTime<x2, x3, v2, v3>::value,
	};
};

template<number x1, number x2, number v1>
struct ExtendedEuclidCompileTime<x1, x2, v1, 1>
{
	enum Variables : number
	{
		value = x2,
	};
};

template<number N>
struct MultiplicativeInverse
{
	enum Variables : number
	{
		value = ExtendedEuclidCompileTime<number(-1), 1, ~N + 1, N>::value,
	};
};

// For even numbers
template<number N, number K, bool IsEven>
struct MultiplicativeInverseEvenImpl
{
	enum Variables : number
	{
		value = MultiplicativeInverseEvenImpl<N / 2, K + 1, ((N / 2) % 2) == 0>::value,
		shift = MultiplicativeInverseEvenImpl<N / 2, K + 1, ((N / 2) % 2) == 0>::shift,
	};
};

template<number N, number K>
struct MultiplicativeInverseEvenImpl<N, K, false>
{
	enum Variables : number
	{
		value = MultiplicativeInverse<N>::value,
		shift = K,
	};
};

template<number N>
struct MultiplicativeInverseEven
{
	enum Variables : number
	{
		value = MultiplicativeInverseEvenImpl<N, 0, (N % 2) == 0>::value,
		shift = MultiplicativeInverseEvenImpl<N, 0, (N % 2) == 0>::shift,
	};
};

#pragma warning(pop)
