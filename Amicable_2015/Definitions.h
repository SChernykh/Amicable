#pragma once

typedef unsigned char byte;
typedef unsigned long long int number;

#define RUN_PERFORMANCE_TEST 1

// The code can use some conjectures which are true for all known amicable pairs, but are still not proved
// Set USE_CONJECTURES to 0 when doing an exhaustive search
// Set USE_CONJECTURES to 1 when the guarantee to find all amicable pairs in the specified range is not needed
#define USE_CONJECTURES 0

enum SearchLimit : number
{
	// Search up to 10^3, 1 pair
	//SQRT2 = 32,
	//value = 1000,

	// Search up to 10^4, 4 pairs
	//SQRT2 = 100,
	//value = 10000,

	// Search up to 10^5, 8 pairs
	//SQRT2 = 317,
	//value = 100000,

	// Search up to 10^6, 29 pairs
	//SQRT2 = 1000,
	//value = 1000000,

	// Search up to 10^7, 66 pairs
	//SQRT2 = 3163,
	//value = 10000000,

	// Search up to 10^8, 128 pairs
	//SQRT2 = 10000,
	//value = 100000000,

	// Search up to 10^9, 350 pairs
	//SQRT2 = 31623,
	//value = 1000000000,

	// Search up to 10^10, 841 pairs
	//SQRT2 = 100000,
	//value = 10000000000,

	// Search up to 10^11, 1913 pairs
	SQRT2 = 316228,
	value = 100000000000,

	// Search up to 10^12, 4302 pairs
	//SQRT2 = 1000000,
	//value = 1000000000000,

	// Search up to 10^13, 9877 pairs
	//SQRT2 = 3162278,
	//value = 10000000000000,

	// Search up to 10^14, 21855 pairs
	//SQRT2 = 10000000,
	//value = 100000000000000,

	// Search up to 10^15, 47728 pairs
	//SQRT2 = 31622777,
	//value = 1000000000000000,

	// Search up to 10^16, 103673 pairs
	//SQRT2 = 100000000,
	//value = 10000000000000000,

	// Search up to 10^17, 224748 pairs
	//SQRT2 = 316227767,
	//value = 100000000000000000,

	// Search up to 10^18, ~485800 pairs
	//SQRT2 = 1000000000,
	//value = 1000000000000000000,

	// Search up to 10^19, ~1048700 pairs
	//SQRT2 = 3162277661,
	//value = 10000000000000000000,

	// We need to calculate all primes to factorize numbers up to S(n) where n <= SearchLimit::value
	// Since S(n) / n < 2 / 1 * 3 / 2 * 5 / 4 * ... * 47 / 46 = 7.2095926010299534105882364948723... for all 64-bit numbers, S(n) < n * 7.2095926010299534105882364948723...
	// So we need to calculate all primes up to sqrt(n * 7.2095926010299534105882364948723...) = SearchLimit::SQRT2 * 2.6850684536953528869874068825892...
	PrimesUpToSqrtLimitValue = (SQRT2 * 268507) / 100000 + 1,

	// The first prime gap larger than 256 is 436273009 - 436273291
	// So we don't need a multiplier up until then
	ShiftMultiplier = (PrimesUpToSqrtLimitValue < 436273009) ? 1 : 2,

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

	// Linear search starts when p >= SearchLimit::LinearLimit
	// It iterates over 1 * p, 2 * p, 3 * p, ... until it reaches SearchLimit::value
	// It assumes that:
	// 1) If k is abundant then k * p is also abundant and vice versa
	// 2) If k if deficient then k * p if also deficient and vice versa
	// We can guarantee these assumptions only if p / k > 2 which means linear limit must be > sqrt(limit * 2)
	// But we don't want to store too many values for linear search, so we only store values up to 10000000
#pragma warning(push)
#pragma warning(disable : 4296)
	LinearLimit = max(value / 10000000, (SearchLimit::SQRT2 * 141422) / 100000 + 1),
#pragma warning(pop)
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
// The first actually missed known pair is (8428612252940935965, 10094556152778199011) with pair sum 1.004142971339777187789188417355 * 2^64
static_assert(SearchLimit::value <= 7000000000000000000ULL, "Search limit is too large, some pairs can be missed");

// Linear search starts when p >= SearchLimit::LinearLimit
// It iterates over 1 * p, 2 * p, 3 * p, ... until it reaches SearchLimit::value
// It assumes that:
// 1) If k is abundant then k * p is also abundant and vice versa
// 2) If k if deficient then k * p if also deficient and vice versa
// We can guarantee these assumptions only if p / k > 2 which means linear limit must be > sqrt(limit * 2)
static_assert(SearchLimit::LinearLimit > (SearchLimit::SQRT2 * 1.41422) / 100000 + 1, "Linear search limit is too low");

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

enum { CompileTimePrimesCount = 64 };

enum { InlinePrimesInSearch = 10 };

#pragma warning(push)
#pragma warning(disable : 4307)
template<int index, number M, number N, bool finish>
struct MaxSearchDepthHelper
{
	enum
	{
		value = MaxSearchDepthHelper<index + 1, M * CompileTimePrimes<index + 1>::value, N, (N / CompileTimePrimes<index + 1>::value) < M>::value,
	};
};

template<int index, number M, number N>
struct MaxSearchDepthHelper<index, M, N, true>
{
	enum
	{
		value = index - 1 - InlinePrimesInSearch,
	};
};

template<number N>
struct MaxSearchDepth
{
	enum
	{
		value = MaxSearchDepthHelper<InlinePrimesInSearch, 1, N, false>::value,
	};
};
#pragma warning(pop)

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

number MaximumSumOfDivisors3(const number a, const number p, const number q);

#define IF_CONSTEXPR(X) \
	__pragma(warning(push)) \
	__pragma(warning(disable:4127)) \
	static_assert((X) || !(X), "Error: "#X" is not a constexpr"); \
	if (X) \
	__pragma(warning(pop))

#if RUN_PERFORMANCE_TEST
#define COUT(X)
#else
#define COUT(X)
#endif

extern unsigned __int64(__fastcall *udiv128)(unsigned __int64 numhi, unsigned __int64 numlo, unsigned __int64 den, unsigned __int64* rem);
