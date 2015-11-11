#pragma once

enum
{
	// There are exactly 192725 primes below 2^(64/3)
	// We can use this table for factorization when p^3 <= N < 2^64
	ReciprocalsTableSize = 192725,
};

// Reciprocals are calculated using algorithm published in http://www.agner.org/optimize/optimizing_assembly.pdf (section 16.9 "Integer division by a constant")
#pragma pack(push, 1)
struct SReciprocal
{
	number reciprocal;
	unsigned char increment;
	unsigned char shift;

	NOINLINE void Init(const number aDivisor)
	{
		unsigned long bitIndex;
		_BitScanReverse64(&bitIndex, aDivisor);
		shift = static_cast<unsigned char>(bitIndex);

		number remainder;
		number quotient = udiv128(number(1) << shift, 0, aDivisor, &remainder);
		if (remainder * 2 < aDivisor)
		{
			reciprocal = quotient;
			increment = 1;
		}
		else
		{
			reciprocal = quotient + 1;
			increment = 0;
		}
	}

	FORCEINLINE bool Divide(const number n, const number divisor, number& q) const
	{
		number highProduct;
		_umul128(n + increment, reciprocal, &highProduct);
		q = highProduct >> shift;
		return q * divisor == n;
	}

	FORCEINLINE number DivideNoRemainder(const number n) const
	{
		number highProduct;
		_umul128(n + increment, reciprocal, &highProduct);
		return highProduct >> shift;
	}
};
#pragma pack(pop)

extern SReciprocal privPrimeReciprocals[ReciprocalsTableSize];
#define PrimeReciprocals ((const SReciprocal* const)(privPrimeReciprocals))

enum
{
	ShiftTableSize = 129982215,
};

extern std::vector<byte> PrimesUpToSqrtLimit;
extern unsigned int privPrimesUpToSqrtLimitSorted[ShiftTableSize];
extern number privPrimesUpToSqrtLimitSortedCount;
extern byte privNextPrimeShifts[ShiftTableSize];
extern std::vector<std::pair<unsigned int, unsigned int>> privLinearSearchData;
#define NextPrimeShifts ((const byte* const)(privNextPrimeShifts))
#define PrimesUpToSqrtLimitSorted ((const unsigned int* const)(privPrimesUpToSqrtLimitSorted))
#define PrimesUpToSqrtLimitSortedCount ((const unsigned int)(privPrimesUpToSqrtLimitSortedCount))
#define LinearSearchData ((const std::vector<std::pair<unsigned int, unsigned int>>&)(privLinearSearchData))

struct SumEstimateData
{
	number P;
	number Q;
};

enum
{
	SumEstimatesSize = 16,
	IS_NUM_ELIGIBLE_BEGIN = 16,
};

extern const SumEstimateData* privSumEstimates[SumEstimatesSize];

#define SumEstimates ((const SumEstimateData * const * const)(privSumEstimates))

struct SmallFactorNumbers
{
	std::vector<std::pair<number, number>> myNumbers;
};

extern SmallFactorNumbers g_SmallFactorNumbersData;
extern number g_MaxSumRatios[MaxSearchDepth<SearchLimit::value>::value + 1];

template<number Index>
struct MaxSumMDiv2Index
{
	template<number PrimeProduct>
	FORCEINLINE static number get(number D, number smallestFactorInM1)
	{
		enum { p = CompileTimePrimes<Index + InlinePrimesInSearch + 1>::value };

		if ((D > SearchLimit::value / (PrimeProduct * p)) || (smallestFactorInM1 <= p))
			return Index;

		return MaxSumMDiv2Index<Index + 1>::get<PrimeProduct * p>(D, smallestFactorInM1);
	}
};

template<>
struct MaxSumMDiv2Index<MaxSearchDepth<SearchLimit::value>::value>
{
	template<number PrimeProduct>
	FORCEINLINE static number get(number, number)
	{
		return MaxSearchDepth<SearchLimit::value>::value;
	}
};

FORCEINLINE number GetMaxSumMDiv2(number D, number sumD, number smallestFactorInM1)
{
	number h;
	_umul128(sumD, g_MaxSumRatios[MaxSumMDiv2Index<0>::get<1>(D, smallestFactorInM1)], &h);
	return h;
}

FORCEINLINE number GCD(number a, number b)
{
	while (b)
	{
		const number prev_a = a;
		a = b;
		b = prev_a % b;
	}
	return a;
};

void PrimeTablesInit(bool isSubmit);
number CalculatePrimes(number aLowerBound, number anUpperBound, std::vector<byte>& anOutPrimes);

// Works only for numbers up to SearchLimit::PrimesUpToSqrtLimitValue
bool IsPrime(number n);
