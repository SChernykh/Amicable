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
extern number privPrimesUpToSqrtLimitSortedCount;
extern byte privNextPrimeShifts[ShiftTableSize];
extern std::vector<std::pair<unsigned int, unsigned int>> privLinearSearchData;
extern number privPrimeInverses[CompileTimePrimesCount * 2];

#define NextPrimeShifts ((const byte* const)(privNextPrimeShifts))
#define PrimesUpToSqrtLimitSortedCount ((const unsigned int)(privPrimesUpToSqrtLimitSortedCount))
#define LinearSearchData ((const std::vector<std::pair<unsigned int, unsigned int>>&)(privLinearSearchData))

#define PrimeInverses ((const number*)(privPrimeInverses))

FORCEINLINE number CalculateInverse(number n)
{
	number x1 = number(-1);
	number x2 = 1;
	number v1 = ~n + 1;
	number v2 = n;
	do
	{
		const number q = v1 / v2;
		const number x3 = x1 - q * x2;
		const number v3 = v1 - q * v2;
		x1 = x2;
		x2 = x3;
		v1 = v2;
		v2 = v3;
	} while (v2 > 1);
	return x2;
}

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
extern number g_MaxSumRatios[64];

FORCEINLINE number GetMaxSumMDiv2(number D, number sumD)
{
	unsigned long index;
	_BitScanReverse64(&index, D);

	number h;
	_umul128(sumD, g_MaxSumRatios[index], &h);
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

template<number Index>
FORCEINLINE bool IsValidCandidate(const number D, const number sumD, number gcd, number sumGCD)
{
	enum { p = CompileTimePrimes<Index>::value };

	number qa = D * MultiplicativeInverse<p>::value;
	if (qa <= number(-1) / p)
	{
		number qb = sumD * MultiplicativeInverse<p>::value;
		if (qb <= number(-1) / p)
		{
			number n = p;
			number curSum = p + 1;
			number a = qa;
			number b = qb;
			for (;;)
			{
				qa = a * MultiplicativeInverse<p>::value;
				if (qa > number(-1) / p)
					break;
				qb = b * MultiplicativeInverse<p>::value;
				if (qb > number(-1) / p)
					break;
				n *= p;
				curSum += n;
				a = qa;
				b = qb;
			}
			gcd *= n;
			sumGCD *= curSum;
		}
	}
	return IsValidCandidate<Index + 1>(D, sumD, gcd, sumGCD);
}

template<> FORCEINLINE bool IsValidCandidate<InlinePrimesInSearch + 1>(const number D, const number sumD, number gcd, number sumGCD)
{
	// "sumGCD * (sumD - D) < gcd * sumD" must be true
	number sum1[2];
	sum1[0] = _umul128(sumGCD, sumD - D, &sum1[1]);

	number sum2[2];
	sum2[0] = _umul128(gcd, sumD, &sum2[1]);

	return ((sum1[1] < sum2[1]) || ((sum1[1] == sum2[1]) && (sum1[0] < sum2[0])));
}

FORCEINLINE bool IsValidCandidate(const number D, const number sumD)
{
	// If one of numbers is odd, then it's very unlikely for it to not be a candidate
	// Such numbers exist, but they're very rare, so it's better to skip factorization and exit
	// My tests show that this check speeds up the search
	const number allBitsCombined = D | sumD;
	if (allBitsCombined & 1)
		return true;

	unsigned long index;
	_BitScanForward64(&index, allBitsCombined);

	const number gcd = number(1) << index;
	const number sumGCD = gcd * 2 - 1;
	return IsValidCandidate<1>(D, sumD, gcd, sumGCD);
}

void PrimeTablesInit(bool isSubmit);
number CalculatePrimes(number aLowerBound, number anUpperBound, std::vector<byte>& anOutPrimes);

// Works only for numbers up to SearchLimit::PrimesUpToSqrtLimitValue
bool IsPrime(number n);

class PrimesUpToSqrtLimitIterator
{
public:
	explicit PrimesUpToSqrtLimitIterator(number aStartPrime)
		: mySieveChunk(0xfafd7bbef7ffffffULL & ~number(3))
		, mySieveData((const number*)(PrimesUpToSqrtLimit.data()))
		, myPossiblePrimesForModuloPtr(NumbersCoprimeToModulo)
		, myModuloIndex(0)
		, myBitIndexShift(0)
		, myCurrentPrime(1)
	{
		while (myCurrentPrime < aStartPrime)
			operator++();
	}

	FORCEINLINE number Get() const { return myCurrentPrime; }

	FORCEINLINE PrimesUpToSqrtLimitIterator& operator++()
	{
		if (myCurrentPrime < 11)
		{
			myCurrentPrime = (0xB0705320UL >> (myCurrentPrime * 4)) & 15;
			return *this;
		}

		while (!mySieveChunk)
		{
			mySieveChunk = *(++mySieveData);

			const number NewValues = (PrimeTableParameters::Modulo / 2) | (number(PrimeTableParameters::Modulo / 2) << 16) | (number(PrimeTableParameters::Modulo) << 32) |
				(16 << 8) | (32 << 24) | (number(0) << 40);

			myModuloIndex += ((NewValues >> myBitIndexShift) & 255) * 2;
			myBitIndexShift = (NewValues >> (myBitIndexShift + 8)) & 255;

			myPossiblePrimesForModuloPtr = NumbersCoprimeToModulo + myBitIndexShift;
		}

		unsigned long bitIndex;
		_BitScanForward64(&bitIndex, mySieveChunk);
		mySieveChunk &= (mySieveChunk - 1);

		myCurrentPrime = myModuloIndex + myPossiblePrimesForModuloPtr[bitIndex];

		return *this;
	}

private:
	number mySieveChunk;
	const number* mySieveData;
	const unsigned int* myPossiblePrimesForModuloPtr;
	number myModuloIndex;
	number myBitIndexShift;
	number myCurrentPrime;
};
