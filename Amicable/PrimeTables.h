#pragma once

void PrimeTablesInit();
number CalculatePrimes(number aLowerBound, number anUpperBound, std::vector<byte>& anOutPrimes);
bool IsPrime(number n);

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

extern CACHE_ALIGNED SReciprocal privPrimeReciprocals[ReciprocalsTableSize];
#define PrimeReciprocals ((const SReciprocal* const)(privPrimeReciprocals))

#pragma pack(push, 1)
struct LinearSearchDataEntry
{
	LinearSearchDataEntry() {}

	LinearSearchDataEntry(unsigned int _value, unsigned int _sum, unsigned char is_over_abundant_mask)
		: value(_value)
		, sum(_sum)
		, is_not_over_abundant_mask(static_cast<unsigned char>(~is_over_abundant_mask))
	{
	}

	unsigned int value;
	unsigned int sum;
	unsigned char is_not_over_abundant_mask;
};
#pragma pack(pop)

extern std::vector<byte> MainPrimeTable;
extern byte bitOffset[PrimeTableParameters::Modulo];
extern byte* privNextPrimeShifts;
extern std::vector<LinearSearchDataEntry> privCandidatesData;
extern CACHE_ALIGNED unsigned char privCandidatesDataMask[5 * 7 * 11];
extern std::pair<number, number>* privPrimeInverses;
extern CACHE_ALIGNED std::pair<number, number> privPrimeInverses2[CompileTimePrimesCount];

#define NextPrimeShifts ((const byte* const)(privNextPrimeShifts))
#define CandidatesData ((const std::vector<LinearSearchDataEntry>&)(privCandidatesData))
#define CandidatesDataMask ((const unsigned char*)(privCandidatesDataMask))

#define PrimeInverses ((const std::pair<number, number>*)(privPrimeInverses))
#define PrimeInverses2 ((const std::pair<number, number>*)(privPrimeInverses2))

FORCEINLINE number GetLinearSearchDataRemainder(const number n)
{
	number highProduct;
	static_assert(ARRAYSIZE(privCandidatesDataMask) == 5 * 7 * 11, "!!! Recalculate these constants (12265886968492584971ULL and 8) if ARRAYSIZE(privLinearSearchDataIndex) changes !!!");
	_umul128(n, 12265886968492584971ULL, &highProduct);
	return n - (highProduct >> 8) * ARRAYSIZE(privCandidatesDataMask);
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
extern CACHE_ALIGNED number privSumEstimatesBeginP[SumEstimatesSize];
extern CACHE_ALIGNED number privSumEstimatesBeginQ[SumEstimatesSize];

#define SumEstimates ((const SumEstimateData * const * const)(privSumEstimates))
#define SumEstimatesBeginP (((const number* const)(privSumEstimatesBeginP)))
#define SumEstimatesBeginQ (((const number* const)(privSumEstimatesBeginQ)))

extern CACHE_ALIGNED byte privIsNotOverAbundantMod385[128 * 512];

#define IsNotOverAbundantMod385 ((const byte*)(privIsNotOverAbundantMod385))

FORCEINLINE number GCD(number a, number b)
{
	if (a == 0) return b;
	if (b == 0) return a;

	unsigned long shift;
	_BitScanForward64(&shift, a | b);

	unsigned long index_a;
	_BitScanForward64(&index_a, a);
	a >>= index_a;

	do
	{
		unsigned long index_b;
		_BitScanForward64(&index_b, b);
		b >>= index_b;

		const number a1 = a;
		const number b1 = b;
		a = (a1 > b1) ? b1 : a1;
		b = (a1 > b1) ? a1 : b1;

		b -= a;
	} while (b);

	return (a << shift);
}

class PrimeIterator
{
public:
	explicit FORCEINLINE PrimeIterator(const number* aSieveData = reinterpret_cast<const number*>(MainPrimeTable.data()))
		: mySieveChunk(*aSieveData & ~number(1))
		, mySieveData(aSieveData)
		, myPossiblePrimesForModuloPtr(NumbersCoprimeToModulo)
		, myModuloIndex(0)
		, myBitIndexShift(0)
		, myCurrentPrime(2)
	{
	}

	explicit PrimeIterator(number aStartNumber, const number* aSieveData = (const number*)(MainPrimeTable.data()))
		: mySieveData(aSieveData)
		, myPossiblePrimesForModuloPtr(NumbersCoprimeToModulo)
	{
		if (aStartNumber < 2)
		{
			aStartNumber = 2;
		}

		const number chunkIndex = ((aStartNumber / PrimeTableParameters::Modulo) * (PrimeTableParameters::NumOffsets / Byte::Bits)) / sizeof(number);
		mySieveData = aSieveData + chunkIndex;
		mySieveChunk = *mySieveData;
		myModuloIndex = (chunkIndex / 3) * PrimeTableParameters::Modulo * 4 + (chunkIndex % 3) * PrimeTableParameters::Modulo;
		myBitIndexShift = (chunkIndex % 3) * 16;
		myPossiblePrimesForModuloPtr = NumbersCoprimeToModulo + myBitIndexShift;

		myCurrentPrime = static_cast<number>(-1);
		for (;;)
		{
			operator++();
			if (myCurrentPrime >= aStartNumber)
			{
				break;
			}
		}
	}

	FORCEINLINE number Get() const { return myCurrentPrime; }

	FORCEINLINE PrimeIterator& operator++()
	{
		if (myCurrentPrime < 7)
		{
			myCurrentPrime = (0x705320UL >> (myCurrentPrime * 4)) & 15;
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

struct Factor
{
	number p;
	unsigned int k;
	int index;
	number p_inv;
	number q_max;
};

FORCEINLINE byte OverAbundant(const Factor* f, int last_factor_index, const number value, const number sum, number sum_for_gcd_coeff)
{
	number g = 1;
	number sum_g = 1;
	number sum_for_gcd = sum * sum_for_gcd_coeff;

	const Factor* last_factor = f + last_factor_index;

	if (f->p == 2)
	{
		DWORD power_of_2;
		_BitScanForward64(&power_of_2, sum_for_gcd);
		const DWORD k = static_cast<DWORD>(f->k);
		if (power_of_2 > k)
		{
			power_of_2 = k;
		}
		sum_for_gcd >>= power_of_2;
		g <<= power_of_2;
		sum_g <<= power_of_2;
		sum_g = sum_g * 2 - 1;
		++f;
	}

	while (f <= last_factor)
	{
		const number prev_sum_g = sum_g;
		for (unsigned int j = 0; j < f->k; ++j)
		{
			const number q = sum_for_gcd * f->p_inv;
			if (q > f->q_max)
			{
				break;
			}
			sum_for_gcd = q;
			g *= f->p;
			sum_g = sum_g * f->p + prev_sum_g;
		}
		++f;
	}

	number n1[2];
	number n2[2];
	n1[0] = _umul128(sum_g - g, sum - value, &n1[1]);
	n2[0] = _umul128(g, value, &n2[1]);
#if _MSC_VER >= 1900
	return _subborrow_u64(_subborrow_u64(1, n2[0], n1[0], &n2[0]), n2[1], n1[1], &n2[1]);
#else
	return static_cast<byte>((n2[1] < n1[1]) || ((n2[1] == n1[1]) && (n2[0] <= n1[0])));
#endif
}

NOINLINE byte OverAbundantNoInline(const Factor* f, int last_factor_index, const number value, const number sum, number sum_for_gcd_coeff);
