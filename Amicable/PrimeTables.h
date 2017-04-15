#pragma once

void PrimeTablesInit(num64 startPrime, num64 primeLimit, const char* stopAt);
bool IsPrime(num64 n);

enum
{
	// There are exactly 192725 primes below 2^(64/3)
	// We can use this table for factorization when p^3 <= N < 2^64
	ReciprocalsTableSize = 192725,

	// There are exactly 325161 primes below 10^(20/3)
	// We can use this table for factorization when p^3 <= N < 10^20
	ReciprocalsTableSize128 = 325161,
};

// Reciprocals are calculated using algorithm published in http://www.agner.org/optimize/optimizing_assembly.pdf (section 16.9 "Integer division by a constant")
#pragma pack(push, 1)
struct SReciprocal
{
	num64 reciprocal;
	unsigned char increment;
	unsigned char shift;

	FORCEINLINE void Init(const num64 aDivisor)
	{
		unsigned long bitIndex;
		_BitScanReverse64(&bitIndex, aDivisor);
		shift = static_cast<unsigned char>(bitIndex);

		num64 remainder;
		num64 quotient = udiv128(num64(1) << shift, 0, aDivisor, &remainder);
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

	FORCEINLINE bool Divide(const num64 n, const num64 divisor, num64& q) const
	{
		num64 highProduct;
		_umul128(n + increment, reciprocal, &highProduct);
		q = highProduct >> shift;
		return q * divisor == n;
	}

	FORCEINLINE num64 DivideNoRemainder(const num64 n) const
	{
		num64 highProduct;
		_umul128(n + increment, reciprocal, &highProduct);
		return highProduct >> shift;
	}
};
#pragma pack(pop)

extern CACHE_ALIGNED SReciprocal privPrimeReciprocals[ReciprocalsTableSize];
#define PrimeReciprocals ((const SReciprocal* const)(privPrimeReciprocals))

#pragma pack(push, 1)
struct AmicableCandidate
{
	AmicableCandidate() {}
	AmicableCandidate(num64 _value, num64 _sum, unsigned char _is_over_abundant_mask);

	num64 value;
	num64 sum : 56;
	num64 is_over_abundant_mask : 8;
};
#pragma pack(pop)

struct InverseData128
{
	unsigned int shift;
	unsigned int shift_bits;
	num128 inverse;
	num128 max_value;
};

extern CACHE_ALIGNED byte MainPrimeTable[404061114];
extern byte bitOffset[PrimeTableParameters::Modulo];
extern CACHE_ALIGNED byte privNextPrimeShifts[ReciprocalsTableSize];
extern std::vector<AmicableCandidate> privCandidatesData;
extern CACHE_ALIGNED unsigned char privCandidatesDataMask[5 * 7 * 11];
extern CACHE_ALIGNED std::pair<num64, num64> privPrimeInverses[CompileTimePrimesCount];
extern CACHE_ALIGNED std::pair<num64, num64> privPrimeInverses2[CompileTimePrimesCount];
extern CACHE_ALIGNED num64 privPrimeInverses3[ReciprocalsTableSize];
extern CACHE_ALIGNED num64 privPrimeInverses4[ReciprocalsTableSize];

extern CACHE_ALIGNED std::pair<num128, num128> privPrimeInverses128[ReciprocalsTableSize128];
extern CACHE_ALIGNED InverseData128* privPowersOfP_128DivisibilityData[ReciprocalsTableSize128];

#define NextPrimeShifts ((const byte* const)(privNextPrimeShifts))
#define CandidatesData ((const std::vector<AmicableCandidate>&)(privCandidatesData))
#define CandidatesDataMask ((const unsigned char*)(privCandidatesDataMask))

#define PrimeInverses ((const std::pair<num64, num64>*)(privPrimeInverses))
#define PrimeInverses2 ((const std::pair<num64, num64>*)(privPrimeInverses2))
#define PrimeInverses3 ((const num64*)(privPrimeInverses3))
#define PrimeInverses4 ((const num64*)(privPrimeInverses4))

#define PrimeInverses128 ((const std::pair<num128, num128>*)(privPrimeInverses128))
#define PowersOfP_128DivisibilityData ((const InverseData128* const*)(privPowersOfP_128DivisibilityData))

FORCEINLINE num64 Mod385(const num64 n)
{
	num64 highProduct;
	static_assert(ARRAYSIZE(privCandidatesDataMask) == 5 * 7 * 11, "!!! Recalculate these constants (12265886968492584971ULL and 8) if ARRAYSIZE(privLinearSearchDataIndex) changes !!!");
	_umul128(n, 12265886968492584971ULL, &highProduct);
	return n - (highProduct >> 8) * ARRAYSIZE(privCandidatesDataMask);
}

struct SumEstimateData
{
	num64 P;
	num64 Q;
};

enum
{
	SumEstimatesSize = 16,
	IS_NUM_ELIGIBLE_BEGIN = 16,
};

extern const SumEstimateData* privSumEstimates[SumEstimatesSize];
extern CACHE_ALIGNED num64 privSumEstimatesBeginP[SumEstimatesSize];
extern CACHE_ALIGNED num64 privSumEstimatesBeginQ[SumEstimatesSize];

#define SumEstimates ((const SumEstimateData * const * const)(privSumEstimates))
#define SumEstimatesBeginP (((const num64* const)(privSumEstimatesBeginP)))
#define SumEstimatesBeginQ (((const num64* const)(privSumEstimatesBeginQ)))

FORCEINLINE num64 GCD(num64 a, num64 b)
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

		const num64 a1 = a;
		const num64 b1 = b;
		a = (a1 > b1) ? b1 : a1;
		b = (a1 > b1) ? a1 : b1;

		b -= a;
	} while (b);

	return (a << shift);
}

class PrimeIterator
{
public:
	explicit FORCEINLINE PrimeIterator(const num64* aSieveData = reinterpret_cast<const num64*>(MainPrimeTable))
		: mySieveChunk(0xfafd7bbef7fffffeULL)
		, mySieveData(aSieveData)
		, myPossiblePrimesForModuloPtr(NumbersCoprimeToModulo)
		, myModuloIndex(0)
		, myBitIndexShift(0)
		, myCurrentPrime(2)
	{
	}

	explicit PrimeIterator(num64 aStartNumber, const num64* aSieveData = (const num64*)(MainPrimeTable))
		: mySieveData(aSieveData)
		, myPossiblePrimesForModuloPtr(NumbersCoprimeToModulo)
	{
		if (aStartNumber < 2)
		{
			aStartNumber = 2;
		}

		const num64 chunkIndex = ((aStartNumber / PrimeTableParameters::Modulo) * (PrimeTableParameters::NumOffsets / ByteParams::Bits)) / sizeof(num64);
		mySieveData = aSieveData + chunkIndex;
		mySieveChunk = *mySieveData;
		myModuloIndex = (chunkIndex / 3) * PrimeTableParameters::Modulo * 4 + (chunkIndex % 3) * PrimeTableParameters::Modulo;
		myBitIndexShift = (chunkIndex % 3) * 16;
		myPossiblePrimesForModuloPtr = NumbersCoprimeToModulo + myBitIndexShift;

		myCurrentPrime = static_cast<num64>(-1);
		for (;;)
		{
			operator++();
			if (myCurrentPrime >= aStartNumber)
			{
				break;
			}
		}
	}

	FORCEINLINE num64 Get() const { return myCurrentPrime; }

	FORCEINLINE num64 GetNext() const
	{
		if (myCurrentPrime < 7)
		{
			return (0x705320UL >> (myCurrentPrime * 4)) & 15;
		}

		num64 sieveChunk = mySieveChunk;
		if (sieveChunk)
		{
			unsigned long bitIndex;
			_BitScanForward64(&bitIndex, sieveChunk);
			return myModuloIndex + myPossiblePrimesForModuloPtr[bitIndex];
		}

		const num64* sieveData = mySieveData;
		num64 moduloIndex = myModuloIndex;
		num64 bitIndexShift = myBitIndexShift;
		const unsigned int* possiblePrimesForModuloPtr = myPossiblePrimesForModuloPtr;
		const num64 NewValues = (PrimeTableParameters::Modulo / 2) | (num64(PrimeTableParameters::Modulo / 2) << 16) | (num64(PrimeTableParameters::Modulo) << 32) | (16 << 8) | (32 << 24) | (num64(0) << 40);
		do
		{
			sieveChunk = *(++sieveData);
			moduloIndex += ((NewValues >> bitIndexShift) & 255) * 2;
			bitIndexShift = (NewValues >> (bitIndexShift + 8)) & 255;
			possiblePrimesForModuloPtr = NumbersCoprimeToModulo + bitIndexShift;
		} while (!sieveChunk);

		unsigned long bitIndex;
		_BitScanForward64(&bitIndex, sieveChunk);
		return moduloIndex + possiblePrimesForModuloPtr[bitIndex];
	}

	FORCEINLINE PrimeIterator& operator++()
	{
		if (myCurrentPrime < 7)
		{
			myCurrentPrime = (0x705320UL >> (myCurrentPrime * 4)) & 15;
			return *this;
		}

		if (mySieveChunk)
		{
			unsigned long bitIndex;
			_BitScanForward64(&bitIndex, mySieveChunk);
			mySieveChunk &= (mySieveChunk - 1);
			myCurrentPrime = myModuloIndex + myPossiblePrimesForModuloPtr[bitIndex];
			return *this;
		}

		const num64 NewValues = (PrimeTableParameters::Modulo / 2) | (num64(PrimeTableParameters::Modulo / 2) << 16) | (num64(PrimeTableParameters::Modulo) << 32) | (16 << 8) | (32 << 24) | (num64(0) << 40);
		do
		{
			mySieveChunk = *(++mySieveData);
			myModuloIndex += ((NewValues >> myBitIndexShift) & 255) * 2;
			myBitIndexShift = (NewValues >> (myBitIndexShift + 8)) & 255;
			myPossiblePrimesForModuloPtr = NumbersCoprimeToModulo + myBitIndexShift;
		} while (!mySieveChunk);

		unsigned long bitIndex;
		_BitScanForward64(&bitIndex, mySieveChunk);
		mySieveChunk &= (mySieveChunk - 1);
		myCurrentPrime = myModuloIndex + myPossiblePrimesForModuloPtr[bitIndex];
		return *this;
	}

private:
	num64 mySieveChunk;
	const num64* mySieveData;
	const unsigned int* myPossiblePrimesForModuloPtr;
	num64 myModuloIndex;
	num64 myBitIndexShift;
	num64 myCurrentPrime;
};

struct Factor
{
	PrimeIterator p;
	unsigned int k;
	int index;
	num64 p_inv;
	num64 q_max;
};

template<num64 sum_coeff_max_factor>
FORCEINLINE byte OverAbundant(const Factor* f, int last_factor_index, const num64 value, const num64 sum, const num64 sum_coeff)
{
	num64 g = 1;
	num64 sum_g = 1;
	num64 sum_for_gcd = sum;

	const Factor* last_factor = f + last_factor_index;

	if (f->p.Get() == 2)
	{
		DWORD power_of_2;
		_BitScanForward64(&power_of_2, sum_for_gcd);

		IF_CONSTEXPR(sum_coeff_max_factor > 1)
		{
			DWORD power_of_2_sum_coeff;
			_BitScanForward64(&power_of_2_sum_coeff, sum_coeff);
			power_of_2 += power_of_2_sum_coeff;
		}

		const DWORD k = static_cast<DWORD>(f->k);
		if (power_of_2 > k)
		{
			power_of_2 = k;
		}
		g <<= power_of_2;
		sum_g <<= power_of_2;
		sum_g = sum_g * 2 - 1;
		++f;
	}

	while (f <= last_factor)
	{
		const num64 prev_sum_g = sum_g;
		for (unsigned int j = 0; j < f->k; ++j)
		{
			const num64 q = sum_for_gcd * f->p_inv;
			if (q > f->q_max)
			{
				IF_CONSTEXPR(sum_coeff_max_factor > 2)
				{
					if ((f->p.Get() <= sum_coeff_max_factor) && (sum_coeff * f->p_inv <= f->q_max))
					{
						g *= f->p.Get();
						sum_g = sum_g * f->p.Get() + prev_sum_g;
					}
				}
				break;
			}
			sum_for_gcd = q;
			g *= f->p.Get();
			sum_g = sum_g * f->p.Get() + prev_sum_g;
		}
		++f;
	}

	num64 n1[2];
	num64 n2[2];
	n1[0] = _umul128(sum_g - g, sum - value, &n1[1]);
	n2[0] = _umul128(g, value, &n2[1]);
	return leq128(n2[0], n2[1], n1[0], n1[1]);
}

FORCEINLINE byte whole_branch_deficient(const num128& Limit, num128 value, num128 sum, const Factor* f)
{
	if (sum - value >= value)
	{
		return false;
	}

	num128 sum1 = sum;
	num128 value1 = value;
	PrimeIterator it = f->p;
	for (;;)
	{
		++it;
		const num64 p = it.Get();
		if (value1 * p >= Limit)
		{
			break;
		}

		// sigma(p^k) / p^k =
		// (p^(k+1) - 1) / (p^k * (p - 1)) = 
		// (p - p^-k) / (p-1) <
		// p / (p-1)
		value1 *= p - 1;
		sum1 *= p;
	}
	sum1 -= value1;

	return (sum1 < value1);
}
