#pragma once

void PrimeTablesInit(num64 startPrime, num64 primeLimit, const char* stopAt);
bool IsPrime(num64 n);

enum PrimeTablesParams : num64
{
	// There are exactly 325161 primes below 10^(20/3)
	// We can use this table for factorization when p^3 <= N < 10^20
	// Set it to 325184 because it's divisible by 32
	ReciprocalsTableSize128 = 325184,

	PowersOfP_128DivisibilityData_count = 990107,

	MainPrimeTableSize = 404061114,
};

static_assert(ReciprocalsTableSize128 % 32 == 0, "ReciprocalsTableSize128 must be divisible by 32");

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

extern CACHE_ALIGNED SReciprocal privPrimeReciprocals[ReciprocalsTableSize128];
#define PrimeReciprocals ((const SReciprocal* const)(privPrimeReciprocals))

struct AmicableCandidate
{
	AmicableCandidate() {}
	AmicableCandidate::AmicableCandidate(num64 _value, num64 _sum) : value(_value), sum(_sum) {}

	num64 value;
	num64 sum;
};

struct InverseData128
{
	unsigned int shift;
	unsigned int shift_bits;
	num128 inverse;
};

// Can store primes up to 2^37
struct PrimeCompactData
{
	num64 base : 37;
	num64 offsets : 27;
};

static_assert(sizeof(PrimeCompactData) == sizeof(num64), "PrimeCompactData has invalid size");

extern byte bitOffset[PrimeTableParameters::Modulo];
extern unsigned int PrimesCompactAllocationSize;
extern PrimeCompactData* privPrimesCompact;
extern unsigned int NumPrimes;
extern std::vector<AmicableCandidate> privCandidatesData;
extern CACHE_ALIGNED unsigned char privCandidatesDataMask[5 * 7 * 11];
extern CACHE_ALIGNED std::pair<num64, num64> privPrimeInverses[ReciprocalsTableSize128];

extern CACHE_ALIGNED num128 privPrimeInverses128[ReciprocalsTableSize128];
extern CACHE_ALIGNED num128 privPowersOf2_128DivisibilityData[128];
extern InverseData128* privPowersOfP_128DivisibilityData_base;
extern CACHE_ALIGNED InverseData128* privPowersOfP_128DivisibilityData[ReciprocalsTableSize128];

extern CACHE_ALIGNED num64 privSumEstimates128[ReciprocalsTableSize128 / 16];

#define PrimesCompact ((const PrimeCompactData* const)(privPrimesCompact))
#define CandidatesData ((const std::vector<AmicableCandidate>&)(privCandidatesData))
#define CandidatesDataMask ((const unsigned char*)(privCandidatesDataMask))

#define PrimeInverses ((const std::pair<num64, num64>*)(privPrimeInverses))

#define PrimeInverses128 ((const num128*)(privPrimeInverses128))
#define PowersOf2_128DivisibilityData ((const num128*)(privPowersOf2_128DivisibilityData))
#define PowersOfP_128DivisibilityData_base ((const InverseData128*)(privPowersOfP_128DivisibilityData_base))
#define PowersOfP_128DivisibilityData ((const InverseData128* const*)(privPowersOfP_128DivisibilityData))

#define SumEstimates128 ((const num64*)(privSumEstimates128))

// Can be zero if search limit is <= 10^20
#define SumEstimates128Shift 0

template<typename T>
FORCEINLINE num64 GetNthPrime(T n)
{
	const PrimeCompactData& data = privPrimesCompact[n >> 2];
	return data.base + ((data.offsets >> ((~n & 3) * 9)) & 511);
}

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
	SumEstimatesSize2 = 8192,
	SumEstimatesSize2_GPU = SumEstimatesSize2 / 16,
	IS_NUM_ELIGIBLE_BEGIN = 16,
};

extern CACHE_ALIGNED std::pair<num64, num64> PQ[SumEstimatesSize][SumEstimatesSize2];
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
	FORCEINLINE PrimeIterator() : myIndex(0), myCurrentPrime(2) {}

	explicit NOINLINE PrimeIterator(num64 aStartNumber)
	{
		if (aStartNumber <= 2)
		{
			myIndex = 0;
			myCurrentPrime = 2;
			return;
		}

		unsigned int a = 0;
		unsigned int b = NumPrimes;
		do
		{
			const unsigned int c = (a + b) >> 1;
			if (GetNthPrime(c) >= aStartNumber)
			{
				b = c;
			}
			else
			{
				a = c + 1;
			}
		} while (a < b);

		myIndex = b;
		myCurrentPrime = GetNthPrime(b);
	}

	FORCEINLINE num64 Get() const { return myCurrentPrime; }
	FORCEINLINE num64 GetNext() const { return GetNthPrime(myIndex + 1); }

	FORCEINLINE PrimeIterator& operator++()
	{
		++myIndex;
		myCurrentPrime = GetNthPrime(myIndex);
		return *this;
	}

private:
	unsigned int myIndex;
	num64 myCurrentPrime;
};

struct Factor
{
	PrimeIterator p;
	unsigned int k;
	int index;
	num128 p_inv128;
};

template<num64 sum_coeff_max_factor>
FORCEINLINE byte OverAbundant(const Factor* f, int last_factor_index, const num128 value, const num128 sum, const num64 sum_coeff)
{
	num128 g = 1;
	num128 sum_g = 1;
	num128 sum_for_gcd = sum;

	const Factor* last_factor = f + last_factor_index;

	if (f->p.Get() == 2)
	{
		unsigned long power_of_2;
		if (LowWord(sum_for_gcd))
		{
			_BitScanForward64(&power_of_2, LowWord(sum_for_gcd));
		}
		else
		{
			_BitScanForward64(&power_of_2, HighWord(sum_for_gcd));
			power_of_2 += 64;
		}

		IF_CONSTEXPR(sum_coeff_max_factor > 1)
		{
			unsigned long power_of_2_sum_coeff;
			_BitScanForward64(&power_of_2_sum_coeff, sum_coeff);
			power_of_2 += power_of_2_sum_coeff;
		}

		const unsigned long k = f->k;
		if (power_of_2 > k)
		{
			power_of_2 = k;
		}

		if (power_of_2 < 64)
		{
			g = num64(1) << power_of_2;
			sum_g = num64(1) << power_of_2;
		}
		else
		{
			g = CombineNum128(0, num64(1) << (power_of_2 - 64));
			sum_g = CombineNum128(0, num64(1) << (power_of_2 - 64));
		}
		sum_g += sum_g - 1;
		++f;
	}

	while (f <= last_factor)
	{
		const num128 prev_sum_g = sum_g;
		for (unsigned int j = 0; j < f->k; ++j)
		{
			const num128 q = sum_for_gcd * f->p_inv128;
			if (q > sum_for_gcd)
			{
				static_assert(sum_coeff_max_factor <= 2, "The following code was commented because it's never used in OpenCL version");
				//IF_CONSTEXPR(sum_coeff_max_factor > 2)
				//{
				//	if ((f->p.Get() <= sum_coeff_max_factor) && (sum_coeff * f->p_inv <= f->q_max))
				//	{
				//		g *= f->p.Get();
				//		sum_g = sum_g * f->p.Get() + prev_sum_g;
				//	}
				//}
				break;
			}
			sum_for_gcd = q;
			g *= f->p.Get();
			sum_g = sum_g * f->p.Get() + prev_sum_g;
		}
		++f;
	}

	return (sum - value) * (sum_g - g) >= value * g;
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
