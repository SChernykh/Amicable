#pragma once

void PrimeTablesInit(bool doLargePrimes = true);
bool IsPrime(number n);

enum
{
	// There are exactly 192725 primes below 2^(64/3)
	// We can use this table for factorization when p^3 <= N < 2^64
	// Use 192768 because it's divisible by 256
	ReciprocalsTableSize = 192768,
};

// Reciprocals are calculated using algorithm published in http://www.agner.org/optimize/optimizing_assembly.pdf (section 16.9 "Integer division by a constant")
#pragma pack(push, 1)
struct SReciprocal
{
	number reciprocal;
	unsigned char increment;
	unsigned char shift;

	FORCEINLINE void Init(const number aDivisor)
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
struct AmicableCandidate
{
	AmicableCandidate() {}
	AmicableCandidate(number _value, number _sum, unsigned char _is_over_abundant_mask);

	unsigned int value;
	unsigned int sum;
	unsigned char is_over_abundant_mask;
};
#pragma pack(pop)

struct uint2
{
	unsigned int x;
	unsigned int y;
};

extern byte bitOffset[PrimeTableParameters::Modulo];
extern unsigned int PrimesCompactAllocationSize;
extern uint2* privPrimesCompact;
extern unsigned int NumPrimes;
extern std::vector<AmicableCandidate> privCandidatesData;
extern CACHE_ALIGNED unsigned char privCandidatesDataMask[5 * 7 * 11];
extern CACHE_ALIGNED std::pair<number, number> privPrimeInverses[ReciprocalsTableSize];

#define PrimesCompact ((const uint2* const)(privPrimesCompact))
#define CandidatesData ((const std::vector<AmicableCandidate>&)(privCandidatesData))
#define CandidatesDataMask ((const unsigned char*)(privCandidatesDataMask))

#define PrimeInverses ((const std::pair<number, number>*)(privPrimeInverses))

template<typename T>
FORCEINLINE number GetNthPrime(T n)
{
	uint2 data = privPrimesCompact[n >> 2];
	data.x += (data.y >> (30 - (n & 3) * 10)) & 1023;
	return (static_cast<number>(data.x) << 1) + 1;
}

FORCEINLINE number Mod385(const number n)
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
	SumEstimatesSize2 = 8192,
	SumEstimatesSize2_GPU = SumEstimatesSize2 / 16,
	IS_NUM_ELIGIBLE_BEGIN = 16,
};

extern CACHE_ALIGNED std::pair<number, number> PQ[SumEstimatesSize][SumEstimatesSize2];
extern const SumEstimateData* privSumEstimates[SumEstimatesSize];
extern CACHE_ALIGNED number privSumEstimatesBeginP[SumEstimatesSize];
extern CACHE_ALIGNED number privSumEstimatesBeginQ[SumEstimatesSize];

#define SumEstimates ((const SumEstimateData * const * const)(privSumEstimates))
#define SumEstimatesBeginP (((const number* const)(privSumEstimatesBeginP)))
#define SumEstimatesBeginQ (((const number* const)(privSumEstimatesBeginQ)))

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
	FORCEINLINE PrimeIterator() : myIndex(0), myCurrentPrime(2) {}

	explicit NOINLINE PrimeIterator(number aStartNumber)
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

	FORCEINLINE number Get() const
	{
		return myCurrentPrime;
	}

	FORCEINLINE PrimeIterator& operator++()
	{
		++myIndex;
		myCurrentPrime = GetNthPrime(myIndex);
		return *this;
	}

private:
	unsigned int myIndex;
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

template<number sum_coeff_max_factor>
FORCEINLINE byte OverAbundant(const Factor* f, int last_factor_index, const number value, const number sum, const number sum_coeff)
{
	number g = 1;
	number sum_g = 1;
	number sum_for_gcd = sum;

	const Factor* last_factor = f + last_factor_index;

	if (f->p == 2)
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
		const number prev_sum_g = sum_g;
		for (unsigned int j = 0; j < f->k; ++j)
		{
			const number q = sum_for_gcd * f->p_inv;
			if (q > f->q_max)
			{
				IF_CONSTEXPR(sum_coeff_max_factor > 2)
				{
					if ((f->p <= sum_coeff_max_factor) && (sum_coeff * f->p_inv <= f->q_max))
					{
						g *= f->p;
						sum_g = sum_g * f->p + prev_sum_g;
					}
				}
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
	return leq128(n2[0], n2[1], n1[0], n1[1]);
}

NOINLINE byte OverAbundantNoInline(const Factor* f, int last_factor_index, const number value, const number sum, number sum_for_gcd_coeff);

template<number Limit>
FORCEINLINE bool whole_branch_deficient(number value, number sum, const Factor* f)
{
	if (sum - value >= value)
	{
		return false;
	}

	number sumHi = 0;
	number value1 = value;
	number p = f->p;
	int index = f->index;
	for (;;)
	{
		++index;
		p = GetNthPrime(index);
		number h;
		value = _umul128(value, p, &h);
		if ((value >= Limit) || h)
		{
			break;
		}

		// sigma(p^k) / p^k =
		// (p^(k+1) - 1) / (p^k * (p - 1)) = 
		// (p - p^-k) / (p-1) <
		// p / (p-1)
		value1 *= p - 1;
		sum = _umul128(sum, p, &sumHi);
	}
	sub128(sum, sumHi, value1, 0, &sum, &sumHi);

	return ((sum < value1) && !sumHi);
}
