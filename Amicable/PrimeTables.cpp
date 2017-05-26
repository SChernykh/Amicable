#include "stdafx.h"
#include "PrimeTables.h"
#include <algorithm>
#include "sprp64.h"
#include "primesieve.hpp"

CACHE_ALIGNED SReciprocal privPrimeReciprocals[ReciprocalsTableSize];
CACHE_ALIGNED byte* privNextPrimeShifts = (byte*) AllocateSystemMemory(ReciprocalsTableSize, false);
const SumEstimateData* privSumEstimates[SumEstimatesSize];
CACHE_ALIGNED num64 privSumEstimatesBeginP[SumEstimatesSize];
CACHE_ALIGNED num64 privSumEstimatesBeginQ[SumEstimatesSize];
std::vector<AmicableCandidate> privCandidatesData;
CACHE_ALIGNED unsigned char privCandidatesDataMask[5 * 7 * 11];
CACHE_ALIGNED std::pair<num64, num64> privPrimeInverses[CompileTimePrimesCount];
CACHE_ALIGNED std::pair<num64, num64> privPrimeInverses2[CompileTimePrimesCount];
CACHE_ALIGNED num64 privPrimeInverses3[ReciprocalsTableSize];
CACHE_ALIGNED num64 privPrimeInverses4[ReciprocalsTableSize];

CACHE_ALIGNED std::pair<num128, num128> privPrimeInverses128[ReciprocalsTableSize128];
CACHE_ALIGNED std::pair<num128, num128> privPowersOf2_128DivisibilityData[128];
CACHE_ALIGNED InverseData128* privPowersOfP_128DivisibilityData[ReciprocalsTableSize128];
CACHE_ALIGNED num64 privSumEstimates128[ReciprocalsTableSize128 / 16];

byte* MainPrimeTable = (byte*) AllocateSystemMemory(MainPrimeTableSize, false);

byte bitOffset[PrimeTableParameters::Modulo];
num64 bitMask[PrimeTableParameters::Modulo];

struct MainPrimeTableInitializer
{
	MainPrimeTableInitializer() : nPrimes(0), prev_p(0) {}

	FORCEINLINE void operator()(num64 p)
	{
		if ((nPrimes > 0) && (nPrimes < ReciprocalsTableSize))
		{
			privPrimeReciprocals[nPrimes].Init(p);
			privNextPrimeShifts[nPrimes - 1] = static_cast<byte>((p - prev_p) / ShiftMultiplier);
		}

		prev_p = p;
		++nPrimes;

		if (p >= 11)
		{
			const num64 bit = bitOffset[p % PrimeTableParameters::Modulo];
			const num64 k = (p / PrimeTableParameters::Modulo) * PrimeTableParameters::NumOffsets + bit;
			MainPrimeTable[k / ByteParams::Bits] |= (1 << (k % ByteParams::Bits));
		}
	}

	num64 nPrimes;
	num64 prev_p;
};

static NOINLINE num64 CalculateMainPrimeTable()
{
	const num64 upperBound = ((SearchLimit::MainPrimeTableBound / PrimeTableParameters::Modulo) + 16) * PrimeTableParameters::Modulo;
	const size_t arraySize = static_cast<size_t>((upperBound + PrimeTableParameters::Modulo) / PrimeTableParameters::Modulo * (PrimeTableParameters::NumOffsets / ByteParams::Bits));
	if (arraySize > MainPrimeTableSize)
	{
		std::cerr << "MainPrimeTable is too small: it should be at least " << arraySize << " elements" << std::endl;
		abort();
	}
	MainPrimeTable[0] = 1;

	MainPrimeTableInitializer p;
	primesieve::PrimeSieve sieve;
	sieve.sieveTemplated(0, upperBound, p);

	return p.nPrimes;
}

bool IsPrime(num64 n)
{
	if (n >= SearchLimit::MainPrimeTableBound)
	{
		return efficient_mr64(n);
	}

	if (n <= 63)
	{
		const num64 mask = 
			(num64(1) << 2) |
			(num64(1) << 3) |
			(num64(1) << 5) |
			(num64(1) << 7) |
			(num64(1) << 11) |
			(num64(1) << 13) |
			(num64(1) << 17) |
			(num64(1) << 19) |
			(num64(1) << 23) |
			(num64(1) << 29) |
			(num64(1) << 31) |
			(num64(1) << 37) |
			(num64(1) << 41) |
			(num64(1) << 43) |
			(num64(1) << 47) |
			(num64(1) << 53) |
			(num64(1) << 59) |
			(num64(1) << 61);

		return (mask & (num64(1) << n)) != 0;
	}

	const num64 bit = bitOffset[n % PrimeTableParameters::Modulo];
	const num64 k = (n / PrimeTableParameters::Modulo) * PrimeTableParameters::NumOffsets + bit;
	if (bit >= PrimeTableParameters::NumOffsets)
		return false;

	return (MainPrimeTable[k / ByteParams::Bits] & (1 << (k % ByteParams::Bits))) != 0;
}

struct NumberAndSumOfDivisors
{
	NumberAndSumOfDivisors() : N(1), sumN(1), ratio(0.0) {}

	num128 N;
	num128 sumN;
	double ratio;
};

NOINLINE void GetSuperAbundantNumber(PrimeIterator p, const num64 maxPower, const num128 maxN, NumberAndSumOfDivisors cur, NumberAndSumOfDivisors& best)
{
	++p;

	const num128 startSum = cur.sumN;
	for (num64 k = 1; k <= maxPower; ++k)
	{
		const num128 nextN = cur.N * p.Get();
		if (nextN >= maxN)
		{
			cur.ratio = Num128ToDouble(cur.sumN - cur.N) / Num128ToDouble(cur.N);
			if (cur.ratio > best.ratio)
			{
				best = cur;
			}
			return;
		}

		cur.N = nextN;
		cur.sumN = cur.sumN * p.Get() + startSum;

		GetSuperAbundantNumber(p, k, maxN, cur, best);
	}
}

NOINLINE num64 GetMaxSumRatio(const PrimeIterator& p, const num128 limit)
{
	NumberAndSumOfDivisors cur;
	NumberAndSumOfDivisors result;
	GetSuperAbundantNumber(p, num64(-1), limit, cur, result);

	const num128 q = result.sumN / result.N;
	num128 r = result.sumN - result.N * q;

	if (q > 1)
	{
		return num64(-1);
	}

	if (HighWord(r) == 0)
	{
		return LowWord(CombineNum128(0, LowWord(r)) / result.N) + 1;
	}
	else
	{
		unsigned long index;
		_BitScanReverse64(&index, HighWord(r));
		if (index < 12)
		{
			r <<= 63 - index;
			return (LowWord(r / result.N) + 1) << (index + 1);
		}
		else
		{
			// If index >= 12, we get less than 52 bits of precision from integer division
			// So it makes more sense to use floating point instead of integer division
			return static_cast<num64>(Num128ToDouble(r) / Num128ToDouble(result.N) * 18446744073709551616.0) + 1;
		}
	}
}

AmicableCandidate::AmicableCandidate(num64 _value, num64 _sum, unsigned char _is_over_abundant_mask)
	: value(_value)
	, sum(_sum - _value * 2)
	, is_over_abundant_mask(_is_over_abundant_mask)
{
}

static num64 g_MaxPrime;
static num64 g_LargestCandidate;

NOINLINE void SearchCandidates(Factor* factors, const num64 value, const num64 sum, int depth)
{
	if (sum - value >= value)
	{
		unsigned char is_over_abundant_mask = 0;
		is_over_abundant_mask |= OverAbundant<5>(factors, depth - 1, value, sum, 2 * 5) << 1;
		is_over_abundant_mask |= OverAbundant<7>(factors, depth - 1, value, sum, 2 * 7) << 2;
		is_over_abundant_mask |= OverAbundant<11>(factors, depth - 1, value, sum, 2 * 11) << 4;
		is_over_abundant_mask |= (((is_over_abundant_mask & 0x06) || OverAbundant<7>(factors, depth - 1, value, sum, 2 * 5 * 7)) ? byte(1) : byte(0)) << 3;
		is_over_abundant_mask |= (((is_over_abundant_mask & 0x12) || OverAbundant<11>(factors, depth - 1, value, sum, 2 * 5 * 11)) ? byte(1) : byte(0)) << 5;
		is_over_abundant_mask |= (((is_over_abundant_mask & 0x14) || OverAbundant<11>(factors, depth - 1, value, sum, 2 * 7 * 11)) ? byte(1) : byte(0)) << 6;
		is_over_abundant_mask |= (((is_over_abundant_mask & 0x7E) || OverAbundant<11>(factors, depth - 1, value, sum, 2 * 5 * 7 * 11)) ? byte(1) : byte(0)) << 7;

		privCandidatesData.emplace_back(value, sum, is_over_abundant_mask);
	}

	int start_i = (depth == 0) ? 0 : (factors[depth - 1].index + 1);

	Factor& f = factors[depth];
	if (depth > 0)
	{
		f.p = factors[depth - 1].p;
		++f.p;
	}
	else
	{
		f.p = PrimeIterator();
	}

	// A check to ensure that m is not divisible by 6
	if (depth == 1)
	{
		// factors[0].p is 2
		// factors[1].p is 3
		// change factors[1].p to 5
		if (start_i == 1)
		{
			++f.p;
			start_i = 2;
		}
	}

	// Check only 2, 3, 5 as the smallest prime factor because the smallest abundant num64 coprime to 2*3*5 is ~2*10^25
	const unsigned int max_prime = static_cast<unsigned int>((depth > 0) ? (g_MaxPrime + 1) : 7);
	for (f.index = start_i; f.p.Get() < max_prime; ++f.index, ++f.p)
	{
		num64 h;
		num64 next_value = _umul128(value, f.p.Get(), &h);
		if ((next_value > g_LargestCandidate) || h)
		{
			return;
		}
		num64 next_sum = sum * (f.p.Get() + 1);

		f.k = 1;
		PRAGMA_WARNING(suppress : 4146)
		f.p_inv = -modular_inverse64(f.p.Get());
		f.q_max = num64(-1) / f.p.Get();
		f.p_inv128 = -modular_inverse128(f.p.Get());
		f.q_max128 = NUM128_MAX / f.p.Get();

		for (;;)
		{
			if (next_sum - next_value >= next_value)
			{
				if (OverAbundant<2>(factors, depth, next_value, next_sum, 2))
				{
					goto next;
				}
			}
			else if (whole_branch_deficient(g_LargestCandidate, next_value, next_sum, &f))
			{
				goto next;
			}

			SearchCandidates(factors, next_value, next_sum, depth + 1);
			next:

			next_value = _umul128(next_value, f.p.Get(), &h);
			if ((next_value > g_LargestCandidate) || h)
			{
				break;
			}
			next_sum = next_sum * f.p.Get() + sum;
			++f.k;
		}
	}
}

NOINLINE void GenerateCandidates()
{
	privCandidatesData.reserve(std::min<num64>(1591361426, g_LargestCandidate / 30));

	Factor factors[MaxPrimeFactors];
	SearchCandidates(factors, 1, 1, 0);

	std::sort(privCandidatesData.begin(), privCandidatesData.end(), [](const AmicableCandidate& a, const AmicableCandidate& b){ return a.value < b.value; });
}

void PrimeTablesInit(num64 startPrime, num64 primeLimit, const char* stopAt)
{
	// Make sure all floating point calculations round up
	ForceRoundUpFloatingPoint();

	memset(bitOffset, -1, sizeof(bitOffset));
	for (byte b = 0; b < PrimeTableParameters::NumOffsets; ++b)
	{
		bitOffset[NumbersCoprimeToModulo[b]] = b;
		bitMask[NumbersCoprimeToModulo[b]] = ~(1ULL << b);
	}

	for (unsigned int i = 0; i < 385; ++i)
	{
		unsigned int index = 0;
		if (i * MultiplicativeInverse<5>::value <= num64(-1) / 5) index += 1;
		if (i * MultiplicativeInverse<7>::value <= num64(-1) / 7) index += 2;
		if (i * MultiplicativeInverse<11>::value <= num64(-1) / 11) index += 4;
		privCandidatesDataMask[i] = static_cast<byte>(1 << index);
	}

	CalculateMainPrimeTable();

	PrimeIterator it(3);
	for (num64 index = 1; index < ARRAYSIZE(privPrimeInverses3); ++it, ++index)
	{
		const num64 p = it.Get();
		const num64 p_max = num64(-1) / p;

		PRAGMA_WARNING(suppress : 4146)
		const num64 p_inv = -modular_inverse64(p);

		if (p_inv * p != 1)
		{
			std::cerr << "modular_inverse64 failed for p = " << p << std::endl;
			abort();
		}

		if (index < CompileTimePrimesCount)
		{
			privPrimeInverses[index].first = p_inv;
			privPrimeInverses[index].second = p_max;
			privPrimeInverses2[index].first = p_inv;
			privPrimeInverses2[index].second = p_max;
		}
		privPrimeInverses3[index] = p_inv;
		privPrimeInverses4[index] = p_inv;
	}

	num128 curPowerOf2 = 4;
	for (num64 i = 1; i < 128; ++i, curPowerOf2 += curPowerOf2)
	{
		const num128 value = curPowerOf2 - 1;
		privPowersOf2_128DivisibilityData[i].first = -modular_inverse128(value);
		privPowersOf2_128DivisibilityData[i].second = NUM128_MAX / value;
	}

	const size_t inverse_ptr128_size = 990059;
	InverseData128* inverse_ptr128_buf = reinterpret_cast<InverseData128*>(AllocateSystemMemory(inverse_ptr128_size * sizeof(InverseData128), false));
	InverseData128* inverse_ptr128 = inverse_ptr128_buf;
	InverseData128* inverse_ptr128_end = inverse_ptr128_buf + inverse_ptr128_size;

	it = PrimeIterator(3);
	for (num64 index = 0; index < ReciprocalsTableSize128; ++it, ++index)
	{
		const num64 p = it.Get();
		const num128 p_inv = -modular_inverse128(p);
		if (p_inv * p != 1)
		{
			std::cerr << "modular_inverse128 failed for p = " << p << std::endl;
			abort();
		}
		privPrimeInverses128[index].first = p_inv;
		privPrimeInverses128[index].second = NUM128_MAX / p;

		privPowersOfP_128DivisibilityData[index] = inverse_ptr128;
		num128 p1 = p;
		num128 sum1 = p + 1;
		const num128 sum_limit = SearchLimit::value * 3;
		do
		{
			if (LowWord(sum1) == 0)
			{
				std::cerr << "shift is too large for p = " << p << std::endl;
				abort();
			}

			if (inverse_ptr128 < inverse_ptr128_end)
			{
				unsigned long k;
				_BitScanForward64(&k, LowWord(sum1));
				inverse_ptr128->shift = k;
				inverse_ptr128->shift_bits = (1U << k) - 1;

				const num128 value = sum1 >> k;
				inverse_ptr128->inverse = -modular_inverse128(value);
				inverse_ptr128->max_value = NUM128_MAX / value;
				if (inverse_ptr128->inverse * value != 1)
				{
					std::cerr << "modular_inverse128 failed for value = " << value << std::endl;
					abort();
				}
			}

			p1 *= p;
			sum1 += p1;
			++inverse_ptr128;
		} while (sum1 <= sum_limit);

		if (((index + 1) % 16) == 0)
		{
			const num64 r = (GetMaxSumRatio(it, SearchLimit::value) + ((num64(1) << SumEstimates128Shift) - 1)) >> SumEstimates128Shift;

			// Check if "SearchLimit::value * r" overflows
			if (((SearchLimit::value * r) % r) != 0)
			{
				std::cerr << "Increase SumEstimates128Shift" << std::endl;
				abort();
			}

			privSumEstimates128[index / 16] = r;
		}
	}

	{
		const num64 p = it.Get();
		if ((num128(p) * p) * p < SearchLimit::value)
		{
			std::cerr << "Increase ReciprocalsTableSize128: it must be at least the number of primes below cube root of SearchLimit::value" << std::endl;
			abort();
		}
	}

	if (inverse_ptr128 > inverse_ptr128_end)
	{
		std::cerr << "inverse_ptr128_size must be at least " << (inverse_ptr128 - inverse_ptr128_buf) << std::endl;
		abort();
	}

	// Gather data for linear search and do preliminary filtering
	// All filters combined leave 971348 numbers out of the first 1000000
	// It's a great speed-up compared to the recursive search
	if ((startPrime && primeLimit) || !stopAt)
	{
		g_LargestCandidate = LowWord(SearchLimit::value / std::max<num64>(SearchLimit::LinearLimit, startPrime));
		g_MaxPrime = g_LargestCandidate / 4;
		GenerateCandidates();
	}

	// PQ corresponds to tables P and Q in lemma 2.1 from
	// "Computation of All the Amicable Pairs Below 10^10 By H.J.J.te Riele": http://www.ams.org/journals/mcom/1986-47-175/S0025-5718-1986-0842142-3/S0025-5718-1986-0842142-3.pdf
	// The only difference is that we calculate exact (hence better) upper bounds for S(m)/m instead of inexact estimates
	const num64 maxI = 16384;
	std::vector<std::pair<num64, num64>> PQ[SumEstimatesSize];
	for (num64 j = 0; j < SumEstimatesSize; ++j)
	{
		PQ[j].resize(maxI);
	}

	PrimeIterator prevP(1);
	PrimeIterator p(2);
	num64 PQ_size = maxI;
	auto MultiplyWithSaturation = [](const num64 a, const num64 b)
	{
		num64 h;
		const num64 result = _umul128(a, b, &h);
		return h ? num64(-1) : result;
	};
	for (num64 i = 0; (i < maxI) && (p.Get() <= std::max<num64>(SearchLimit::MainPrimeTableBound, CompileTimePrimes<CompileTimePrimesCount>::value)); ++i, ++p)
	{
		num64 j = 1;
		PrimeIterator q(p);
		++q;

		PQ[0][i].first = p.Get();
		PQ[0][i].second = GetMaxSumRatio(prevP, MultiplyWithSaturation(p.Get(), q.Get()));

		for (; (j < SumEstimatesSize) && (q.Get() <= std::max<num64>(SearchLimit::MainPrimeTableBound, CompileTimePrimes<CompileTimePrimesCount>::value)); ++j)
		{
			num64 highProductP;
			const num64 mulP = _umul128(PQ[j - 1][i].first, q.Get(), &highProductP);
			++q;
			if (highProductP)
			{
				const std::pair<num64, num64> k(num64(-1), PQ[j - 1][i].second);
				for (; j < SumEstimatesSize; ++j)
					PQ[j][i] = k;
				break;
			}
			else
			{
				PQ[j][i] = std::pair<num64, num64>(mulP, GetMaxSumRatio(prevP, MultiplyWithSaturation(mulP, q.Get())));
			}
		}
		if (p.Get() > 65536)
		{
			PQ_size = i;
			break;
		}
		prevP = p;
	}

	for (num64 i = 0; i < PQ_size; ++i)
	{
		for (num64 j = 0; j < SumEstimatesSize; ++j)
		{
			if (PQ[j][i].first != num64(-1))
			{
				--PQ[j][i].first;
			}
		}
	}

	SumEstimateData* data = new SumEstimateData[SumEstimatesSize * ((PQ_size - IS_NUM_ELIGIBLE_BEGIN + 7) / 8)];
	for (num64 j = 0; j < SumEstimatesSize; ++j)
	{
		privSumEstimates[j] = data - (IS_NUM_ELIGIBLE_BEGIN / 8);
		for (num64 i = IS_NUM_ELIGIBLE_BEGIN; i < PQ_size; i += 8)
		{
			data->P = PQ[j][i].first;
			data->Q = PQ[j][i].second;
			++data;
		}
	}
	for (num64 j = 0; j < SumEstimatesSize; ++j)
	{
		privSumEstimatesBeginP[j] = (j + 1 < SumEstimatesSize) ? privSumEstimates[j + 1][IS_NUM_ELIGIBLE_BEGIN / 8].P : num64(-1);
		privSumEstimatesBeginQ[j] = privSumEstimates[j][IS_NUM_ELIGIBLE_BEGIN / 8].Q;
	}
}
