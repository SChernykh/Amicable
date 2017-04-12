#include "stdafx.h"
#include "PrimeTables.h"
#include <algorithm>
#include "sprp64.h"
#include "primesieve.hpp"

CACHE_ALIGNED SReciprocal privPrimeReciprocals[ReciprocalsTableSize];
byte* privNextPrimeShifts = nullptr;
const SumEstimateData* privSumEstimates[SumEstimatesSize];
CACHE_ALIGNED num64 privSumEstimatesBeginP[SumEstimatesSize];
CACHE_ALIGNED num64 privSumEstimatesBeginQ[SumEstimatesSize];
std::vector<AmicableCandidate> privCandidatesData;
CACHE_ALIGNED unsigned char privCandidatesDataMask[5 * 7 * 11];
CACHE_ALIGNED std::pair<num64, num64> privPrimeInverses[CompileTimePrimesCount];
CACHE_ALIGNED std::pair<num64, num64> privPrimeInverses2[CompileTimePrimesCount];
CACHE_ALIGNED num64 privPrimeInverses3[ReciprocalsTableSize];
CACHE_ALIGNED num64 privPrimeInverses4[ReciprocalsTableSize];

byte* MainPrimeTable = nullptr;
byte bitOffset[PrimeTableParameters::Modulo];
num64 bitMask[PrimeTableParameters::Modulo];

struct MainPrimeTableInitializer
{
	MainPrimeTableInitializer() : nPrimes(0), prev_p(0) {}

	FORCEINLINE void operator()(num64 p)
	{
		if (nPrimes > 0)
		{
			if (nPrimes < ReciprocalsTableSize)
			{
				privPrimeReciprocals[nPrimes].Init(p);
			}

			privNextPrimeShifts[nPrimes * 2 - 2] = static_cast<byte>((p - prev_p) / ShiftMultiplier);
			privNextPrimeShifts[nPrimes * 2 - 1] = privCandidatesDataMask[Mod385(prev_p + 1)];
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

static num64 CalculateMainPrimeTable()
{
	// https://en.wikipedia.org/wiki/Prime_gap#Numerical_results
	// Since we operate in the range 1..2^64, a gap = PrimeTableParameters::Modulo * 9 = 1890 is enough
	const num64 upperBound = ((SearchLimit::MainPrimeTableBound / PrimeTableParameters::Modulo) + 10) * PrimeTableParameters::Modulo;
	const size_t arraySize = static_cast<size_t>((upperBound + PrimeTableParameters::Modulo) / PrimeTableParameters::Modulo * (PrimeTableParameters::NumOffsets / ByteParams::Bits));
	MainPrimeTable = reinterpret_cast<byte*>(AllocateSystemMemory(arraySize, false));
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
	NumberAndSumOfDivisors() : N(1), sumLow(1), sumHigh(0), q(1), r(0) {}

	num64 N;
	num64 sumLow;
	num64 sumHigh;
	num64 q;
	num64 r;
};

NOINLINE void GetSuperAbundantNumber(PrimeIterator p, const num64 maxPower, const num64 maxN, NumberAndSumOfDivisors cur, NumberAndSumOfDivisors& best)
{
	++p;

	num64 curSumLow = cur.sumLow;
	num64 curSumHigh = cur.sumHigh;
	num64 h;
	for (num64 k = 1; k <= maxPower; ++k)
	{
		const num64 nextN = _umul128(cur.N, p.Get(), &h);
		if (h || (nextN >= maxN))
		{
			// cur.sum / cur.N > best.sum / best.N
			cur.q = udiv128(cur.sumHigh, cur.sumLow, cur.N, &cur.r);

			if (cur.q > best.q)
			{
				best = cur;
			}
			else if (cur.q == best.q)
			{
				// cur.r / cur.N > best.r / best.N
				// cur.r * best.N > best.r * cur.N
				num64 a[2];
				a[0] = _umul128(cur.r, best.N, &a[1]);

				num64 b[2];
				b[0] = _umul128(best.r, cur.N, &b[1]);

				if ((a[1] > b[1]) || ((a[1] == b[1]) && (a[0] > b[0])))
					best = cur;
			}
			return;
		}

		cur.N = nextN;

		curSumLow = _umul128(curSumLow, p.Get(), &h);
		curSumHigh = curSumHigh * p.Get() + h;

		num64 carry = 0;
		if (cur.sumLow > ~curSumLow)
			carry = 1;
		cur.sumLow += curSumLow;
		cur.sumHigh += curSumHigh + carry;

		GetSuperAbundantNumber(p, k, maxN, cur, best);
	}
}

static std::vector<std::pair<num64, num64>> locTmpFactorization;

NOINLINE num64 GetMaxSumRatio(const PrimeIterator& p, const num64 limit, num64* numberWihMaxSumRatio = nullptr)
{
	NumberAndSumOfDivisors cur;
	NumberAndSumOfDivisors result;
	GetSuperAbundantNumber(p, num64(-1), limit, cur, result);

	num64 r;
	num64 q = udiv128(result.sumHigh, result.sumLow, result.N, &r);
	if (numberWihMaxSumRatio)
		*numberWihMaxSumRatio = result.N;
	if (q > 1)
		return num64(-1);

	return udiv128(r, 0, result.N, &r) + 1;
}

AmicableCandidate::AmicableCandidate(num64 _value, num64 _sum, unsigned char _is_over_abundant_mask)
	: value(_value)
	, sum(_sum - _value * 2)
	, is_over_abundant_mask(_is_over_abundant_mask)
{
}

#pragma pack(push, 1)
struct PrimeData
{
	PrimeData(unsigned int _p, num64 _p_inv, num64 _q_max) : p(_p), p_inv(_p_inv), q_max(_q_max) {}

	unsigned int p;
	num64 p_inv;
	num64 q_max;
};
#pragma pack(pop)

static std::vector<PrimeData> g_PrimeData;
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
	f.p = (depth == 0) ? 2U : g_PrimeData[static_cast<unsigned int>(factors[depth - 1].index + 1)].p;

	// A check to ensure that m is not divisible by 6
	if (depth == 1)
	{
		// factors[0].p is 2
		// factors[1].p is 3
		// change factors[1].p to 5
		if (start_i == 1)
		{
			f.p = 5;
		}
		if (start_i == 1)
		{
			start_i = 2;
		}
	}

	// Check only 2, 3, 5 as the smallest prime factor because the smallest abundant num64 coprime to 2*3*5 is ~2*10^25
	const unsigned int max_prime = static_cast<unsigned int>((depth > 0) ? (g_MaxPrime + 1) : 7);
	for (f.index = start_i; f.p < max_prime; ++f.index, f.p = (static_cast<unsigned int>(f.index) < g_PrimeData.size()) ? g_PrimeData[static_cast<unsigned int>(f.index)].p : max_prime)
	{
		num64 h;
		num64 next_value = _umul128(value, f.p, &h);
		if ((next_value > g_LargestCandidate) || h)
		{
			return;
		}
		num64 next_sum = sum * (f.p + 1);

		f.k = 1;
		f.p_inv = g_PrimeData[static_cast<unsigned int>(f.index)].p_inv;
		f.q_max = g_PrimeData[static_cast<unsigned int>(f.index)].q_max;

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

			next_value = _umul128(next_value, f.p, &h);
			if ((next_value > g_LargestCandidate) || h)
			{
				break;
			}
			next_sum = next_sum * f.p + sum;
			++f.k;
		}
	}
}

NOINLINE void GenerateCandidates()
{
	privCandidatesData.reserve(std::min<num64>(178832709, g_LargestCandidate / 30));
	{
		const num64 primeDataCount = 87348706;
		g_PrimeData.reserve(primeDataCount);
		g_PrimeData.emplace_back(2, 0, 0);
		for (num64 p = 3, index = 1; index < primeDataCount; p += NextPrimeShifts[index * 2] * ShiftMultiplier, ++index)
		{
			PRAGMA_WARNING(suppress : 4146)
			g_PrimeData.emplace_back(static_cast<unsigned int>(p), -modular_inverse64(p), num64(-1) / p);
			if (p > g_MaxPrime)
			{
				break;
			}
		}

		Factor factors[16];
		SearchCandidates(factors, 1, 1, 0);

		std::vector<PrimeData> tmp;
		g_PrimeData.swap(tmp);
	}
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

	const double nPrimesBound = static_cast<double>((static_cast<num64>(SearchLimit::MainPrimeTableBound) < 100000) ? 100000 : static_cast<num64>(SearchLimit::MainPrimeTableBound));
	const num64 nPrimesEstimate = static_cast<num64>(nPrimesBound / (log(nPrimesBound) - 1.1));
	privNextPrimeShifts = reinterpret_cast<byte*>(AllocateSystemMemory((nPrimesEstimate + (nPrimesEstimate & 1)) * 2, false));

	for (unsigned int i = 0; i < 385; ++i)
	{
		unsigned int index = 0;
		if (i * MultiplicativeInverse<5>::value <= num64(-1) / 5) index += 1;
		if (i * MultiplicativeInverse<7>::value <= num64(-1) / 7) index += 2;
		if (i * MultiplicativeInverse<11>::value <= num64(-1) / 11) index += 4;
		privCandidatesDataMask[i] = static_cast<byte>(1 << index);
	}

	CalculateMainPrimeTable();

	for (num64 p = 3, index = 1; index < ARRAYSIZE(privPrimeInverses3); p += NextPrimeShifts[index * 2] * ShiftMultiplier, ++index)
	{
		const num64 p_max = num64(-1) / p;

		PRAGMA_WARNING(suppress : 4146)
		const num64 p_inv = -modular_inverse64(p);

		if (p_inv * p != 1)
		{
			std::cerr << "modular_inverse64 failed";
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

	// Gather data for linear search and do preliminary filtering
	// All filters combined leave 971348 numbers out of the first 1000000
	// It's a great speed-up compared to the recursive search
	if ((startPrime && primeLimit) || !stopAt)
	{
		g_LargestCandidate = (SearchLimit::value / std::max<num64>(SearchLimit::LinearLimit, startPrime)).lo;
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
