#include "stdafx.h"
#include "PrimeTables.h"
#include <algorithm>
#include "sprp64.h"
#include "primesieve.hpp"

CACHE_ALIGNED SReciprocal privPrimeReciprocals[ReciprocalsTableSize];
unsigned int PrimesCompactAllocationSize = 0;
uint2* privPrimesCompact = nullptr;
unsigned int NumPrimes = 0;
CACHE_ALIGNED std::pair<number, number> PQ[SumEstimatesSize][SumEstimatesSize2];
const SumEstimateData* privSumEstimates[SumEstimatesSize];
CACHE_ALIGNED number privSumEstimatesBeginP[SumEstimatesSize];
CACHE_ALIGNED number privSumEstimatesBeginQ[SumEstimatesSize];
std::vector<AmicableCandidate> privCandidatesData;
CACHE_ALIGNED unsigned char privCandidatesDataMask[5 * 7 * 11];
CACHE_ALIGNED std::pair<number, number> privPrimeInverses[ReciprocalsTableSize];

byte bitOffset[PrimeTableParameters::Modulo];
number bitMask[PrimeTableParameters::Modulo];

static FORCEINLINE void SetNthPrime(unsigned int n, number p)
{
	p >>= 1;
	if ((n & 3) == 0)
	{
		privPrimesCompact[n >> 2].x = static_cast<unsigned int>(p);
	}
	else
	{
		p -= privPrimesCompact[n >> 2].x;
		privPrimesCompact[n >> 2].y |= p << (30 - (n & 3) * 10);
	}
}

struct MainPrimeTableInitializer
{
	FORCEINLINE void operator()(number p)
	{
		SetNthPrime(NumPrimes, p);

		if (NumPrimes < ReciprocalsTableSize)
		{
			InitReciprocals(p);
		}

		++NumPrimes;
	}

	NOINLINE void InitReciprocals(number p)
	{
		if (p < 3)
		{
			return;
		}

		privPrimeReciprocals[NumPrimes].Init(p);

		const number p_max = number(-1) / p;

		PRAGMA_WARNING(suppress : 4146)
		const number p_inv = -modular_inverse64(p);

		if (p_inv * p != 1)
		{
			std::cerr << "modular_inverse64 failed";
			abort();
		}

		privPrimeInverses[NumPrimes].first = p_inv;
		privPrimeInverses[NumPrimes].second = p_max;
	}
};

static NOINLINE void CalculateMainPrimeTable()
{
	// https://en.wikipedia.org/wiki/Prime_gap#Numerical_results
	// Since we operate in the range 1..2^64, a gap = PrimeTableParameters::Modulo * 9 = 1890 is enough
	const number upperBound = ((SearchLimit::MainPrimeTableBound / PrimeTableParameters::Modulo) + 10) * PrimeTableParameters::Modulo;
	MainPrimeTableInitializer p;
	primesieve::PrimeSieve sieve;
	sieve.sieveTemplated(0, upperBound, p);
}

bool IsPrime(number n)
{
	if (n <= 63)
	{
		const number mask = 
			(number(1) << 2) |
			(number(1) << 3) |
			(number(1) << 5) |
			(number(1) << 7) |
			(number(1) << 11) |
			(number(1) << 13) |
			(number(1) << 17) |
			(number(1) << 19) |
			(number(1) << 23) |
			(number(1) << 29) |
			(number(1) << 31) |
			(number(1) << 37) |
			(number(1) << 41) |
			(number(1) << 43) |
			(number(1) << 47) |
			(number(1) << 53) |
			(number(1) << 59) |
			(number(1) << 61);

		return (mask & (number(1) << n)) != 0;
	}

	return efficient_mr64(n);
}

struct NumberAndSumOfDivisors
{
	NumberAndSumOfDivisors() : N(1), sumLow(1), sumHigh(0), q(1), r(0) {}

	number N;
	number sumLow;
	number sumHigh;
	number q;
	number r;
};

NOINLINE void GetSuperAbundantNumber(PrimeIterator p, const number maxPower, const number maxN, NumberAndSumOfDivisors cur, NumberAndSumOfDivisors& best)
{
	++p;

	number curSumLow = cur.sumLow;
	number curSumHigh = cur.sumHigh;
	number h;
	for (number k = 1; k <= maxPower; ++k)
	{
		const number nextN = _umul128(cur.N, p.Get(), &h);
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
				number a[2];
				a[0] = _umul128(cur.r, best.N, &a[1]);

				number b[2];
				b[0] = _umul128(best.r, cur.N, &b[1]);

				if ((a[1] > b[1]) || ((a[1] == b[1]) && (a[0] > b[0])))
					best = cur;
			}
			return;
		}

		cur.N = nextN;

		curSumLow = _umul128(curSumLow, p.Get(), &h);
		curSumHigh = curSumHigh * p.Get() + h;

		number carry = 0;
		if (cur.sumLow > ~curSumLow)
			carry = 1;
		cur.sumLow += curSumLow;
		cur.sumHigh += curSumHigh + carry;

		GetSuperAbundantNumber(p, k, maxN, cur, best);
	}
}

NOINLINE number GetMaxSumRatio(const PrimeIterator& p, const number limit, number* numberWihMaxSumRatio = nullptr)
{
	NumberAndSumOfDivisors cur;
	NumberAndSumOfDivisors result;
	GetSuperAbundantNumber(p, number(-1), limit, cur, result);

	number r;
	number q = udiv128(result.sumHigh, result.sumLow, result.N, &r);
	if (numberWihMaxSumRatio)
		*numberWihMaxSumRatio = result.N;
	if (q > 1)
		return number(-1);

	return udiv128(r, 0, result.N, &r) + 1;
}

AmicableCandidate::AmicableCandidate(number _value, number _sum)
	: value(static_cast<unsigned int>(_value))
	, sum(static_cast<unsigned int>(_sum - _value * 2))
{
	if (_sum - _value * 2 > UINT_MAX)
	{
		std::cerr << "sigma(" << _value << ") = " << _sum << " is too high" << std::endl;
		abort();
	}
}

#pragma pack(push, 1)
struct PrimeData
{
	PrimeData(unsigned int _p, number _p_inv, number _q_max) : p(_p), p_inv(_p_inv), q_max(_q_max) {}

	unsigned int p;
	number p_inv;
	number q_max;
};
#pragma pack(pop)

static std::vector<PrimeData> g_PrimeData;
static number g_MaxPrime = SearchLimit::LinearLimit / 4;
static number g_LargestCandidate = SearchLimit::value / SearchLimit::LinearLimit;

NOINLINE void SearchCandidates(Factor* factors, const number value, const number sum, int depth)
{
	if (sum - value >= value)
	{
		privCandidatesData.emplace_back(value, sum);
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

	// Check only 2, 3, 5 as the smallest prime factor because the smallest abundant number coprime to 2*3*5 is ~2*10^25
	const unsigned int max_prime = static_cast<unsigned int>((depth > 0) ? (g_MaxPrime + 1) : 7);
	for (f.index = start_i; f.p < max_prime; ++f.index, f.p = (static_cast<unsigned int>(f.index) < g_PrimeData.size()) ? g_PrimeData[static_cast<unsigned int>(f.index)].p : max_prime)
	{
		number h;
		number next_value = _umul128(value, f.p, &h);
		if ((next_value > g_LargestCandidate) || h)
		{
			return;
		}
		number next_sum = sum * (f.p + 1);

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

void GenerateCandidates()
{
	privCandidatesData.reserve(std::min<number>(77432320, g_LargestCandidate / 30));
	{
		const number primeDataCount = 16441820;
		g_PrimeData.reserve(primeDataCount);
		for (number index = 0; index < primeDataCount; ++index)
		{
			const number p = GetNthPrime(index);
			PRAGMA_WARNING(suppress : 4146)
			g_PrimeData.emplace_back(static_cast<unsigned int>(p), -modular_inverse64(p), number(-1) / p);
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

void PrimeTablesInit(number startPrime, number primeLimit, const char* stopAt)
{
	// Make sure all floating point calculations round up
	ForceRoundUpFloatingPoint();

	memset(bitOffset, -1, sizeof(bitOffset));
	for (byte b = 0; b < PrimeTableParameters::NumOffsets; ++b)
	{
		bitOffset[NumbersCoprimeToModulo[b]] = b;
		bitMask[NumbersCoprimeToModulo[b]] = ~(1ULL << b);
	}

	const double nPrimesBound = static_cast<double>((static_cast<number>(SearchLimit::MainPrimeTableBound) < 100000) ? 100000 : static_cast<number>(SearchLimit::MainPrimeTableBound));
	unsigned int numPrimesEstimate = static_cast<unsigned int>(nPrimesBound / (log(nPrimesBound) - 1.1));
	PrimesCompactAllocationSize = numPrimesEstimate * 2;
	PrimesCompactAllocationSize = ((PrimesCompactAllocationSize / 4096) + 1) * 4096;
	privPrimesCompact = reinterpret_cast<uint2*>(AllocateSystemMemory(PrimesCompactAllocationSize, false));

	for (unsigned int i = 0; i < 385; ++i)
	{
		unsigned int index = 0;
		if (i * MultiplicativeInverse<5>::value <= number(-1) / 5) index += 1;
		if (i * MultiplicativeInverse<7>::value <= number(-1) / 7) index += 2;
		if (i * MultiplicativeInverse<11>::value <= number(-1) / 11) index += 4;
		privCandidatesDataMask[i] = static_cast<byte>(1 << index);
	}

	CalculateMainPrimeTable();

	// Gather data for linear search and do preliminary filtering
	// All filters combined leave 971348 numbers out of the first 1000000
	// It's a great speed-up compared to the recursive search
	if ((startPrime && primeLimit) || !stopAt)
	{
		g_LargestCandidate = SearchLimit::value / std::max<number>(SearchLimit::LinearLimit, startPrime);
		g_MaxPrime = g_LargestCandidate / 4;
		GenerateCandidates();
	}

	// PQ corresponds to tables P and Q in lemma 2.1 from
	// "Computation of All the Amicable Pairs Below 10^10 By H.J.J.te Riele": http://www.ams.org/journals/mcom/1986-47-175/S0025-5718-1986-0842142-3/S0025-5718-1986-0842142-3.pdf
	// The only difference is that we calculate exact (hence better) upper bounds for S(m)/m instead of inexact estimates

	PrimeIterator prevP(1);
	PrimeIterator p(2);
	number PQ_size = SumEstimatesSize2;
	auto MultiplyWithSaturation = [](const number a, const number b)
	{
		number h;
		const number result = _umul128(a, b, &h);
		return h ? number(-1) : result;
	};
	for (number i = 0; (i < SumEstimatesSize2) && (p.Get() <= std::max<number>(SearchLimit::MainPrimeTableBound, CompileTimePrimes<CompileTimePrimesCount>::value)); ++i, ++p)
	{
		number j = 1;
		PrimeIterator q(p);
		++q;

		PQ[0][i].first = p.Get();
		PQ[0][i].second = GetMaxSumRatio(prevP, MultiplyWithSaturation(p.Get(), q.Get()));

		for (; (j < SumEstimatesSize) && (q.Get() <= std::max<number>(SearchLimit::MainPrimeTableBound, CompileTimePrimes<CompileTimePrimesCount>::value)); ++j)
		{
			number highProductP;
			const number mulP = _umul128(PQ[j - 1][i].first, q.Get(), &highProductP);
			++q;
			if (highProductP)
			{
				const std::pair<number, number> k(number(-1), PQ[j - 1][i].second);
				for (; j < SumEstimatesSize; ++j)
					PQ[j][i] = k;
				break;
			}
			else
			{
				PQ[j][i] = std::pair<number, number>(mulP, GetMaxSumRatio(prevP, MultiplyWithSaturation(mulP, q.Get())));
			}
		}
		if (p.Get() > 65536)
		{
			PQ_size = i;
			break;
		}
		prevP = p;
	}

	for (number i = 0; i < PQ_size; ++i)
	{
		for (number j = 0; j < SumEstimatesSize; ++j)
		{
			if (PQ[j][i].first != number(-1))
			{
				--PQ[j][i].first;
			}
		}
	}

	SumEstimateData* data = new SumEstimateData[SumEstimatesSize * (PQ_size - IS_NUM_ELIGIBLE_BEGIN)];
	for (number j = 0; j < SumEstimatesSize; ++j)
	{
		privSumEstimates[j] = data - IS_NUM_ELIGIBLE_BEGIN;
		for (number i = IS_NUM_ELIGIBLE_BEGIN; i < PQ_size; ++i)
		{
			data->P = PQ[j][i].first;
			data->Q = PQ[j][i].second;
			++data;
		}
	}
	for (number j = 0; j < SumEstimatesSize; ++j)
	{
		privSumEstimatesBeginP[j] = (j + 1 < SumEstimatesSize) ? privSumEstimates[j + 1][IS_NUM_ELIGIBLE_BEGIN].P : number(-1);
		privSumEstimatesBeginQ[j] = privSumEstimates[j][IS_NUM_ELIGIBLE_BEGIN].Q;
	}
}
