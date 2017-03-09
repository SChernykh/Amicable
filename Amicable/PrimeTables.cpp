#include "stdafx.h"
#include "PrimeTables.h"
#include <algorithm>
#include "sprp64.h"
#include "primesieve.hpp"

CACHE_ALIGNED SReciprocal privPrimeReciprocals[ReciprocalsTableSize];
byte* privNextPrimeShifts = nullptr;
const SumEstimateData* privSumEstimates[SumEstimatesSize];
CACHE_ALIGNED number privSumEstimatesBeginP[SumEstimatesSize];
CACHE_ALIGNED number privSumEstimatesBeginQ[SumEstimatesSize];
std::vector<AmicableCandidate> privCandidatesData;
CACHE_ALIGNED unsigned char privCandidatesDataMask[5 * 7 * 11];
CACHE_ALIGNED std::pair<number, number> privPrimeInverses[CompileTimePrimesCount];
CACHE_ALIGNED std::pair<number, number> privPrimeInverses2[CompileTimePrimesCount];

byte* MainPrimeTable = nullptr;
byte bitOffset[PrimeTableParameters::Modulo];
number bitMask[PrimeTableParameters::Modulo];

struct MainPrimeTableInitializer
{
	MainPrimeTableInitializer() : nPrimes(0), prev_p(0) {}

	FORCEINLINE void operator()(number p)
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
			const number bit = bitOffset[p % PrimeTableParameters::Modulo];
			const number k = (p / PrimeTableParameters::Modulo) * PrimeTableParameters::NumOffsets + bit;
			MainPrimeTable[k / Byte::Bits] |= (1 << (k % Byte::Bits));
		}
	}

	number nPrimes;
	number prev_p;
};

static number CalculateMainPrimeTable()
{
	// https://en.wikipedia.org/wiki/Prime_gap#Numerical_results
	// Since we operate in the range 1..2^64, a gap = PrimeTableParameters::Modulo * 9 = 1890 is enough
	const number upperBound = ((SearchLimit::MainPrimeTableBound / PrimeTableParameters::Modulo) + 10) * PrimeTableParameters::Modulo;
	const size_t arraySize = static_cast<size_t>((upperBound + PrimeTableParameters::Modulo) / PrimeTableParameters::Modulo * (PrimeTableParameters::NumOffsets / Byte::Bits));
	MainPrimeTable = reinterpret_cast<byte*>(AllocateSystemMemory(arraySize, false));
	MainPrimeTable[0] = 1;

	MainPrimeTableInitializer p;
	primesieve::PrimeSieve sieve;
	sieve.sieveTemplated(0, upperBound, p);

	return p.nPrimes;
}

bool IsPrime(number n)
{
	if (n >= SearchLimit::MainPrimeTableBound)
	{
		return efficient_mr64(n);
	}

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

	const number bit = bitOffset[n % PrimeTableParameters::Modulo];
	const number k = (n / PrimeTableParameters::Modulo) * PrimeTableParameters::NumOffsets + bit;
	if (bit >= PrimeTableParameters::NumOffsets)
		return false;

	return (MainPrimeTable[k / Byte::Bits] & (1 << (k % Byte::Bits))) != 0;
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

static std::vector<std::pair<number, number>> locTmpFactorization;

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

AmicableCandidate::AmicableCandidate(number _value, number _sum, unsigned char _is_over_abundant_mask)
	: value(static_cast<unsigned int>(_value))
	, sum(static_cast<unsigned int>(_sum - _value * 2))
	, is_over_abundant_mask(static_cast<unsigned char>(_is_over_abundant_mask))
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

NOINLINE void SearchCandidates(Factor* factors, const number value, const number sum, int depth)
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

	// Check only 2, 3, 5 as the smallest prime factor because the smallest abundant number coprime to 2*3*5 is ~2*10^25
	const unsigned int max_prime = static_cast<unsigned int>((depth > 0) ? ((SearchLimit::LinearLimit / 20) + 1) : 7);
	for (f.index = start_i; f.p < max_prime; ++f.index, f.p = g_PrimeData[static_cast<unsigned int>(f.index)].p)
	{
		number h;
		number next_value = _umul128(value, f.p, &h);
		if ((next_value >= SearchLimit::value / SearchLimit::LinearLimit) || h)
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
			else if (whole_branch_deficient<SearchLimit::value / SearchLimit::LinearLimit>(next_value, next_sum, &f))
			{
				goto next;
			}

			SearchCandidates(factors, next_value, next_sum, depth + 1);
			next:

			next_value = _umul128(next_value, f.p, &h);
			if ((next_value >= SearchLimit::value / SearchLimit::LinearLimit) || h)
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
	privCandidatesData.reserve(Min<77432115, SearchLimit::value / SearchLimit::LinearLimit / 30>::value);
	{
		const number primeDataCount = 16441820;
		g_PrimeData.reserve(primeDataCount);
		g_PrimeData.emplace_back(2, 0, 0);
		for (number p = 3, index = 1; index < primeDataCount; p += NextPrimeShifts[index * 2] * ShiftMultiplier, ++index)
		{
			PRAGMA_WARNING(suppress : 4146)
			g_PrimeData.emplace_back(static_cast<unsigned int>(p), -modular_inverse64(p), number(-1) / p);
		}

		Factor factors[16];
		SearchCandidates(factors, 1, 1, 0);

		std::vector<PrimeData> tmp;
		g_PrimeData.swap(tmp);
	}
	std::sort(privCandidatesData.begin(), privCandidatesData.end(), [](const AmicableCandidate& a, const AmicableCandidate& b){ return a.value < b.value; });
}

void PrimeTablesInit(bool doLargePrimes)
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
	const number nPrimesEstimate = static_cast<number>(nPrimesBound / (log(nPrimesBound) - 1.1));
	privNextPrimeShifts = reinterpret_cast<byte*>(AllocateSystemMemory((nPrimesEstimate + (nPrimesEstimate & 1)) * 2, false));

	for (unsigned int i = 0; i < 385; ++i)
	{
		unsigned int index = 0;
		if (i * MultiplicativeInverse<5>::value <= number(-1) / 5) index += 1;
		if (i * MultiplicativeInverse<7>::value <= number(-1) / 7) index += 2;
		if (i * MultiplicativeInverse<11>::value <= number(-1) / 11) index += 4;
		privCandidatesDataMask[i] = static_cast<byte>(1 << index);
	}

	CalculateMainPrimeTable();

	for (number p = 3, index = 1; index < CompileTimePrimesCount; p += NextPrimeShifts[index * 2] * ShiftMultiplier, ++index)
	{
		const number p_max = number(-1) / p;

		PRAGMA_WARNING(suppress : 4146)
		const number p_inv = -modular_inverse64(p);

		if (p_inv * p != 1)
		{
			std::cerr << "modular_inverse64 failed";
			abort();
		}

		privPrimeInverses[index].first = p_inv;
		privPrimeInverses[index].second = p_max;
		privPrimeInverses2[index].first = p_inv;
		privPrimeInverses2[index].second = p_max;
	}

	// Gather data for linear search and do preliminary filtering
	// All filters combined leave 971348 numbers out of the first 1000000
	// It's a great speed-up compared to the recursive search
	if (doLargePrimes)
	{
		GenerateCandidates();
	}

	// PQ corresponds to tables P and Q in lemma 2.1 from
	// "Computation of All the Amicable Pairs Below 10^10 By H.J.J.te Riele": http://www.ams.org/journals/mcom/1986-47-175/S0025-5718-1986-0842142-3/S0025-5718-1986-0842142-3.pdf
	// The only difference is that we calculate exact (hence better) upper bounds for S(m)/m instead of inexact estimates
	const number maxI = 16384;
	std::vector<std::pair<number, number>> PQ[SumEstimatesSize];
	for (number j = 0; j < SumEstimatesSize; ++j)
	{
		PQ[j].resize(maxI);
	}

	PrimeIterator prevP(1);
	PrimeIterator p(2);
	number PQ_size = maxI;
	auto MultiplyWithSaturation = [](const number a, const number b)
	{
		number h;
		const number result = _umul128(a, b, &h);
		return h ? number(-1) : result;
	};
	for (number i = 0; (i < maxI) && (p.Get() <= std::max<number>(SearchLimit::MainPrimeTableBound, CompileTimePrimes<CompileTimePrimesCount>::value)); ++i, ++p)
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
