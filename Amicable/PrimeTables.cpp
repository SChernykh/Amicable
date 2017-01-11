#include "stdafx.h"
#include "PrimeTables.h"
#include <algorithm>
#include "sprp64.h"

CACHE_ALIGNED SReciprocal privPrimeReciprocals[ReciprocalsTableSize];
byte* privNextPrimeShifts = nullptr;
const SumEstimateData* privSumEstimates[SumEstimatesSize];
CACHE_ALIGNED number privSumEstimatesBeginP[SumEstimatesSize];
CACHE_ALIGNED number privSumEstimatesBeginQ[SumEstimatesSize];
std::vector<AmicableCandidate> privCandidatesData;
CACHE_ALIGNED unsigned char privCandidatesDataMask[5 * 7 * 11];
std::pair<number, number>* privPrimeInverses = nullptr;
CACHE_ALIGNED std::pair<number, number> privPrimeInverses2[CompileTimePrimesCount];

std::vector<byte> MainPrimeTable;
byte bitOffset[PrimeTableParameters::Modulo];
number bitMask[PrimeTableParameters::Modulo];

number CalculatePrimes(number aLowerBound, number anUpperBound, std::vector<byte>& anOutPrimes)
{
	aLowerBound -= (aLowerBound % PrimeTableParameters::Modulo);

	static THREAD_LOCAL unsigned char* locGetRemainderCode = 0;
	unsigned char* GetRemainderCode = locGetRemainderCode;
	if (!GetRemainderCode)
	{
		GetRemainderCode = locGetRemainderCode = (unsigned char*) AllocateSystemMemory(16, true);
		unsigned char code[] =
		{
			0xB8,0,0,0,0, // mov eax, 0
			0xBA,0,0,0,0, // mov edx, 0
			0xF7,0xF1,    // div ecx
			0x8B,0xC2,    // mov eax, edx
			0xC3,         // ret
		};
		memcpy(GetRemainderCode, code, sizeof(code));
	}
	*(unsigned int*)(GetRemainderCode + 1) = static_cast<unsigned int>(aLowerBound);
	*(unsigned int*)(GetRemainderCode + 6) = static_cast<unsigned int>(aLowerBound >> 32);
	number (*GetRemainder)(number) = (number(*)(number))((void*)(GetRemainderCode));

	// https://en.wikipedia.org/wiki/Prime_gap#Numerical_results
	// Since we operate in the range 1..2^64, a gap = PrimeTableParameters::Modulo * 9 = 1890 is enough
	anUpperBound = ((anUpperBound / PrimeTableParameters::Modulo) + 10) * PrimeTableParameters::Modulo;
	const size_t arraySize = static_cast<size_t>((anUpperBound - aLowerBound + PrimeTableParameters::Modulo) / PrimeTableParameters::Modulo * (PrimeTableParameters::NumOffsets / Byte::Bits));
	anOutPrimes.resize(arraySize);
	byte* primes = anOutPrimes.data();
	memset(primes, -1, arraySize);
	number d = 11;
	number sqr_d = 121;

	const number* sieveData = (const number*)(MainPrimeTable.data()) + 1;
	number moduloIndex = 0;
	number bitIndexShift = 0;
	// Precomputed first 64 bits of MainPrimeTable
	// They're not available when we're only starting to calculate MainPrimeTable
	// So I just hard-coded them here
	number curSieveChunk = 0xfafd7bbef7ffffffULL & ~number(3);
	const unsigned int* PossiblePrimesForModuloPtr = NumbersCoprimeToModulo;

	const number rangeSize = anUpperBound - aLowerBound;
	const number GetRemainderBound = aLowerBound >> 32;
	while (sqr_d <= anUpperBound)
	{
		number k = sqr_d;
		if (k < aLowerBound)
		{
			if (d * 2 > GetRemainderBound)
				break;
			number k1 = aLowerBound % (d * 2);
			number correction = 0;
			if (k1 > d)
				correction = d * 2;
			k = d + correction - k1;
		}
		else
			k -= aLowerBound;

		for (; k <= rangeSize; k += d * 2)
		{
			const number mask = bitMask[k % PrimeTableParameters::Modulo];
			if (mask)
			{
				number* chunk = reinterpret_cast<number*>(primes + (k / PrimeTableParameters::Modulo) * (PrimeTableParameters::NumOffsets / Byte::Bits));
				*chunk &= mask;
			}
		}

		while (!curSieveChunk)
		{
			curSieveChunk = *(sieveData++);

			const number NextValuesModuloIndex = (PrimeTableParameters::Modulo / 2) | (number(PrimeTableParameters::Modulo / 2) << 16) | (number(PrimeTableParameters::Modulo) << 32);
			const number NextValuesBitIndexShift = 16 | (32 << 16) | (number(0) << 32);

			moduloIndex += ((NextValuesModuloIndex >> bitIndexShift) & 255) * 2;
			bitIndexShift = (NextValuesBitIndexShift >> bitIndexShift) & 255;

			PossiblePrimesForModuloPtr = NumbersCoprimeToModulo + bitIndexShift;
		}

		unsigned long bitIndex;
		_BitScanForward64(&bitIndex, curSieveChunk);
		curSieveChunk &= (curSieveChunk - 1);

		d = moduloIndex + PossiblePrimesForModuloPtr[bitIndex];
		sqr_d = d * d;
	}

	while (sqr_d <= anUpperBound)
	{
		number k = sqr_d;
		if (k < aLowerBound)
		{
			const number k1 = GetRemainder(d * 2);
			number correction = 0;
			if (k1 > d)
				correction = d * 2;
			k = d + correction - k1;
		}
		else
			k -= aLowerBound;

		for (; k <= rangeSize; k += d * 2)
		{
			const number mask = bitMask[k % PrimeTableParameters::Modulo];
			if (mask)
			{
				number* chunk = reinterpret_cast<number*>(primes + (k / PrimeTableParameters::Modulo) * (PrimeTableParameters::NumOffsets / Byte::Bits));
				*chunk &= mask;
			}
		}

		while (!curSieveChunk)
		{
			curSieveChunk = *(sieveData++);

			const number NextValuesModuloIndex = (PrimeTableParameters::Modulo / 2) | (number(PrimeTableParameters::Modulo / 2) << 16) | (number(PrimeTableParameters::Modulo) << 32);
			const number NextValuesBitIndexShift = 16 | (32 << 16) | (number(0) << 32);

			moduloIndex += ((NextValuesModuloIndex >> bitIndexShift) & 255) * 2;
			bitIndexShift = (NextValuesBitIndexShift >> bitIndexShift) & 255;

			PossiblePrimesForModuloPtr = NumbersCoprimeToModulo + bitIndexShift;
		}

		unsigned long bitIndex;
		_BitScanForward64(&bitIndex, curSieveChunk);
		curSieveChunk &= (curSieveChunk - 1);

		d = moduloIndex + PossiblePrimesForModuloPtr[bitIndex];
		sqr_d = d * d;
	}
	return aLowerBound;
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

template<number sum_coeff_max_factor>
NOINLINE byte OverAbundantNoInline(const Factor* f, int last_factor_index, const number value, const number sum, number sum_for_gcd_coeff)
{
	return OverAbundant<sum_coeff_max_factor>(f, last_factor_index, value, sum, sum_for_gcd_coeff);
}

template<int depth>
FORCEINLINE void SearchCandidates(Factor* factors, const number value, const number sum)
{
	if ((sum - value >= value) && !OverAbundantNoInline<2>(factors, depth - 1, value, sum, 2))
	{
		unsigned char is_over_abundant_mask = 0;
		is_over_abundant_mask |= OverAbundantNoInline<5>(factors, depth - 1, value, sum, 2 * 5) << 1;
		is_over_abundant_mask |= OverAbundantNoInline<7>(factors, depth - 1, value, sum, 2 * 7) << 2;
		is_over_abundant_mask |= OverAbundantNoInline<11>(factors, depth - 1, value, sum, 2 * 11) << 4;
		is_over_abundant_mask |= (((is_over_abundant_mask & 0x06) || OverAbundantNoInline<7>(factors, depth - 1, value, sum, 2 * 5 * 7)) ? byte(1) : byte(0)) << 3;
		is_over_abundant_mask |= (((is_over_abundant_mask & 0x12) || OverAbundantNoInline<11>(factors, depth - 1, value, sum, 2 * 5 * 11)) ? byte(1) : byte(0)) << 5;
		is_over_abundant_mask |= (((is_over_abundant_mask & 0x14) || OverAbundantNoInline<11>(factors, depth - 1, value, sum, 2 * 7 * 11)) ? byte(1) : byte(0)) << 6;
		is_over_abundant_mask |= (((is_over_abundant_mask & 0x7E) || OverAbundantNoInline<11>(factors, depth - 1, value, sum, 2 * 5 * 7 * 11)) ? byte(1) : byte(0)) << 7;

		privCandidatesData.emplace_back(AmicableCandidate(static_cast<unsigned int>(value), static_cast<unsigned int>(sum), is_over_abundant_mask));
	}

	int start_i = (depth == 0) ? 0 : (factors[depth - 1].index + 1);

	Factor& f = factors[depth];
	f.p = (depth == 0) ? 2 : (factors[depth - 1].p + NextPrimeShifts[factors[depth - 1].index * 2] * ShiftMultiplier);

	// A check to ensure that m is not divisible by 6
	IF_CONSTEXPR(depth == 1)
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
	for (f.index = start_i; f.p < max_prime; f.p += NextPrimeShifts[f.index * 2] * ShiftMultiplier, ++f.index)
	{
		number h;
		number next_value = _umul128(value, f.p, &h);
		if ((next_value >= SearchLimit::value / SearchLimit::LinearLimit) || h)
		{
			return;
		}
		number next_sum = sum * (f.p + 1);

		f.k = 1;
		f.p_inv = PrimeInverses[f.index].first;
		f.q_max = PrimeInverses[f.index].second;

		for (;;)
		{
			SearchCandidates<depth + 1>(factors, next_value, next_sum);

			next_value = _umul128(next_value, f.p, &h);
			if ((next_value >= SearchLimit::value / SearchLimit::LinearLimit) || h)
			{
				break;
			}
			next_sum = next_sum * f.p + sum;
			++f.k;
		}

		// Workaround for shifting from 2 to 3 because NextPrimeShifts[0] is 0
		IF_CONSTEXPR(depth == 0)
		{
			if (f.p == 2)
			{
				f.p = 3;
			}
		}
	}
}

template<> FORCEINLINE void SearchCandidates<10>(Factor*, const number, const number) {}

void GenerateCandidates()
{
	privCandidatesData.reserve(SearchLimit::value / SearchLimit::LinearLimit / 30);

	Factor factors[16];
	SearchCandidates<0>(factors, 1, 1);

	std::sort(privCandidatesData.begin(), privCandidatesData.end(), [](const AmicableCandidate& a, const AmicableCandidate& b){ return a.value < b.value; });
}

void PrimeTablesInit(bool doSmallPrimes, bool doLargePrimes)
{
	memset(bitOffset, -1, sizeof(bitOffset));
	for (byte b = 0; b < PrimeTableParameters::NumOffsets; ++b)
	{
		bitOffset[NumbersCoprimeToModulo[b]] = b;
		bitMask[NumbersCoprimeToModulo[b]] = ~(1ULL << b);
	}

	CalculatePrimes(0, SearchLimit::MainPrimeTableBound, MainPrimeTable);

	// Make sure all floating point calculations round up
	ForceRoundUpFloatingPoint();

	number nPrimes = 0;
	number nPrimeInverses = 0;

	number curPrimeInversesBound = SearchLimit::PrimeInversesBound;
	if (!doSmallPrimes)
	{
		curPrimeInversesBound = SearchLimit::LinearLimit / 20;
	}

	if (IsPopcntAvailable())
	{
		number* data = reinterpret_cast<number*>(MainPrimeTable.data());
		const number primeInversesBoundBytes = ((curPrimeInversesBound / PrimeTableParameters::Modulo) + 1) * (PrimeTableParameters::NumOffsets / Byte::Bits) - 1;

		number* e = reinterpret_cast<number*>(MainPrimeTable.data() + primeInversesBoundBytes);
		for (; data <= e; ++data)
		{
			nPrimes += __popcnt64(*data);
		}
		nPrimeInverses = nPrimes;

		e = reinterpret_cast<number*>(MainPrimeTable.data() + MainPrimeTable.size() - 1);
		for (; data <= e; ++data)
		{
			nPrimes += __popcnt64(*data);
		}
	}
	else
	{
		PrimeIterator it(2);
		for (; it.Get() <= curPrimeInversesBound; ++it)
		{
			++nPrimes;
		}
		nPrimeInverses = nPrimes;

		for (; it.Get() <= SearchLimit::MainPrimeTableBound; ++it)
		{
			++nPrimes;
		}
	}

	if (nPrimeInverses < CompileTimePrimesCount)
	{
		nPrimeInverses = CompileTimePrimesCount;
	}

	privPrimeInverses = reinterpret_cast<std::pair<number, number>*>(AllocateSystemMemory((((nPrimeInverses + 3) / 4) * 4) * sizeof(std::pair<number, number>) + (nPrimes + (nPrimes & 1)) * 2, false));
	privNextPrimeShifts = reinterpret_cast<byte*>(privPrimeInverses + (((nPrimeInverses + 3) / 4) * 4));

	for (unsigned int i = 0; i < 385; ++i)
	{
		unsigned int index = 0;
		if (i * MultiplicativeInverse<5>::value <= number(-1) / 5) index += 1;
		if (i * MultiplicativeInverse<7>::value <= number(-1) / 7) index += 2;
		if (i * MultiplicativeInverse<11>::value <= number(-1) / 11) index += 4;
		privCandidatesDataMask[i] = static_cast<byte>(1 << index);
	}

	nPrimes = 0;
	for (PrimeIterator it(2); it.Get() <= SearchLimit::MainPrimeTableBound;)
	{
		const number p = it.Get();

		if ((p > 2) && (nPrimes < ReciprocalsTableSize))
		{
			privPrimeReciprocals[nPrimes].Init(p);
		}

		++it;
		const number q = it.Get();
		if (q - p >= 256 * ShiftMultiplier)
		{
			std::cerr << "Primes gap is 512 or greater, decrease MainPrimeTableBound";
			abort();
		}
		privNextPrimeShifts[nPrimes * 2] = static_cast<byte>((q - p) / ShiftMultiplier);
		privNextPrimeShifts[nPrimes * 2 + 1] = privCandidatesDataMask[Mod385(p + 1)];
		++nPrimes;
	}

	for (number p = 3, index = 1; p <= curPrimeInversesBound; p += privNextPrimeShifts[index * 2] * ShiftMultiplier, ++index)
	{
		PRAGMA_WARNING(suppress : 4146)
		const number p_inv = -modular_inverse64(p);
		const number p_max = number(-1) / p;
		if (p_inv * p != 1)
		{
			std::cerr << "modular_inverse64 failed";
			abort();
		}

		privPrimeInverses[index].first = p_inv;
		privPrimeInverses[index].second = p_max;
		if (index < CompileTimePrimesCount)
		{
			privPrimeInverses2[index].first = p_inv;
			privPrimeInverses2[index].second = p_max;
		}
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
