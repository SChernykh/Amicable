#include "stdafx.h"
#include "PrimeTables.h"
#include "Factorize.h"
#include <algorithm>
#include <process.h>
#include "sprp64.h"

SReciprocal privPrimeReciprocals[ReciprocalsTableSize];
number privPrimesUpToSqrtLimitSortedCount;
byte privNextPrimeShifts[ShiftTableSize];
const SumEstimateData* privSumEstimates[SumEstimatesSize];
std::vector<std::pair<unsigned int, unsigned int>> privLinearSearchData;
number privPrimeInverses[CompileTimePrimesCount * 2];

std::vector<byte> PrimesUpToSqrtLimit;
byte bitOffset[PrimeTableParameters::Modulo];
number bitMask[PrimeTableParameters::Modulo];

number CalculatePrimes(number aLowerBound, number anUpperBound, std::vector<byte>& anOutPrimes)
{
	aLowerBound -= (aLowerBound % PrimeTableParameters::Modulo);

	static __declspec(thread) unsigned char* locGetRemainderCode = 0;
	unsigned char* GetRemainderCode = locGetRemainderCode;
	if (!GetRemainderCode)
	{
		GetRemainderCode = locGetRemainderCode = (unsigned char*) VirtualAlloc(0, 16, MEM_COMMIT, PAGE_EXECUTE_READWRITE);
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
	unsigned __int64(__fastcall *GetRemainder)(unsigned __int64) = (unsigned __int64(__fastcall*)(unsigned __int64))((void*)(GetRemainderCode));

	// https://en.wikipedia.org/wiki/Prime_gap#Numerical_results
	// Since we operate in the range 1..2^64, a gap = PrimeTableParameters::Modulo * 9 = 1890 is enough
	anUpperBound = ((anUpperBound / PrimeTableParameters::Modulo) + 10) * PrimeTableParameters::Modulo;
	const size_t arraySize = static_cast<size_t>((anUpperBound - aLowerBound + PrimeTableParameters::Modulo) / PrimeTableParameters::Modulo * (PrimeTableParameters::NumOffsets / Byte::Bits));
	anOutPrimes.resize(arraySize);
	byte* primes = anOutPrimes.data();
	memset(primes, -1, arraySize);
	number d = 11;
	number sqr_d = 121;

	const number* sieveData = (const number*)(PrimesUpToSqrtLimit.data()) + 1;
	number moduloIndex = 0;
	number bitIndexShift = 0;
	// Precomputed first 64 bits of PrimesUpToSqrtLimit
	// They're not available when we're only starting to calculate PrimesUpToSqrtLimit table
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

			const number NewValues =	(PrimeTableParameters::Modulo / 2)	| (number(PrimeTableParameters::Modulo / 2) << 16)	| (number(PrimeTableParameters::Modulo) << 32) |
										(16 << 8)							| (32 << 24)										| (number(0) << 40);

			moduloIndex += ((NewValues >> bitIndexShift) & 255) * 2;
			bitIndexShift = (NewValues >> (bitIndexShift + 8)) & 255;

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

			const number NewValues =	(PrimeTableParameters::Modulo / 2)	| (number(PrimeTableParameters::Modulo / 2) << 16)	| (number(PrimeTableParameters::Modulo) << 32) |
										(16 << 8)							| (32 << 24)										| (number(0) << 40);

			moduloIndex += ((NewValues >> bitIndexShift) & 255) * 2;
			bitIndexShift = (NewValues >> (bitIndexShift + 8)) & 255;

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
	if (n >= SearchLimit::PrimesUpToSqrtLimitValue)
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

	return (PrimesUpToSqrtLimit[k / Byte::Bits] & (1 << (k % Byte::Bits))) != 0;
}

SmallFactorNumbers g_SmallFactorNumbersData;

template<number Index>
FORCEINLINE void CalculateSmallFactorNumbersInternal(number D, number sumD, number smallestFactorInM1)
{
	enum { p = CompileTimePrimes<Index>::value };
	number curSum = sumD;
	for (;;)
	{
		CalculateSmallFactorNumbersInternal<Index + 1>(D, curSum, smallestFactorInM1);
		if (D > SearchLimit::value / p)
			break;
		IF_CONSTEXPR(p == 3)
		{
			// Smaller number in an amicable pair must not be divisible by 6
			if (D % 2 == 0)
				return;
		}
		D *= p;
		sumD *= p;
		curSum += sumD;
	}
}

template<>
FORCEINLINE void CalculateSmallFactorNumbersInternal<0>(number D, number sumD, number smallestFactorInM1)
{
	const bool doDivision = (D > 1);
	SReciprocal recD, recSumD;
	if (doDivision)
		recD.Init(D);

	const bool isSumDPowerOf2 = ((sumD & (sumD - 1)) == 0);
	unsigned long sumDShift = 0;
	if (doDivision)
	{
		if (isSumDPowerOf2)
			_BitScanForward64(&sumDShift, sumD);
		else
			recSumD.Init(sumD);
	}

	g_SmallFactorNumbersData.myNumbers.clear();
	number curSum = sumD;
	for (;;)
	{
		CalculateSmallFactorNumbersInternal<1>(D, curSum, smallestFactorInM1);
		if (D > SearchLimit::value / 2)
			break;
		D *= 2;
		sumD *= 2;
		curSum += sumD;
	}

	std::sort(g_SmallFactorNumbersData.myNumbers.begin(), g_SmallFactorNumbersData.myNumbers.end(), [](const std::pair<number, number>& a, const std::pair<number, number>& b) { return a.first < b.first; });
	if (doDivision)
	{
		for (std::pair<number, number>& n : g_SmallFactorNumbersData.myNumbers)
		{
			n.first = recD.DivideNoRemainder(n.first);
			if (isSumDPowerOf2)
				n.second >>= sumDShift;
			else
				n.second = recSumD.DivideNoRemainder(n.second);
		}
	}
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

NOINLINE void GetSuperAbundantNumber(PrimesUpToSqrtLimitIterator p, const number maxPower, const number maxN, NumberAndSumOfDivisors cur, NumberAndSumOfDivisors& best)
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

template<>
NOINLINE void CalculateSmallFactorNumbersInternal<InlinePrimesInSearch + 1>(number D, number sumD, number /*smallestFactorInM1*/)
{
	// Check for overflow when calculating S(D)
	if (sumD < D)
		__debugbreak();

#if USE_CONJECTURES
	// If a number is divisible by 2^3, then it must not be divisible by 5 or 7
	// This is not proved, but true for all known amicable pairs
	if ((D % 8) == 0)
		if (((D % 5) == 0) || ((D % 7) == 0))
			return;
#endif

	// Numbers <= SearchLimit::value are represented as N = D*M where D=2^k0 * ... * 31^k10 * M1
	// D is given, M can be anything satisfying with these restrictions:
	// 1) M is coprime to D
	// 2) All prime factors in M are > 31 and < than all prime factors in M1
	// 3) D * M <= SearchLimit::value
	//
	// If N is a smaller member of an amicable pair, then S(N) = S(D*M) = S(D)*S(M) must be > 2 * N
	// S(D) * S(M) > 2 * N
	// S(D) * S(M) > 2 * D * M
	// S(D) / D > 2 * M / S(M) for at least one M such that D * M <= SearchLimit::value
	// S(D) / D > 2 * min(M / S(M))
	// Since all FP calculations round up here, we must transform it to this: f(...) > C
	// where f(...) are a series of calculations which will be rounded up
	// Thus we neutralize rounding errors
	// Calculating minimum is not correct when FP calculations round up
	// So we also change minimum to maximum
	//
	// S(D) > D * 2 * min(M / S(M))
	// S(D) > D * 2 / max(S(M) / M)
	// S(D) * max(S(M) / M) > D * 2
	// S(D) * max(S(M) / M) / 2 > D
	bool is_abundant = sumD - D > D;
	if (!is_abundant)
	{
		PrimesUpToSqrtLimitIterator p(CompileTimePrimes<InlinePrimesInSearch>::value);
		NumberAndSumOfDivisors cur;
		NumberAndSumOfDivisors result;
		GetSuperAbundantNumber(p, number(-1), (SearchLimit::value / D) + 1, cur, result);
		const number D1 = D * result.N;
		const number sumD1 = sumD * result.sumLow;
		is_abundant = sumD1 - D1 > D1;
	}
	if (is_abundant)
	{
		// The partial factorization must be deficient in order for the bigger number to be deficient too
		const unsigned int gcd = static_cast<unsigned int>(GCD(D, sumD));
		number sumGCD;
		locTmpFactorization.clear();
		Factorize(gcd, sumGCD, locTmpFactorization);
		if (sumGCD - gcd < gcd)
		{
			// It must be not just deficient, but "deficient enough"
			// "sumGCD * (sumD - D) < gcd * sumD" must be true
			// See comments for privLinearSearchData below
			number sum1[2];
			sum1[0] = _umul128(sumGCD, sumD - D, &sum1[1]);

			number sum2[2];
			sum2[0] = _umul128(gcd, sumD, &sum2[1]);

			if ((sum1[1] < sum2[1]) || ((sum1[1] == sum2[1]) && (sum1[0] < sum2[0])))
				g_SmallFactorNumbersData.myNumbers.push_back(std::make_pair(D, sumD));
		}
	}
}

NOINLINE number GetMaxSumRatio(const PrimesUpToSqrtLimitIterator& p, const number limit, number* numberWihMaxSumRatio = nullptr)
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

number g_MaxSumRatios[64];

void PrimeTablesInit(bool isSubmit)
{
	memset(bitOffset, -1, sizeof(bitOffset));
	for (byte b = 0; b < PrimeTableParameters::NumOffsets; ++b)
	{
		bitOffset[NumbersCoprimeToModulo[b]] = b;
		bitMask[NumbersCoprimeToModulo[b]] = ~(1ULL << b);
	}

	CalculatePrimes(0, PrimesUpToSqrtLimitValue, PrimesUpToSqrtLimit);

	// Make sure all floating point calculations round up
	_control87(_RC_UP, _MCW_RC);

	{
		PrimesUpToSqrtLimitIterator p(CompileTimePrimes<InlinePrimesInSearch>::value);
		number limit = SearchLimit::value;
		while (limit > CompileTimePrimes<InlinePrimesInSearch + 1>::value)
		{
			NumberAndSumOfDivisors cur;
			NumberAndSumOfDivisors result;
			GetSuperAbundantNumber(p, number(-1), limit, cur, result);
			const number D = SearchLimit::value / result.N;
			long index;
			_BitScanReverse64(reinterpret_cast<unsigned long*>(&index), D);
			number r;
			const number ratio = udiv128(result.sumLow, 0, result.N * 2, &r) + 1;
			while ((index >= 0) && !g_MaxSumRatios[index])
			{
				g_MaxSumRatios[index] = ratio;
				--index;
			}
			limit = result.N;
		}
		for (long index = ARRAYSIZE(g_MaxSumRatios) - 1; (index >= 0) && !g_MaxSumRatios[index]; --index)
		{
			g_MaxSumRatios[index] = number(1) << 63;
		}
	}

	number nPrimes = 0;
	for (PrimesUpToSqrtLimitIterator it(2); it.Get() <= PrimesUpToSqrtLimitValue;)
	{
		const number p = it.Get();

		if ((p > 2) && (nPrimes < ReciprocalsTableSize))
			privPrimeReciprocals[nPrimes].Init(p);

		++it;
		const number q = it.Get();
		if (q - p >= 256 * SearchLimit::ShiftMultiplier)
		{
			__debugbreak();
		}
		if (nPrimes < ARRAYSIZE(privNextPrimeShifts))
			privNextPrimeShifts[nPrimes] = static_cast<byte>((q - p) / SearchLimit::ShiftMultiplier);
		++nPrimes;
	}
	if (nPrimes > ARRAYSIZE(privNextPrimeShifts))
	{
		std::cerr << "ShiftTableSize must be >= " << nPrimes;
		__debugbreak();
	}

	CalculateSmallFactorNumbersInternal<0>(1, 1, number(-1));

	{
		privPrimesUpToSqrtLimitSortedCount = nPrimes - InlinePrimesInSearch - 1;

		if (isSubmit)
			return;

		// Gather data for linear search and do preliminary filtering
		// All filters combined leave 971348 (971686 when USE_CONJECTURES is on) numbers out of the first 1000000
		// It's a great speed-up compared to the recursive search
		std::vector<std::pair<number, number>> f;
		for (unsigned int i = 1; i <= SearchLimit::value / SearchLimit::LinearLimit; ++i)
		{
			// Smaller number in an amicable pair must not be divisible by 6
			if ((i % 6) == 0)
				continue;

#if USE_CONJECTURES
			// If a number is divisible by 2^3, then it must not be divisible by 5 or 7
			// This is not proved, but true for all known amicable pairs
			if ((i % 8) == 0)
				if (((i % 5) == 0) || ((i % 7) == 0))
					continue;
#endif

			number sumI;
			f.clear();
			Factorize(i, sumI, f);
			// Smaller number in an amicable pair must not be deficient
			if (sumI >= i * 2)
			{
				// Linear search goes through numbers of a form "k * p" where p > k * 2
				// We collect all values of k here
				// Since we know k, we can calculate partial factorization of bigger number in an amicable pair
				// n = S(m) - m = S(k * p) - k * p = S(k) * (p + 1) - k * p
				// So factorization will contain GCD(S(k) * (p + 1), k * p)
				// We know k, S(k) and we also know that p + 1 is even, hence sumI * 2
				const unsigned int gcd = static_cast<unsigned int>(GCD(i, sumI * 2));
				number sumGCD;
				f.clear();
				Factorize(gcd, sumGCD, f);

				// This partial factorization must be deficient in order for the bigger number to be deficient too
				// In fact, we can use even stronger restrictions:
				// If (M,N) is an amicable pair, M<N, then S(M)=S(N)=M+N
				//
				// S(M)/M-1=N/M
				// S(N)/N-1=M/N
				//
				// let M=m*M1, N=n*N1, gcd(m,M1)=1, gcd(n,N1)=1, M1 > 1, N1 > 1
				// S(M)/M - 1 > S(m)/m - 1
				// S(N)/N - 1 > S(n)/n - 1
				// let S(m)/m - 1 = x
				// S(M)/M - 1 = y > x
				// S(N)/N - 1 = 1/y < 1/x
				// if S(n)/n - 1 >= 1/x then S(N)/N - 1 > S(n)/n - 1 >= 1/x, and (M,N) can't be an amicable pair
				// so S(n)/n must be < 1 + 1/x
				//
				// x = sumI / i - 1
				// x = (sumI - i) / i
				// 1 / x = i / (sumI - i)
				//
				// sumGCD / gcd < 1 + i / (sumI - i)
				// sumGCD < gcd * (1 + i / (sumI - i))
				// sumGCD * (sumI - i) < gcd * (sumI - i) * (1 + i / (sumI - i))
				// sumGCD * (sumI - i) < gcd * ((sumI - i) + i)
				// sumGCD * (sumI - i) < gcd * sumI
				if (sumGCD * (sumI - i) < gcd * sumI)
					privLinearSearchData.push_back(std::pair<unsigned int, unsigned int>(i, static_cast<unsigned int>(sumI)));
			}
		}
	}

	// PQ corresponds to tables P and Q in lemma 2.1 from
	// "Computation of All the Amicable Pairs Below 10^10 By H.J.J.te Riele": http://www.ams.org/journals/mcom/1986-47-175/S0025-5718-1986-0842142-3/S0025-5718-1986-0842142-3.pdf
	// The only difference is that we calculate exact (hence better) upper bounds for S(m)/m instead of inexact estimates
	const number maxI = 16384;
	std::vector<std::pair<number, number>> PQ[SumEstimatesSize];
	for (number j = 0; j < SumEstimatesSize; ++j)
		PQ[j].resize(maxI);

	PrimesUpToSqrtLimitIterator prevP(1);
	PrimesUpToSqrtLimitIterator p(2);
	number PQ_size = maxI;
	auto MultiplyWithSaturation = [](const number a, const number b)
	{
		number h;
		const number result = _umul128(a, b, &h);
		return h ? number(-1) : result;
	};
	for (number i = 0; (i < maxI) && (p.Get() <= max(SearchLimit::PrimesUpToSqrtLimitValue, CompileTimePrimes<CompileTimePrimesCount>::value)); ++i, ++p)
	{
		number j = 1;
		PrimesUpToSqrtLimitIterator q(p);
		++q;

		PQ[0][i].first = p.Get();
		PQ[0][i].second = GetMaxSumRatio(prevP, MultiplyWithSaturation(p.Get(), q.Get()));

		for (; (j < SumEstimatesSize) && (q.Get() <= max(SearchLimit::PrimesUpToSqrtLimitValue, CompileTimePrimes<CompileTimePrimesCount>::value)); ++j)
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
		for (number j = 0; j < SumEstimatesSize; ++j)
			if (PQ[j][i].first != number(-1))
				--PQ[j][i].first;

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

	number index = 0;
	for (PrimesUpToSqrtLimitIterator it(2); index < CompileTimePrimesCount; ++it, ++index)
	{
		if (index > 0)
		{
			privPrimeInverses[index * 2] = CalculateInverse(it.Get());
			privPrimeInverses[index * 2 + 1] = number(-1) / it.Get();
		}
	}

	// cppcheck-suppress memleak
}
