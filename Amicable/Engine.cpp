#include "stdafx.h"
#include "PrimeTables.h"
#include "Engine.h"
#include "RangeGen.h"
#include "RangeSieve.h"
#if AMICABLE_X64
#include <immintrin.h>
#endif
#include <limits>
#include <atomic>

//#pragma optimize("", off)
//#undef FORCEINLINE
//#define FORCEINLINE NOINLINE

PRAGMA_WARNING(push, 1)
PRAGMA_WARNING(disable : 4091 4917 4625 4626 5026 5027)
#include <boinc_api.h>
PRAGMA_WARNING(pop)

// "a" has at most two factors (a = p1 * p2 or a = p^2 or a is prime)
// We know that p1 >= p and p2 <= q
//
// We have 3 cases:
// 1) a is prime, S(a) = a + 1
// 2) a = p * q, where q = a / p
//		S(a) = (p + 1) * (a / p + 1) = a + p + a / p + 1
//
//		Let's find the maximum possible value for S(a):
//
//		f(x) = a + x + a / x + 1 -> max, 2 < p <= x <= a / p
//		f'(x) = 1 - a / x^2
//		we have an optimum of f(x) when x = sqrt(a)
//		f''(x) = 2 * a / x^3 > 0, so it's a minimum
//		the minimum is reached when x = sqrt(a)
//		the maximum is reached when x = p or x = a / p
//
//		x1 = p:		f(x1) = a + p + a / p + 1
//		x2 = a / p:	f(x2) = a + a / p + p + 1
//
//		So the maximum possible value for S(a) is f(x1) = f(x2) = a + p + a / p + 1
// 3) a = p^2, S(a) = 1 + p + p^2 = 1 + p + a
//
// Value from case (2) is always greater than value from case (3):
// a + p + a / p + 1 - (1 + p + a) = a + p + a / p + 1 - 1 - p - a = a / p > 0
//
// So, S(a) <= a + p + a / p + 1 = a + p + q + 1
FORCEINLINE num64 MaximumSumOfDivisors2(const num64 a, const num64 p, const num64 q)
{
	num64 m1;
	if (AddAndDetectOverflow(a, p + q + 1, &m1))
	{
		return num64(-1);
	}
	return m1;
}

// a = pqr, p0 <= p < q < r, p0^4 > a, S(a) = (p+1)(q+1)(r+1)
// F(p,q) = (p + 1)(q + 1)(a / pq + 1) = (pq + p + q + 1)(a / pq + 1) = a + a / q + a / p + a / pq + pq + p + q + 1
//
// F(p,q) = a + a / p + a / q + a / pq + pq + p + q + 1 -> max
// p >= p0
// p <= q <= sqrt(a) / p
//
// F'p(p,q) = -a / p^2 - a / p^2q + q + 1
// F'q(p,q) = -a / q^2 - a / pq^2 + p + 1
//
// F''p(p,q) = 2a / p^3 + 2a / p^3q > 0
// F''pq(p,q) = a / p^2q^2 + 1 > 0
//
// F''q(p,q) = 2a / q^3 + 2a / pq^3 > 0
// F''qp(p,q) = a / p^2q^2 + 1 > 0
//
// All second partial derivatives are always positive, therefore F(p, q) doesn't have maximum inside the given area (p >= p0, p <= q <= sqrt(a) / p)
//
// 1) p=p0: F(p0,q) = F0(q) = a + a / p0 + a / q + a / p0q + p0q + p0 + q + 1
//		F0(q) = a + a / p0 + a / q + a / p0q + p0q + p0 + q + 1
//		F0(q) = q * (p0 + 1) + (a + a / p0 + p0 + 1) + (a + a / p0) / q
//
//		F0'(q) = p0 + 1 - (a + a / p0) / q^2
//		F0''(q) = 2 * (a + a / p0) / q^3 > 0
//		therefore F0(q) has a maximum when q = p0 or q = sqrt(a) / p0
//
//		F0(p0) = a + a / p0 + a / p0 + a / p0^2 + p0^2 + p0 + p0 + 1
//		F0(p0) = a + 2*a / p0 + a / p0^2 + p0^2 + 2*p0 + 1
//		F0(p0) = a + 2 * (a / p0 + p0) + 1 + p0^2 + a / p0^2
//		since p0^4 > a then p0^2 > a / p0^2
//		a / p0^2 < p0^2
//		F0(p0) < a + 2 * (a / p0 + p0) + 1 + p0^2 + p0^2
//		F0(p0) < a + 2 * (a / p0 + p0) + 1 + 2 * p0^2
//		F0(p0) < a + 2 * (a / p0 + p0 + p0^2) + 1
//		F0(p0) < a + (a / p0 + p0^2 + p0) * 2 + 1
//		^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//
//		F0(sqrt(a) / p0) = a + a / p0 + a / (sqrt(a) / p0) + a / (p0 * sqrt(a) / p0) + p0 * (sqrt(a) / p0) + p0 + sqrt(a) / p0 + 1
//		F0(sqrt(a) / p0) = a + a / p0 + a / sqrt(a) * p0 + a / sqrt(a) + sqrt(a) + p0 + sqrt(a) / p0 + 1
//		F0(sqrt(a) / p0) = a + a / p0 + sqrt(a) * p0 + sqrt(a) + sqrt(a) + p0 + sqrt(a) / p0 + 1
//		F0(sqrt(a) / p0) = a + a / p0 + sqrt(a) * (p0 + 2 + 1 / p0) + p0 + 1
//		since p0^4 > a then p0^2 > sqrt(a)
//		sqrt(a) < p0^2
//		F0(sqrt(a) / p0) < a + a / p0 + p0^2 * (p0 + 2 + 1 / p0) + p0 + 1
//		F0(sqrt(a) / p0) < a + a / p0 + p0^3 + 2 * p0^2 + p0 * 2 + 1
//		F0(sqrt(a) / p0) < a + a / p0 + p0^3 + (p0^2 + p0) * 2 + 1
//		^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//
//		F0(p0) - F0(sqrt(a) / p0) = a / p0 - p0^3
//		since p0^4 > a then p0^3 > a / p0
//		therefore F0(p0) - F0(sqrt(a) / p0) < 0
//		F0(p0) < F0(sqrt(a) / p0)
//
// 2) q=p: F(p,p) = F1(p) = a + a / p + a / p + a / p^2 + p^2 + p + p + 1
//		F1(p) = a + a / p + a / p + a / p^2 + p^2 + p + p + 1
//		F1(p) = a + 2 * a / p + a / p^2 + p^2 + 2 * p + 1
//		F1(p) = p^2 + 2 * p + 1 + a * (1 + 2 / p + 1 / p^2)
//
//		F1'(p) = 2 * p + 2 + a * (-2 / p^2 - 2 / p^3)
//		F1''(p) = 2 + a * (4 / p^3 + 6 / p^4) > 0
//		therefore F1(p) has a maximum when p = p0 or p = sqrt(a) / p0
//
//		F1(p0) = p0^2 + 2 * p0 + 1 + a * (1 + 2 / p0 + 1 / p0^2)
//		F1(p0) = p0^2 + 2 * p0 + 1 + a + 2 * a / p0 + a / p0^2
//		since p0^4 > a then p0^2 > a / p0^2
//		a / p0^2 < p0^2
//		F1(p0) < p0^2 + 2 * p0 + 1 + a + 2 * a / p0 + p0^2
//		F1(p0) < a + 2 * a / p0 + 2 * p0^2 + 2 * p0 + 1
//		F1(p0) < a + (a / p0 + p0^2 + p0) * 2 + 1 = F0(p0)
//
//		F1(sqrt(a) / p0) = (sqrt(a) / p0)^2 + 2 * (sqrt(a) / p0) + 1 + a * (1 + 2 / (sqrt(a) / p0) + 1 / (sqrt(a) / p0)^2)
//		F1(sqrt(a) / p0) = a / p0^2 + 2 * (sqrt(a) / p0) + 1 + a * (1 + 2 / sqrt(a) * p0 + 1 / (a / p0^2))
//		F1(sqrt(a) / p0) = a / p0^2 + 2 * p0 * sqrt(a) + 1 + a * (1 + 2 * p0 / sqrt(a) + p0^2 / a)
//		F1(sqrt(a) / p0) = a / p0^2 + 2 * p0 * sqrt(a) + 1 + a + 2 * a * p0 / sqrt(a) + p0^2
//		F1(sqrt(a) / p0) = a / p0^2 + 2 * p0 * sqrt(a) + 1 + a + 2 * sqrt(a) * p0 + p0^2
//		F1(sqrt(a) / p0) = a / p0^2 + 4 * p0 * sqrt(a) + 1 + a + p0^2
//		a / p0^2 < p0^2
//		sqrt(a) < p0^2
//		F1(sqrt(a) / p0) < p0^2 + 4 * p0 * p0^2 + 1 + a + p0^2
//		F1(sqrt(a) / p0) < a + 2 * p0^2 + 4 * p0^3 + 1
//		F1(sqrt(a) / p0) < a + 4 * p0^3 + 2 * p0^2 + 1
//		^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//
//		F1(sqrt(a) / p0) - F0(sqrt(a) / p0) = (4 * p0^3) - (a / p0 + p0^3 + p0 * 2)
//		p0^4 > a then p0^3 > a / p0
//		(4 * p0^3) - (a / p0 + p0^3 + p0 * 2) > (4 * p0^3) - (p0^3 + p0^3 + p0 * 2)
//		(4 * p0^3) - (p0^3 * 2 + p0 * 2)
//		(4 * p0^3) - p0^3 * 2 - p0 * 2
//		p0^3 * 2 - p0 * 2
//		2 * p0 * (p0^2 - p0)
//		2 * p0^2 * (p0 - 1) > 0
//
//		F1(sqrt(a) / p0) > F0(sqrt(a) / p0)
//
// 3) q = sqrt(a) / p: F(p, sqrt(a) / p) = F2(p) = a + a / p + a / (sqrt(a) / p) + a / (p * sqrt(a) / p) + p * (sqrt(a) / p) + p + sqrt(a) / p + 1
//		F2(p) = a + a / p + a / (sqrt(a) / p) + a / (p * sqrt(a) / p) + p * (sqrt(a) / p) + p + sqrt(a) / p + 1
//		F2(p) = a + a / p + a / sqrt(a) * p + a / sqrt(a) + sqrt(a) + p + sqrt(a) / p + 1
//		F2(p) = a + a / p + sqrt(a) * p + sqrt(a) * 2 + p + sqrt(a) / p + 1
//		F2(p) = a + a / p + sqrt(a) * (p + 2 + 1 / p) + p + 1
//
//		F2'(p) = -a / p^2 + sqrt(a) * (1 - 1 / p^2) + 1
//		F2''(p) = 2 * a / p^3 + sqrt(a) * 2 / p^3 > 0
//		therefore F2(p) has a maximum when p = p0 or q = p
//		q = p
//		sqrt(a) / p = p
//		sqrt(a) = p^2
//		p^4 = a
//		since p >= p0 and p0^4 > a then we only have to check p = p0:
//
//		F2(p0) = a + a / p0 + sqrt(a) * (p0 + 2 + 1 / p0) + p0 + 1
//		p0^4 > a, p0^2 > sqrt(a), sqrt(a) < p0^2
//		F2(p0) < a + a / p0 + p0^2 * (p0 + 2 + 1 / p0) + p0 + 1
//		F2(p0) < a + a / p0 + p0^3 + 2 * p0^2 + p0 + p0 + 1
//		F2(p0) < a + a / p0 + p0^3 + 2 * p0^2 + 2 * p0 + 1
//		F2(p0) < a + a / p0 + p0^3 + 2 * (p0^2 + p0) + 1
//		^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//
//		F2(p0) - F1(sqrt(a) / p0) = a / p0 + 2 * p0 - 3 * p0^3
//		p0^4 > a, p0^3 > a / p0, a / p0 < p0^3
//		a / p0 + 2 * p0 - 3 * p0^3 < p0^3 + 2 * p0 - 3 * p0^3
//		p0^3 + 2 * p0 - 3 * p0^3
//		2 * p0 - 2 * p0^3
//		2 * p0 * (1 - p0^2) < 0
//		which means that F2(p0) < F1(sqrt(a) / p0)
//
//		The result: the largest estimate is
//		F1(sqrt(a) / p0) = F(sqrt(a) / p0, sqrt(a) / p0) = (sqrt(a) / p0 + 1)(sqrt(a) / p0 + 1)(a / (a / p0^2) + 1)
//		(sqrt(a) / p0 + 1)^2 * (a / a * p0^2 + 1)
//		(sqrt(a) / p0 + 1)^2 * (p0^2 + 1)
//		(a / p0^2 + 2 * sqrt(a) / p0 + 1) * (p0^2 + 1)
//		(a / p0^2 + 2 * sqrt(a) / p0 + 1) * p0^2 + (a / p0^2 + 2 * sqrt(a) / p0 + 1)
//		a + 2 * p0 * sqrt(a) + p0^2 + a / p0^2 + 2 * sqrt(a) / p0 + 1
//		sqrt(a) < p0^2
//		a + 2 * p0 * p0^2 + p0^2 + a / p0^2 + 2 * p0^2 / p0 + 1
//		a + 4 * p0^3 + p0^2 + a / p0^2 + 1
//		a / p0^2 < p0^2
//		a + 4 * p0^3 + 2 * p0^2 + 1
FORCEINLINE num64 MaximumSumOfDivisors3(const num64 a, const num64 p0, const num64 a_div_p0)
{
	// F1(sqrt(a) / p0) is the analytically deduced sup(S(a)) - see comments above
	//
	// I use F0(p0) instead: I analyzed all cases when "a" has 3 factors
	// and concluded that F0(p0) will always be correct - see comments in CheckMaximumSumOfDivisors3()
	// Moreover, comments in CheckMaximumSumOfDivisors3 show that F0(p0) is always strictly greater than S(a)
	// So we can just remove "+ 1".

	//  F0(p0)
	num64 m1;
	if (AddAndDetectOverflow(a, (a_div_p0 + p0 * p0 + p0) * 2 /*+ 1*/, &m1))
	{
		return num64(-1);
	}
	return m1;
	
	// F0(sqrt(a) / p0)
	//const num64 p0_2 = p0 * p0;
	//return a + a_div_p0 + p0_2 * p0 + (p0_2 + p0) * 2 + 1;

	// F1(sqrt(a) / p0)
	//return a + 4 * p0 * p0 * p0 + 2 * p0 * p0 + 1;
}

NOINLINE num64 MaximumSumOfDivisors3NoInline(const num64 a, const num64 p0, const num64 a_div_p0)
{
	return MaximumSumOfDivisors3(a, p0, a_div_p0);
}

// Lemma 2.1 from http://www.ams.org/journals/mcom/1986-47-175/S0025-5718-1986-0842142-3/S0025-5718-1986-0842142-3.pdf
FORCEINLINE num64 GetCoeffForMaximumSumOfDivisorsN(const num64 m, const num64 i, num64& j)
{
	for (; ; --j)
	{
		const SumEstimateData* data = SumEstimates[j] + (i >> 3);
		if (m > data->P)
			return data->Q;
	}
}

FORCEINLINE num64 MaximumSumOfDivisorsN(const num64 m, const num64 i, num64& j)
{
	num64 highProduct;
	_umul128(m, GetCoeffForMaximumSumOfDivisorsN(m, i, j), &highProduct);

	num64 m1;
	if (AddAndDetectOverflow(m, highProduct, &m1))
	{
		return num64(-1);
	}

	return m1;
}

template<int numPrimesCheckedSoFar>
FORCEINLINE bool IsNumEligible(const num64 a, const num64 sumA, const num64 targetSum, num64& j)
{
	num64 highProduct;
	const num64 maxPossibleSum = _umul128(sumA, MaximumSumOfDivisorsN(a, numPrimesCheckedSoFar, j), &highProduct);
	return (maxPossibleSum >= targetSum) || highProduct;
}

#define M(X) template<> FORCEINLINE bool IsNumEligible<X>(const num64, const num64, const num64, num64&) { return true; }
	M(0)M(1)M(2)M(3)M(4)M(5)M(6)M(7)M(8)M(9)
	M(10)M(11)M(12)M(13)M(14)M(15)
#undef M

template<> FORCEINLINE bool IsNumEligible<IS_NUM_ELIGIBLE_BEGIN>(const num64 a, const num64 sumA, const num64 targetSum, num64& j)
{
	for (num64 k = 0; ; ++k)
	{
		if (a <= SumEstimatesBeginP[k])
		{
			num64 highProduct;
			_umul128(a, SumEstimatesBeginQ[k], &highProduct);
			j = k;

			num64 a1;
			if (AddAndDetectOverflow(a, highProduct, &a1))
			{
				return true;
			}

			const num64 maxPossibleSum = _umul128(sumA, a1, &highProduct);
			return (maxPossibleSum >= targetSum) || highProduct;
		}
	}
}

template<int PrimeIndex>
FORCEINLINE bool CheckDivisibility(num64& a, num64& sumA, const num64 targetSum, num64& n2TargetSum, num64& j)
{
	enum {p = CompileTimePrimes<PrimeIndex>::value};

	IF_CONSTEXPR(((PrimeIndex & 7) == 0) && (PrimeIndex != 56) && (PrimeIndex != 72) && (PrimeIndex != 88))
	{
		if (!IsNumEligible<PrimeIndex>(a, sumA, targetSum, j))
			return false;
	}

	// M / N = (M * value) mod 2^64
	// This means that (M * value) must be <= (num64(-1) / N)
	// Otherwise it's not divisible by N
	num64 q = a * PrimeInverses[PrimeIndex].first;
	if (UNLIKELY(q <= PrimeInverses[PrimeIndex].second))
	{
		const num64 prevSumA = sumA;
		num64 nextSumA = prevSumA * (p + 1);
		num64 curPower = 0;
		for (;;)
		{
			a = q;
			q *= PrimeInverses2[PrimeIndex].first;
			if (q > PrimeInverses2[PrimeIndex].second)
			{
				break;
			}
			nextSumA = nextSumA * p + prevSumA;
			++curPower;
		}
		sumA = nextSumA;

		// If a num64 is fully factored, then exit immediately
		if (a < CompileTimePrimes<PrimeIndex + 1>::value * CompileTimePrimes<PrimeIndex + 1>::value)
			return true;

		// found new prime factor, let's check that N is not abundant yet
		num64 h;
		const num64 minSum = _umul128(sumA, a + 1, &h);
		if ((minSum > targetSum) || h)
			return false;

		IF_CONSTEXPR((PrimeIndex > IS_NUM_ELIGIBLE_BEGIN * 2) && ((PrimeIndex & 7) != 0) && (((PrimeIndex + 1) & 7) != 0))
		{
			if (!IsNumEligible<PrimeIndex + 1>(a, sumA, targetSum, j))
				return false;
		}

		// and then check that targetSum is divisible by partial sum
		enum Sums : num64
		{
			s1 = num64(p) + 1,
			s2 = num64(p) * p + s1,
			s3 = num64(p) * p * p + s2,
			s4 = num64(p) * p * p * p + s3,
			s5 = num64(p) * p * p * p * p + s4,
		};

		static CACHE_ALIGNED const num64 InverseS[5][3] = {
			{MultiplicativeInverseEven<s1>::value, MultiplicativeInverseEven<s1>::shift, num64(-1) / s1},
			{MultiplicativeInverse<s2>::value, 0, num64(-1) / s2},
			{MultiplicativeInverseEven<s3>::value, MultiplicativeInverseEven<s3>::shift, num64(-1) / s3},
			{MultiplicativeInverse<s4>::value, 0, num64(-1) / s4},
			{MultiplicativeInverseEven<s5>::value, MultiplicativeInverseEven<s5>::shift, num64(-1) / s5},
		};

		if (curPower < ARRAYSIZE(InverseS))
		{
			n2TargetSum = _rotr64(n2TargetSum * InverseS[curPower][0], static_cast<int>(InverseS[curPower][1]));
			if (n2TargetSum > InverseS[curPower][2])
				return false;
		}
	}

	return CheckDivisibility<PrimeIndex + 1>(a, sumA, targetSum, n2TargetSum, j);
}

template<> FORCEINLINE bool CheckDivisibility<CompileTimePrimesCount>(num64&, num64&, const num64, num64&, num64&) { return true; }

static THREAD_LOCAL unsigned int NumFoundPairs;

void SetNumFoundPairsInThisThread(unsigned int value)
{
	NumFoundPairs = value;
}

unsigned int GetNumFoundPairsInThisThread()
{
	return NumFoundPairs;
}

bool g_PrintNumbers = true;
FILE* g_outputFile = nullptr;

std::vector<num128>* g_FoundNumbersLog = nullptr;

NOINLINE void NumberFound(const num128 n1)
{
	++NumFoundPairs;
	if (g_FoundNumbersLog)
	{
		g_FoundNumbersLog->push_back(n1);
	}
	if (g_PrintNumbers)
	{
		char buf[40];
		fprintf(g_outputFile, "%s\n", itoa128(n1, buf, sizeof(buf)));
		fflush(g_outputFile);
	}
}

// Fast integer cube root. Returns num64 m such that m^3 <= n < (m+1)^3 for all n > 0
FORCEINLINE num64 IntegerCbrt(const num64 n)
{
	if (n < 8)
	{
		return 1;
	}

	unsigned long index;
	_BitScanReverse64(&index, n);

	index = (index * ((65536 / 3) + 1)) >> 16;

	num64 result = num64(1) << index;
	for (num64 cur_bit = result >> 1; cur_bit > 0; cur_bit >>= 1)
	{
		const num64 k = result | cur_bit;
		if ((k <= 2642245) && (k * k * k <= n))
		{
			result = k;
		}
	}
	return result;
}

FORCEINLINE num64 Root4(const num64 n)
{
	return static_cast<num64>(static_cast<__int64>(sqrt(sqrt(n))));
}

// The trial division of CheckPairInternal after its first CompileTimePrimesCount primes, from the prime p with index
// numPrimesCheckedSoFar (as in PrimeInverses3: 1 is 3): n2 has no prime factors below p, sum is sigma of the part of the number
// divided out so far, n2TargetSum = targetSum / sum (exact), indexForMaximumSumOfDivisorsN as set by IsNumEligible
FORCEINLINE void CheckPairInternalTail(const num128& n1, const num64 targetSum, num64 n2TargetSum, num64 n2, num64 sum, num64 p, num64 numPrimesCheckedSoFar, num64 indexForMaximumSumOfDivisorsN)
{
	num64 n2_sqrt4 = Root4(n2);

	while (p <= n2_sqrt4)
	{
		num64 h;
		num64 q = n2 * PrimeInverses3[numPrimesCheckedSoFar];
		_umul128(q, p, &h);
		if (UNLIKELY(h == 0))
		{
			num64 n = p;
			num64 curSum = p + 1;
			n2 = q;

			while (p <= n2)
			{
				q = n2 * PrimeInverses4[numPrimesCheckedSoFar];
				_umul128(q, p, &h);
				if (h != 0)
					break;
				n *= p;
				curSum += n;
				n2 = q;
			}
			sum *= curSum;

			if (n2 == 1)
			{
				if (sum == targetSum)
					NumberFound(n1);
				return;
			}

			// found new prime factor, let's check that N is not abundant yet
			if (sum * (n2 + 1) > targetSum)
				return;

			// and then check that targetSum is still divisible by sum
			n2TargetSum = targetSum / sum;
			if ((targetSum % sum) != 0)
				return;
			n2_sqrt4 = Root4(n2);
		}

		p += static_cast<unsigned int>(NextPrimeShifts[numPrimesCheckedSoFar]);

		++numPrimesCheckedSoFar;
		if (((numPrimesCheckedSoFar & 7) == 0) && (MaximumSumOfDivisorsN(n2, numPrimesCheckedSoFar, indexForMaximumSumOfDivisorsN) < n2TargetSum))
			return;
	}

	// Here p^4 > n2, so n2 can have at most 3 factors
	if (n2TargetSum & 1)
	{
		// if sum of divisors is odd, then it must be a perfect square
		if (!IsPerfectSquareCandidate(n2))
			return;

		p = static_cast<num64>(static_cast<__int64>(sqrt(n2)));
		const num64 p2 = p * p;
		if (p2 != n2)
			return;

		// a perfect square which has at most 3 factors can be only of the form "n=p^2"
		// so we can just check this case and exit
		if (p2 + p + 1 == n2TargetSum)
			NumberFound(n1);

		return;
	}

	// Check the case when n2 is prime (it has 1 factor)
	if (n2 + 1 == n2TargetSum)
	{
		if (IsPrime(n2))
		{
			NumberFound(n1);
		}
		return;
	}

	// n2 must have 2 or 3 factors in order to have such n2TargetSum
	// if n2 has 3 factors, then minimal possible sum of divisors is n2 + n2^(2/3) + n2^(1/3) + 1
	// if target sum is less than this, then n2 must not have 3 factors in order to have such n2TargetSum
	// So we only run the loop up to cube root of n2 only if target sum is >= than this
	num64 n2_cbrt = IntegerCbrt(n2);
	if (n2TargetSum > n2 + n2_cbrt * (n2_cbrt + 1))
	{
		while (p <= n2_cbrt)
		{
			num64 h;
			num64 q = n2 * PrimeInverses3[numPrimesCheckedSoFar];
			_umul128(q, p, &h);
			if (UNLIKELY(h == 0))
			{
				num64 n = p;
				num64 curSum = p + 1;
				n2 = q;

				while (p <= n2)
				{
					q = n2 * PrimeInverses4[numPrimesCheckedSoFar];
					_umul128(q, p, &h);
					if (h != 0)
						break;
					n *= p;
					curSum += n;
					n2 = q;
				}
				sum *= curSum;

				if (n2 == 1)
				{
					if (sum == targetSum)
						NumberFound(n1);
					return;
				}

				// found new prime factor, let's check that N is not abundant yet
				if (sum * (n2 + 1) > targetSum)
					return;

				// and then check that targetSum is still divisible by sum
				n2TargetSum = targetSum / sum;
				if ((targetSum % sum) != 0)
					return;

				n2_cbrt = IntegerCbrt(n2);
			}

			if ((numPrimesCheckedSoFar & 15) == 0)
			{
				const num64 q1 = PrimeReciprocals[numPrimesCheckedSoFar].DivideNoRemainder(n2);
				if (MaximumSumOfDivisors3(n2, p, q1) < n2TargetSum)
					return;
			}

			p += static_cast<unsigned int>(NextPrimeShifts[numPrimesCheckedSoFar]);
			++numPrimesCheckedSoFar;
		}
	}

	// Here p^3 > n2, so n2 can have at most 2 factors
	while (p * p <= n2)
	{
		num64 q;
		if (numPrimesCheckedSoFar < ReciprocalsTableSize)
		{
			if (PrimeReciprocals[numPrimesCheckedSoFar].Divide(n2, p, q))
			{
				num64 n = p;
				num64 curSum = p + 1;
				n2 = q;
				while (p <= n2)
				{
					if (!PrimeReciprocals[numPrimesCheckedSoFar].Divide(n2, p, q))
						break;
					n *= p;
					curSum += n;
					n2 = q;
				}
				sum *= curSum;
				break;
			}
			++numPrimesCheckedSoFar;
		}
		else
		{
			q = n2 / p;
			num64 r = n2 % p;
			if (r == 0)
			{
				num64 n = p;
				num64 curSum = p + 1;
				n2 = q;

				while (p <= n2)
				{
					q = n2 / p;
					r = n2 % p;
					if (r != 0)
						break;
					n *= p;
					curSum += n;
					n2 = q;
				}
				sum *= curSum;
				break;
			}
		}

		if (MaximumSumOfDivisors2(n2, p, q) < n2TargetSum)
			return;

		// If n2 + 1 != n2TargetSum then n2 must be composite to satisfy S(N) = target sum
		if (n2 + 1 != n2TargetSum)
		{
			if ((n2TargetSum & 1) == 0)
			{
				// 1) Check the case n2 = p*q, p < q
				// B=(n2TargetSum-n2-1)/2
				// D=B^2-n2
				// p=B-sqrt(D)
				// q=B+sqrt(D)
				const num64 B = (n2TargetSum - n2 - 1) >> 1;
				num128 D = num128(B) * B;
				if (D > n2)
				{
					D -= n2;
					if (IsPerfectSquareCandidate(LowWord(D)))
					{
						const num64 sqrt_D = IntegerSquareRoot(D);
						if (num128(sqrt_D) * sqrt_D == D)
						{
							p = B - sqrt_D;
							q = B + sqrt_D;
							if ((p * q == n2) && IsPrime(p) && IsPrime(q))
							{
								NumberFound(n1);
							}
							return;
						}
					}
				}
			}
			else
			{
				// 2) Check the case n2 = p^2
				p = n2TargetSum - n2 - 1;
				if ((p * p == n2) && IsPrime(p))
				{
					NumberFound(n1);
				}
			}
		}
		else
		{
			// n2 + 1 == n2TargetSum, n2 must be a prime number
			if (IsPrime(n2))
			{
				NumberFound(n1);
			}
		}

		// There are no other cases, so just exit now
		return;
	}

	if (n2 > 1)
		sum *= n2 + 1;

	if (sum == targetSum)
		NumberFound(n1);
}

FORCEINLINE void CheckPairInternal(const num128& n1, const num64 targetSum, num64 n2TargetSum, num64 n2, num64 sum)
{
	num64 indexForMaximumSumOfDivisorsN;
	if (!CheckDivisibility<1>(n2, sum, targetSum, n2TargetSum, indexForMaximumSumOfDivisorsN))
		return;

	if (n2 < CompileTimePrimes<CompileTimePrimesCount>::value * CompileTimePrimes<CompileTimePrimesCount>::value)
	{
		if (n2 > 1)
			sum *= n2 + 1;

		if (sum == targetSum)
			NumberFound(n1);

		return;
	}

	const num64 numPrimesCheckedSoFar = CompileTimePrimesCount;

	if (!IsNumEligible<CompileTimePrimesCount>(n2, sum, targetSum, indexForMaximumSumOfDivisorsN))
		return;

	// We must have sum == targetSum in the end
	// sum is constructed by multiplying it by (p^(e+1) - 1) for each found prime factor p^e
	// p^(e+1) - 1 is an integer num64
	//
	// no matter how we split the sum, in "partial_sum1 * partial_sum2 == targetSum" both partial sums must be integer
	//
	// so if targetSum is not divisible by sum, we'll have this: sum * x == targetSum where x is not integer
	// which means we can stop right now
	n2TargetSum = targetSum / sum;
	if ((targetSum % sum) != 0)
		return;

	CheckPairInternalTail(n1, targetSum, n2TargetSum, n2, sum, CompileTimePrimes<CompileTimePrimesCount>::value, numPrimesCheckedSoFar, indexForMaximumSumOfDivisorsN);
}

NOINLINE static void CheckPairInternalNoInline(const num128& n1, const num64 targetSum, num64 n2TargetSum, num64 n2, num64 sum)
{
	CheckPairInternal(n1, targetSum, n2TargetSum, n2, sum);
}

// n2 is either a prime, squared prime, or has two distinct prime factors
NOINLINE void FinalCheck128(const num128& n1, const num128& targetSum, const num128& n2)
{
	if (LowWord(targetSum) & 1)
	{
		// targetSum is odd, it can only happen if n2 = p^2
		// targetSum = 1 + p + p^2 = 1 + p + n2, p = targetSum - n2 - 1
		const num128 p = targetSum - n2 - 1;
		if ((p * p == n2) && IsPrime(LowWord(p)))
		{
			NumberFound(n1);
		}
	}
	else
	{
		// targetSum is even, n2 is either p or p * q
		if (targetSum == n2 + 1)
		{
			/* May give some false positives if n2 is composite
			   TODO: check that n2 is prime

			   Example:

			   M = 467915168733863040016 = 2^4*17*199*431*19973*1004210069
			   N = 499109513316120575984 = 2^4*1908062567*16348701097
			   targetSum = 31194344582257536000
			   n2 = 31194344582257535999 = 1908062567*16348701097
			*/
			NumberFound(n1);
		}
		else
		{
			// B=(n2TargetSum-n2-1)/2
			// D=B^2-n2
			// p=B-sqrt(D)
			// q=B+sqrt(D)
			const num128 B = (targetSum - n2 - 1) >> 1;
			num128 D = B * B;
			if (D > n2)
			{
				D -= n2;
				if (IsPerfectSquareCandidate(LowWord(D)))
				{
					const num64 sqrtD = IntegerSquareRoot(D);
					if (num128(sqrtD) * sqrtD == D)
					{
						const num128 p = B - sqrtD;
						const num128 q = B + sqrtD;
						if ((p * q == n2) && IsPrime(LowWord(p)) && IsPrime(LowWord(q)))
						{
							NumberFound(n1);
						}
					}
				}
			}
		}
	}
}

FORCEINLINE bool InitialCheck128(const num128 n1, num128& targetSum, num128& n2)
{
	n2 = targetSum - n1;

	unsigned long powerOf2;
	if (LowWord(n2))
	{
		_BitScanForward64(&powerOf2, LowWord(n2));
		n2 >>= powerOf2;
	}
	else
	{
		_BitScanForward64(&powerOf2, HighWord(n2));
		n2 = HighWord(n2) >> powerOf2;
		powerOf2 += 64;
	}

	if (powerOf2 > 0)
	{
		targetSum *= PowersOf2_128DivisibilityData[powerOf2].first;
		if (targetSum > PowersOf2_128DivisibilityData[powerOf2].second)
		{
			return false;
		}
	}

	// sigma(n2) >= n2 + 1 > targetSum, so it can't be an amicable pair
	// Without this check, such numbers would go through the whole 128-bit trial division loop below if n2 has no small factors
	if (n2 >= targetSum)
	{
		return false;
	}

	if (HighWord(targetSum) == 0)
	{
		return true;
	}

	const num128* prime_inverse = PrimeInverses128;
	const num64* max_sum_ratio = SumEstimates128;
	do
	{
		for (const num128* e = prime_inverse + 16; prime_inverse < e; ++prime_inverse)
		{
			num128 q = n2 * (*prime_inverse);
			if (q > n2)
			{
				continue;
			}

			const InverseData128* data = PowersOfP_128DivisibilityData[prime_inverse - PrimeInverses128];
			for (;;)
			{
				n2 = q;
				q *= *prime_inverse;
				if (HighWord(q) > HighWord(n2))
				{
					break;
				}
				++data;
			}

			if (data->shift)
			{
				if (LowWord(targetSum) & data->shift_bits)
				{
					return false;
				}
				targetSum >>= data->shift;
			}

			targetSum *= data->inverse;
			if (targetSum > data->max_value)
			{
				return false;
			}

			if (n2 >= targetSum)
			{
				return false;
			}

			if (HighWord(targetSum) == 0)
			{
				return true;
			}
		}

#if SumEstimates128Shift == 0
		if (n2 + HighWord(n2 * (*max_sum_ratio)) < targetSum)
#else
		if (n2 + ((n2 * (*max_sum_ratio)) >> (64 - SumEstimates128Shift)) < targetSum)
#endif
		{
			return false;
		}
		++max_sum_ratio;

	} while (prime_inverse < PrimeInverses128 + ReciprocalsTableSize128);

	FinalCheck128(n1, targetSum, n2);
	return false;
}

FORCEINLINE void CheckPair128(const num128 n1, num128 targetSum)
{
	num128 n2;
	if (InitialCheck128(n1, targetSum, n2))
	{
		CheckPairInternalNoInline(n1, LowWord(targetSum), LowWord(targetSum), LowWord(n2), 1);
	}
}

NOINLINE void CheckPair128NoInline(const num128 n1, const num128 targetSum)
{
	CheckPair128(n1, targetSum);
}

// N = r^2 s (DecideAtMost2Factors): the sign of f(r) - D for f(r) = (r + 1)(r + N / r^2), which is the sign of (r + 1)(r^3 + N) - D r^2
// (r < 2^27, N < 2^78, D < 2^54: the products are below 2^120)
FORCEINLINE int R2S_Sign(const num64 r, const num128& N, const num64 D)
{
	const num128 b = ((num128(r * r) * r) + N) * (r + 1);
	const num128 c = num128(D) * (r * r);
	return (b > c) ? 1 : ((b < c) ? -1 : 0);
}

// f'(r) > 0 for f(r) = (r + 1)(r + N / r^2): r^3 f'(r) = 2 r^4 + r^3 - N r - 2 N = r^3 (2 r + 1) - N (r + 2)
FORCEINLINE bool R2S_Rising(const num64 r, const num128& N)
{
	return ((num128(r * r) * r) * (2 * r + 1) > N * (r + 2)) != 0;
}

// Numbers with few prime factors (the same as DecideAtMost2Factors in the GPU version, kernel.cl). D = targetSum - N - 1 must be
// sigma(N) - N - 1, the sum of the divisors of N between 1 and N. If N has 4 or more prime factors, D >= N / p >= N^(3/4) (p: the
// smallest prime factor), with 3 different ones D > r1 r2 + r1 r3 + r2 r3 >= 3 N^(2/3) (AM-GM), N = r^2 s (r != s) has
// D = (r + 1)(r + s) > r^2 + N / r >= 1.8899 N^(2/3) and N = r^3 has D = r + r^2. So if D < 3 floor(N^(1/3))^2 (and N >= 2^20:
// N^(3/4) > 3 N^(2/3)), N can only be the other member of a pair if it is prime (D = 0), p^2 (D = p), p * q (D = p + q, p < q),
// r^3 (D = r + r^2) or r^2 s (D = (r + 1)(r + s)) with primes p, q, r, s: that is decided here without trial divisions (every factor
// that a report assumes prime is tested). N > 1, targetSum > N. Returns true if decided (n1 reported if it is amicable).
FORCEINLINE bool DecideAtMost2Factors(const num128& n1, const num128& targetSum, const num128& N)
{
	const num128 D128 = targetSum - N - 1;
	if (D128 == 0)
	{
		if (HighWord(N) != 0)
		{
			return false; // (no primality test for N >= 2^64 here)
		}
		if (IsPrime(LowWord(N)))
		{
			NumberFound(n1);
		}
		return true;
	}
	if ((HighWord(D128) != 0) || (N < (num64(1) << 20)))
	{
		return false;
	}
	const num64 D = LowWord(D128);

	// a floating point test first: most numbers have D >> N^(2/3) (with the 1% margin it only skips numbers with D >= 3 r^2)
	const double rf = cbrt(Num128ToDouble(N));
	if (static_cast<double>(D) > rf * rf * 3.03)
	{
		return false;
	}

	// r = floor(N^(1/3)) (N < 2^78: r < 2^26)
	num64 r = static_cast<num64>(rf);
	while ((num128(r * r) * r) > N)
	{
		--r;
	}
	while ((num128((r + 1) * (r + 1)) * (r + 1)) <= N)
	{
		++r;
	}
	const num64 r2 = r * r;
	if (D >= r2 * 3)
	{
		return false;
	}

	// N = r^3: D = r + r^2 (a cube has no other of these forms)
	if ((num128(r2) * r) == N)
	{
		if ((D == r + r2) && IsPrime(r))
		{
			NumberFound(n1);
		}
		return true;
	}

	// N = p^2 (D = p) or N = p * q (D = p + q, p < q): if N has such a pair of divisors, it is not r^2 s with D = (r + 1)(r + s)
	{
		const num128 DD = num128(D) * D;
		if (DD == N)
		{
			if (IsPrime(D))
			{
				NumberFound(n1);
			}
			return true;
		}
		// (q - p)^2 = D^2 - 4N
		const num128 N4 = N << 2;
		if (DD > N4)
		{
			const num128 E = DD - N4;
			const num64 sq = IntegerSquareRoot(E);
			if ((num128(sq) * sq) == E)
			{
				const num64 p = (D - sq) >> 1;
				const num64 q = (D + sq) >> 1;
				if (IsPrime(p) && IsPrime(q))
				{
					NumberFound(n1);
				}
				return true;
			}
		}
	}

	// N = r^2 s (D >= 1.8899 N^(2/3)): D = f(r) = (r + 1)(r + N / r^2). f is convex for r > 0, with its minimum near (N / 2)^(1/3), and
	// D < 3 N^(2/3) puts both solutions in (0.34, 1.54) N^(1/3). Exact bisections (R2S_Sign: the sign of f(r) - D, R2S_Rising: f'(r) > 0)
	// find the integers next to them.
	if (D * 25 >= r2 * 47)
	{
		const num64 lo = std::max<num64>((r * 34) / 100, 2);
		const num64 hi = (r * 154) / 100 + 2;
		// the minimum: the first x in [lo, hi] with f'(x) > 0
		num64 a = lo;
		num64 b = hi;
		while (a < b)
		{
			const num64 m = (a + b) >> 1;
			if (R2S_Rising(m, N))
			{
				b = m;
			}
			else
			{
				a = m + 1;
			}
		}
		const num64 r0 = a;
		num64 candidates[2];
		// f decreases in [lo, r0 - 1] (r0 - 1 <= the real minimum): the last x with f(x) >= D
		a = lo;
		b = (r0 > lo) ? (r0 - 1) : lo;
		while (a < b)
		{
			const num64 m = (a + b + 1) >> 1;
			if (R2S_Sign(m, N, D) >= 0)
			{
				a = m;
			}
			else
			{
				b = m - 1;
			}
		}
		candidates[0] = a;
		// f increases in [r0, hi]: the first x with f(x) >= D
		a = r0;
		b = hi;
		while (a < b)
		{
			const num64 m = (a + b) >> 1;
			if (R2S_Sign(m, N, D) >= 0)
			{
				b = m;
			}
			else
			{
				a = m + 1;
			}
		}
		candidates[1] = a;
		for (int j = 0; j < 2; ++j)
		{
			const num64 x = candidates[j];
			if (R2S_Sign(x, N, D) != 0)
			{
				continue;
			}
			// D = (x + 1)(x + s), N = x^2 s
			const num64 s = D / (x + 1) - x;
			if ((num128(x + 1) * (x + s) == D) && (num128(x * x) * s == N) && (s != x))
			{
				if (IsPrime(x) && IsPrime(s))
				{
					NumberFound(n1);
				}
				return true;
			}
		}
	}

	// none of the forms fits
	return true;
}

// ---- trial division of CheckCofactor in doubles ----
// n < 2^78 = n2 2^52 + n1 2^26 + n0 (26-bit limbs) is congruent to x = n2 (2^52 mod p) + n1 (2^26 mod p) + n0 < 2^52 + 2^26 (p < 2^25), an
// exact double; r = x - floor(x * (1 / p)) p is exact and in (-p, 2p) (the floor is off by at most one), and p | n <=> |r| is 0 or p.
// The tables hold the first TD_SIZE primes of PrimeInverses128 (index 0: 3). 8 primes per AVX-512 instruction, 4 with AVX2 + FMA,
// instead of a 128-bit product per prime (k4: 90% of the time was this loop).
enum : unsigned int { TD_SIZE = ReciprocalsTableSize128 }; // (all primes of the 128-bit tables: up to 2.15e7 < 2^25)
static double* g_tdC2 = nullptr;
static double* g_tdC1 = nullptr;
static double* g_tdP = nullptr;
static double* g_tdInvP = nullptr;
static int g_tdLevel = 0; // 0: scalar, 1: AVX2 + FMA, 2: AVX-512

#if defined(_MSC_VER)
#define ENGINE_TARGET_AVX2
#define ENGINE_TARGET_AVX512
#else
#define ENGINE_TARGET_AVX2 __attribute__((target("avx2,fma")))
#define ENGINE_TARGET_AVX512 __attribute__((target("avx512f")))
#endif

void CheckCofactorInit(int simdLevel)
{
	if (g_tdP)
	{
		return;
	}
	g_tdC2 = static_cast<double*>(AllocateSystemMemory(TD_SIZE * sizeof(double) * 4, false));
	g_tdC1 = g_tdC2 + TD_SIZE;
	g_tdP = g_tdC1 + TD_SIZE;
	g_tdInvP = g_tdP + TD_SIZE;
	for (unsigned int i = 0; i < TD_SIZE; ++i)
	{
		// the prime: the inverse of its 128-bit inverse (mod 2^64)
		const num64 inv = LowWord(PrimeInverses128[i]);
		num64 p = inv;
		for (int k = 0; k < 5; ++k)
		{
			p *= 2 - inv * p;
		}
		const num64 c1 = (num64(1) << 26) % p;
		const num64 c2 = (c1 * c1) % p;
		g_tdC2[i] = static_cast<double>(c2);
		g_tdC1[i] = static_cast<double>(c1);
		g_tdP[i] = static_cast<double>(p);
		g_tdInvP[i] = 1.0 / static_cast<double>(p);
	}
	if ((g_tdP[TD_SIZE - 1] >= 33554432.0) || (simdLevel < 0))
	{
		std::cerr << "CheckCofactorInit: bad table" << std::endl;
		abort();
	}
	int level = simdLevel;
#if AMICABLE_X64
	if (level == 1)
	{
		// AVX2 without FMA: scalar
#if defined(_MSC_VER)
		int r[4];
		__cpuid(r, 1);
		if (((r[2] >> 12) & 1) == 0)
		{
			level = 0;
		}
#else
		if (!__builtin_cpu_supports("fma"))
		{
			level = 0;
		}
#endif
	}
#else
	// (the portable trial division)
	level = 0;
#endif
	g_tdLevel = level;
}

static FORCEINLINE unsigned int FindDivisorScalar(const num128& n, unsigned int i0, unsigned int i1)
{
	for (unsigned int i = i0; i < i1; ++i)
	{
		if (n * PrimeInverses128[i] <= n)
		{
			return i;
		}
	}
	return i1;
}

#if AMICABLE_X64
// The first index j in [i0, i1) (i1 <= TD_SIZE) whose prime divides n (n < 2^78), or i1
static ENGINE_TARGET_AVX512 unsigned int FindDivisorAVX512(const num128& n, unsigned int i0, unsigned int i1)
{
	const num64 lo = LowWord(n);
	const num64 hi = HighWord(n);
	const __m512d n0 = _mm512_set1_pd(static_cast<double>(lo & 0x3FFFFFF));
	const __m512d n1 = _mm512_set1_pd(static_cast<double>((lo >> 26) & 0x3FFFFFF));
	const __m512d n2 = _mm512_set1_pd(static_cast<double>((lo >> 52) | (hi << 12)));
	const __m512d zero = _mm512_setzero_pd();
	unsigned int i = i0 & ~7U;
	unsigned int first = 0xFFU << (i0 - i);
	for (; i < i1; i += 8)
	{
		const __m512d p = _mm512_load_pd(g_tdP + i);
		const __m512d x = _mm512_fmadd_pd(n2, _mm512_load_pd(g_tdC2 + i), _mm512_fmadd_pd(n1, _mm512_load_pd(g_tdC1 + i), n0));
		const __m512d q = _mm512_roundscale_pd(_mm512_mul_pd(x, _mm512_load_pd(g_tdInvP + i)), _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		const __m512d r = _mm512_abs_pd(_mm512_fnmadd_pd(q, p, x));
		const unsigned int k = (static_cast<unsigned int>(_mm512_cmp_pd_mask(r, zero, _CMP_EQ_OQ)) | static_cast<unsigned int>(_mm512_cmp_pd_mask(r, p, _CMP_EQ_OQ))) & first;
		first = 0xFF;
		if (k)
		{
			unsigned long t;
			_BitScanForward(&t, k);
			_mm256_zeroupper();
			return ((i + t) < i1) ? (i + t) : i1;
		}
	}
	_mm256_zeroupper();
	return i1;
}

static ENGINE_TARGET_AVX2 unsigned int FindDivisorAVX2(const num128& n, unsigned int i0, unsigned int i1)
{
	const num64 lo = LowWord(n);
	const num64 hi = HighWord(n);
	const __m256d n0 = _mm256_set1_pd(static_cast<double>(lo & 0x3FFFFFF));
	const __m256d n1 = _mm256_set1_pd(static_cast<double>((lo >> 26) & 0x3FFFFFF));
	const __m256d n2 = _mm256_set1_pd(static_cast<double>((lo >> 52) | (hi << 12)));
	const __m256d zero = _mm256_setzero_pd();
	const __m256d signBit = _mm256_set1_pd(-0.0);
	unsigned int i = i0 & ~3U;
	unsigned int first = 0xFU << (i0 - i);
	for (; i < i1; i += 4)
	{
		const __m256d p = _mm256_load_pd(g_tdP + i);
		const __m256d x = _mm256_fmadd_pd(n2, _mm256_load_pd(g_tdC2 + i), _mm256_fmadd_pd(n1, _mm256_load_pd(g_tdC1 + i), n0));
		const __m256d q = _mm256_round_pd(_mm256_mul_pd(x, _mm256_load_pd(g_tdInvP + i)), _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		const __m256d r = _mm256_andnot_pd(signBit, _mm256_fnmadd_pd(q, p, x));
		const unsigned int k = static_cast<unsigned int>(_mm256_movemask_pd(_mm256_or_pd(_mm256_cmp_pd(r, zero, _CMP_EQ_OQ), _mm256_cmp_pd(r, p, _CMP_EQ_OQ)))) & first;
		first = 0xF;
		if (k)
		{
			unsigned long t;
			_BitScanForward(&t, k);
			_mm256_zeroupper();
			return ((i + t) < i1) ? (i + t) : i1;
		}
	}
	_mm256_zeroupper();
	return i1;
}

#endif

// The first index j in [i0, i1) whose prime (PrimeInverses128) divides n, or i1
static FORCEINLINE unsigned int FindDivisor(const num128& n, unsigned int i0, unsigned int i1)
{
#if AMICABLE_X64
	if ((g_tdLevel == 0) || (HighWord(n) >= (num64(1) << 14)) || (i0 >= TD_SIZE))
	{
		return FindDivisorScalar(n, i0, i1);
	}
	const unsigned int e = (i1 < TD_SIZE) ? i1 : static_cast<unsigned int>(TD_SIZE);
	const unsigned int j = (g_tdLevel == 2) ? FindDivisorAVX512(n, i0, e) : FindDivisorAVX2(n, i0, e);
	if ((j < e) || (e == i1))
	{
		return j;
	}
	return FindDivisorScalar(n, e, i1);
#else
	return FindDivisorScalar(n, i0, i1);
#endif
}

// The trial division of InitialCheck128 stops after the block of 16 primes b if n2 + n2 * SumEstimates128[b] < targetSum (the
// bound for the factors >= the prime 16 (b + 1); it decreases with b): the first such block end after i (binary search), or
// ReciprocalsTableSize128
static FORCEINLINE unsigned int SumEstimateStop(const num128& n2, const num128& targetSum, unsigned int i)
{
	unsigned int lo = i >> 4;
	unsigned int hi = ReciprocalsTableSize128 / 16;
	while (lo < hi)
	{
		const unsigned int mid = (lo + hi) >> 1;
		if (n2 + ((n2 * SumEstimates128[mid]) >> (64 - SumEstimates128Shift)) < targetSum)
		{
			hi = mid;
		}
		else
		{
			lo = mid + 1;
		}
	}
	return (lo < ReciprocalsTableSize128 / 16) ? (16 * (lo + 1)) : static_cast<unsigned int>(ReciprocalsTableSize128);
}

// The first index >= i of the double tables whose prime is > x (x below the largest prime of the tables): galloping search (the
// index is usually close to i), then bisection
static FORCEINLINE unsigned int FirstPrimeAbove(const double x, unsigned int i)
{
	unsigned int lo = i;
	unsigned int step = 1;
	while ((lo + step < TD_SIZE - 1) && (g_tdP[lo + step] <= x))
	{
		lo += step;
		step <<= 1;
	}
	unsigned int hi = (lo + step < TD_SIZE - 1) ? (lo + step) : (TD_SIZE - 1);
	while (lo < hi)
	{
		const unsigned int mid = (lo + hi) >> 1;
		if (g_tdP[mid] > x)
		{
			hi = mid;
		}
		else
		{
			lo = mid + 1;
		}
	}
	return lo;
}

// If N is the other member of a pair with D = targetSum - N - 1 > 0, every prime factor of N is > L = N / D (N / r is a divisor
// between 1 and N), N has at most K = floor(ln N / ln L) of them, and the smallest one is <= K L + K + 2 (GPU version: DecideByBand).
// bandHi: the index of the first prime above that bound (in the double tables), ReciprocalsTableSize128 if there is no bound (D = 0,
// L < 2, or beyond the tables); bandLo: the index of the first prime above L (i if L is below the prime with index i)
// (Lout, Eout: L and the bound, Eout = infinity if there is none)
// BandBounds: L and that bound E alone (false if there is none)
static FORCEINLINE bool BandBounds(const num128& n2, const num128& targetSum, double& Lout, double& Eout)
{
	const num128 D = targetSum - n2 - 1;
	if (D == 0)
	{
		return false;
	}
	const double N = Num128ToDouble(n2);
	const double L = N / Num128ToDouble(D);
	if (!(L >= 2.0) || (L >= 9.0e18))
	{
		return false;
	}
	// K <= bitlen(N) / floor(log2 L) (bitlen(N) > log2 N): an upper bound without logarithms (a larger K only widens the band)
	unsigned long t;
	unsigned int bitsN;
	if (HighWord(n2))
	{
		_BitScanReverse64(&t, HighWord(n2));
		bitsN = static_cast<unsigned int>(t) + 65;
	}
	else
	{
		_BitScanReverse64(&t, LowWord(n2));
		bitsN = static_cast<unsigned int>(t) + 1;
	}
	_BitScanReverse64(&t, static_cast<num64>(L));
	unsigned int K = bitsN / static_cast<unsigned int>(t);
	if (K < 2)
	{
		K = 2;
	}
	Lout = L;
	Eout = (K * L + K + 2.0) * (1.0 + 1e-9) + 2.0;
	return true;
}

static FORCEINLINE void BandRange(const num128& n2, const num128& targetSum, unsigned int i, unsigned int& bandLo, unsigned int& bandHi, double& Lout, double& Eout)
{
	bandLo = i;
	bandHi = ReciprocalsTableSize128;
	Lout = 0.0;
	Eout = std::numeric_limits<double>::infinity();
	double L, E;
	if ((g_tdP == nullptr) || (i >= TD_SIZE) || !BandBounds(n2, targetSum, L, E))
	{
		return;
	}
	Lout = L;
	Eout = E;
	if (E < g_tdP[i])
	{
		// (most numbers: no prime of the band is left)
		bandHi = i;
		return;
	}
	if (!(E < g_tdP[TD_SIZE - 1]))
	{
		return;
	}
	bandHi = FirstPrimeAbove(E, i);
	// (L slightly too small only makes the band longer)
	const double Llo = L * (1.0 - 1e-9) - 1.0;
	if (Llo > g_tdP[i])
	{
		bandLo = FirstPrimeAbove(Llo, i);
	}
}

// The first index >= i of the double tables whose prime p has p^3 > n2 (exact), ReciprocalsTableSize128 if it is beyond them
static FORCEINLINE unsigned int CubeStop(const num128& n2, unsigned int i)
{
	if ((g_tdP == nullptr) || (i >= TD_SIZE))
	{
		return ReciprocalsTableSize128;
	}
	const num64 pmax = static_cast<num64>(g_tdP[TD_SIZE - 1]);
	if ((num128(pmax * pmax) * pmax) <= n2)
	{
		return ReciprocalsTableSize128;
	}
	unsigned int lo = i;
	unsigned int hi = TD_SIZE - 1;
	while (lo < hi)
	{
		const unsigned int mid = (lo + hi) >> 1;
		const num64 p = static_cast<num64>(g_tdP[mid]);
		if ((num128(p * p) * p) > n2)
		{
			hi = mid;
		}
		else
		{
			lo = mid + 1;
		}
	}
	return lo;
}

// ---- band sieve (GPU version: SieveBand) ----
static const signed char BAND_CLS30[30] = { -1, 0, -1, -1, -1, -1, -1, 1, -1, -1, -1, 2, -1, 3, -1, -1, -1, 4, -1, 5, -1, -1, -1, 6, -1, -1, -1, -1, -1, 7 };
// The multiples d0 + k p prime to 30 (p prime to 30): d0 + k p == p (m + k) (mod 30) with m = d0 p^-1 mod 30, so k runs through
// u - m for the units u = 1, 7, 11, 13, 17, 19, 23, 29 (+ 30 t) above m. BAND_INV30: u^-1 mod 30; BAND_FIRST30_J / _K: the index of
// the first unit >= m and its distance from m; BAND_GAP30: the gaps between the units.
static const unsigned char BAND_INV30[30] = { 0, 1, 0, 0, 0, 0, 0, 13, 0, 0, 0, 11, 0, 7, 0, 0, 0, 23, 0, 19, 0, 0, 0, 17, 0, 0, 0, 0, 0, 29 };
static const unsigned char BAND_FIRST30_J[30] = { 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7 };
static const unsigned char BAND_FIRST30_K[30] = { 1, 0, 5, 4, 3, 2, 1, 0, 3, 2, 1, 0, 1, 0, 3, 2, 1, 0, 1, 0, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0 };
static const unsigned char BAND_GAP30[8] = { 6, 4, 2, 4, 2, 4, 6, 2 };
static THREAD_LOCAL num64 t_bandStats[3];

struct BandRootCache
{
	enum : unsigned int { NO_ROOT = 0xFFFFFFFFU };
	num128 A1 = 0;
	num128 B1 = 0;
	unsigned int iLo = 0;
	unsigned int iHi = 0;
	std::vector<unsigned int> roots; // (of the primes iLo .. iHi - 1)
	std::vector<unsigned int> offsets;
};
static THREAD_LOCAL BandRootCache* t_bandRoots = nullptr;
#ifdef RANGE_SIEVE_VERIFY_BAND
std::atomic<num64> g_bandVerifyCount(0);
std::atomic<num64> g_bandVerifyErrors(0);
#endif

double CofactorTableMaxPrime()
{
	return g_tdP ? g_tdP[TD_SIZE - 1] : 0.0;
}

double CofactorTablePrimeCount(double lo, double hi)
{
	if ((g_tdP == nullptr) || !(hi > lo))
	{
		return 0.0;
	}
	const unsigned int iLo = (lo < g_tdP[0]) ? 0 : FirstPrimeAbove(lo, 0);
	const unsigned int iHi = (hi < g_tdP[0]) ? 0 : FirstPrimeAbove(hi, iLo);
	return static_cast<double>(iHi - iLo);
}

void CofactorBandStats(num64* out3)
{
	for (int k = 0; k < 3; ++k)
	{
		out3[k] = t_bandStats[k];
		t_bandStats[k] = 0;
	}
}

// x mod p (exact) of the prime with index i of the double tables, x < 2^78 (as in FindDivisor)
static FORCEINLINE num64 ModTable(const num128& x, unsigned int i)
{
	const num64 lo = LowWord(x);
	const num64 hi = HighWord(x);
	const double y = static_cast<double>((lo >> 52) | (hi << 12)) * g_tdC2[i] + static_cast<double>((lo >> 26) & 0x3FFFFFF) * g_tdC1[i] + static_cast<double>(lo & 0x3FFFFFF);
	const double p = g_tdP[i];
	double r = y - floor(y * g_tdInvP[i]) * p;
	if (r < 0.0)
	{
		r += p;
	}
	else if (r >= p)
	{
		r -= p;
	}
	return static_cast<num64>(r);
}

// a^-1 mod m (0 < a < m < 2^32, gcd(a, m) = 1)
static FORCEINLINE num64 InvMod32Engine(num64 a, num64 m)
{
	long long t = 0, nt = 1;
	unsigned int r = static_cast<unsigned int>(m), nr = static_cast<unsigned int>(a);
	while (nr)
	{
		const unsigned int q = r / nr;
		const long long tmp = t - static_cast<long long>(q) * nt;
		t = nt;
		nt = tmp;
		const unsigned int tmp2 = r - q * nr;
		r = nr;
		nr = tmp2;
	}
	return static_cast<num64>((t < 0) ? (t + static_cast<long long>(m)) : t);
}

// ---- roots of the band primes: -B1 / A1 mod p, NO_ROOT if p | A1 (gcd(A1, B1) = 1), A1, B1 < 2^78 ----
// SIMD: a^-1 = a^(p - 2) mod p (Fermat, right-to-left binary powering) in doubles, for primes p >= 1024. The residues
// h - round(h / p) p (exact with FMA) are within 0.504 p (round may be off by one near a half-integer: the error of h * (1 / p) is
// below 2^-8 for h < 2^53), within p/2 + 1 for the products h < 2^49 of the powering: exact for p < 2^25.
static void BandRootsScalar(const num128& A1, const num128& B1, unsigned int i0, unsigned int i1, unsigned int* out)
{
	for (unsigned int i = i0; i < i1; ++i)
	{
		const num64 p = static_cast<num64>(g_tdP[i]);
		const num64 a = ModTable(A1, i);
		unsigned int root = BandRootCache::NO_ROOT;
		if (a != 0)
		{
			const num64 b = ModTable(B1, i);
			root = static_cast<unsigned int>(((b ? (p - b) : 0) * InvMod32Engine(a, p)) % p);
		}
		out[i - i0] = root;
	}
}

// (root - P0) mod p, NO_ROOT stays
static void BandOffsetsScalar(const num64 P0, unsigned int i0, unsigned int i1, const unsigned int* roots, unsigned int* out)
{
	const num128 P0w = P0;
	for (unsigned int i = i0; i < i1; ++i)
	{
		const unsigned int root = roots[i - i0];
		if (root == BandRootCache::NO_ROOT)
		{
			out[i - i0] = root;
			continue;
		}
		const num64 p = static_cast<num64>(g_tdP[i]);
		const num64 P0p = ModTable(P0w, i);
		out[i - i0] = static_cast<unsigned int>((root >= P0p) ? (root - P0p) : (root + p - P0p));
	}
}

#if AMICABLE_X64
#define BAND_RN (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)

// x y mod p within p/2 + 1 (|x|, |y| < 0.504 p)
static ENGINE_TARGET_AVX512 FORCEINLINE __m512d BandMulMod512(const __m512d x, const __m512d y, const __m512d p, const __m512d ip)
{
	const __m512d h = _mm512_mul_pd(x, y);
	return _mm512_fnmadd_pd(_mm512_roundscale_pd(_mm512_mul_pd(h, ip), BAND_RN), p, h);
}

// x mod p within 0.504 p, x = x2 2^52 + x1 2^26 + x0 < 2^78
static ENGINE_TARGET_AVX512 FORCEINLINE __m512d BandReduce512(const __m512d x2, const __m512d x1, const __m512d x0, const unsigned int i, const __m512d p, const __m512d ip)
{
	const __m512d x = _mm512_fmadd_pd(x2, _mm512_loadu_pd(g_tdC2 + i), _mm512_fmadd_pd(x1, _mm512_loadu_pd(g_tdC1 + i), x0));
	return _mm512_fnmadd_pd(_mm512_roundscale_pd(_mm512_mul_pd(x, ip), BAND_RN), p, x);
}

// (4 independent blocks of 8 primes per step of the powering: its latency chain is ~25 multiplications long)
#define BAND_R512_INIT(j) \
	const __m512d p##j = _mm512_loadu_pd(g_tdP + i + 8 * j); \
	const __m512d ip##j = _mm512_loadu_pd(g_tdInvP + i + 8 * j); \
	const __m512d a##j = BandReduce512(A2, A1v, A0, i + 8 * j, p##j, ip##j); \
	const __m512i e##j = _mm512_cvtepi32_epi64(_mm512_cvttpd_epi32(_mm512_sub_pd(p##j, two))); \
	__m512d res##j = one; \
	__m512d base##j = a##j;
#define BAND_R512_STEP(j) \
	res##j = _mm512_mask_mov_pd(res##j, _mm512_test_epi64_mask(e##j, bit), BandMulMod512(res##j, base##j, p##j, ip##j)); \
	base##j = BandMulMod512(base##j, base##j, p##j, ip##j);
#define BAND_R512_DONE(j) \
	{ \
		const __m512d b##j = BandReduce512(B2, B1v, B0, i + 8 * j, p##j, ip##j); \
		__m512d r##j = BandMulMod512(_mm512_sub_pd(zero, b##j), res##j, p##j, ip##j); \
		r##j = _mm512_mask_add_pd(r##j, _mm512_cmp_pd_mask(r##j, zero, _CMP_LT_OQ), r##j, p##j); \
		r##j = _mm512_mask_mov_pd(r##j, _mm512_cmp_pd_mask(a##j, zero, _CMP_EQ_OQ), noRoot); \
		_mm256_storeu_si256(reinterpret_cast<__m256i*>(out + (i - i0) + 8 * j), _mm512_cvttpd_epi32(r##j)); \
	}

static ENGINE_TARGET_AVX512 void BandRootsAVX512(const num128& A1, const num128& B1, unsigned int i0, unsigned int i1, unsigned int* out)
{
	const __m512d A0 = _mm512_set1_pd(static_cast<double>(LowWord(A1) & 0x3FFFFFF));
	const __m512d A1v = _mm512_set1_pd(static_cast<double>((LowWord(A1) >> 26) & 0x3FFFFFF));
	const __m512d A2 = _mm512_set1_pd(static_cast<double>((LowWord(A1) >> 52) | (HighWord(A1) << 12)));
	const __m512d B0 = _mm512_set1_pd(static_cast<double>(LowWord(B1) & 0x3FFFFFF));
	const __m512d B1v = _mm512_set1_pd(static_cast<double>((LowWord(B1) >> 26) & 0x3FFFFFF));
	const __m512d B2 = _mm512_set1_pd(static_cast<double>((LowWord(B1) >> 52) | (HighWord(B1) << 12)));
	const __m512d zero = _mm512_setzero_pd();
	const __m512d one = _mm512_set1_pd(1.0);
	const __m512d two = _mm512_set1_pd(2.0);
	const __m512d noRoot = _mm512_set1_pd(-1.0); // (converts to NO_ROOT)
	unsigned int i = i0;
	if (i1 - i0 >= 8)
	{
		// the bits of the largest exponent
		unsigned long t;
		_BitScanReverse64(&t, static_cast<num64>(g_tdP[i1 - 1]) - 2);
		const unsigned int nbits = static_cast<unsigned int>(t) + 1;
		for (; i + 32 <= i1; i += 32)
		{
			BAND_R512_INIT(0)
			BAND_R512_INIT(1)
			BAND_R512_INIT(2)
			BAND_R512_INIT(3)
			__m512i bit = _mm512_set1_epi64(1);
			for (unsigned int k = 0; k < nbits; ++k)
			{
				BAND_R512_STEP(0)
				BAND_R512_STEP(1)
				BAND_R512_STEP(2)
				BAND_R512_STEP(3)
				bit = _mm512_add_epi64(bit, bit);
			}
			BAND_R512_DONE(0)
			BAND_R512_DONE(1)
			BAND_R512_DONE(2)
			BAND_R512_DONE(3)
		}
		for (; i + 8 <= i1; i += 8)
		{
			BAND_R512_INIT(0)
			__m512i bit = _mm512_set1_epi64(1);
			for (unsigned int k = 0; k < nbits; ++k)
			{
				BAND_R512_STEP(0)
				bit = _mm512_add_epi64(bit, bit);
			}
			BAND_R512_DONE(0)
		}
	}
	_mm256_zeroupper();
	BandRootsScalar(A1, B1, i, i1, out + (i - i0));
}

static ENGINE_TARGET_AVX512 void BandOffsetsAVX512(const num64 P0, unsigned int i0, unsigned int i1, const unsigned int* roots, unsigned int* out)
{
	const __m512d P0v = _mm512_set1_pd(static_cast<double>(P0 & 0x3FFFFFF));
	const __m512d P1v = _mm512_set1_pd(static_cast<double>((P0 >> 26) & 0x3FFFFFF));
	const __m512d P2v = _mm512_set1_pd(static_cast<double>(P0 >> 52));
	const __m512d zero = _mm512_setzero_pd();
	const __m512d noRoot = _mm512_set1_pd(-1.0);
	unsigned int i = i0;
	for (; i + 8 <= i1; i += 8)
	{
		const __m512d p = _mm512_loadu_pd(g_tdP + i);
		const __m512d ip = _mm512_loadu_pd(g_tdInvP + i);
		// (NO_ROOT: -1)
		const __m512d root = _mm512_cvtepi32_pd(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(roots + (i - i0))));
		__m512d d = _mm512_sub_pd(root, BandReduce512(P2v, P1v, P0v, i, p, ip));
		d = _mm512_mask_add_pd(d, _mm512_cmp_pd_mask(d, zero, _CMP_LT_OQ), d, p);
		d = _mm512_mask_sub_pd(d, _mm512_cmp_pd_mask(d, p, _CMP_GE_OQ), d, p);
		d = _mm512_mask_mov_pd(d, _mm512_cmp_pd_mask(root, zero, _CMP_LT_OQ), noRoot);
		_mm256_storeu_si256(reinterpret_cast<__m256i*>(out + (i - i0)), _mm512_cvttpd_epi32(d));
	}
	_mm256_zeroupper();
	BandOffsetsScalar(P0, i, i1, roots + (i - i0), out + (i - i0));
}

static ENGINE_TARGET_AVX2 FORCEINLINE __m256d BandMulMod256(const __m256d x, const __m256d y, const __m256d p, const __m256d ip)
{
	const __m256d h = _mm256_mul_pd(x, y);
	return _mm256_fnmadd_pd(_mm256_round_pd(_mm256_mul_pd(h, ip), BAND_RN), p, h);
}

static ENGINE_TARGET_AVX2 FORCEINLINE __m256d BandReduce256(const __m256d x2, const __m256d x1, const __m256d x0, const unsigned int i, const __m256d p, const __m256d ip)
{
	const __m256d x = _mm256_fmadd_pd(x2, _mm256_loadu_pd(g_tdC2 + i), _mm256_fmadd_pd(x1, _mm256_loadu_pd(g_tdC1 + i), x0));
	return _mm256_fnmadd_pd(_mm256_round_pd(_mm256_mul_pd(x, ip), BAND_RN), p, x);
}

#define BAND_R256_INIT(j) \
	const __m256d p##j = _mm256_loadu_pd(g_tdP + i + 4 * j); \
	const __m256d ip##j = _mm256_loadu_pd(g_tdInvP + i + 4 * j); \
	const __m256d a##j = BandReduce256(A2, A1v, A0, i + 4 * j, p##j, ip##j); \
	const __m256i e##j = _mm256_cvtepi32_epi64(_mm256_cvttpd_epi32(_mm256_sub_pd(p##j, two))); \
	__m256d res##j = one; \
	__m256d base##j = a##j;
#define BAND_R256_STEP(j) \
	res##j = _mm256_blendv_pd(res##j, BandMulMod256(res##j, base##j, p##j, ip##j), _mm256_castsi256_pd(_mm256_cmpeq_epi64(_mm256_and_si256(e##j, bit), bit))); \
	base##j = BandMulMod256(base##j, base##j, p##j, ip##j);
#define BAND_R256_DONE(j) \
	{ \
		const __m256d b##j = BandReduce256(B2, B1v, B0, i + 4 * j, p##j, ip##j); \
		__m256d r##j = BandMulMod256(_mm256_sub_pd(zero, b##j), res##j, p##j, ip##j); \
		r##j = _mm256_add_pd(r##j, _mm256_and_pd(_mm256_cmp_pd(r##j, zero, _CMP_LT_OQ), p##j)); \
		r##j = _mm256_blendv_pd(r##j, noRoot, _mm256_cmp_pd(a##j, zero, _CMP_EQ_OQ)); \
		_mm_storeu_si128(reinterpret_cast<__m128i*>(out + (i - i0) + 4 * j), _mm256_cvttpd_epi32(r##j)); \
	}

static ENGINE_TARGET_AVX2 void BandRootsAVX2(const num128& A1, const num128& B1, unsigned int i0, unsigned int i1, unsigned int* out)
{
	const __m256d A0 = _mm256_set1_pd(static_cast<double>(LowWord(A1) & 0x3FFFFFF));
	const __m256d A1v = _mm256_set1_pd(static_cast<double>((LowWord(A1) >> 26) & 0x3FFFFFF));
	const __m256d A2 = _mm256_set1_pd(static_cast<double>((LowWord(A1) >> 52) | (HighWord(A1) << 12)));
	const __m256d B0 = _mm256_set1_pd(static_cast<double>(LowWord(B1) & 0x3FFFFFF));
	const __m256d B1v = _mm256_set1_pd(static_cast<double>((LowWord(B1) >> 26) & 0x3FFFFFF));
	const __m256d B2 = _mm256_set1_pd(static_cast<double>((LowWord(B1) >> 52) | (HighWord(B1) << 12)));
	const __m256d zero = _mm256_setzero_pd();
	const __m256d one = _mm256_set1_pd(1.0);
	const __m256d two = _mm256_set1_pd(2.0);
	const __m256d noRoot = _mm256_set1_pd(-1.0);
	unsigned int i = i0;
	if (i1 - i0 >= 4)
	{
		unsigned long t;
		_BitScanReverse64(&t, static_cast<num64>(g_tdP[i1 - 1]) - 2);
		const unsigned int nbits = static_cast<unsigned int>(t) + 1;
		for (; i + 16 <= i1; i += 16)
		{
			BAND_R256_INIT(0)
			BAND_R256_INIT(1)
			BAND_R256_INIT(2)
			BAND_R256_INIT(3)
			__m256i bit = _mm256_set1_epi64x(1);
			for (unsigned int k = 0; k < nbits; ++k)
			{
				BAND_R256_STEP(0)
				BAND_R256_STEP(1)
				BAND_R256_STEP(2)
				BAND_R256_STEP(3)
				bit = _mm256_add_epi64(bit, bit);
			}
			BAND_R256_DONE(0)
			BAND_R256_DONE(1)
			BAND_R256_DONE(2)
			BAND_R256_DONE(3)
		}
		for (; i + 4 <= i1; i += 4)
		{
			BAND_R256_INIT(0)
			__m256i bit = _mm256_set1_epi64x(1);
			for (unsigned int k = 0; k < nbits; ++k)
			{
				BAND_R256_STEP(0)
				bit = _mm256_add_epi64(bit, bit);
			}
			BAND_R256_DONE(0)
		}
	}
	_mm256_zeroupper();
	BandRootsScalar(A1, B1, i, i1, out + (i - i0));
}

static ENGINE_TARGET_AVX2 void BandOffsetsAVX2(const num64 P0, unsigned int i0, unsigned int i1, const unsigned int* roots, unsigned int* out)
{
	const __m256d P0v = _mm256_set1_pd(static_cast<double>(P0 & 0x3FFFFFF));
	const __m256d P1v = _mm256_set1_pd(static_cast<double>((P0 >> 26) & 0x3FFFFFF));
	const __m256d P2v = _mm256_set1_pd(static_cast<double>(P0 >> 52));
	const __m256d zero = _mm256_setzero_pd();
	const __m256d noRoot = _mm256_set1_pd(-1.0);
	unsigned int i = i0;
	for (; i + 4 <= i1; i += 4)
	{
		const __m256d p = _mm256_loadu_pd(g_tdP + i);
		const __m256d ip = _mm256_loadu_pd(g_tdInvP + i);
		const __m256d root = _mm256_cvtepi32_pd(_mm_loadu_si128(reinterpret_cast<const __m128i*>(roots + (i - i0))));
		__m256d d = _mm256_sub_pd(root, BandReduce256(P2v, P1v, P0v, i, p, ip));
		d = _mm256_add_pd(d, _mm256_and_pd(_mm256_cmp_pd(d, zero, _CMP_LT_OQ), p));
		d = _mm256_sub_pd(d, _mm256_and_pd(_mm256_cmp_pd(d, p, _CMP_GE_OQ), p));
		d = _mm256_blendv_pd(d, noRoot, _mm256_cmp_pd(root, zero, _CMP_LT_OQ));
		_mm_storeu_si128(reinterpret_cast<__m128i*>(out + (i - i0)), _mm256_cvttpd_epi32(d));
	}
	_mm256_zeroupper();
	BandOffsetsScalar(P0, i, i1, roots + (i - i0), out + (i - i0));
}
#endif

bool CofactorBandSetup(CofactorBand& band)
{
	if ((g_tdP == nullptr) || !(band.hi < g_tdP[TD_SIZE - 1]) || !(band.hi > band.lo))
	{
		return false;
	}
	const unsigned int iLo = (band.lo < g_tdP[0]) ? 0 : FirstPrimeAbove(band.lo, 0);
	const unsigned int iHi = FirstPrimeAbove(band.hi, iLo);
	for (unsigned int k = 0; k <= BAND_SUBBANDS; ++k)
	{
		band.sub[k] = iLo + static_cast<unsigned int>((static_cast<num64>(iHi - iLo) * k + BAND_SUBBANDS - 1) / BAND_SUBBANDS);
	}
	return true;
}

unsigned int CofactorBandSubBandOf(const CofactorBand& band, num64 prime)
{
	if ((g_tdP == nullptr) || (band.sub[0] >= band.sub[BAND_SUBBANDS]) || (static_cast<double>(prime) < g_tdP[band.sub[0]]) || !(static_cast<double>(prime) <= g_tdP[band.sub[BAND_SUBBANDS] - 1]))
	{
		return 0;
	}
	const unsigned int i = FirstPrimeAbove(static_cast<double>(prime) - 0.5, band.sub[0]);
	if (static_cast<num64>(g_tdP[i]) != prime)
	{
		return 0;
	}
	unsigned int sb = 1;
	while (i >= band.sub[sb])
	{
		++sb;
	}
	return sb;
}

unsigned int BandSieve(const num128& A1, const num128& B1, num64 P0, num64 p_hi, CofactorBand& band, unsigned char* subBands)
{
	if ((HighWord(A1) >= (num64(1) << 14)) || (HighWord(B1) >= (num64(1) << 14)) || (p_hi - P0 >= (num64(1) << 30)) || !CofactorBandSetup(band))
	{
		return 0;
	}
	const unsigned int iLo = band.sub[0];
	const unsigned int iHi = band.sub[BAND_SUBBANDS];
#if AMICABLE_X64
	const int level = (g_tdP[iLo] >= 1024.0) ? g_tdLevel : 0;
#endif
	// the roots of the band primes: the same for all chunks of a range, kept per thread (also for a part of the cached primes:
	// the band of a chunk moves a little with p_hi)
	if (!t_bandRoots)
	{
		t_bandRoots = new BandRootCache();
	}
	BandRootCache& rc = *t_bandRoots;
	if (!((rc.A1 == A1) && (rc.B1 == B1) && (rc.iLo <= iLo) && (iHi <= rc.iHi)))
	{
		rc.A1 = A1;
		rc.B1 = B1;
		rc.iLo = iLo;
		rc.iHi = iHi;
		rc.roots.resize(iHi - iLo);
#if AMICABLE_X64
		if (level == 2)
		{
			BandRootsAVX512(A1, B1, iLo, iHi, rc.roots.data());
		}
		else if (level == 1)
		{
			BandRootsAVX2(A1, B1, iLo, iHi, rc.roots.data());
		}
		else
#endif
		{
			BandRootsScalar(A1, B1, iLo, iHi, rc.roots.data());
		}
	}
	// the first position >= P0 of each prime (x == root (mod p))
	rc.offsets.resize(iHi - iLo);
	unsigned int* off = rc.offsets.data();
	const unsigned int* roots = rc.roots.data() + (iLo - rc.iLo);
#if AMICABLE_X64
	if (level == 2)
	{
		BandOffsetsAVX512(P0, iLo, iHi, roots, off);
	}
	else if (level == 1)
	{
		BandOffsetsAVX2(P0, iLo, iHi, roots, off);
	}
	else
#endif
	{
		BandOffsetsScalar(P0, iLo, iHi, roots, off);
	}

	// Decreasing sub-bands and every hit stores its sub-band: the last store to a position is from its smallest band prime (no load:
	// the 4-bit version had to test the nibble first). NO_ROOT > span: p divides A1, not B1, never A1 p + B1. Only the positions
	// prime to 30 have a byte (P0 is a multiple of 210): the wheel steps over the others. (32 bits: span < 2^30, p < 2^25)
	const unsigned int span = static_cast<unsigned int>(p_hi - P0);
	for (unsigned int sb = BAND_SUBBANDS; sb >= 1; --sb)
	{
		const unsigned char v = static_cast<unsigned char>(sb);
		for (unsigned int i = band.sub[sb - 1]; i < band.sub[sb]; ++i)
		{
			const unsigned int d0 = off[i - iLo];
			if (d0 > span)
			{
				continue;
			}
			const unsigned int p = static_cast<unsigned int>(g_tdP[i]);
			const unsigned int m = ((d0 % 30) * BAND_INV30[p % 30]) % 30;
			unsigned int j = BAND_FIRST30_J[m];
			unsigned int d = d0 + BAND_FIRST30_K[m] * p;
			while (d <= span)
			{
				// position 8 (d / 30) + class
				const unsigned int q30 = d / 30;
				subBands[q30 * 8 + static_cast<unsigned int>(BAND_CLS30[d - q30 * 30])] = v;
				d += BAND_GAP30[j] * p;
				j = (j + 1) & 7;
			}
		}
	}
	return iHi - iLo;
}

bool CofactorBandBounds(const num128& n2, const num128& targetSum, double& L, double& E)
{
	return BandBounds(n2, targetSum, L, E);
}

num64 g_testBandP = 0;
num64 g_testBandPasses = 0;

// The sub-band of the number of the range at p in the band sieve (0: no prime of the band divides A1 p + B1); n2: its cofactor
// (RANGE_SIEVE_VERIFY_BAND)
static FORCEINLINE unsigned int BandSubBand(const CofactorBand& band, const num64 p, const num128& n2)
{
	const num64 d = p - band.P0;
	const num64 q30 = d / 30;
	const num64 k = q30 * 8 + static_cast<num64>(BAND_CLS30[d - q30 * 30]);
	const unsigned int sb = band.subBands[k];
#ifdef RANGE_SIEVE_VERIFY_BAND
	// (every 1024th lookup: the scalar trial division by all primes of a band takes up to ~10 ms)
	static THREAD_LOCAL unsigned int t_verifyTick = 0;
	if (((t_verifyTick++) & 1023) == 0)
	{
		// the sub-band of the smallest prime of the band that divides n2 (all prime factors of A1 p + B1 above Q are in n2, and g
		// has no prime in the band)
		const unsigned int jt = FindDivisorScalar(n2, band.sub[0], band.sub[BAND_SUBBANDS]);
		unsigned int st = 0;
		if (jt < band.sub[BAND_SUBBANDS])
		{
			st = 1;
			while (jt >= band.sub[st])
			{
				++st;
			}
		}
		if (st != sb)
		{
			g_bandVerifyErrors += 1;
		}
		g_bandVerifyCount += 1;
	}
#else
	(void)n2;
#endif
	return sb;
}

// FinalCheck128 without the report: true if n1 is amicable (n2 is prime, p^2 or p q)
static bool FinalCheck128Decide(const num128& targetSum, const num128& n2)
{
	if (LowWord(targetSum) & 1)
	{
		const num128 p = targetSum - n2 - 1;
		return (p * p == n2) && IsPrime(LowWord(p));
	}
	if (targetSum == n2 + 1)
	{
		// (as in FinalCheck128: no test of the primality of n2 >= 2^64)
		return (HighWord(n2) != 0) || IsPrime(LowWord(n2));
	}
	const num128 B = (targetSum - n2 - 1) >> 1;
	num128 D = B * B;
	if (D > n2)
	{
		D -= n2;
		if (IsPerfectSquareCandidate(LowWord(D)))
		{
			const num64 sqrtD = IntegerSquareRoot(D);
			if (num128(sqrtD) * sqrtD == D)
			{
				const num128 p = B - sqrtD;
				const num128 q = B + sqrtD;
				return (p * q == n2) && IsPrime(LowWord(p)) && IsPrime(LowWord(q));
			}
		}
	}
	return false;
}

NOINLINE void CheckCofactor(const num128 n1, num128 targetSum, num128 n2, const unsigned int startIndex, const CofactorBand* band, num64 p, const num128& targetOrig)
{
	if (n2 == 1)
	{
		if (targetSum == 1)
		{
			NumberFound(n1);
		}
		return;
	}

	// sigma(n2) >= n2 + 1
	if (n2 >= targetSum)
	{
		return;
	}

	// The bound of the smallest prime factor first (as at the start of the loop below): most numbers have no prime left below it, and
	// most numbers of the near-degenerate ranges have no prime of the band (their D ~ 3 N^(2/3) takes DecideAtMost2Factors a cube
	// root and often a square root)
	if ((g_tdP != nullptr) && (startIndex < TD_SIZE))
	{
		double L, E;
		if (BandBounds(n2, targetSum, L, E))
		{
			if (E < g_tdP[startIndex])
			{
				return;
			}
			if (band && (E <= band->hi) && (L >= band->L0))
			{
				if (BandSubBand(*band, p, n2) == 0)
				{
					++t_bandStats[0];
					return;
				}
				if (p == g_testBandP)
				{
					++g_testBandPasses;
				}
			}
		}
	}

	if (DecideAtMost2Factors(n1, targetSum, n2))
	{
		return;
	}

	// Trial division from startIndex up to the first of: the band bound for the smallest prime factor (most numbers stop right away),
	// the block where the sum estimate of InitialCheck128 rejects n2, the cube root of n2 (then n2 is prime, p^2 or p q with primes
	// p, q above it: FinalCheck128). The same loop for 128-bit and 64-bit numbers (before: the 64-bit ones went through the scalar
	// trial division of CheckPairInternal up to n2^(1/3), ~0.4 ms for each near-degenerate number of t1).
	unsigned int i = startIndex;
	bool bandChecked = (band == nullptr);
	// (a factor taken from the band skips the primes of (Q, lo]: a factor there means that n1 is not amicable, so only a report has
	// to be checked again, from scratch)
	bool verifyReports = false;
	auto Report = [&]()
	{
		if (verifyReports)
		{
			CheckPair128NoInline(n1, targetOrig);
		}
		else
		{
			NumberFound(n1);
		}
	};
	for (;;)
	{
		unsigned int bandLo, bandStop;
		double Lb, Eb;
		BandRange(n2, targetSum, i, bandLo, bandStop, Lb, Eb);
		if (bandStop <= i)
		{
			return;
		}

		unsigned int j = ReciprocalsTableSize128;
		if (!bandChecked)
		{
			// (once: the sub-band is about A1 p + B1, whose prime factors > Q were all still in n2 at the first check)
			bandChecked = true;
			if ((Eb <= band->hi) && (Lb >= band->L0))
			{
				const unsigned int sb = BandSubBand(*band, p, n2);
				if (sb == 0)
				{
					// no prime of the band divides N
					++t_bandStats[0];
					return;
				}
				++t_bandStats[1];
				const unsigned int s0 = (band->sub[sb - 1] > i) ? band->sub[sb - 1] : i;
				const unsigned int jb = (s0 < band->sub[sb]) ? FindDivisor(n2, s0, band->sub[sb]) : band->sub[sb];
				if (jb < band->sub[sb])
				{
					if (g_tdP[jb] < Lb)
					{
						// a prime factor below L
						return;
					}
					j = jb;
					verifyReports = true;
				}
			}
			else
			{
				++t_bandStats[2];
			}
		}

		if (j == ReciprocalsTableSize128)
		{
			// (a band that ends soon: the other stops would only end the scan a little earlier)
			unsigned int sumStop = ReciprocalsTableSize128;
			unsigned int cubeStop = ReciprocalsTableSize128;
			unsigned int stop = bandStop;
			if (bandStop - i > 256)
			{
				sumStop = SumEstimateStop(n2, targetSum, i);
				cubeStop = CubeStop(n2, i);
				if (sumStop < stop)
				{
					stop = sumStop;
				}
				if (cubeStop < stop)
				{
					stop = cubeStop;
				}
			}
			// the primes up to L can't divide n2 if it is the other member of a pair: the scan starts above L, and the primes it skipped
			// are only tested when it finds a factor or before a conclusion that needs all of them (the cube root, the end of the tables)
			const unsigned int from = (bandLo < stop) ? bandLo : stop;
			j = (from < stop) ? FindDivisor(n2, from, stop) : stop;
			if ((from > i) && ((j < stop) || (stop != bandStop) || (bandStop >= ReciprocalsTableSize128)))
			{
				if (FindDivisor(n2, i, from) < from)
				{
					// a prime factor up to L
					return;
				}
			}
			if (j >= stop)
			{
				if ((stop == bandStop) && (bandStop < ReciprocalsTableSize128))
				{
					// no prime factor in the band
					return;
				}
				if ((stop == sumStop) && (sumStop < ReciprocalsTableSize128))
				{
					return;
				}
				// no prime factor below the cube root of n2 (or of the search limit): at most 2 are left
				if (FinalCheck128Decide(targetSum, n2))
				{
					Report();
				}
				return;
			}
		}

		// the prime with index j divides n2
		const num128 inv = PrimeInverses128[j];
		num128 q = n2 * inv;
		const InverseData128* data = PowersOfP_128DivisibilityData[j];
		for (;;)
		{
			n2 = q;
			q *= inv;
			if (HighWord(q) > HighWord(n2))
			{
				break;
			}
			++data;
		}

		if (data->shift)
		{
			if (LowWord(targetSum) & data->shift_bits)
			{
				return;
			}
			targetSum >>= data->shift;
		}

		targetSum *= data->inverse;
		if (targetSum > data->max_value)
		{
			return;
		}

		if (n2 == 1)
		{
			if (targetSum == 1)
			{
				Report();
			}
			return;
		}

		if (n2 >= targetSum)
		{
			return;
		}

		i = j + 1;
	}
}

// The mask of the numbers p + 1 (mod 385) that can't give an amicable number M0 p for an abundant M0 (r.sum - r.value >= r.value):
// SearchRange (and the test of all known pairs)
unsigned int RangeOverAbundantMask(const RangeData& r)
{
	unsigned int is_over_abundant_mask = 0;
	if (r.sum - r.value >= r.value)
	{
		const Factor* f = r.factors;
		const int k = r.last_factor_index;
		is_over_abundant_mask |= OverAbundant<5>(f, k, r.value, r.sum, 2 * 5) << 1;
		is_over_abundant_mask |= OverAbundant<7>(f, k, r.value, r.sum, 2 * 7) << 2;
		is_over_abundant_mask |= OverAbundant<11>(f, k, r.value, r.sum, 2 * 11) << 4;
		is_over_abundant_mask |= (((is_over_abundant_mask & 0x06) || OverAbundant<7>(f, k, r.value, r.sum, 2 * 5 * 7)) ? byte(1) : byte(0)) << 3;
		is_over_abundant_mask |= (((is_over_abundant_mask & 0x12) || OverAbundant<11>(f, k, r.value, r.sum, 2 * 5 * 11)) ? byte(1) : byte(0)) << 5;
		is_over_abundant_mask |= (((is_over_abundant_mask & 0x14) || OverAbundant<11>(f, k, r.value, r.sum, 2 * 7 * 11)) ? byte(1) : byte(0)) << 6;
		is_over_abundant_mask |= (((is_over_abundant_mask & 0x7E) || OverAbundant<11>(f, k, r.value, r.sum, 2 * 5 * 7 * 11)) ? byte(1) : byte(0)) << 7;
	}
	return is_over_abundant_mask;
}

NOINLINE num64 SearchRange(const RangeData& r)
{
	num128 prime_limit = (SearchLimit::value - 1) / r.value;
	if (prime_limit > SearchLimit::MainPrimeTableBound)
	{
		prime_limit = SearchLimit::MainPrimeTableBound;
	}

	if (r.sum - r.value < r.value)
	{
		const num128 prime_limit2 = (r.sum - 1) / (r.value * 2 - r.sum);
		if (prime_limit > prime_limit2)
		{
			prime_limit = prime_limit2;
		}
	}
	const unsigned int is_over_abundant_mask = RangeOverAbundantMask(r);

	// Long ranges are much faster with the arithmetic-progression sieve
	if (g_UseRangeSieve && (LowWord(prime_limit) >= r.start_prime))
	{
		const num64 num_primes = RangeSieveCountPrimes(r.start_prime, LowWord(prime_limit));
		if (num_primes >= g_RangeSieveMinPrimes)
		{
			// (long ranges are checked in chunks, shared with idle threads)
			if (!RangeGen::ShareRange(r, r.start_prime, LowWord(prime_limit), is_over_abundant_mask, num_primes))
			{
				SearchRangeSieve(r, r.start_prime, LowWord(prime_limit), is_over_abundant_mask, num_primes);
			}
			return num_primes;
		}
	}

	PrimeIterator q(r.start_prime);
	num64 result = 0;
	const num128 m = r.value;
	const num128 sum_m = r.sum;
	while (q.Get() <= LowWord(prime_limit))
	{
		const num64 prev_q = q.Get();
		++result;
		++q;
		if ((CandidatesDataMask[Mod385(prev_q + 1)] & is_over_abundant_mask) == 0)
		{
			num128 n1 = m * prev_q;
			num128 targetSum = sum_m * (prev_q + 1);
			num128 n2;
			if (InitialCheck128(n1, targetSum, n2))
			{
				CheckPairInternal(n1, LowWord(targetSum), LowWord(targetSum), LowWord(n2), 1);
			}
		}
	}

	return result;
}

NOINLINE num64 SearchRangeSquared(const RangeData& r)
{
	num64 prime_limit = static_cast<num64>(sqrt(Num128ToDouble(SearchLimit::value / r.value))) + 1;
	if (r.sum - r.value < r.value)
	{
		// r.sum * (p^2 + p + 1) = r.value * p^2 * 2
		// (r.sum - r.value * 2) * p^2 + r.sum * (p + 1) = 0
		// (r.value * 2 - r.sum) * p^2 - r.sum * (p + 1) = 0
		// (r.value * 2 - r.sum) * p^2 - r.sum * p - r.sum = 0
		// (r.value * 2 - r.sum) / r.sum * p^2 - p - 1 = 0
		// (r.value * 2 / r.sum - 1) * p^2 - p - 1 = 0
		const double A = Num128ToDouble(r.value * 2 - r.sum) / Num128ToDouble(r.sum);
		const num64 prime_limit2 = static_cast<num64>((1.0 + sqrt(1.0 + 4.0 * A)) / (2.0 * A)) + 1;
		if (prime_limit > prime_limit2)
		{
			prime_limit = prime_limit2;
		}
	}

	PrimeIterator q(r.start_prime);
	num64 result = 0;
	const num128 m = r.value;
	const num128 sum_m = r.sum;
	while (q.Get() < prime_limit)
	{
		const num64 q2 = q.Get() * q.Get();
		CheckPair128NoInline(m * q2, sum_m * (q2 + q.Get() + 1));
		++q;
		++result;
	}

	return result;
}

NOINLINE num64 SearchRangeCubed(const RangeData& r)
{
	PrimeIterator q(r.start_prime);
	num64 result = 0;
	const num128 m = r.value;
	const num128 sum_m = r.sum;
	for (;;)
	{
		const num128 q2 = q.Get() * q.Get();
		const num128 q3 = q2 * q.Get();
		const num128 value = m * q3;
		if (value >= SearchLimit::value)
		{
			break;
		}
		const num128 sum = sum_m * (q3 + q2 + q.Get() + 1);
		if (value >= sum - value)
		{
			break;
		}
		CheckPair128NoInline(value, sum);
		++q;
		++result;
	}
	return result;
}

#include <primesieve/PrimeGenerator.hpp>

namespace primesieve
{
	class PrimeGeneratorLargePrimes final : public PrimeGenerator
	{
	public:
		PrimeGeneratorLargePrimes(PrimeSieve& ps, const PreSieve& preSieve) : PrimeGenerator(ps, preSieve) {}

		FORCEINLINE void Init(const num64 rangeBegin)
		{
			auto it = std::lower_bound(CandidatesData.begin(), CandidatesData.end(), rangeBegin, [](const AmicableCandidate& candidate, num64 k)
			{
				num64 m_hi;
				const num64 m_lo = _umul128(k, candidate.GetValue(), &m_hi);
				return (SearchLimit::value > CombineNum128(m_lo, m_hi));
			});
			last_candidate = CandidatesData.data() + (it - CandidatesData.begin()) - 1;
		}

		NOINLINE virtual void generatePrimes(const byte_t* sieve, uint64_t sieveSize) override
		{
			uint64_t base = getSegmentLow();
			for (uint64_t i = 0; i < sieveSize; i += 8, base += NUMBERS_PER_BYTE * 8)
			{
				uint64_t bits = littleendian_cast<uint64_t>(&sieve[i]); 
				while (bits != 0)
				{
					const num64 curPrime = nextPrime(&bits, base);

					const unsigned int mask = CandidatesDataMask[Mod385(curPrime + 1)];
					while (last_candidate >= CandidatesData.data())
					{
						num64 m_hi;
						const num64 m_lo = _umul128(curPrime, last_candidate->GetValue(), &m_hi);
						if (SearchLimit::value > CombineNum128(m_lo, m_hi))
						{
							break;
						}
						--last_candidate;
					}

					for (const AmicableCandidate* candidate = CandidatesData.data(); candidate <= last_candidate; ++candidate)
					{
						if ((candidate->is_over_abundant_mask & mask) == 0)
						{
							CheckPair128(num128(candidate->GetValue()) * curPrime, num128(candidate->GetSum()) * (curPrime + 1));
						}
					}
				}
			}
		}

	private:
		const AmicableCandidate* last_candidate;

		DISALLOW_COPY_AND_ASSIGN(PrimeGeneratorLargePrimes);
	};
}

NOINLINE void SearchLargePrimes(volatile num64* SharedCounterForSearch, const num64 StartPrime, const num64 PrimeLimit, num64 &sharedCounterValue)
{
	primesieve::PrimeSieve s;
	sharedCounterValue = _InterlockedIncrement(SharedCounterForSearch) - 1;

	while (sharedCounterValue < RangeGen::LargePrimesSplitSize)
	{
		const num64 curRangeBegin = LowWord((num128(PrimeLimit - StartPrime + 1) * sharedCounterValue) >> RangeGen::LargePrimesSplitShift) + StartPrime;
		const num64 curRangeEnd = LowWord((num128(PrimeLimit - StartPrime + 1) * (sharedCounterValue + 1)) >> RangeGen::LargePrimesSplitShift) + StartPrime - 1;

		if (curRangeBegin <= curRangeEnd)
		{
			s.setStart(curRangeBegin);
			s.setStop(curRangeEnd);

			primesieve::PreSieve preSieve(curRangeBegin, curRangeEnd);
			primesieve::PrimeGeneratorLargePrimes primeGen(s, preSieve);
			primeGen.Init(curRangeBegin);

			PrimeIterator p(3);
			while (p.Get() <= preSieve.getMaxPrime())
			{
				++p;
			}
			while (p.Get() <= primeGen.getSqrtStop())
			{
				primeGen.addSievingPrime(static_cast<primesieve::uint_t>(p.Get()));
				++p;
			}

			primeGen.sieve();
		}

		const num64 numFinished = _InterlockedIncrement(SharedCounterForSearch + 1);
		const double f = numFinished / static_cast<double>(RangeGen::LargePrimesSplitSize);
		boinc_fraction_done((f > 1.0) ? 1.0 : f);

		sharedCounterValue = _InterlockedIncrement(SharedCounterForSearch) - 1;
	}
}
