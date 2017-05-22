#include "stdafx.h"
#include "PrimeTables.h"
#include "Engine.h"
#include "RangeGen.h"

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
	return a + p + q + 1;
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
	return a + (a_div_p0 + p0 * p0 + p0) * 2 /*+ 1*/;
	
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

NOINLINE void NumberFound(const num128 n1)
{
	++NumFoundPairs;
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

FORCEINLINE void CheckPairInternal(const num128 n1, const num64 targetSum, num64 n2TargetSum, num64 n2, num64 sum)
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

	num64 numPrimesCheckedSoFar = CompileTimePrimesCount;

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
	num64 n2_sqrt4 = Root4(n2);

	num64 p = CompileTimePrimes<CompileTimePrimesCount>::value;
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

		p += static_cast<unsigned int>(NextPrimeShifts[numPrimesCheckedSoFar] * ShiftMultiplier);

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

			p += static_cast<unsigned int>(NextPrimeShifts[numPrimesCheckedSoFar] * ShiftMultiplier);
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

NOINLINE static void CheckPairInternalNoInline(const num128 n1, const num64 targetSum, num64 n2TargetSum, num64 n2, num64 sum)
{
	CheckPairInternal(n1, targetSum, n2TargetSum, n2, sum);
}

// n2 is either a prime, squared prime, or has two distinct prime factors
NOINLINE void FinalCheck128(const num128 n1, const num128& targetSum, const num128& n2)
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
			// May give some false positives if n2 is composite
			// TODO: check that n2 is prime
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

	if (HighWord(targetSum) == 0)
	{
		return true;
	}

	const std::pair<num128, num128>* prime_inverse = PrimeInverses128;
	const num64* max_sum_ratio = SumEstimates128;
	do
	{
		for (const std::pair<num128, num128>* e = prime_inverse + 16; prime_inverse < e; ++prime_inverse)
		{
			num128 q = n2 * prime_inverse->first;
			if (q > prime_inverse->second)
			{
				continue;
			}

			const InverseData128* data = PowersOfP_128DivisibilityData[prime_inverse - PrimeInverses128];
			for (;;)
			{
				n2 = q;
				q *= prime_inverse->first;
				if (q > prime_inverse->second)
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

		if (n2 + HighWord(n2 * (*max_sum_ratio)) < targetSum)
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

NOINLINE num64 SearchRange(const RangeData& r)
{
	num128 prime_limit = (SearchLimit::value - 1) / r.value;
	if (prime_limit > SearchLimit::MainPrimeTableBound)
	{
		prime_limit = SearchLimit::MainPrimeTableBound;
	}

	unsigned int is_over_abundant_mask = 0;
	if (r.sum - r.value < r.value)
	{
		const num128 prime_limit2 = (r.sum - 1) / (r.value * 2 - r.sum);
		if (prime_limit > prime_limit2)
		{
			prime_limit = prime_limit2;
		}
	}
	else
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
			CheckPair128NoInline(m * prev_q, sum_m * (prev_q + 1));
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
		CheckPair128NoInline(value, sum_m * (q3 + q2 + q.Get() + 1));
		++q;
		++result;
	}
	return result;
}

#include <primesieve/SieveOfEratosthenes-inline.hpp>

namespace primesieve
{
	class PrimeFinderLargePrimes final : public PrimeFinder
	{
	public:
		PrimeFinderLargePrimes(PrimeSieve& ps, const PreSieve& preSieve) : PrimeFinder(ps, preSieve) {}

		FORCEINLINE void Init(const num64 rangeBegin)
		{
			auto it = std::lower_bound(CandidatesData.begin(), CandidatesData.end(), rangeBegin, [](const AmicableCandidate& candidate, num64 k)
			{
				num64 h;
				const num64 m = _umul128(k, candidate.value, &h);
				return ((SearchLimit::value > m) && (h == 0));
			});
			last_candidate = CandidatesData.data() + (it - CandidatesData.begin()) - 1;
		}

		NOINLINE virtual void segmentFinished(const byte_t* sieve, uint_t sieveSize)
		{
			uint64_t base = getSegmentLow();
			for (uint_t i = 0; i < sieveSize; i += 8, base += NUMBERS_PER_BYTE * 8)
			{
				uint64_t bits = littleendian_cast<uint64_t>(&sieve[i]); 
				while (bits != 0)
				{
					const num64 curPrime = getNextPrime(&bits, base);

					const unsigned int mask = CandidatesDataMask[Mod385(curPrime + 1)];
					while (last_candidate >= CandidatesData.data())
					{
						num64 h;
						const num64 m = _umul128(curPrime, last_candidate->value, &h);
						if ((SearchLimit::value > m) && (h == 0))
						{
							break;
						}
						--last_candidate;
					}

					for (const AmicableCandidate* candidate = CandidatesData.data(); candidate <= last_candidate; ++candidate)
					{
						if ((candidate->is_over_abundant_mask & mask) == 0)
						{
							CheckPair128NoInline(num128(candidate->value) * curPrime, num128(candidate->sum + candidate->value * 2) * (curPrime + 1));
						}
					}
				}
			}
		}

	private:
		const AmicableCandidate* last_candidate;

		DISALLOW_COPY_AND_ASSIGN(PrimeFinderLargePrimes);
	};
}

NOINLINE void SearchLargePrimes(volatile num64* SharedCounterForSearch, const num64 StartPrime, const num64 PrimeLimit, num64 &sharedCounterValue)
{
	primesieve::PrimeSieve s;
	sharedCounterValue = _InterlockedIncrement(SharedCounterForSearch) - 1;

	while (sharedCounterValue < RangeGen::LargePrimesSplitSize)
	{
		const num64 curRangeBegin = StartPrime + (PrimeLimit - StartPrime + 1) * sharedCounterValue / RangeGen::LargePrimesSplitSize;
		const num64 curRangeEnd = StartPrime + (PrimeLimit - StartPrime + 1) * (sharedCounterValue + 1) / RangeGen::LargePrimesSplitSize - 1;

		if (curRangeBegin <= curRangeEnd)
		{
			s.setStart(curRangeBegin);
			s.setStop(curRangeEnd);

			primesieve::PreSieve preSieve(curRangeBegin, curRangeEnd);
			primesieve::PrimeFinderLargePrimes finder(s, preSieve);
			finder.Init(curRangeBegin);

			unsigned int p = 3;
			const byte* shift = NextPrimeShifts + 1;
			while (p <= preSieve.getLimit())
			{
				p += static_cast<unsigned int>(*shift * ShiftMultiplier);
				++shift;
			}
			while (p <= finder.getSqrtStop())
			{
				finder.addSievingPrime(p);
				p += static_cast<unsigned int>(*shift * ShiftMultiplier);
				++shift;
			}

			finder.sieve();
		}

		const num64 numFinished = _InterlockedIncrement(SharedCounterForSearch + 1);
		const double f = numFinished / static_cast<double>(RangeGen::LargePrimesSplitSize);
		boinc_fraction_done((f > 1.0) ? 1.0 : f);

		sharedCounterValue = _InterlockedIncrement(SharedCounterForSearch) - 1;
	}
}
