#include "stdafx.h"
#include "PrimeTables.h"
#include "Engine.h"

static number EvenOddCount = 0;

NOINLINE void CheckEvenOdd(const number n1, const number n2)
{
	++EvenOddCount;
	number sqrt_n2 = static_cast<number>(sqrt(static_cast<double>(n2)));
	number n3 = sqrt_n2 * sqrt_n2;
	while (n3 < n2)
	{
		n3 += sqrt_n2 * 2 + 1;
		++sqrt_n2;
	}
	if (n3 == n2)
		CheckPair(n1, n1 + n2);
}

FORCEINLINE void CheckPowersOf2EvenOdd(const number aValue, const number aSumOfDivisors)
{
	number nextValue = aValue * 2;
	number n = aSumOfDivisors * 2;
	number sumN = aSumOfDivisors * 3;
	for (;;)
	{
		if (nextValue < sumN - nextValue)
			CheckEvenOdd(nextValue, sumN - nextValue);
		if (nextValue > SearchLimit::value / 2)
			return;
		nextValue *= 2;
		n *= 2;
		sumN += n;
	}
}

FORCEINLINE void CheckPowersOf3EvenOdd(const number aValue, const number aSumOfDivisors)
{
	if (aValue > SearchLimit::value / 9)
		return;
	number nextValue = aValue * 9;
	number n = aSumOfDivisors * 9;
	number sumN = aSumOfDivisors * 13;
	for (;;)
	{
		CheckPowersOf2EvenOdd(nextValue, sumN);
		if (nextValue > SearchLimit::value / 9)
			return;
		nextValue *= 9;
		n *= 3;
		sumN += n;
		n *= 3;
		sumN += n;
	}
}

template<int NumDistinctOddPrimeFactors>
FORCEINLINE void SearchEvenOdd(const number aValue, const number aSumOfDivisors, const number aPreviousPrimeFactor)
{
	CheckPowersOf3EvenOdd(aValue, aSumOfDivisors);
	number curPrime = 5;
	const byte* shift = NextPrimeShifts + 2;
	while (curPrime < aPreviousPrimeFactor)
	{
		const number p2 = curPrime * curPrime;
		number highProductNextValue;
		number nextValue = _umul128(aValue, p2, &highProductNextValue);
		if ((nextValue > SearchLimit::value) || highProductNextValue)
			break;
		number n = aSumOfDivisors * curPrime;
		number sumN = aSumOfDivisors + n;
		n *= curPrime;
		sumN += n;
		for (;;)
		{
			SearchEvenOdd<NumDistinctOddPrimeFactors + 1>(nextValue, sumN, curPrime);
			nextValue = _umul128(nextValue, p2, &highProductNextValue);
			if ((nextValue > SearchLimit::value) || highProductNextValue)
				break;
			n *= curPrime;
			sumN += n;
			n *= curPrime;
			sumN += n;
		}
		curPrime += *(shift++) * SearchLimit::ShiftMultiplier;
	}
}

template<> FORCEINLINE void SearchEvenOdd<(MaxSearchDepth<SearchLimit::value>::value + InlinePrimesInSearch) / 2>(const number aValue, const number aSumOfDivisors, const number)
{
	CheckPowersOf2EvenOdd(aValue, aSumOfDivisors);
}

void RunSearchEvenOdd()
{
	_control87(_RC_DOWN, _MCW_RC);
	SearchEvenOdd<0>(1, 1, SearchLimit::SQRT2);
	std::cout << EvenOddCount << " even-odd pairs checked" << std::endl;
}
