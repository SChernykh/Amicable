#pragma once

#include "PrimeTables.h"

template<int PrimeIndex>
FORCEINLINE void FactorizeInternal(number a, number& sumA, std::vector<std::pair<number, number>>& aFactorization)
{
	enum {p = CompileTimePrimes<PrimeIndex>::value};

	if (p * p > a)
	{
		if (a > 1)
		{
			sumA *= (a + 1);
			if (aFactorization.empty() || (aFactorization.back().first != a))
				aFactorization.push_back(std::pair<number, number>(a, 1));
			else
				++aFactorization.back().second;
		}
		return;
	}

	number q = a / p;
	if (q * p == a)
	{
		std::pair<number, number> factor(p, 1);

		number n = p;
		number curSum = p + 1;
		a = q;

		while (p <= a)
		{
			q = a / p;
			if (q * p != a)
				break;
			++factor.second;
			n *= p;
			curSum += n;
			a = q;
		}
		sumA *= curSum;
		aFactorization.push_back(factor);
		if (IsPrime(a))
		{
			aFactorization.push_back(std::make_pair(a, 1));
			sumA *= a + 1;
			return;
		}
	}

	FactorizeInternal<PrimeIndex + 1>(a, sumA, aFactorization);
}

FORCEINLINE bool DivideFast(const number aNumPrimesCheckedSoFar, const number a, const number p, number& q)
{
	if (aNumPrimesCheckedSoFar < ReciprocalsTableSize)
		return PrimeReciprocals[aNumPrimesCheckedSoFar].Divide(a, p, q);

	q = a / p;
	return ((a % p) == 0);
}

template<> FORCEINLINE void FactorizeInternal<CompileTimePrimesCount>(number a, number& sumA, std::vector<std::pair<number, number>>& aFactorization)
{
	number numPrimesCheckedSoFar = CompileTimePrimesCount;
	const byte* shift = NextPrimeShifts + CompileTimePrimesCount - 1;
	number p = static_cast<number>(CompileTimePrimes<CompileTimePrimesCount - 1>::value) + *(shift++) * SearchLimit::ShiftMultiplier;
	while (p * p <= a)
	{
		number q;
		if (DivideFast(numPrimesCheckedSoFar, a, p, q))
		{
			std::pair<number, number> factor(p, 1);

			number n = p;
			number curSum = p + 1;
			a = q;

			while (p <= a)
			{
				if (!DivideFast(numPrimesCheckedSoFar, a, p, q))
					break;
				++factor.second;
				n *= p;
				curSum += n;
				a = q;
			}
			sumA *= curSum;
			aFactorization.push_back(factor);
			if (IsPrime(a))
			{
				aFactorization.push_back(std::make_pair(a, 1));
				sumA *= a + 1;
				return;
			}
		}

		p = (p <= SearchLimit::PrimesUpToSqrtLimitValue) ? (p + *(shift++) * SearchLimit::ShiftMultiplier) : (p + 2);
		++numPrimesCheckedSoFar;
	}
	if (a > 1)
	{
		sumA *= (a + 1);
		if (aFactorization.empty() || (aFactorization.back().first != a))
			aFactorization.push_back(std::pair<number, number>(a, 1));
		else
			++aFactorization.back().second;
	}
}

void Factorize(number a, number& sumA, std::vector<std::pair<number, number>>& aFactorization);
