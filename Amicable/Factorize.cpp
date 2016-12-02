#include "stdafx.h"
#include "Factorize.h"

NOINLINE void Factorize(number a, number& sumA, std::vector<std::pair<number, number>>& aFactorization)
{
	sumA = 1;
	unsigned long bitIndex;
	_BitScanForward64(&bitIndex, a);
	if (bitIndex > 0)
	{
		a >>= bitIndex;
		sumA = (number(1) << (bitIndex + 1)) - 1;
		aFactorization.push_back(std::make_pair(2, bitIndex));
		if (IsPrime(a))
		{
			aFactorization.push_back(std::make_pair(a, 1));
			sumA *= a + 1;
			return;
		}
	}
	FactorizeInternal<1>(a, sumA, aFactorization);
}
