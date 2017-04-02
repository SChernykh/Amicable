#pragma once

#include "CompileTimeFunctions.h"

enum
{
	ShiftMultiplier = 2,
	MaxPrimeFactors = 16,
};

#ifndef DYNAMIC_SEARCH_LIMIT
#define DYNAMIC_SEARCH_LIMIT 0
#endif

#if DYNAMIC_SEARCH_LIMIT

namespace SearchLimit
{
	extern number value;
	extern number LinearLimit;
	extern number MainPrimeTableBound;
	extern number PrimeInversesBound;
	extern number SafeLimit;
}

#else

enum SearchLimit : number
{
	// Search up to 10^3, 1 pair
	//value = 1000,

	// Search up to 10^4, 5 pairs
	//value = 10000,

	// Search up to 10^5, 13 pairs
	//value = 100000,

	// Search up to 10^6, 42 pairs
	//value = 1000000,

	// Search up to 10^7, 108 pairs
	//value = 10000000,

	// Search up to 10^8, 236 pairs
	//value = 100000000,

	// Search up to 10^9, 586 pairs
	//value = 1000000000,

	// Search up to 10^10, 1427 pairs
	//value = 10000000000,

	// Search up to 10^11, 3340 pairs
	//value = 100000000000,

	// Search up to 10^12, 7642 pairs
	//value = 1000000000000,

	// Search up to 10^13, 17519 pairs
	//value = 10000000000000,

	// Search up to 10^14, 39374 pairs
	//value = 100000000000000,

	// Search up to 10^15, 87102 pairs
	//value = 1000000000000000,

	// Search up to 10^16, 190775 pairs
	//value = 10000000000000000,

	// Search up to 10^17, 415523 pairs
	//value = 100000000000000000,

	// Search up to 10^18, 901312 pairs
	//value = 1000000000000000000,

	// Search up to 10^19, ~1955600 pairs
	//value = 10000000000000000000ULL,

	// Search up to 2^64, ~2401900 pairs
	value = 18446744073709551615ULL,

	// Linear search starts when p >= SearchLimit::LinearLimit
	// It iterates over 1 * p, 2 * p, 3 * p, ... until it reaches SearchLimit::value
	// It assumes that:
	// 1) If k is abundant then k * p is also abundant and vice versa
	// 2) If k if deficient then k * p if also deficient and vice versa
	// We can guarantee these assumptions only if p / k > 2 which means linear limit must be > sqrt(limit * 2)
	// But we don't want to store too many values for linear search, so we only store values up to 10000000
	LinearLimit = ((CompileTimeSQRT<SearchLimit::value>::value + 1) * 1414213563) / 1000000000 + 1,

	MainPrimeTableBound = Max<LinearLimit, 1000>::value,
	PrimeInversesBound = Max<CompileTimeSQRT<SearchLimit::value / 4>::value, CompileTimePrimes<CompileTimePrimesCount>::value>::value,

	// Safe upper bound for the largest prime factor
	//
	// Lemma: m = k * p (p is the largest prime factor) can't be a smaller member of an amicable pair when k < 20
	//
	// Proof.
	//
	// n = S(m) - m must be > m, so let's check that for all k < 20 it's false:
	//
	// 1) For all p > 31:
	//
	// S(1 * p) - p = (p + 1) - p = 1 < p
	// S(2 * p) - 2 * p = 3 * (p + 1) - 2 * p = p + 3 < 2 * p when p > 3
	// S(3 * p) - 3 * p = 4 * (p + 1) - 3 * p = p + 4 < 3 * p when p > 2
	// S(4 * p) - 4 * p = 7 * (p + 1) - 4 * p = 3 * p + 7 < 4 * p when p > 7
	// S(5 * p) - 5 * p = 6 * (p + 1) - 5 * p = p + 6 < 5 * p when p > 1
	// S(6 * p) - 6 * p = 12 * (p + 1) - 6 * p = 6 * p + 12 > 6 * p
	// But! 6 * p = 0 (mod 6) and S(6 * p) = 12 * (p + 1) - even number, so it can't be an amicable number
	//
	// S(7 * p) - 7 * p = 8 * (p + 1) - 7 * p = p + 8 < 7 * p when p > 1
	// S(8 * p) - 8 * p = 15 * (p + 1) - 8 * p = 7 * p + 15 < 8 * p when p > 15
	// S(9 * p) - 9 * p = 13 * (p + 1) - 9 * p = 4 * p + 13 < 9 * p when p > 2
	// S(10 * p) - 10 * p = 18 * (p + 1) - 10 * p = 8 * p + 18 < 10 * p when p > 9
	// S(11 * p) - 11 * p = 12 * (p + 1) - 11 * p = p + 12 < 11 * p when p > 1
	// S(12 * p) - 12 * p = 28 * (p + 1) - 12 * p = 16 * p + 28 > 12 * p
	// But! 12 * p = 0 (mod 6) and S(12 * p) = 28 * (p + 1) - even number, so it can't be an amicable number 
	//
	// S(13 * p) - 13 * p = 14 * (p + 1) - 13 * p = p + 14 < 13 * p when p > 1
	// S(14 * p) - 14 * p = 24 * (p + 1) - 14 * p = 10 * p + 24 < 14 * p when p > 6
	// S(15 * p) - 15 * p = 24 * (p + 1) - 15 * p = 9 * p + 24 < 15 * p when p > 4
	// S(16 * p) - 16 * p = 31 * (p + 1) - 16 * p = 15 * p + 31 < 16 * p when p > 31
	// S(17 * p) - 17 * p = 18 * (p + 1) - 17 * p = p + 18 < 17 * p when p > 1
	// S(18 * p) - 18 * p = 39 * (p + 1) - 18 * p = 21 * p + 39 > 18 * p
	// But! 18 * p = 0 (mod 6) and S(18 * p) = 39 * (p + 1) - even number because (p + 1) is even, so it can't be an amicable number 
	//
	// S(19 * p) - 19 * p = 20 * (p + 1) - 19 * p = p + 20 < 19 * p when p > 2
	// S(20 * p) - 20 * p = 42 * (p + 1) - 20 * p = 22 * p + 42 > 20 * p, so we must stop here
	//
	// 2) Amicable pairs where smaller member has the form k * p, k < 20, p <= 31
	// They must be <= 19 * 31 = 589, so it can only be 220 = 20 * 11, so k = 20 in this case (m = k * p)
	//
	// Lemma is proved.
	//
	// We can use SearchLimit::value / 20 as a safe upper bound
	// because if largest prime factor p > SearchLimit::value / 20 then we'll only have numbers of the form k * p where k < 20
	// SearchLimit::value / 20 is also the best possible bound because 220 = 20 * 11
	SafeLimit = value / 20,
};

#endif

enum PrimeTableParameters
{
	Modulo = 210,
	NumOffsets = 48,
};

static_assert(PrimeTableParameters::Modulo < 256, "Modulo is too large");

const unsigned int NumbersCoprimeToModulo[NumOffsets * 2] = {
	  1,  11,  13,  17,  19,  23,  29,  31,
	 37,  41,  43,  47,  53,  59,  61,  67,
	 71,  73,  79,  83,  89,  97, 101, 103,
	107, 109, 113, 121, 127, 131, 137, 139,
	143, 149, 151, 157, 163, 167, 169, 173,
	179, 181, 187, 191, 193, 197, 199, 209,

	  1 + 210,  11 + 210,  13 + 210,  17 + 210,  19 + 210,  23 + 210,  29 + 210,  31 + 210,
	 37 + 210,  41 + 210,  43 + 210,  47 + 210,  53 + 210,  59 + 210,  61 + 210,  67 + 210,
	 71 + 210,  73 + 210,  79 + 210,  83 + 210,  89 + 210,  97 + 210, 101 + 210, 103 + 210,
	107 + 210, 109 + 210, 113 + 210, 121 + 210, 127 + 210, 131 + 210, 137 + 210, 139 + 210,
	143 + 210, 149 + 210, 151 + 210, 157 + 210, 163 + 210, 167 + 210, 169 + 210, 173 + 210,
	179 + 210, 181 + 210, 187 + 210, 191 + 210, 193 + 210, 197 + 210, 199 + 210, 209 + 210,
};

enum ByteParams
{
	Bits = 8,
};

FORCEINLINE bool IsPerfectSquareCandidate(const number n)
{
	// If N is a perfect square, then N mod 64 must be one of the following: 0, 1, 4, 9, 16, 17, 25, 33, 36, 41, 49, 57
	enum Modulo64SquareCheck : number
	{
		value =
		(number(1) << 0) |
		(number(1) << 1) |
		(number(1) << 4) |
		(number(1) << 9) |
		(number(1) << 16) |
		(number(1) << 17) |
		(number(1) << 25) |
		(number(1) << 33) |
		(number(1) << 36) |
		(number(1) << 41) |
		(number(1) << 49) |
		(number(1) << 57)
	};
	return (((number(1) << n) & Modulo64SquareCheck::value) != 0);
}

FORCEINLINE number StrToNumber(const char* s)
{
	number result = 0;
	for (;;)
	{
		const char c = *(s++);
		if ((c < '0') || (c > '9'))
		{
			break;
		}
		result = result * 10 + (c - '0');
	}
	return result;
}

void atoi128(const char* s, number &numlo, number &numhi);
