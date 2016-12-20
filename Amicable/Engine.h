#pragma once

struct RangeData;

void SearchRange(const RangeData& r);
void SearchRangeSquared(const RangeData& r);
void SearchRangeCubed(const RangeData& r);

void SearchLargePrimes(volatile number* SharedCounterForSearch, const number StartPrime, const number PrimeLimit);

void CheckPairNoInline(const number n1, const number targetSum);
void CheckPair128NoInline(const number n1, number targetSumLow, number targetSumHigh);
number MaximumSumOfDivisors3NoInline(const number a, const number p0, const number a_div_p0);
	
void CheckDivisibilityBench();

void SetNumFoundPairsInThisThread(unsigned int value);
unsigned int GetNumFoundPairsInThisThread();

extern bool g_PrintNumbers;
