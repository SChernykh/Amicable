#pragma once

struct RangeData;

number SearchRange(const RangeData& r);
number SearchRangeSquared(const RangeData& r);
number SearchRangeCubed(const RangeData& r);

void SearchLargePrimes(volatile number* SharedCounterForSearch, const number StartPrime, const number PrimeLimit, number &sharedCounterValue);

void CheckPairNoInline(const number n1, const number targetSum);
void CheckPair128NoInline(const number n1, number targetSumLow, number targetSumHigh);
number MaximumSumOfDivisors3NoInline(const number a, const number p0, const number a_div_p0);
	
void SetNumFoundPairsInThisThread(unsigned int value);
unsigned int GetNumFoundPairsInThisThread();

extern bool g_PrintNumbers;
extern FILE* g_outputFile;
