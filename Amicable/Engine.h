#pragma once

struct RangeData;

num64 SearchRange(const RangeData& r);
num64 SearchRangeSquared(const RangeData& r);
num64 SearchRangeCubed(const RangeData& r);

void SearchLargePrimes(volatile num64* SharedCounterForSearch, const num64 StartPrime, const num64 PrimeLimit, num64 &sharedCounterValue);

void CheckPairNoInline(const num64 n1, const num64 targetSum);
void CheckPair128NoInline(const num128 n1, const num128 targetSum);
num64 MaximumSumOfDivisors3NoInline(const num64 a, const num64 p0, const num64 a_div_p0);
	
void SetNumFoundPairsInThisThread(unsigned int value);
unsigned int GetNumFoundPairsInThisThread();

extern bool g_PrintNumbers;
extern FILE* g_outputFile;
