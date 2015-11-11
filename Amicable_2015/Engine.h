#pragma once

extern number TotalCPUcycles;
extern __declspec(thread) number NumFoundPairs;

void CheckPair(const number n1, const number targetSum);
number RunSearch(number numThreadsOverride = 0);
void RunSearchEvenOdd();
