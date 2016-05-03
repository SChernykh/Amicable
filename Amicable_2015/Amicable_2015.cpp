#include "stdafx.h"
#include "PrimeTables.h"
#include "Tests.h"
#include "Engine.h"

int main(int argc, char* argv[])
{
	PrimeTablesInit(false);

	if ((argc >= 2) && (strcmp(argv[1], "/test") == 0))
	{
		if (!TestCheckPair())
		{
			std::cerr << "TestCheckPair() failed" << std::endl;
			return 1;
		}
		if (!TestLinearSearchData())
		{
			std::cerr << "TestLinearSearchData() failed" << std::endl;
			return 2;
		}
		if (!TestSmallFactorNumbersData())
		{
			std::cerr << "TestSmallFactorNumbersData() failed" << std::endl;
			return 3;
		}
		if (!TestIsValidCandidate())
		{
			std::cerr << "TestIsValidCandidate() failed" << std::endl;
			return 4;
		}
		if (!TestGetMaxSumMDiv2())
		{
			std::cerr << "TestGetMaxSumMDiv2() failed" << std::endl;
			return 5;
		}
		if (!TestMaximumSumOfDivisors3())
		{
			std::cerr << "MaximumSumOfDivisors3 test failed" << std::endl;
			return 6;
		}
		std::cout << "All tests passed" << std::endl;
		return 0;
	}

	number numThreads = 0;
	if ((argc == 3) && (strcmp(argv[1], "/threads") == 0))
	{
		numThreads = static_cast<number>(atoi(argv[2]));
		if (numThreads > 64)
			numThreads = 64;
	}

#if RUN_PERFORMANCE_TEST
	for (;;)
#endif
	{
		LARGE_INTEGER f, t1, t2;
		QueryPerformanceFrequency(&f);
		QueryPerformanceCounter(&t1);
		const number n = RunSearch(numThreads);
		QueryPerformanceCounter(&t2);
		printf("completed in %f seconds\n%I64u pairs found\n%e CPU cycles\n\n",
				static_cast<double>(t2.QuadPart - t1.QuadPart) / f.QuadPart,
				n,
				static_cast<double>(TotalCPUcycles)
		);
	}
}
