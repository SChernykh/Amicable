#include "stdafx.h"
#include "PrimeTables.h"
#include "Engine.h"
#include "RangeGen.h"
#include "PGO_Helper.h"
#include "Tests.h"

int main(int argc, char* argv[])
{
	PrimeTablesInit();

	if ((argc > 1) && (strcmp(argv[1], "/instrument") == 0))
	{
		// Quickly generate profile data for all hot (most used) code paths
		ProfileGuidedOptimization_Instrument();
		return 0;
	}

	if ((argc > 1) && (strcmp(argv[1], "/test") == 0))
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
		if (!TestMaximumSumOfDivisors3())
		{
			std::cerr << "MaximumSumOfDivisors3 test failed" << std::endl;
			return 3;
		}
		if (!TestPrimeSieve())
		{
			std::cerr << "primesieve test failed" << std::endl;
			return 4;
		}
		std::cout << "All tests passed" << std::endl;
		return 0;
	}

	number numThreads = 0;
	if ((argc == 3) && (strcmp(argv[1], "/threads") == 0))
	{
		numThreads = static_cast<number>(atoi(argv[2]));
	}

#if RUN_PERFORMANCE_TEST
	for (;;)
#endif
	{
		LARGE_INTEGER f, t1, t2;
		QueryPerformanceFrequency(&f);
		QueryPerformanceCounter(&t1);
		RangeGen::Run(numThreads);
		QueryPerformanceCounter(&t2);
		printf("completed in %f seconds\n%u pairs found\n%e CPU cycles\n\n",
			static_cast<double>(t2.QuadPart - t1.QuadPart) / f.QuadPart,
			NumFoundPairs,
			static_cast<double>(RangeGen::cpu_cycles)
		);
	}
#if !RUN_PERFORMANCE_TEST
	return 0;
#endif
}
