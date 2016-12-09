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
		std::cout << "Testing CheckPair()...";
		flush(std::cout);
		if (!TestCheckPair())
		{
			std::cerr << "CheckPair() test failed" << std::endl;
			return 1;
		}
		std::cout << "done" << std::endl << "Testing amicable candidates...";
		flush(std::cout);
		if (!TestAmicableCandidates())
		{
			std::cerr << "Amicable candidates test failed" << std::endl;
			return 2;
		}
		std::cout << "done" << std::endl << "Testing MaximumSumOfDivisors3()...";
		flush(std::cout);
		if (!TestMaximumSumOfDivisors3())
		{
			std::cerr << "MaximumSumOfDivisors3() test failed" << std::endl;
			return 3;
		}
		std::cout << "done" << std::endl << "Testing primesieve...";
		flush(std::cout);
		if (!TestPrimeSieve())
		{
			std::cerr << "primesieve test failed" << std::endl;
			return 4;
		}
		std::cout << "done" << std::endl << "All tests passed" << std::endl;
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
		Timer t;
		RangeGen::Run(numThreads);
		const double dt = t.getElapsedTime();

		printf("completed in %f seconds\n%u pairs found\n%e CPU cycles\n\n",
			dt,
			GetNumFoundPairsInThisThread(),
			static_cast<double>(RangeGen::cpu_cycles)
		);
	}
#if !RUN_PERFORMANCE_TEST
	return 0;
#endif
}
