#include "stdafx.h"
#include "PrimeTables.h"
#include "Engine.h"
#include "RangeGen.h"
#include "PGO_Helper.h"
#include "Tests.h"

int main(int argc, char* argv[])
{
#if DYNAMIC_SEARCH_LIMIT
	if (argc < 2)
	{
		std::cerr << "You must specify search limit as a first command-line argument" << std::endl;
		return 0;
	}
	SearchLimit::value = static_cast<number>(StrToNumber(argv[1]));
	if (SearchLimit::value < 1000)
	{
		SearchLimit::value = 1000;
	}
	SearchLimit::LinearLimit = static_cast<number>(sqrt(SearchLimit::value * 2.0)) + 1;
	SearchLimit::MainPrimeTableBound = std::max<number>(SearchLimit::LinearLimit, 1000);
	SearchLimit::PrimeInversesBound = std::max<number>(static_cast<number>(sqrt(SearchLimit::value / 4)), CompileTimePrimes<CompileTimePrimesCount>::value);
	SearchLimit::SafeLimit = SearchLimit::value / 20;
#endif

	PrimeTablesInit();

	number numThreads = 0;
	char* startFrom = nullptr;
	char* stopAt = nullptr;
	unsigned int largestPrimePower = 1;
	number startPrime = 0;
	number primeLimit = 0;

	// Parse command line parameters
	for (int i = 1; i < argc; ++i)
	{
		if (strcmp(argv[i], "/bench") == 0)
		{
			g_PrintNumbers = false;
		}

		if (strcmp(argv[i], "/instrument") == 0)
		{
			// Quickly generate profile data for all hot (most used) code paths
			ProfileGuidedOptimization_Instrument();
			return 0;
		}

		if (strcmp(argv[i], "/test") == 0)
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

		if ((strcmp(argv[i], "/threads") == 0) && (i + 1 < argc))
		{
			numThreads = static_cast<number>(atoi(argv[++i]));
		}

		if ((strcmp(argv[i], "/from") == 0) && (i + 1 < argc))
		{
			startFrom = argv[++i];
		}

		if ((strcmp(argv[i], "/to") == 0) && (i + 1 < argc))
		{
			stopAt = argv[++i];
		}

		if (((strcmp(argv[i], "/largest_prime_power") == 0) || (strcmp(argv[i], "/lpp") == 0)) && (i + 1 < argc))
		{
			largestPrimePower = static_cast<unsigned int>(atoi(argv[++i]));
			if (largestPrimePower < 1)
			{
				largestPrimePower = 1;
			}
			if (largestPrimePower > 3)
			{
				largestPrimePower = 3;
			}
		}

		if (((strcmp(argv[i], "/large_primes_range") == 0) || (strcmp(argv[i], "/lpr") == 0)) && (i + 2 < argc))
		{
			startPrime = static_cast<unsigned int>(StrToNumber(argv[++i]));
			primeLimit = static_cast<unsigned int>(StrToNumber(argv[++i]));
		}
	}

	do
	{
		Timer t;
		RangeGen::Run(numThreads, startFrom, stopAt, largestPrimePower, startPrime, primeLimit);
		const double dt = t.getElapsedTime();
		printf("completed in %f seconds\n%u pairs found\n\n", dt, GetNumFoundPairsInThisThread());
	} while (!g_PrintNumbers);

	return 0;
}
