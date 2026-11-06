#include "stdafx.h"
#include "PrimeTables.h"
#include "Engine.h"
#include "RangeGen.h"
#include "PGO_Helper.h"
#include "Tests.h"
#include "RangeSieve.h"

double boinc_process_cpu_time_correction = 0.0;

PRAGMA_WARNING(push, 1)
PRAGMA_WARNING(disable : 4091 4917 4625 4626 5026 5027)
#include <boinc_api.h>
#include <diagnostics.h>
PRAGMA_WARNING(pop)

#if defined(_MSC_VER) && (_MSC_VER >= 1400)
void AppInvalidParameterHandler(const wchar_t* /*expression*/, const wchar_t* /*function*/, const wchar_t* /*file*/, unsigned int /*line*/, uintptr_t /*pReserved*/)
{
	DebugBreak();
}
#endif

int main(int argc, char* argv[])
{
#if defined(_MSC_VER) && (_MSC_VER >= 1400)
	_set_invalid_parameter_handler(AppInvalidParameterHandler);
#endif

	if (!TestNum128Division())
	{
		std::cerr << "128 bit division test failed" << std::endl;
		return 1;
	}

	BOINC_OPTIONS options;
	boinc_options_defaults(options);
	options.multi_thread = true;

	boinc_init_diagnostics(BOINC_DIAG_REDIRECTSTDERR | BOINC_DIAG_TRACETOSTDERR);

	const int boinc_init_result = boinc_init_options(&options);
	if (boinc_init_result)
	{
		char buf[256];
		fprintf(stderr, "%s boinc_init returned %d\n", boinc_msg_prefix(buf, sizeof(buf)), boinc_init_result);
		exit(boinc_init_result);
	}

	boinc_fraction_done(0.0);

	char* startFrom = nullptr;
	char* stopAt = nullptr;
	unsigned int largestPrimePower = 1;
	num64 startPrime = 0;
	num64 primeLimit = 0;

	// Parse command line parameters
	for (int i = 1; i < argc; ++i)
	{
		if (strcmp(argv[i], "/bench") == 0)
		{
			g_PrintNumbers = false;
		}

		if (strcmp(argv[i], "/nosieve") == 0)
		{
			g_UseRangeSieve = false;
		}

		if (strcmp(argv[i], "/instrument") == 0)
		{
			// Quickly generate profile data for all hot (most used) code paths
			PrimeTablesInit(LowWord(SearchLimit::value / 1000000), SearchLimit::SafeLimit, nullptr);
			ProfileGuidedOptimization_Instrument();
			return 0;
		}

		if (strcmp(argv[i], "/test") == 0)
		{
			PrimeTablesInit(SearchLimit::LinearLimit, SearchLimit::SafeLimit, nullptr);
			g_PrintNumbers = false;

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
			std::cout << "done" << std::endl << "Testing range sieve...";
			flush(std::cout);
			if (!TestRangeSieve())
			{
				std::cerr << "Range sieve test failed" << std::endl;
				return 5;
			}
			std::cout << "done" << std::endl << "All tests passed" << std::endl;
			return 0;
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
			startPrime = StrToNumber(argv[++i]);
			primeLimit = StrToNumber(argv[++i]);
		}

		if ((strcmp(argv[i], "/task_size") == 0) && (i + 1 < argc))
		{
			RangeGen::total_numbers_to_check = atof(argv[++i]);
		}

		// testing: the vector code of the range sieve scan (0 = SSE2, 1 = AVX2, 2 = AVX-512; at most what the CPU supports)
		if ((strcmp(argv[i], "/sieve_scan") == 0) && (i + 1 < argc))
		{
			g_RangeSieveScanLevel = atoi(argv[++i]);
		}

		// testing: the band sieve of the near-degenerate ranges (0 = off)
		if ((strcmp(argv[i], "/sieve_band") == 0) && (i + 1 < argc))
		{
			g_RangeSieveBand = atoi(argv[++i]);
		}
		if ((strcmp(argv[i], "/sieve_band_min") == 0) && (i + 1 < argc))
		{
			g_RangeSieveBandMinSpanLog2 = atoi(argv[++i]);
		}

		// testing: the layout of the range sieve (30 or 210)
		if ((strcmp(argv[i], "/sieve_layout") == 0) && (i + 1 < argc))
		{
			g_RangeSieveLayout = atoi(argv[++i]);
		}

		// testing: force the depth tier of the range sieve (0..5)
		if ((strcmp(argv[i], "/sieve_tier") == 0) && (i + 1 < argc))
		{
			g_RangeSieveTierOverride = atoi(argv[++i]);
		}

		// testing: the deepest tier for each class of the window height W (8 comma-separated values)
		if ((strcmp(argv[i], "/sieve_wtiers") == 0) && (i + 1 < argc))
		{
			const char* s = argv[++i];
			for (unsigned int k = 0; (k < 8) && *s; ++k)
			{
				g_RangeSieveMaxTierByW[k] = atoi(s);
				while (*s && (*s != ','))
				{
					++s;
				}
				if (*s == ',')
				{
					++s;
				}
			}
		}
	}

	PrimeTablesInit(startPrime, primeLimit, stopAt);

	APP_INIT_DATA aid;
	boinc_get_init_data(aid);

	std::string resolved_name;
	const int boinc_resolve_result = boinc_resolve_filename_s("output.txt", resolved_name);
	if (boinc_resolve_result)
	{
		char buf[256];
		fprintf(stderr, "%s boinc_resolve_filename returned %d\n", boinc_msg_prefix(buf, sizeof(buf)), boinc_resolve_result);
		exit(boinc_resolve_result);
	}

	g_outputFile = boinc_fopen(resolved_name.c_str(), "ab+");
	RangeGen::Run(static_cast<num64>(ceil(aid.ncpus)), startFrom, stopAt, largestPrimePower, startPrime, primeLimit);
	RangeSievePrintStats();
#ifdef RANGE_SIEVE_VERIFY
	std::cerr << "Range sieve verification: " << g_RangeSieveVerifyCount << " sums checked, " << g_RangeSieveVerifyErrors << " mismatches" << std::endl;
#endif

	// Now sort all numbers found so far
	rewind(g_outputFile);
	std::vector<num128> found_numbers;
	while (!feof(g_outputFile))
	{
		char buf[32];
		if (!fgets(buf, sizeof(buf), g_outputFile))
		{
			break;
		}
		const num128 m = atoi128(buf);
		if (m != 0)
		{
			found_numbers.push_back(m);
		}
	}
	fclose(g_outputFile);

	std::sort(found_numbers.begin(), found_numbers.end());

	// Remove all duplicates
	found_numbers.erase(std::unique(found_numbers.begin(), found_numbers.end()), found_numbers.end());

	// And write remaining sorted numbers back to the output file
	g_outputFile = boinc_fopen(resolved_name.c_str(), "wb");
	for (num128 m : found_numbers)
	{
		char buf[40];
		fprintf(g_outputFile, "%s\n", itoa128(m, buf, sizeof(buf)));
	}
	fclose(g_outputFile);

	boinc_finish(0);
	return 0;
}
