#include "stdafx.h"
#include "PrimeTables.h"
#include "RangeGen.h"
#include "Amicable_OpenCL.h"

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

// Largest prime that can occur in the factorization of M0 for the ranges between startFrom and stopAt (RangeGen order):
// all of them start with the common prefix of the two factorizations, and every other prime factor f satisfies
// prefix * f * (next prime after f) < SearchLimit::value. The table must also hold every prime of the two factorizations
// themselves (and the next prime after it): the primes of the prefix, and those that RangeGen parses although no range has them
// (/from 2^3*29*31*367*49606367 /to 2^3*29*31*367*49606367*49606387 read past the end of the table)
static num64 GetRangesPrimeBound(const char* startFrom, const char* stopAt)
{
	auto parse = [](const char* s, std::vector<std::pair<num64, num64>>& factors)
	{
		while (*s)
		{
			num64 p = 0;
			while ((*s >= '0') && (*s <= '9'))
			{
				p = p * 10 + static_cast<num64>(*s - '0');
				++s;
			}
			num64 k = 1;
			if (*s == '^')
			{
				++s;
				k = 0;
				while ((*s >= '0') && (*s <= '9'))
				{
					k = k * 10 + static_cast<num64>(*s - '0');
					++s;
				}
			}
			if (p < 2)
			{
				return false;
			}
			factors.emplace_back(p, k);
			if (*s == '*')
			{
				++s;
			}
			else if (*s)
			{
				return false;
			}
		}
		return true;
	};

	std::vector<std::pair<num64, num64>> a, b;
	if (!parse(startFrom, a) || !parse(stopAt, b))
	{
		return 0;
	}

	num128 prefix = 1;
	for (size_t i = 0; (i < a.size()) && (i < b.size()) && (a[i] == b[i]); ++i)
	{
		for (num64 k = 0; k < a[i].second; ++k)
		{
			prefix *= a[i].first;
		}
	}

	num64 bound = IntegerSquareRoot((SearchLimit::value - 1) / prefix) + 1000000;
	for (const std::vector<std::pair<num64, num64>>* v : { &a, &b })
	{
		for (const std::pair<num64, num64>& f : *v)
		{
			if (bound < f.first + 1000000)
			{
				bound = f.first + 1000000;
			}
		}
	}
	return (bound < SearchLimit::MainPrimeTableBound) ? bound : 0;
}

static void ParseList(const char* s, unsigned int* values, unsigned int maxCount, unsigned int& count)
{
	count = 0;
	while (*s && (count < maxCount))
	{
		values[count++] = static_cast<unsigned int>(StrToNumber(s));
		while ((*s >= '0') && (*s <= '9'))
		{
			++s;
		}
		if (*s == ',')
		{
			++s;
		}
		else
		{
			break;
		}
	}
}

int main(int argc, char* argv[])
{
#if defined(_MSC_VER) && (_MSC_VER >= 1400)
	_set_invalid_parameter_handler(AppInvalidParameterHandler);
#endif

	BOINC_OPTIONS options;
	boinc_options_defaults(options);
	options.normal_thread_priority = true;
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
	bool testLargePrimes = false;

	// Parse command line parameters
	for (int i = 1; i < argc; ++i)
	{
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

		// The original algorithm (trial division of every number) instead of the sieve
		if (strcmp(argv[i], "/nosieve") == 0)
		{
			g_UseSieve = false;
		}

		// Self-tests of the sieve, optionally followed by the sampling stride of the known pairs (1 = all pairs)
		// "/test_lp" tests the large primes stage
		if ((strcmp(argv[i], "/test") == 0) || (strcmp(argv[i], "/test_lp") == 0) || (strcmp(argv[i], "/test_reduction") == 0))
		{
			g_RunTests = true;
			testLargePrimes = (strcmp(argv[i], "/test_lp") == 0);
			g_TestReductionOnly = (strcmp(argv[i], "/test_reduction") == 0);
			if ((i + 1 < argc) && (argv[i + 1][0] >= '0') && (argv[i + 1][0] <= '9'))
			{
				g_TestStride = static_cast<unsigned int>(StrToNumber(argv[++i]));
			}
		}

		// Profiling of the host part of the sieve (no GPU work, no results)
		// Without the necessary conditions of SieveReduction (for comparison)
		if (strcmp(argv[i], "/noreduce") == 0)
		{
			g_SieveReduce = false;
		}
		// Without the divisibility kill of the wheel primes (for comparison)
		if (strcmp(argv[i], "/nokill") == 0)
		{
			g_SieveKill = false;
		}
		if ((strcmp(argv[i], "/analyze") == 0) && (i + 1 < argc))
		{
			g_AnalyzeRate = atof(argv[++i]);
			g_SieveDryRun = true;
		}
		if ((strcmp(argv[i], "/explore") == 0) && (i + 1 < argc))
		{
			g_ExploreNumbers = static_cast<num64>(atof(argv[++i]));
		}
		if (strcmp(argv[i], "/sieve_dry_run") == 0)
		{
			g_SieveDryRun = true;
		}

		// Positions per batch of ranges (log2): every batch ends with a wait on the host and a possible checkpoint
		if ((strcmp(argv[i], "/sieve_batch") == 0) && (i + 1 < argc))
		{
			g_SieveBatchShift = static_cast<unsigned int>(std::min<num64>(40, std::max<num64>(16, StrToNumber(argv[++i]))));
		}

		// Sieve batches on the GPU at once (1-4: the GPU runs the kernels of one batch while the host waits for another)
		if ((strcmp(argv[i], "/sieve_pipes") == 0) && (i + 1 < argc))
		{
			g_SievePipes = static_cast<int>(std::min<num64>(4, std::max<num64>(1, StrToNumber(argv[++i]))));
		}

		// Helper threads that prepare the ranges for the sieve (0-8, 0: all on the main thread)
		if ((strcmp(argv[i], "/sieve_threads") == 0) && (i + 1 < argc))
		{
			g_SieveThreads = static_cast<int>(std::min<num64>(8, StrToNumber(argv[++i])));
		}

		// GPU time of the sieve kernels (OpenCL profiling)
		if (strcmp(argv[i], "/sieve_profile") == 0)
		{
			g_SieveProfile = true;
		}

		// Sieve tuning
		if ((strcmp(argv[i], "/sieve_wg") == 0) && (i + 1 < argc))
		{
			g_SieveConfig.wgSize = static_cast<unsigned int>(StrToNumber(argv[++i]));
		}
		if ((strcmp(argv[i], "/sieve_seg") == 0) && (i + 1 < argc))
		{
			g_SieveConfig.segJ = static_cast<unsigned int>(StrToNumber(argv[++i]));
		}
		if ((strcmp(argv[i], "/sieve_chunk") == 0) && (i + 1 < argc))
		{
			g_SieveConfig.chunkSegments = static_cast<unsigned int>(StrToNumber(argv[++i]));
		}
		if ((strcmp(argv[i], "/sieve_hits") == 0) && (i + 1 < argc))
		{
			g_SieveConfig.hitsPerTask = static_cast<unsigned int>(StrToNumber(argv[++i]));
		}
		if ((strcmp(argv[i], "/sieve_wheel") == 0) && (i + 1 < argc))
		{
			g_SieveConfig.wheelThreshold = static_cast<unsigned int>(StrToNumber(argv[++i]));
		}
		if ((strcmp(argv[i], "/sieve_high_w") == 0) && (i + 1 < argc))
		{
			g_SieveConfig.highTierW = atof(argv[++i]);
		}
		if ((strcmp(argv[i], "/sieve_q") == 0) && (i + 1 < argc))
		{
			ParseList(argv[++i], g_SieveConfig.tierQ, SieveMaxTiers, g_SieveConfig.numTiers);
		}
		if ((strcmp(argv[i], "/sieve_tier_min") == 0) && (i + 1 < argc))
		{
			unsigned int values[SieveMaxTiers] = {};
			unsigned int n = 0;
			ParseList(argv[++i], values, SieveMaxTiers, n);
			for (unsigned int k = 0; k < n; ++k)
			{
				g_SieveConfig.tierMinPrimes[k] = values[k];
			}
		}
	}

	if (g_RunTests)
	{
		// The tests only need small prime tables (and no candidates: the large primes test creates its own)
		startFrom = nullptr;
		stopAt = nullptr;
		startPrime = testLargePrimes ? 10000000000000000000ULL : 0;
		primeLimit = startPrime;
		g_RangesPrimeTableBound = 10000000;
	}
	else if (g_UseSieve && (largestPrimePower == 1) && !(startPrime && primeLimit) && startFrom && stopAt)
	{
		g_RangesPrimeTableBound = GetRangesPrimeBound(startFrom, stopAt);
	}

	std::cerr << "Initializing prime tables..." << std::flush;
	PrimeTablesInit(startPrime, primeLimit, g_RunTests ? "" : stopAt);
	std::cerr << "done" << std::endl << std::flush;

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

	bool cl_result;
	{
		OpenCL cl(aid.project_preferences);
		cl_result = cl.Run(argc, argv, startFrom, stopAt, largestPrimePower, startPrime, primeLimit);
	}

	if (!cl_result)
	{
		boinc_finish(-1);
		return -1;
	}

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
