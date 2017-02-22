#include "stdafx.h"
#include "PrimeTables.h"
#include "RangeGen.h"
#include "OpenCL.h"

PRAGMA_WARNING(push, 1)
PRAGMA_WARNING(disable : 4091 4917)
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

	char* startFrom = nullptr;
	char* stopAt = nullptr;
	unsigned int largestPrimePower = 1;
	number startPrime = 0;
	number primeLimit = 0;

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
	}

	std::cerr << "Initializing prime tables..." << std::flush;
	PrimeTablesInit((startPrime && primeLimit) || !stopAt);
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
		OpenCL cl;
		cl_result = cl.Run(argc, argv, startFrom, stopAt, largestPrimePower);
	}

	if (!cl_result)
	{
		boinc_finish(-1);
		return -1;
	}

	// Now sort all numbers found so far
	rewind(g_outputFile);
	std::vector<number> found_numbers;
	while (!feof(g_outputFile))
	{
		char buf[32];
		if (!fgets(buf, sizeof(buf), g_outputFile))
		{
			break;
		}
		const number m = StrToNumber(buf);
		if (m)
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
	for (number m : found_numbers)
	{
		fprintf(g_outputFile, "%llu\n", m);
	}
	fclose(g_outputFile);

	boinc_finish(0);
	return 0;
}
