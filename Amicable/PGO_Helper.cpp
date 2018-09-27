#include "stdafx.h"
#include "PrimeTables.h"
#include "Engine.h"
#include "RangeGen.h"
#include "PGO_Helper.h"
#include "sprp64.h"
#include <setjmp.h>
#include <signal.h>

static THREAD_LOCAL bool locIsWorkerThread = false;
static THREAD_LOCAL jmp_buf locWorkerThreadJumpBuffer;

static RangeData locRange1;
static RangeData locRange2;

static NOINLINE void ProfileGuidedOptimization_Instrument_WorkerThread(num64 data)
{
	locIsWorkerThread = true;

	Timer t;

	PRAGMA_WARNING(suppress : 4611)
	if (setjmp(locWorkerThreadJumpBuffer))
	{
		const double dt = t.getElapsedTime();
		printf("Thread %llu finished in %.3f seconds\n", data, dt);
		return;
	}

	TRY
	{
		RangeData r;
		num64 k = 0;
		num64 sharedCounterValue = 0;

		switch (data)
		{
		case 0:
			SearchRange(locRange1);
			break;

		case 1:
			SearchRange(locRange2);
			break;

		case 2:
			r.start_prime = 7;
			r.index_start_prime = 3;
			r.value = 20;
			r.sum = 42;
			SearchRangeSquared(r);
			break;

		case 3:
			r.start_prime = 7;
			r.index_start_prime = 3;
			r.value = 20;
			r.sum = 42;
			SearchRangeCubed(r);
			break;

		case 4:
			SearchLargePrimes(&k, 1000000000000000, 1000000000000000 + 1000000, sharedCounterValue);
			break;

		case 5:
			SearchLargePrimes(&k, SearchLimit::SafeLimit / 50, (SearchLimit::SafeLimit / 50) + 1000000000, sharedCounterValue);
			break;

		case 6:
			for (unsigned int lpr = 1; lpr <= 3; ++lpr)
			{
				RangeData tmp = {};
				RangeGen::Init(nullptr, nullptr, nullptr, nullptr, lpr);
				for (int i = 0; i < 25000; ++i)
				{
					if (!RangeGen::Iterate(tmp))
					{
						break;
					}
				}
			}
			break;

		}
	}
	EXCEPT(EXCEPTION_EXECUTE_HANDLER)
	{
	}

	longjmp(locWorkerThreadJumpBuffer, 1);
}

void sigsegv_handler(int)
{
	if (locIsWorkerThread)
	{
		longjmp(locWorkerThreadJumpBuffer, 1);
	}
}

NOINLINE void ProfileGuidedOptimization_Instrument()
{
	printf("Collecting profile data...\n");

	signal(SIGSEGV, sigsegv_handler);

	g_PrintNumbers = false;

	Factor stopAtFactors[MaxPrimeFactors + 1] = {};

	char startFrom[32];
	char stopAt[32];

	strcpy_s(startFrom, "2^2*11*19*37*487*29023");
	strcpy_s(stopAt, "2^2*11*19*37*499*6323");
	RangeGen::Init(startFrom, stopAt, &locRange1, stopAtFactors, 1);

	strcpy_s(startFrom, "2^3*11*67*71*73*1097");
	strcpy_s(stopAt, "2^3*11*67*71*83*251*587");
	RangeGen::Init(startFrom, stopAt, &locRange2, stopAtFactors, 1);

	RangeGen::Init(nullptr, nullptr, nullptr, nullptr, 1);

	const num64 NumWorkerThreads = 7;

	std::vector<std::thread> threads;
	threads.reserve(NumWorkerThreads);
	for (num64 i = 0; i < NumWorkerThreads; ++i)
	{
		threads.emplace_back(std::thread(ProfileGuidedOptimization_Instrument_WorkerThread, i));
	}

	Sleep(1000);

	printf("Terminating threads which are still active...\n");

	// Just disable access to some crucial data structures to crash running threads.
	// Access violation will be caught and the thread will then finish gracefully, saving all profiling data.
	memset(privPowersOfP_128DivisibilityData, 0, sizeof(privPowersOfP_128DivisibilityData));
	memset(privSumEstimates, 0, sizeof(privSumEstimates));

	DisableAccessToMemory(privNextPrimeShifts, ReciprocalsTableSize);
	DisableAccessToMemory(MainPrimeTable, MainPrimeTableSize);

	for (num64 i = 0; i < NumWorkerThreads; ++i)
	{
		threads[i].join();
	}
}
