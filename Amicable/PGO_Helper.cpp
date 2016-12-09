#include "stdafx.h"
#include "PrimeTables.h"
#include "Engine.h"
#include "RangeGen.h"
#include "PGO_Helper.h"
#include <setjmp.h>
#include <signal.h>

static THREAD_LOCAL bool locIsWorkerThread = false;
static THREAD_LOCAL jmp_buf locWorkerThreadJumpBuffer;

static NOINLINE void ProfileGuidedOptimization_Instrument_WorkerThread(number data)
{
	locIsWorkerThread = true;

	Timer t;

	SUPPRESS_WARNING(4611)
	if (setjmp(locWorkerThreadJumpBuffer))
	{
		const double dt = t.getElapsedTime();
		printf("Thread %llu finished in %.3f seconds\n", data, dt);
		return;
	}

	TRY
	{
		RangeData r;
		r.factors[0].p = 2;
		r.factors[0].k = 8;
		r.factors[0].index = 0;
		r.factors[1].p = CompileTimePrimes<CompileTimePrimesCount>::value;
		r.factors[1].k = 1;
		r.factors[1].index = CompileTimePrimesCount;
		r.factors[1].p_inv = MultiplicativeInverse<CompileTimePrimes<CompileTimePrimesCount>::value>::value;
		r.factors[1].q_max = number(-1) / CompileTimePrimes<CompileTimePrimesCount>::value;
		r.value = 256 * CompileTimePrimes<CompileTimePrimesCount>::value;
		r.sum = 511 * (CompileTimePrimes<CompileTimePrimesCount>::value + 1);
		r.start_prime = CompileTimePrimes<CompileTimePrimesCount + 1>::value;
		r.index_start_prime = CompileTimePrimesCount + 1;
		r.last_factor_index = 1;

		switch (data)
		{
		case 0:
			SearchRange(r);
			break;

		case 1:
			SearchRangeSquared(r);
			break;

		case 2:
			r.start_prime = 7;
			r.index_start_prime = 3;
			r.value = 20;
			r.sum = 42;
			SearchRangeCubed(r);
			break;

		case 3:
			{
				number k = 0;
				SearchLargePrimes(&k, CompileTimeParams::MainPrimeTableBound + 1, CompileTimeParams::SafeLimit);
			}
			break;

		case 4:
			{
				number k = 0;
				SearchLargePrimes(&k, CompileTimeParams::SafeLimit / 5, CompileTimeParams::SafeLimit);
			}
			break;

		case 5:
			{
				RangeData tmp;
				do {} while (RangeGen::Iterate(tmp));
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

	RangeGen::Init();

	const number NumWorkerThreads = 6;

	std::vector<std::thread> threads;
	threads.reserve(NumWorkerThreads);
	for (number i = 0; i < NumWorkerThreads; ++i)
	{
		threads.emplace_back(std::thread(ProfileGuidedOptimization_Instrument_WorkerThread, i));
	}

	Sleep(1000);

	printf("Terminating threads which are still active...\n");

	// Just disable access to some crucial data structures to crash running threads.
	// Access violation will be caught and the thread will then finish gracefully, saving all profiling data.
	DisableAccessToMemory(privPrimeInverses, 65536);

	for (number i = 0; i < threads.size(); ++i)
	{
		threads[i].join();
	}
}
