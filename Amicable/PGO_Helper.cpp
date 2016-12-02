#include "stdafx.h"
#include "PrimeTables.h"
#include "Engine.h"
#include "RangeGen.h"
#include "PGO_Helper.h"
#include <process.h>

static NOINLINE unsigned int __stdcall ProfileGuidedOptimization_Instrument_WorkerThread(void* data)
{
	LARGE_INTEGER freq, t1, t2;
	QueryPerformanceFrequency(&freq);
	QueryPerformanceCounter(&t1);
	__try
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

		switch (number(data))
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
	__except(EXCEPTION_EXECUTE_HANDLER)
	{
	}
	QueryPerformanceCounter(&t2);
	printf("Thread %llu finished in %.3f seconds\n", number(data), static_cast<double>(t2.QuadPart - t1.QuadPart) / freq.QuadPart);
	return 0; 
}

NOINLINE void ProfileGuidedOptimization_Instrument()
{
	printf("Collecting profile data...\n");

	RangeGen::Init();

	HANDLE h[6];
	for (number i = 0; i < ARRAYSIZE(h); ++i)
	{
		h[i] = reinterpret_cast<HANDLE>(_beginthreadex(0, 0, ProfileGuidedOptimization_Instrument_WorkerThread, (void*) i, 0, 0));
	}

	if (WaitForMultipleObjects(ARRAYSIZE(h), h, true, 1000) == WAIT_TIMEOUT)
	{
		printf("Terminating threads which are still active...\n");

		// Just disable access to some crucial data structures to crash running threads.
		// Access violation will be caught and the thread will then finish gracefully, saving all profiling data.
		DWORD oldProtect;
		VirtualProtect(reinterpret_cast<LPVOID>(privPrimeInverses), 65536, PAGE_NOACCESS, &oldProtect);

		WaitForMultipleObjects(ARRAYSIZE(h), h, true, INFINITE);
	}

	for (number i = 0; i < ARRAYSIZE(h); ++i)
	{
		CloseHandle(h[i]);
	}
}
