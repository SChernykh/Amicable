#include "stdafx.h"
#include "PrimeTables.h"
#include "Engine.h"
#include "RangeGen.h"
#include "PGO_Helper.h"
#include "Tests.h"

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

int main()
{
	num64 startPrime = 100000000000;
	num64 primeLimit = 100000010000000;

	PrimeTablesInit(startPrime, primeLimit, nullptr);

	num64 SharedCounterForSearch[2] = {};
	num64 sharedCounterValue = 0;

	SearchLargePrimes(SharedCounterForSearch, startPrime, primeLimit, sharedCounterValue);

	return 0;
}
