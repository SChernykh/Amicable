#include "stdafx.h"
#include "Definitions.h"

#if 1

const num128 SearchLimit::value = atoi128("1000000000000000000000");// 10^21
const num64 SearchLimit::LinearLimit = 44721359550ULL;				// sqrt(value * 2) + 1
const num64 SearchLimit::MainPrimeTableBound = 44721359550ULL;		// sqrt(value * 2) + 1
const num64 SearchLimit::PrimeInversesBound = 15811388301ULL;		// sqrt(value) / 2
const num64 SearchLimit::SafeLimit = num64(-1);						// value / 20

#else

const num128 SearchLimit::value = atoi128("18446744073709551616");	// 2^64
const num64 SearchLimit::LinearLimit = 6074001000ULL;				// sqrt(value * 2) + 1
const num64 SearchLimit::MainPrimeTableBound = 6074001000ULL;		// sqrt(value * 2) + 1
const num64 SearchLimit::PrimeInversesBound = 2147483648ULL;		// sqrt(value) / 2
const num64 SearchLimit::SafeLimit = 922337203685477580;			// value / 20

#endif
