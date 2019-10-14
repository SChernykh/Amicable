#include "stdafx.h"
#include "Definitions.h"

#if 1

const num128 SearchLimit::value = atoi128("1000000000000000000000");	// 10^21
const num64 SearchLimit::LinearLimit = 100000000000ULL;				// 10^11
const num64 SearchLimit::MainPrimeTableBound = 100000000000ULL;		// 10^11
const num64 SearchLimit::RangeGenPrimeBound = 15811388301ULL;		// sqrt(value) / 2
const num64 SearchLimit::SafeLimit = 18446744073709551615ULL;		// value / 20

#else

const num128 SearchLimit::value = atoi128("18446744073709551616");	// 2^64
const num64 SearchLimit::LinearLimit = 6074001000ULL;				// sqrt(value * 2) + 1
const num64 SearchLimit::MainPrimeTableBound = 6074001000ULL;		// sqrt(value * 2) + 1
const num64 SearchLimit::RangeGenPrimeBound = 2147483648ULL;		// sqrt(value) / 2
const num64 SearchLimit::SafeLimit = 922337203685477580;			// value / 20

#endif
