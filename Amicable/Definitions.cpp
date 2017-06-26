#include "stdafx.h"
#include "Definitions.h"

#if 1

const num128 SearchLimit::value = atoi128("100000000000000000000");	// 10^20
const num64 SearchLimit::LinearLimit = 14142135624ULL;				// sqrt(value * 2) + 1
const num64 SearchLimit::MainPrimeTableBound = 14142135624ULL;		// sqrt(value * 2) + 1
const num64 SearchLimit::RangeGenPrimeBound = 5000000000ULL;		// sqrt(value) / 2
const num64 SearchLimit::SafeLimit = 5000000000000000000ULL;		// value / 20

#else

const num128 SearchLimit::value = atoi128("18446744073709551616");	// 2^64
const num64 SearchLimit::LinearLimit = 6074001000ULL;				// sqrt(value * 2) + 1
const num64 SearchLimit::MainPrimeTableBound = 6074001000ULL;		// sqrt(value * 2) + 1
const num64 SearchLimit::RangeGenPrimeBound = 2147483648ULL;		// sqrt(value) / 2
const num64 SearchLimit::SafeLimit = 922337203685477580;			// value / 20

#endif
