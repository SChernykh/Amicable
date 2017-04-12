#include "stdafx.h"
#include "Definitions.h"

const num128 SearchLimit::value = atoi128("100000000000000000000"); // 10^20
const num64 SearchLimit::LinearLimit = 14142135624ULL;				// sqrt(value * 2) + 1
const num64 SearchLimit::MainPrimeTableBound = 14142135624ULL;		// sqrt(value * 2) + 1
const num64 SearchLimit::PrimeInversesBound = 5000000000ULL;		// sqrt(value) / 2
const num64 SearchLimit::SafeLimit = 5000000000000000000ULL;		// value / 20
