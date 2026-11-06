#pragma once

// Platform/compiler-dependent code must be here

#ifdef _MSC_VER
#include "Platform_MSVC_Windows.h"
#elif __GNUG__
#include "Platform_GCC_Linux.h"
#else
static_assert(false, "This compiler is not supported");
#endif

// x86-64 (SSE2 included): the SIMD code paths (AVX2 and AVX-512 chosen at runtime by CPUID). Other architectures (ARM, RISC-V) only
// compile the portable code; AMICABLE_PORTABLE selects it on x86-64 too (to test it).
#if (defined(__x86_64__) || defined(_M_X64)) && !defined(AMICABLE_PORTABLE)
#define AMICABLE_X64 1
#else
#define AMICABLE_X64 0
#endif

// A fast monotonic counter for shares of time (the time stamp counter on x86-64, steady_clock elsewhere)
#if AMICABLE_X64
FORCEINLINE num64 ReadCycleCounter()
{
	return __rdtsc();
}
#else
#include <chrono>
FORCEINLINE num64 ReadCycleCounter()
{
	return static_cast<num64>(std::chrono::steady_clock::now().time_since_epoch().count());
}
#endif
