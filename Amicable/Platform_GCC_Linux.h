#pragma once

#ifndef NOINLINE
#define NOINLINE __attribute__ ((noinline))
#endif

#ifdef FORCEINLINE
#undef FORCEINLINE
#endif

#define FORCEINLINE __attribute__((always_inline)) inline

#define CACHE_ALIGNED __attribute__((aligned(64)))

#define THREAD_LOCAL thread_local

template <typename T, number N> char(*ArraySizeHelper(T(&)[N]))[N];
#define ARRAYSIZE(A) (sizeof(*ArraySizeHelper(A)))

typedef unsigned long DWORD;
typedef long long __int64;

FORCEINLINE void _BitScanReverse64(unsigned long* index, unsigned long long mask)
{
	*index = static_cast<unsigned long>(63 - __builtin_clzll(mask));
}

FORCEINLINE void _BitScanForward64(unsigned long* index, unsigned long long mask)
{
	*index = static_cast<unsigned long>(__builtin_ctzll(mask));
}

FORCEINLINE number _rotr64(number value, int shift)
{
	return (value >> shift) | (value << (64 - shift));
}

#include <pthread.h>
#include <string.h>
#include <sys/mman.h>
#include <fenv.h>
#include <cpuid.h>
#include <time.h>
#include <alloca.h>

#define CRITICAL_SECTION pthread_mutex_t
#define CRITICAL_SECTION_INITIALIZER PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP
#define InitializeCriticalSection(X) ((void)(X))
#define DeleteCriticalSection(X) ((void)(X))

FORCEINLINE int TryEnterCriticalSection(CRITICAL_SECTION* lock)
{
	return (pthread_mutex_trylock(lock) == 0) ? 1 : 0;
}

FORCEINLINE int EnterCriticalSection(CRITICAL_SECTION* lock)
{
	return (pthread_mutex_lock(lock) == 0) ? 1 : 0;
}

FORCEINLINE int LeaveCriticalSection(CRITICAL_SECTION* lock)
{
	return (pthread_mutex_unlock(lock) == 0) ? 1 : 0;
}

#define IF_CONSTEXPR(X) \
	static_assert((X) || !(X), "Error: "#X" is not a constexpr"); \
	if (X)

#define __popcnt64 __builtin_popcountll

#define ASSUME(cond) do { if (!(cond)) __builtin_unreachable(); } while (0)
#define UNLIKELY(cond) __builtin_expect(cond, 0)

#define _InterlockedIncrement(X) (__sync_fetch_and_add(X, 1) + 1)

#define sscanf_s sscanf

#define TRY
#define EXCEPT(X)

class Timer
{
public:
	FORCEINLINE explicit Timer()
	{
		clock_gettime(CLOCK_REALTIME, &t1);
	}

	FORCEINLINE double getElapsedTime() const
	{
		timespec t2;
		clock_gettime(CLOCK_REALTIME, &t2);
		return (t2.tv_sec - t1.tv_sec) + (t2.tv_nsec - t1.tv_nsec) * 1e-9;
	}

private:
	timespec t1;
};

FORCEINLINE void Sleep(DWORD ms)
{
	Timer t0;
	for (;;)
	{
		const double sleepTime = ms / 1000.0 - t0.getElapsedTime();
		if (sleepTime <= 0)
		{
			return;
		}

		timespec t;
		t.tv_sec = static_cast<time_t>(sleepTime);
		t.tv_nsec = static_cast<long>((sleepTime - t.tv_sec) * 1e9);
		nanosleep(&t, nullptr);
	}
}

FORCEINLINE bool IsPopcntAvailable()
{
	unsigned int cpuid_data[4] = {};
	__get_cpuid(1, &cpuid_data[0], &cpuid_data[1], &cpuid_data[2], &cpuid_data[3]);
	return (cpuid_data[2] & (1 << 23)) != 0;
}

FORCEINLINE void* AllocateSystemMemory(number size, bool is_executable)
{
	return mmap(0, size, is_executable ? (PROT_READ | PROT_WRITE | PROT_EXEC) : (PROT_READ | PROT_WRITE), MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
}

FORCEINLINE void DisableAccessToMemory(void* ptr, number size)
{
	mprotect(ptr, size, PROT_NONE);
}

FORCEINLINE void ForceRoundUpFloatingPoint()
{
	fesetround(FE_UPWARD);
}

// All code that doesn't pass "-Wpedantic" must be below this comment

#pragma GCC system_header

FORCEINLINE number _umul128(number a, number b, number* h)
{
	const unsigned __int128 result = static_cast<unsigned __int128>(a) * b;
	*h = static_cast<number>(result >> 64);
	return static_cast<number>(result);
}

FORCEINLINE number udiv128(number numhi, number numlo, number den, number* rem)
{
	const unsigned __int128 n = (static_cast<unsigned __int128>(numhi) << 64) + numlo;

	const number result = n / den;
	*rem = n % den;

	return result;
}

FORCEINLINE number mulmod64(number a, number b, number n)
{
	return (static_cast<unsigned __int128>(a) * b) % n;
}

FORCEINLINE void add128(number a_lo, number a_hi, number b_lo, number b_hi, number* result_lo, number* result_hi)
{
	const unsigned __int128 result = ((static_cast<unsigned __int128>(a_hi) << 64) + a_lo) + ((static_cast<unsigned __int128>(b_hi) << 64) + b_lo);
	*result_lo = static_cast<number>(result);
	*result_hi = static_cast<number>(result >> 64);
}

FORCEINLINE void sub128(number a_lo, number a_hi, number b_lo, number b_hi, number* result_lo, number* result_hi)
{
	const unsigned __int128 result = ((static_cast<unsigned __int128>(a_hi) << 64) + a_lo) - ((static_cast<unsigned __int128>(b_hi) << 64) + b_lo);
	*result_lo = static_cast<number>(result);
	*result_hi = static_cast<number>(result >> 64);
}

FORCEINLINE byte leq128(number a_lo, number a_hi, number b_lo, number b_hi)
{
	return (static_cast<__int128>(((static_cast<unsigned __int128>(a_hi) << 64) + a_lo) - ((static_cast<unsigned __int128>(b_hi) << 64) + b_lo)) <= 0) ? 1 : 0;
}

FORCEINLINE void shr128(number& lo, number& hi, unsigned char count)
{
	unsigned __int128 t = (static_cast<unsigned __int128>(hi) << 64) + lo;
	t >>= count;
	lo = static_cast<number>(t);
	hi = static_cast<number>(t >> 64);
}
