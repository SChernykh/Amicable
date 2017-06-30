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

#define sprintf_s sprintf

#define PAUSE __builtin_ia32_pause

template <typename T, num64 N> char(*ArraySizeHelper(T(&)[N]))[N];
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

FORCEINLINE num64 _rotr64(num64 value, int shift)
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
#include <sys/resource.h>
#include <unistd.h>
#include <limits.h>
#include <semaphore.h>

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

FORCEINLINE num64 SetHighestTimerResolution() { return 0; }
FORCEINLINE void SetTimerResoluion(const num64) {}
FORCEINLINE void HiResSleep(const double ms)
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

class Semaphore
{
public:
	FORCEINLINE Semaphore() { sem_init(&mySemaphore, 0, 0); }
	FORCEINLINE ~Semaphore() { sem_destroy(&mySemaphore); }

	FORCEINLINE bool Signal() { return (sem_post(&mySemaphore) == 0); }
	FORCEINLINE bool Wait() { return (sem_wait(&mySemaphore) == 0); }

private:
	sem_t mySemaphore;
};

FORCEINLINE void* AllocateSystemMemory(num64 size, bool is_executable)
{
	return mmap(0, size, is_executable ? (PROT_READ | PROT_WRITE | PROT_EXEC) : (PROT_READ | PROT_WRITE), MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
}

FORCEINLINE void DisableAccessToMemory(void* ptr, num64 size)
{
	mprotect(ptr, size, PROT_NONE);
}

FORCEINLINE void ForceRoundUpFloatingPoint()
{
	fesetround(FE_UPWARD);
}

FORCEINLINE bool SetLowPriority()
{
	errno = 0;
	const int result = nice(NZERO - 1);
	return ((result != -1) || (errno == 0));
}

// All code that doesn't pass "-Wpedantic" must be below this comment

#pragma GCC system_header

typedef unsigned __int128 num128;

FORCEINLINE num64 _umul128(num64 a, num64 b, num64* h)
{
	const unsigned __int128 result = static_cast<unsigned __int128>(a) * b;
	*h = static_cast<num64>(result >> 64);
	return static_cast<num64>(result);
}

FORCEINLINE num64 udiv128(num64 numhi, num64 numlo, num64 den, num64* rem)
{
	const unsigned __int128 n = (static_cast<unsigned __int128>(numhi) << 64) + numlo;

	const num64 result = n / den;
	*rem = n % den;

	return result;
}

FORCEINLINE num64 udiv128_noremainder(num64 numlo, num64 numhi, num64 den)
{
	const num64 result = ((static_cast<num128>(numhi) << 64) + numlo) / den;
	return result;
}

FORCEINLINE num64 mulmod64(num64 a, num64 b, num64 n)
{
	return (static_cast<unsigned __int128>(a) * b) % n;
}

FORCEINLINE void add128(num64 a_lo, num64 a_hi, num64 b_lo, num64 b_hi, num64* result_lo, num64* result_hi)
{
	const unsigned __int128 result = ((static_cast<unsigned __int128>(a_hi) << 64) + a_lo) + ((static_cast<unsigned __int128>(b_hi) << 64) + b_lo);
	*result_lo = static_cast<num64>(result);
	*result_hi = static_cast<num64>(result >> 64);
}

FORCEINLINE void sub128(num64 a_lo, num64 a_hi, num64 b_lo, num64 b_hi, num64* result_lo, num64* result_hi)
{
	const unsigned __int128 result = ((static_cast<unsigned __int128>(a_hi) << 64) + a_lo) - ((static_cast<unsigned __int128>(b_hi) << 64) + b_lo);
	*result_lo = static_cast<num64>(result);
	*result_hi = static_cast<num64>(result >> 64);
}

FORCEINLINE byte leq128(num64 a_lo, num64 a_hi, num64 b_lo, num64 b_hi)
{
	return (static_cast<__int128>(((static_cast<unsigned __int128>(a_hi) << 64) + a_lo) - ((static_cast<unsigned __int128>(b_hi) << 64) + b_lo)) <= 0) ? 1 : 0;
}

FORCEINLINE byte less128(num64 a_lo, num64 a_hi, num64 b_lo, num64 b_hi)
{
	return (static_cast<__int128>(((static_cast<num128>(a_hi) << 64) + a_lo) - ((static_cast<num128>(b_hi) << 64) + b_lo)) < 0) ? 1 : 0;
}

FORCEINLINE void shl128(num64& lo, num64& hi, unsigned char count)
{
	num128 t = (static_cast<num128>(hi) << 64) + lo;
	t <<= count;
	lo = static_cast<num64>(t);
	hi = static_cast<num64>(t >> 64);
}

FORCEINLINE void shr128(num64& lo, num64& hi, unsigned char count)
{
	unsigned __int128 t = (static_cast<unsigned __int128>(hi) << 64) + lo;
	t >>= count;
	lo = static_cast<num64>(t);
	hi = static_cast<num64>(t >> 64);
}

#include "num128.h"

FORCEINLINE num64 LowWord(const num128& a) { return static_cast<num64>(a); }
FORCEINLINE num64 HighWord(const num128& a) { return static_cast<num64>(a >> 64); }

FORCEINLINE num128 CombineNum128(num64 _lo, num64 _hi)
{
	return (static_cast<num128>(_hi) << 64) + _lo;
}

FORCEINLINE double Num128ToDouble(const num128 a) { return static_cast<double>(a); }
