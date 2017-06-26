#pragma once

PRAGMA_WARNING(push, 1)
PRAGMA_WARNING(disable : 4668)

// Target Windows XP or better
#define WINVER 0x0501
#define _WIN32_WINNT 0x0501
#define _WIN32_WINDOWS 0x0501

#include <SDKDDKVer.h>
#include <Windows.h>
#include <tchar.h>
#include <intrin.h>
PRAGMA_WARNING(pop)

#undef min
#undef max

#ifndef NOINLINE
#define NOINLINE __declspec(noinline)
#endif

#ifdef FORCEINLINE
#undef FORCEINLINE
#endif

#define FORCEINLINE __forceinline

#define ASSUME(cond) __assume(cond)
#define UNLIKELY(cond) cond

#define CACHE_ALIGNED __declspec(align(64))

#define THREAD_LOCAL __declspec(thread)

#define PAUSE YieldProcessor

#define CRITICAL_SECTION_INITIALIZER {}

#define IF_CONSTEXPR(X) \
	PRAGMA_WARNING(suppress : 4127) \
	static_assert((X) || !(X), "Error: "#X" is not a constexpr"); \
	if (X)

#define TRY __try
#define EXCEPT(X) __except(X)

class Timer
{
public:
	FORCEINLINE explicit Timer()
	{
		QueryPerformanceFrequency(&f);
		QueryPerformanceCounter(&t1);
	}

	FORCEINLINE double getElapsedTime() const
	{
		LARGE_INTEGER t2;
		QueryPerformanceCounter(&t2);
		return static_cast<double>(t2.QuadPart - t1.QuadPart) / f.QuadPart;
	}

private:
	LARGE_INTEGER f, t1;
};

num64 SetHighestTimerResolution();
void SetTimerResoluion(const num64 res);
void HiResSleep(const double ms);

class Semaphore
{
public:
	FORCEINLINE Semaphore() : mySemaphore(CreateSemaphore(nullptr, 0, 1, nullptr)) {}
	FORCEINLINE ~Semaphore() { CloseHandle(mySemaphore); }

	FORCEINLINE bool Signal() { return (ReleaseSemaphore(mySemaphore, 1, nullptr) != 0); }
	FORCEINLINE bool Wait() { return (WaitForSingleObject(mySemaphore, INFINITE) == WAIT_OBJECT_0);  }

private:
	HANDLE mySemaphore;
};

extern num64(*udiv128)(num64 numhi, num64 numlo, num64 den, num64* rem);
extern num64(*udiv128_noremainder)(num64 numlo, num64 numhi, num64 den);
extern num64(*mulmod64)(num64 a, num64 b, num64 n);

FORCEINLINE void add128(num64 a_lo, num64 a_hi, num64 b_lo, num64 b_hi, num64* result_lo, num64* result_hi)
{
#if _MSC_VER >= 1900
	_addcarry_u64(_addcarry_u64(0, a_lo, b_lo, result_lo), a_hi, b_hi, result_hi);
#else
	const int carry = ((a_lo + b_lo < a_lo) ? 1 : 0);
	*result_lo = a_lo + b_lo;
	*result_hi = a_hi + b_hi + carry;
#endif
}

FORCEINLINE void sub128(num64 a_lo, num64 a_hi, num64 b_lo, num64 b_hi, num64* result_lo, num64* result_hi)
{
#if _MSC_VER >= 1900
	_subborrow_u64(_subborrow_u64(0, a_lo, b_lo, result_lo), a_hi, b_hi, result_hi);
#else
	const int carry = ((a_lo < b_lo) ? 1 : 0);
	*result_lo = a_lo - b_lo;
	*result_hi = a_hi - b_hi - carry;
#endif
}

FORCEINLINE byte leq128(num64 a_lo, num64 a_hi, num64 b_lo, num64 b_hi)
{
#if _MSC_VER >= 1900
	num64 t[2];
	return _subborrow_u64(_subborrow_u64(1, a_lo, b_lo, &t[0]), a_hi, b_hi, &t[1]);
#else
	return static_cast<byte>((a_hi < b_hi) || ((a_hi == b_hi) && (a_lo <= b_lo)));
#endif
}

FORCEINLINE byte less128(num64 a_lo, num64 a_hi, num64 b_lo, num64 b_hi)
{
#if _MSC_VER >= 1900
	return _subborrow_u64(_subborrow_u64(0, a_lo, b_lo, nullptr), a_hi, b_hi, nullptr);
#else
	return static_cast<byte>((a_hi < b_hi) || ((a_hi == b_hi) && (a_lo < b_lo)));
#endif
}

FORCEINLINE void shl128(num64& lo, num64& hi, unsigned char count)
{
	hi = __shiftleft128(lo, hi, count);
	lo <<= count;
}

FORCEINLINE void shr128(num64& lo, num64& hi, unsigned char count)
{
	lo = __shiftright128(lo, hi, count);
	hi >>= count;
}

FORCEINLINE void* AllocateSystemMemory(num64 size, bool is_executable)
{
	return VirtualAlloc(nullptr, size, MEM_COMMIT, static_cast<DWORD>(is_executable ? PAGE_EXECUTE_READWRITE : PAGE_READWRITE));
}

FORCEINLINE void DisableAccessToMemory(void* ptr, num64 size)
{
	DWORD oldProtect;
	VirtualProtect(ptr, size, PAGE_NOACCESS, &oldProtect);
}

FORCEINLINE void ForceRoundUpFloatingPoint()
{
	_control87(_RC_UP, _MCW_RC);
}

FORCEINLINE bool SetLowPriority()
{
	return SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_IDLE) ? true : false;
}

#include "num128.h"

FORCEINLINE num64 LowWord(const num128& a) { return a.lo; }
FORCEINLINE num64 HighWord(const num128& a) { return a.hi; }

FORCEINLINE num128 CombineNum128(num64 _lo, num64 _hi)
{
	num128 result;
	result.lo = _lo;
	result.hi = _hi;
	return result;
}

FORCEINLINE double Num128ToDouble(const num128 a) { return a.hi * 18446744073709551616.0 + a.lo; }
