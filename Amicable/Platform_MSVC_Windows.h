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

extern number(*udiv128)(number numhi, number numlo, number den, number* rem);
extern number(*mulmod64)(number a, number b, number n);

FORCEINLINE void add128(number a_lo, number a_hi, number b_lo, number b_hi, number* result_lo, number* result_hi)
{
#if _MSC_VER >= 1900
	_addcarry_u64(_addcarry_u64(0, a_lo, b_lo, result_lo), a_hi, b_hi, result_hi);
#else
	const int carry = ((a_lo + b_lo < a_lo) ? 1 : 0);
	*result_lo = a_lo + b_lo;
	*result_hi = a_hi + b_hi + carry;
#endif
}

FORCEINLINE void sub128(number a_lo, number a_hi, number b_lo, number b_hi, number* result_lo, number* result_hi)
{
#if _MSC_VER >= 1900
	_subborrow_u64(_subborrow_u64(0, a_lo, b_lo, result_lo), a_hi, b_hi, result_hi);
#else
	const int carry = ((a_lo < b_lo) ? 1 : 0);
	*result_lo = a_lo - b_lo;
	*result_hi = a_hi - b_hi - carry;
#endif
}

FORCEINLINE byte leq128(number a_lo, number a_hi, number b_lo, number b_hi)
{
#if _MSC_VER >= 1900
	number t[2];
	return _subborrow_u64(_subborrow_u64(1, a_lo, b_lo, &t[0]), a_hi, b_hi, &t[1]);
#else
	return static_cast<byte>((a_hi < b_hi) || ((a_hi == b_hi) && (a_lo <= b_lo)));
#endif
}

FORCEINLINE void shr128(number& lo, number& hi, unsigned char count)
{
	lo = __shiftright128(lo, hi, count);
	hi >>= count;
}

FORCEINLINE void* AllocateSystemMemory(number size, bool is_executable)
{
	return VirtualAlloc(nullptr, size, MEM_COMMIT, static_cast<DWORD>(is_executable ? PAGE_EXECUTE_READWRITE : PAGE_READWRITE));
}

FORCEINLINE void DisableAccessToMemory(void* ptr, number size)
{
	DWORD oldProtect;
	VirtualProtect(ptr, size, PAGE_NOACCESS, &oldProtect);
}

FORCEINLINE void ForceRoundUpFloatingPoint()
{
	_control87(_RC_UP, _MCW_RC);
}

FORCEINLINE int GetCurrentPriority()
{
	return GetThreadPriority(GetCurrentThread());
}

FORCEINLINE int SetCurrentPriority(int priority)
{
	return SetThreadPriority(GetCurrentThread(), priority);
}
