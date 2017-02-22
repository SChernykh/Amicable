#include "stdafx.h"
#include "Platform.h"

// Platform/compiler-dependent code must be here

#ifdef _MSC_VER // MSVC & Windows

#pragma section("udiv128", read, execute)
__declspec(allocate("udiv128"))
const unsigned char udiv128Code[] =
{
	0x48, 0x89, 0xD0, // mov rax,rdx
	0x48, 0x89, 0xCA, // mov rdx,rcx
	0x49, 0xF7, 0xF0, // div r8
	0x49, 0x89, 0x11, // mov [r9],rdx
	0xC3              // ret
};

#pragma section("mulmod64", read, execute)
__declspec(allocate("mulmod64"))
const unsigned char mulmod64_code[] = {
		0x48, 0x89, 0xC8, // mov rax,rcx
		0x48, 0xF7, 0xE2, // mul rdx
		0x49, 0xF7, 0xF0, // div r8
		0x48, 0x89, 0xD0, // mov rax,rdx
		0xC3              // ret
};

typedef NTSTATUS (NTAPI *NtQueryTimerResolutionPtr)(OUT PULONG MinimumResolution, OUT PULONG MaximumResolution, OUT PULONG CurrentResolution);
typedef NTSTATUS (NTAPI *NtSetTimerResolutionPtr)(IN ULONG DesiredResolution, IN BOOLEAN SetResolution, OUT PULONG CurrentResolution);
typedef NTSTATUS (NTAPI *NtDelayExecutionPtr)(IN BOOLEAN Alertable, IN PLARGE_INTEGER DelayInterval);

static NtQueryTimerResolutionPtr NtQueryTimerResolution = nullptr;
static NtSetTimerResolutionPtr NtSetTimerResolution = nullptr;
static NtDelayExecutionPtr NtDelayExecution = nullptr;

static bool HiResTimerAvailable = false;

PRAGMA_WARNING(disable : 4191)
static struct HiResTimerInitializer
{
	HiResTimerInitializer()
	{
		const HMODULE hNTDLL = LoadLibraryA("ntdll.dll");
		NtQueryTimerResolution = reinterpret_cast<NtQueryTimerResolutionPtr>(GetProcAddress(hNTDLL, "NtQueryTimerResolution"));
		NtSetTimerResolution = reinterpret_cast<NtSetTimerResolutionPtr>(GetProcAddress(hNTDLL, "NtSetTimerResolution"));
		NtDelayExecution = reinterpret_cast<NtDelayExecutionPtr>(GetProcAddress(hNTDLL, "NtDelayExecution"));
		FreeLibrary(hNTDLL);

		HiResTimerAvailable = (NtQueryTimerResolution && NtSetTimerResolution && NtDelayExecution);
	}
} locHiResTimerInitializer;
PRAGMA_WARNING(default : 4191)

number SetHighestTimerResolution()
{
	if (HiResTimerAvailable)
	{
		ULONG minRes, maxRes, curRes, oldRes;
		NtQueryTimerResolution(&minRes, &maxRes, &curRes);
		NtSetTimerResolution(maxRes, true, &oldRes);
		return oldRes;
	}
	else
	{
		return 0;
	}
}

void SetTimerResoluion(const number res)
{
	if (HiResTimerAvailable)
	{
		ULONG oldRes;
		NtSetTimerResolution(static_cast<ULONG>(res), true, &oldRes);
	}
}

void HiResSleep(const double ms)
{
	if (HiResTimerAvailable)
	{
		LARGE_INTEGER delay;
		delay.QuadPart = static_cast<LONGLONG>(ms * -1e4);
		NtDelayExecution(false, &delay);
	}
	else
	{
		const DWORD delay = static_cast<DWORD>(ms);
		Sleep(delay ? delay : 1);
	}
}

number (*udiv128)(number numhi, number numlo, number den, number* rem) = (number(*)(number, number, number, number*))((const unsigned char*)udiv128Code);
number (*mulmod64)(number a, number b, number n) = (number (*)(number, number, number))((const unsigned char*) mulmod64_code);

#endif
