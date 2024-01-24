#include "stdafx.h"
#include "Platform.h"

PRAGMA_WARNING(push, 1)
PRAGMA_WARNING(disable : 4091 4917 4625 4626 5026 5027)
#include <boinc_api.h>
#include <diagnostics.h>
PRAGMA_WARNING(pop)

// Platform/compiler-dependent code must be here

#ifdef _MSC_VER // MSVC & Windows

#pragma section("udiv128", read, execute)
__declspec(allocate("udiv128"))
static const unsigned char MachineCode[] =
{
	// udiv128
	0x48, 0x89, 0xD0, // mov rax,rdx
	0x48, 0x89, 0xCA, // mov rdx,rcx
	0x49, 0xF7, 0xF0, // div r8
	0x49, 0x89, 0x11, // mov [r9],rdx
	0xC3,             // ret
	0x90, 0x90, 0x90,

	// udiv128_noremainder
	0x48, 0x89, 0xC8, // mov rax,rcx
	0x49, 0xF7, 0xF0, // div r8
	0xC3,             // ret
	0x90, 0x90, 0x90,
	0x90, 0x90, 0x90,
	0x90, 0x90, 0x90,

	// mulmod64
	0x48, 0x89, 0xC8, // mov rax,rcx
	0x48, 0xF7, 0xE2, // mul rdx
	0x49, 0xF7, 0xF0, // div r8
	0x48, 0x89, 0xD0, // mov rax,rdx
	0xC3,             // ret
	0x90, 0x90, 0x90,
};

num64(*udiv128)(num64 numhi, num64 numlo, num64 den, num64* rem) = (num64(*)(num64, num64, num64, num64*))((const unsigned char*)MachineCode);
num64(*udiv128_noremainder)(num64 numlo, num64 numhi, num64 den) = (num64(*)(num64, num64, num64))(((const unsigned char*)MachineCode) + 16);
num64(*mulmod64)(num64 a, num64 b, num64 n) = (num64(*)(num64, num64, num64))(((const unsigned char*)MachineCode) + 32);

typedef NTSTATUS (NTAPI *NtQueryTimerResolutionPtr)(OUT PULONG MinimumResolution, OUT PULONG MaximumResolution, OUT PULONG CurrentResolution);
typedef NTSTATUS (NTAPI *NtSetTimerResolutionPtr)(IN ULONG DesiredResolution, IN BOOLEAN SetResolution, OUT PULONG CurrentResolution);
typedef NTSTATUS (NTAPI *NtDelayExecutionPtr)(IN BOOLEAN Alertable, IN PLARGE_INTEGER DelayInterval);

static NtQueryTimerResolutionPtr NtQueryTimerResolution_ptr = nullptr;
static NtSetTimerResolutionPtr NtSetTimerResolution_ptr = nullptr;
static NtDelayExecutionPtr NtDelayExecution_ptr = nullptr;

static bool HiResTimerAvailable = false;

PRAGMA_WARNING(disable : 4191)
static struct HiResTimerInitializer
{
	HiResTimerInitializer()
	{
		const HMODULE hNTDLL = LoadLibraryA("ntdll.dll");
		NtQueryTimerResolution_ptr = reinterpret_cast<NtQueryTimerResolutionPtr>(GetProcAddress(hNTDLL, "NtQueryTimerResolution"));
		NtSetTimerResolution_ptr = reinterpret_cast<NtSetTimerResolutionPtr>(GetProcAddress(hNTDLL, "NtSetTimerResolution"));
		NtDelayExecution_ptr = reinterpret_cast<NtDelayExecutionPtr>(GetProcAddress(hNTDLL, "NtDelayExecution"));
		FreeLibrary(hNTDLL);

		HiResTimerAvailable = (NtQueryTimerResolution_ptr && NtSetTimerResolution_ptr && NtDelayExecution_ptr);
	}
} locHiResTimerInitializer;
PRAGMA_WARNING(default : 4191)

num64 SetHighestTimerResolution()
{
	if (HiResTimerAvailable)
	{
		ULONG minRes, maxRes, curRes, oldRes;
		NtQueryTimerResolution_ptr(&minRes, &maxRes, &curRes);
		NtSetTimerResolution_ptr(maxRes, true, &oldRes);
		return oldRes;
	}
	else
	{
		return 0;
	}
}

void SetTimerResoluion(const num64 res)
{
	if (HiResTimerAvailable)
	{
		ULONG oldRes;
		NtSetTimerResolution_ptr(static_cast<ULONG>(res), true, &oldRes);
	}
}

void HiResSleep(const double ms)
{
	if (HiResTimerAvailable)
	{
		LARGE_INTEGER delay;
		delay.QuadPart = static_cast<LONGLONG>(ms * -1e4);
		NtDelayExecution_ptr(false, &delay);
	}
	else
	{
		const DWORD delay = static_cast<DWORD>(ms);
		Sleep(delay ? delay : 1);
	}
}

#elif __GNUG__ // Linux, Mac OS

NOINLINE Semaphore::Semaphore() : mySemaphore(SEM_FAILED)
{
	timespec ts;
	current_utc_time(&ts);
	const double cur_timestamp = ts.tv_sec + ts.tv_nsec / 1e9;

	int error_code = 0;
	for (unsigned int counter = 0; mySemaphore == SEM_FAILED; ++counter)
	{
		if (counter > 1000)
		{
			std::cerr << "Failed to create semaphore \"" << myName << "\", errno " << error_code << std::endl;
			boinc_finish(-1);
		}
		sprintf_s(myName, "/amic_%.6f_%u", cur_timestamp, counter);
		mySemaphore = sem_open(myName, O_CREAT | O_EXCL, S_IRWXU, 0);
		error_code = errno;
	}
}

#else
static_assert(false, "This compiler is not supported");
#endif
