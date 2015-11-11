#include "stdafx.h"

void* GetUDIV128()
{
	const unsigned char udiv128Code[] =
	{
		0x48, 0x89, 0xD0, // mov rax,rdx
		0x48, 0x89, 0xCA, // mov rdx,rcx
		0x49, 0xF7, 0xF0, // div r8
		0x49, 0x89, 0x11, // mov [r9],rdx
		0xC3              // ret
	};
	void* p = VirtualAlloc(0, sizeof(udiv128Code), MEM_COMMIT, PAGE_READWRITE);
	if (p)
	{
		memcpy(p, udiv128Code, sizeof(udiv128Code));
		DWORD k;
		VirtualProtect(p, sizeof(udiv128Code), PAGE_EXECUTE_READ, &k);
	}
	return p;
}

unsigned __int64(__fastcall *udiv128)(unsigned __int64 numhi, unsigned __int64 numlo, unsigned __int64 den, unsigned __int64* rem) = (unsigned __int64(__fastcall *)(unsigned __int64, unsigned __int64, unsigned __int64, unsigned __int64*)) GetUDIV128();
