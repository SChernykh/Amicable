// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#pragma warning(disable : 4324)
#pragma warning(disable : 4350)
#pragma warning(disable : 4505)
#pragma warning(disable : 4514)
#if _MSC_VER >= 1900
#pragma warning(disable : 4592)
#endif
#pragma warning(disable : 4710)
#pragma warning(disable : 4711)
#pragma warning(disable : 4820)

#pragma warning(push, 1)
#pragma warning(disable : 4365)
#pragma warning(disable : 4571)

#include <SDKDDKVer.h>

#include <stdio.h>
#include <tchar.h>
#include <Windows.h>
#include <intrin.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

#pragma warning(pop)

#ifndef NOINLINE
#define NOINLINE __declspec(noinline)
#endif

#ifdef FORCEINLINE
#undef FORCEINLINE
#endif

#define FORCEINLINE __forceinline

#define CACHE_ALIGNED __declspec(align(64))

#undef min
#undef max

#include "Definitions.h"

#include <primesieve/config.hpp>
#include <primesieve/callback_t.hpp>
#include <primesieve/PrimeSieve.hpp>
#include <primesieve/primesieve_error.hpp>
#include <primesieve/Callback.hpp>
#include <primesieve/PreSieve.hpp>
#include <primesieve/PrimeFinder.hpp>
#include <primesieve/PrimeGenerator.hpp>
#include <primesieve/pmath.hpp>
#include <primesieve/littleendian_cast.hpp>
