#pragma once

#ifdef _MSC_VER
#define PRAGMA_WARNING(...) __pragma(warning(__VA_ARGS__))
#define _Pragma(...)
#elif __GNUG__
#define PRAGMA_WARNING(...)
#else
static_assert(false, "This compiler is not supported");
#endif

PRAGMA_WARNING(disable : 4619)

PRAGMA_WARNING(disable : 4324)
PRAGMA_WARNING(disable : 4350)
PRAGMA_WARNING(disable : 4505)
PRAGMA_WARNING(disable : 4514)
PRAGMA_WARNING(disable : 4592)
PRAGMA_WARNING(disable : 4710)
PRAGMA_WARNING(disable : 4711)
PRAGMA_WARNING(disable : 4820)

PRAGMA_WARNING(push, 1)
PRAGMA_WARNING(disable : 4265)
PRAGMA_WARNING(disable : 4365)
PRAGMA_WARNING(disable : 4571)

typedef unsigned char byte;
typedef unsigned long long int num64;

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <thread>

PRAGMA_WARNING(pop)

#include "Platform.h"
#include "Definitions.h"

#include <primesieve/config.hpp>
#include <primesieve/callback_t.hpp>
#include <primesieve/PreSieve.hpp>
#include <primesieve/primesieve_error.hpp>
#include <primesieve/Callback.hpp>
#include <primesieve/littleendian_cast.hpp>
#include <primesieve/PrimeFinder.hpp>
#include <primesieve/PrimeGenerator.hpp>
#include <primesieve/pmath.hpp>
#include <primesieve/PrimeSieve.hpp>
