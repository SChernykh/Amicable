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
PRAGMA_WARNING(disable : 4987)
PRAGMA_WARNING(disable : 5039)
PRAGMA_WARNING(disable : 5045)

PRAGMA_WARNING(push, 1)
PRAGMA_WARNING(disable : 4265)
PRAGMA_WARNING(disable : 4365)
PRAGMA_WARNING(disable : 4548)
PRAGMA_WARNING(disable : 4571)
PRAGMA_WARNING(disable : 4623)
PRAGMA_WARNING(disable : 4625)
PRAGMA_WARNING(disable : 4626)
PRAGMA_WARNING(disable : 4774)
PRAGMA_WARNING(disable : 5026)
PRAGMA_WARNING(disable : 5027)

typedef unsigned char byte;
typedef unsigned long long int num64;

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <thread>
#include <string>
#include <sstream>
#include <list>

#include "Platform.h"
#include "Definitions.h"

#include <primesieve/config.hpp>
#include <primesieve/PreSieve.hpp>
#include <primesieve/primesieve_error.hpp>
#include <primesieve/littleendian_cast.hpp>
#include <primesieve/PrimeGenerator.hpp>
#include <primesieve/pmath.hpp>
#include <primesieve/PrimeSieve.hpp>

PRAGMA_WARNING(pop)
