#pragma once

#pragma warning(disable : 4619)

#pragma warning(disable : 4324)
#pragma warning(disable : 4350)
#pragma warning(disable : 4505)
#pragma warning(disable : 4514)
#pragma warning(disable : 4592)
#pragma warning(disable : 4710)
#pragma warning(disable : 4711)
#pragma warning(disable : 4820)

#pragma warning(push, 1)
#pragma warning(disable : 4265)
#pragma warning(disable : 4365)
#pragma warning(disable : 4571)

typedef unsigned char byte;
typedef unsigned long long int number;

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <thread>

#pragma warning(pop)

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
