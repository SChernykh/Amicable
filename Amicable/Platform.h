#pragma once

// Platform/compiler-dependent code must be here

#ifdef _MSC_VER
#include "Platform_MSVC_Windows.h"
#elif __GNUG__
#include "Platform_GCC_Linux.h"
#else
static_assert(false, "This compiler is not supported");
#endif
