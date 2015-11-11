#pragma once

#pragma warning(disable: 4324) // 'struct_name' : structure was padded due to __declspec(align())
#pragma warning(disable: 4350) // behavior change: 'member1' called instead of 'member2'
#pragma warning(disable: 4481) // nonstandard extension used: override specifier 'keyword'
#pragma warning(disable: 4668) // 'symbol' is not defined as a preprocessor macro, replacing with '0' for 'directives'
#pragma warning(disable: 4710) // 'function' : function not inlined
#pragma warning(disable: 4711) // function 'function' selected for inline expansion
#pragma warning(disable: 4820) // 'bytes' bytes padding added after construct 'member_name'
#pragma warning(disable: 4987) // nonstandard extension used: 'throw (...)'

#include <iostream>
#include <vector>

#include <Windows.h>

#include "Definitions.h"

#define NOINLINE DECLSPEC_NOINLINE
