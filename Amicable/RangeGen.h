#pragma once

// Iterates over numbers of the form m*q where q is the largest prime factor and:
//
// 1) m*q < limit
// 2) m is not divisible by 6
// 3) m is abundant (but not overabundant) OR m is deficient, but m*q is abundant
//
// Iteration order is bottom to top: each iteration changes the largest prime.
//
// So the iteration order will be, for example:
//
// 2*5*7*11^k for k >= 1
// 2*5*7*13^k for k >= 1
// 2*5*7*17^k for k >= 1
// ...
// 2*5*7^2*11^k for k >= 1
//
// and so on
//
// This iteration order allows very effective filtering of deficient and overabundant numbers

struct RangeData
{
	Factor factors[16];
	number value;
	number sum;
	number start_prime;
	unsigned int index_start_prime;
	int last_factor_index;
};

class RangeGen
{
public:
	static void Init();
	static void Run(number numThreads);
	static bool Iterate(RangeData& range);

	static number cpu_cycles;

private:
	FORCEINLINE RangeGen() { InitializeCriticalSection(&lock); }
	FORCEINLINE ~RangeGen() { DeleteCriticalSection(&lock);}

	struct StackFrame
	{
		number value, sum, q0, q, sum_q;
		unsigned int max_prime;
		int start_j;
	};

	template<unsigned int largest_prime_power> static bool Iterate(RangeData& range);
	static unsigned int __stdcall WorkerThread(void*);

private:
	static CRITICAL_SECTION lock;
	static StackFrame search_stack[16];
	static CACHE_ALIGNED Factor factors[16];
	static int search_stack_depth;
	static int prev_search_stack_depth;
	static unsigned int cur_largest_prime_power;
	static volatile number SharedCounterForSearch;

	static RangeGen RangeGen_instance;
};
