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
	Factor factors[MaxPrimeFactors];
	num128 value;
	num128 sum;
	num64 start_prime;
	unsigned int index_start_prime;
	int last_factor_index;
};

class RangeGen
{
public:
	static void Init(char* startFrom, char* stopAt, RangeData* outStartFromRange, Factor* outStopAtFactors, unsigned int largestPrimePower);
	static void Run(num64 numThreads, char* startFrom, char* stopAt, unsigned int largestPrimePower, num64 startPrime, num64 primeLimit);
	static bool Iterate(RangeData& range);

	static double total_numbers_to_check;

	enum
	{
		LargePrimesSplitShift = 10,
		LargePrimesSplitSize = 1 << LargePrimesSplitShift,
	};

	// Long ranges of the first stage (the range sieve of SearchRange: p up to 10^11 for one M0, tens of seconds on one thread) are
	// shared: the thread that gets one checks it in chunks of SharedRangeChunk numbers, and threads without other work take chunks
	// too. A thread takes a chunk only when its read-ahead buffer is empty, and sets its checkpoint state to the shared range
	// under the lock, so the checkpoint (the most lagging thread) never passes a range with unfinished chunks.
	// Returns false if the range is not shared (too short, or too many shared ranges): the caller checks it then.
	static bool ShareRange(const RangeData& r, num64 p_lo, num64 p_hi, unsigned int is_over_abundant_mask, num64 num_primes);
	static const num64 SharedRangeChunk = 20 * 210 * 8192; // 34.4M numbers: whole segments of both layouts of the range sieve

private:
	FORCEINLINE RangeGen() { InitializeCriticalSection(&lock); }
	FORCEINLINE ~RangeGen() { DeleteCriticalSection(&lock);}

	struct StackFrame
	{
		num128 value;
		num128 sum;
	};

	template<unsigned int largest_prime_power> static bool Iterate(RangeData& range);

	struct WorkerThreadState
	{
		RangeData curRange;
		unsigned int curLargestPrimePower;
		num64 total_numbers_checked;
	};

	struct WorkerThreadParams
	{
		const RangeData* rangeToCheckFirst;
		const Factor* stopAtFactors;
		num64 startPrime;
		num64 primeLimit;
		num64 sharedCounterValue;
		unsigned int startLargestPrimePower;
		WorkerThreadState stateToSave;
		volatile bool finished;
	};

	static void CheckpointThread(WorkerThreadParams* params, num64 numWorkerThreads);
	static void WorkerThread(WorkerThreadParams* params);

	struct SharedRange
	{
		RangeData range;
		unsigned int is_over_abundant_mask;
		bool inUse;
		num64 id;
		num64 p_next;
		num64 p_hi;
		double primesPerNumber; // for the sieve tier of a chunk
		int tier;               // sieve tier of the next chunks (-1: by the length and window height), see SearchSharedChunk
	};
	enum { MaxSharedRanges = 64 };

	// With the lock held: takes the next chunk of a shared range if there is one
	static bool TakeSharedChunk(RangeData& r, unsigned int& is_over_abundant_mask, num64& p_lo, num64& p_hi, num64& num_primes, unsigned int& slot, num64& id, int& tier);
	static void SearchSharedChunk(const RangeData& r, unsigned int is_over_abundant_mask, num64 p_lo, num64 p_hi, num64 num_primes, unsigned int slot, num64 id, int tier);

private:
	static CRITICAL_SECTION lock;
	static CACHE_ALIGNED StackFrame search_stack[MaxPrimeFactors];
	static CACHE_ALIGNED Factor factors[MaxPrimeFactors];
	static int search_stack_depth;
	static int prev_search_stack_depth;
	static unsigned int cur_largest_prime_power;
	static volatile num64 SharedCounterForSearch[2];
	static num64 total_numbers_checked;
	static volatile bool allDone;

	// (under the lock)
	static SharedRange sharedRanges[MaxSharedRanges];
	static num64 sharedRangeNextId;
	static num64 numWorkerThreadsTotal;
	static num64 numIdleThreads;

	static RangeGen RangeGen_instance;
};
