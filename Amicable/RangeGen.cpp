#include "stdafx.h"
#include "PrimeTables.h"
#include "Engine.h"
#include "RangeGen.h"

CRITICAL_SECTION RangeGen::lock = CRITICAL_SECTION_INITIALIZER;
RangeGen::StackFrame RangeGen::search_stack[16];
CACHE_ALIGNED Factor RangeGen::factors[16];
int RangeGen::search_stack_depth;
int RangeGen::prev_search_stack_depth;
unsigned int RangeGen::cur_largest_prime_power;
volatile number RangeGen::SharedCounterForSearch;
number RangeGen::cpu_cycles;

RangeGen RangeGen::RangeGen_instance;

template<unsigned int largest_prime_power>
NOINLINE bool RangeGen::Iterate(RangeData& range)
{
	StackFrame* s = search_stack + search_stack_depth;
	Factor* f = factors + search_stack_depth;
	int start_i;

recurse_begin:
	const bool is_return = (search_stack_depth < prev_search_stack_depth);
	prev_search_stack_depth = search_stack_depth;
	if (is_return)
	{
		goto recurse_return;
	}

	start_i = (search_stack_depth == 0) ? 0 : (factors[search_stack_depth - 1].index + 1);

	f->p = (search_stack_depth == 0) ? 2 : (factors[search_stack_depth - 1].p + NextPrimeShifts[factors[search_stack_depth - 1].index] * CompileTimeParams::ShiftMultiplier);

	// A check to ensure that m is not divisible by 6
	if (search_stack_depth == 1)
	{
		// factors[0].p is 2
		// factors[1].p is 3
		// change factors[1].p to 5
		if (start_i == 1)
		{
			start_i = 2;
			f->p = 5;
		}
	}

	// Check only 2, 3, 5 as the smallest prime factor because the smallest abundant number coprime to 2*3*5 is ~2*10^25
	s->max_prime = static_cast<unsigned int>((search_stack_depth > 0) ? (MainPrimeTableBound + 1) : 7);
	for (f->index = start_i; f->p < s->max_prime; f->p += NextPrimeShifts[f->index] * CompileTimeParams::ShiftMultiplier, ++f->index)
	{
		number h;
		s[1].value = _umul128(s->value, f->p, &h);
		if ((s[1].value >= SearchLimit::value) || h)
		{
			--search_stack_depth;
			--s;
			--f;
			goto recurse_begin;
		}
		s[1].sum = s->sum * (f->p + 1);

		f->k = 1;
		f->p_inv = PrimeInverses[f->index].first;
		f->q_max = PrimeInverses[f->index].second;

		for (;;)
		{
			s->start_j = f->index + 1;
			s->q0 = f->p + NextPrimeShifts[f->index] * CompileTimeParams::ShiftMultiplier;

			// A check to ensure that m*q is not divisible by 6
			if (search_stack_depth == 0)
			{
				// factors[0].p is 2
				// q is 3
				// change q to 5
				if (s->start_j == 1)
				{
					s->start_j = 2;
					s->q0 = 5;
				}
			}

			{
				number q = s->q0;
				number sum_q = s->q0 + 1;

				IF_CONSTEXPR(largest_prime_power > 1)
				{
					q = _umul128(q, s->q0, &h);
					if (h)
					{
						break;
					}
					sum_q += q;
				}

				IF_CONSTEXPR(largest_prime_power > 2)
				{
					q = _umul128(q, s->q0, &h);
					if (h)
					{
						break;
					}
					sum_q += q;
				}

				s->q = q;
				s->sum_q = sum_q;

				const number value_to_check = _umul128(s[1].value, s->q, &h);
				if ((value_to_check >= SearchLimit::value) || h)
				{
					if (f->k == 1)
					{
						--search_stack_depth;
						--s;
						--f;
						goto recurse_begin;
					}
					break;
				}

				// Skip overabundant numbers
				const bool is_deficient = (s[1].sum - s[1].value < s[1].value);
				if (is_deficient || !OverAbundant(factors, search_stack_depth, s[1].value, s[1].sum, static_cast<number>((largest_prime_power & 1) ? 2 : 1)))
				{
					if (!is_deficient || (s[1].sum * s->sum_q - value_to_check > value_to_check))
					{
						range.value = s[1].value;
						range.sum = s[1].sum;
						range.start_prime = s->q0;
						range.index_start_prime = static_cast<unsigned int>(s->start_j);
						IF_CONSTEXPR(largest_prime_power == 1)
						{
							memcpy(range.factors, factors, sizeof(Factor) * (search_stack_depth + 1));
							range.last_factor_index = search_stack_depth;
						}
						++search_stack_depth;
						return true;
					}

					++search_stack_depth;
					++s;
					++f;
					goto recurse_begin;
				}
			}

recurse_return:
			s[1].value = _umul128(s[1].value, f->p, &h);
			if ((s[1].value >= SearchLimit::value) || h)
			{
				break;
			}
			s[1].sum = s[1].sum * f->p + s->sum;
			++f->k;
		}

		// Workaround for shifting from 2 to 3 because NextPrimeShifts[0] is 0
		if (search_stack_depth == 0)
		{
			if (f->p == 2)
			{
				f->p = 3;
			}
		}
	}
	--search_stack_depth;
	return false;
}

bool RangeGen::Iterate(RangeData& range)
{
	if (search_stack_depth < 0)
	{
		return false;
	}

	if (cur_largest_prime_power == 1)
	{
		if (Iterate<1>(range))
		{
			return true;
		}
		else
		{
			cur_largest_prime_power = 2;
			search_stack_depth = 0;
		}
	}

	if (cur_largest_prime_power == 2)
	{
		if (Iterate<2>(range))
		{
			return true;
		}
		else
		{
			cur_largest_prime_power = 3;
			search_stack_depth = 0;
		}
	}

	if (cur_largest_prime_power == 3)
	{
		return Iterate<3>(range);
	}

	return false;
}

NOINLINE void RangeGen::Init()
{
	search_stack[0].value = 1;
	search_stack[0].sum = 1;
	search_stack_depth = 0;
	cur_largest_prime_power = 1;
	SharedCounterForSearch = 0;
}

NOINLINE void RangeGen::Run(number numThreads)
{
	Init();

	if (numThreads == 0)
	{
		numThreads = std::thread::hardware_concurrency();
	}

	if ((numThreads == 0) || (numThreads > 2048))
	{
		numThreads = 1;
	}

	std::vector<std::thread> threads;
	unsigned int* num_pairs = reinterpret_cast<unsigned int*>(alloca(sizeof(unsigned int) * numThreads));
	for (number i = 0; i < numThreads; ++i)
	{
		threads.emplace_back(std::thread(WorkerThread, num_pairs + i));
	}

	for (number i = 0; i < numThreads; ++i)
	{
		threads[i].join();
	}

	cpu_cycles = 0;
	unsigned int numFoundPairsTotal = 0;
	for (number i = 0; i < numThreads; ++i)
	{
		//ULONG64 c;
		//QueryThreadCycleTime(h[i], &c);
		//cpu_cycles += c;
		numFoundPairsTotal += num_pairs[i];
	}
	SetNumFoundPairsInThisThread(numFoundPairsTotal);
}

NOINLINE void RangeGen::WorkerThread(unsigned int* result)
{
	SetNumFoundPairsInThisThread(0);
	enum
	{
		RangesReadAheadSize = 4,
	};
	RangeData ranges[RangesReadAheadSize];
	unsigned int range_largest_prime_power[RangesReadAheadSize];
	number range_write_index = 0;
	number range_read_index = 0;
	for (;;)
	{
		bool is_locked = TryEnterCriticalSection(&lock) != 0;
		if (!is_locked && (range_read_index == range_write_index))
		{
			EnterCriticalSection(&lock);
			is_locked = true;
		}

		if (is_locked)
		{
			while (range_write_index - range_read_index < RangesReadAheadSize)
			{
				if (!Iterate(ranges[range_write_index % RangesReadAheadSize]))
				{
					break;
				}
				range_largest_prime_power[range_write_index % RangesReadAheadSize] = cur_largest_prime_power;
				++range_write_index;
			}
			LeaveCriticalSection(&lock);
		}

		if (range_read_index == range_write_index)
		{
			break;
		}

		switch (range_largest_prime_power[range_read_index % RangesReadAheadSize])
		{
		case 1:
			SearchRange(ranges[range_read_index % RangesReadAheadSize]);
			break;
		case 2:
			SearchRangeSquared(ranges[range_read_index % RangesReadAheadSize]);
			break;
		case 3:
			SearchRangeCubed(ranges[range_read_index % RangesReadAheadSize]);
			break;
		}
		++range_read_index;
	}

	SearchLargePrimes(&SharedCounterForSearch, CompileTimeParams::MainPrimeTableBound + 1, CompileTimeParams::SafeLimit);

	*result = GetNumFoundPairsInThisThread();
}
