#include "stdafx.h"
#include "PrimeTables.h"
#include "Engine.h"
#include "RangeGen.h"

CRITICAL_SECTION RangeGen::lock = CRITICAL_SECTION_INITIALIZER;
CACHE_ALIGNED RangeGen::StackFrame RangeGen::search_stack[16];
CACHE_ALIGNED Factor RangeGen::factors[16];
int RangeGen::search_stack_depth;
int RangeGen::prev_search_stack_depth;
unsigned int RangeGen::cur_largest_prime_power;
volatile number RangeGen::SharedCounterForSearch;

RangeGen RangeGen::RangeGen_instance;

template<unsigned int largest_prime_power>
NOINLINE bool RangeGen::Iterate(RangeData& range)
{
	StackFrame* s = search_stack + search_stack_depth;
	Factor* f = factors + search_stack_depth;

	int start_i, start_j;
	number q0, q, sum_q;

recurse_begin:
	const bool is_return = (search_stack_depth < prev_search_stack_depth);
	prev_search_stack_depth = search_stack_depth;
	if (is_return)
	{
		if (search_stack_depth < 0)
		{
			return false;
		}
		goto recurse_return;
	}

	start_i = (search_stack_depth == 0) ? 0 : (factors[search_stack_depth - 1].index + 1);

	f->p = (search_stack_depth == 0) ? 2 : (factors[search_stack_depth - 1].p + NextPrimeShifts[factors[search_stack_depth - 1].index] * ShiftMultiplier);

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

	for (f->index = start_i; f->p <= SearchLimit::PrimeInversesBound; f->p += NextPrimeShifts[f->index] * ShiftMultiplier, ++f->index)
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
			start_j = f->index + 1;
			q0 = f->p + NextPrimeShifts[f->index] * ShiftMultiplier;

			// A check to ensure that m*q is not divisible by 6
			if (search_stack_depth == 0)
			{
				// factors[0].p is 2
				// q is 3
				// change q to 5
				if (start_j == 1)
				{
					start_j = 2;
					q0 = 5;
				}
			}

			{
				q = q0;
				sum_q = q0 + 1;

				IF_CONSTEXPR(largest_prime_power > 1)
				{
					q = _umul128(q, q0, &h);
					if (h)
					{
						break;
					}
					sum_q += q;
				}

				IF_CONSTEXPR(largest_prime_power > 2)
				{
					q = _umul128(q, q0, &h);
					if (h)
					{
						break;
					}
					sum_q += q;
				}

				const number value_to_check = _umul128(s[1].value, q, &h);
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
					if (!is_deficient || (s[1].sum * sum_q - value_to_check > value_to_check))
					{
						range.value = s[1].value;
						range.sum = s[1].sum;
						range.start_prime = q0;
						range.index_start_prime = static_cast<unsigned int>(start_j);
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
			// Check only 2, 3, 5 as the smallest prime factor because the smallest abundant number coprime to 2*3*5 is ~2*10^25
			if (f->p >= 5)
			{
				break;
			}
		}
	}
	if (search_stack_depth > 0)
	{
		--search_stack_depth;
		--s;
		--f;
		goto recurse_begin;
	}
	search_stack_depth = -1;
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
			prev_search_stack_depth = 0;
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
			prev_search_stack_depth = 0;
		}
	}

	if (cur_largest_prime_power == 3)
	{
		return Iterate<3>(range);
	}

	return false;
}

template<typename T>
FORCEINLINE unsigned int ParseFactorization(char* factorization, T callback)
{
	unsigned int numFactors = 0;
	int counter = 0;
	number prev_p = 0;
	number p = 0;
	number p1 = 2;
	int index_p1 = 0;
	unsigned int k = 0;
	for (char* ptr = factorization, *prevPtr = factorization; ; ++ptr)
	{
		const char c = *ptr;
		if ((c == '^') || (c == '*') || (c == '\0'))
		{
			*ptr = '\0';
			++counter;
			if (counter == 1)
			{
				p = static_cast<number>(StrToNumber(prevPtr));
			}
			else if (counter == 2)
			{
				k = static_cast<unsigned int>(atoi(prevPtr));
			}

			if (((counter == 1) && (c != '^')) || (counter == 2))
			{
				if (counter == 1)
				{
					k = 1;
				}
				counter = 0;
				if ((numFactors > 0) && (p <= prev_p))
				{
					std::cerr << "Factorization '" << factorization << "' is incorrect: factors must be in increasing order" << std::endl;
					abort();
				}
				prev_p = p;

				if (p1 < p)
				{
					if (p1 == 2)
					{
						p1 = 3;
						++index_p1;
					}
					while (p1 < p)
					{
						p1 += NextPrimeShifts[index_p1++] * ShiftMultiplier;
					}
				}

				if (p1 != p)
				{
					std::cerr << "Factorization '" << factorization << "' is incorrect: " << p << " is not a prime" << std::endl;
					abort();
				}

				if (k < 1)
				{
					std::cerr << "Factorization '" << factorization << "'is incorrect: prime power for " << p << " must be >= 1" << std::endl;
					abort();
				}

				callback(p, k, index_p1, numFactors);

				++numFactors;

				if (numFactors >= MaxPrimeFactors)
				{
					std::cerr << "Factorization '" << factorization << "' is incorrect: too many prime factors" << std::endl;
					abort();
				}
			}
			*ptr = c;
			prevPtr = ptr + 1;
			if (c == '\0')
			{
				break;
			}
		}
	}

	if (numFactors == 0)
	{
		std::cerr << "Factorization '" << factorization << "' is incorrect: too few prime factors" << std::endl;
		abort();
	}

	return numFactors;
}

NOINLINE void RangeGen::Init(char* startFrom, char* stopAt, RangeData* outStartFromRange, Factor* outStopAtFactors, unsigned int largestPrimePower)
{
	search_stack[0].value = 1;
	search_stack[0].sum = 1;
	search_stack_depth = 0;
	prev_search_stack_depth = 0;
	cur_largest_prime_power = largestPrimePower;
	SharedCounterForSearch = 0;
	memset(factors, -1, sizeof(factors));

	if (startFrom)
	{
		const unsigned int numFactors = ParseFactorization(startFrom,
			[startFrom](number p, unsigned int k, int p_index, unsigned int factor_index)
			{
				Factor& f = factors[factor_index];
				f.p = p;
				f.k = k;
				f.index = p_index;
				f.p_inv = PrimeInverses[p_index].first;
				f.q_max = PrimeInverses[p_index].second;
				if ((f.p > 2) && (f.p * f.p_inv != 1))
				{
					std::cerr << "Internal error: PrimeInverses table is incorrect" << std::endl;
					abort();
				}

				StackFrame& cur_s = search_stack[factor_index];
				StackFrame& next_s = search_stack[factor_index + 1];
				next_s.value = cur_s.value;
				next_s.sum = cur_s.sum;
				for (number i = 0; i < k; ++i)
				{
					number h;
					next_s.value = _umul128(next_s.value, p, &h);
					if ((next_s.value >= SearchLimit::value) || h)
					{
						std::cerr << "Factorization '" << startFrom << "' is incorrect: number is too large" << std::endl;
						abort();
					}
					next_s.sum = next_s.sum * p + cur_s.sum;
				}
			}
		);

		search_stack_depth = static_cast<int>(numFactors);
		prev_search_stack_depth = static_cast<int>(numFactors);

		StackFrame* s = search_stack + search_stack_depth;
		Factor* f = factors + search_stack_depth;
		RangeData& range = *outStartFromRange;
		range.value = 0;

		int start_j = f->index + 1;
		number q0 = f->p + NextPrimeShifts[f->index] * ShiftMultiplier;

		// A check to ensure that m*q is not divisible by 6
		if (search_stack_depth == 0)
		{
			// factors[0].p is 2
			// q is 3
			// change q to 5
			if (start_j == 1)
			{
				start_j = 2;
				q0 = 5;
			}
		}

		number q = q0;
		number sum_q = q0 + 1;

		number h;
		const number value_to_check = _umul128(s[1].value, q, &h);
		if ((value_to_check < SearchLimit::value) && !h)
		{
			// Skip overabundant numbers
			const bool is_deficient = (s[1].sum - s[1].value < s[1].value);
			if (is_deficient || !OverAbundant(factors, static_cast<int>(numFactors) - 1, s[1].value, s[1].sum, 2))
			{
				if (!is_deficient || (s[1].sum * sum_q - value_to_check > value_to_check))
				{
					range.value = s[1].value;
					range.sum = s[1].sum;
					range.start_prime = q0;
					range.index_start_prime = static_cast<unsigned int>(start_j);
					memcpy(range.factors, factors, sizeof(Factor) * numFactors);
					range.last_factor_index = static_cast<int>(numFactors) - 1;
				}
			}
		}
	}

	if (stopAt)
	{
		ParseFactorization(stopAt,
			[outStopAtFactors](number p, unsigned int k, int /*p_index*/, unsigned int factor_index)
			{
				outStopAtFactors[factor_index].p = p;
				outStopAtFactors[factor_index].k = k;
			}
		);
	}
}

NOINLINE void RangeGen::Run(number numThreads, char* startFrom, char* stopAt, unsigned int largestPrimePower, number startPrime, number primeLimit)
{
	RangeData startFromRange;
	Factor stopAtFactors[MaxPrimeFactors + 1] = {};
	Init(startFrom, stopAt, &startFromRange, stopAtFactors, largestPrimePower);

	if (numThreads == 0)
	{
		numThreads = std::thread::hardware_concurrency();
	}

	if ((numThreads == 0) || (numThreads > 2048))
	{
		numThreads = 1;
	}

	std::vector<std::thread> threads;

	Timer t;

	unsigned int* num_pairs = reinterpret_cast<unsigned int*>(alloca(sizeof(unsigned int) * numThreads));
	WorkerThreadParams* params = reinterpret_cast<WorkerThreadParams*>(alloca(sizeof(WorkerThreadParams) * numThreads));
	for (number i = 0; i < numThreads; ++i)
	{
		params[i].result = num_pairs + i;
		params[i].rangeToCheckFirst = (startFrom && startFromRange.value && (i == 0)) ? &startFromRange : nullptr;
		params[i].stopAtFactors = stopAt ? stopAtFactors : nullptr;
		params[i].startPrime = startPrime;
		params[i].primeLimit = primeLimit;
		params[i].startLargestPrimePower = largestPrimePower;
		threads.emplace_back(std::thread(WorkerThread, &params[i]));
	}

	for (number i = 0; i < numThreads; ++i)
	{
		threads[i].join();
	}

	unsigned int numFoundPairsTotal = 0;
	for (number i = 0; i < numThreads; ++i)
	{
		numFoundPairsTotal += num_pairs[i];
	}
	SetNumFoundPairsInThisThread(numFoundPairsTotal);
}

NOINLINE void RangeGen::WorkerThread(WorkerThreadParams* params)
{
	SetNumFoundPairsInThisThread(0);

	if (params->startPrime && params->primeLimit)
	{
		if (params->startPrime < SearchLimit::MainPrimeTableBound + 1)
		{
			params->startPrime = SearchLimit::MainPrimeTableBound + 1;
		}
		if (params->primeLimit > SearchLimit::SafeLimit)
		{
			params->primeLimit = SearchLimit::SafeLimit;
		}
		if (params->primeLimit < params->startPrime)
		{
			params->primeLimit = params->startPrime;
		}
		SearchLargePrimes(&SharedCounterForSearch, params->startPrime, params->primeLimit);
	}
	else
	{
		if (params->rangeToCheckFirst)
		{
			SearchRange(*params->rangeToCheckFirst);
		}

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

					if (params->stopAtFactors)
					{
						if (cur_largest_prime_power != params->startLargestPrimePower)
						{
							search_stack_depth = -1;
							break;
						}

						const Factor* f1 = params->stopAtFactors;
						const Factor* f2 = factors;
						while ((f1->p == f2->p) && (f1->k == f2->k))
						{
							++f1;
							++f2;
						}
						if (((f1->p == 0) && (f1->k == 0)) || (f1->p < f2->p) || ((f1->p == f2->p) && (f1->k < f2->k)))
						{
							search_stack_depth = -1;
							break;
						}
					}

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

		if (!params->stopAtFactors)
		{
			SearchLargePrimes(&SharedCounterForSearch, SearchLimit::MainPrimeTableBound + 1, SearchLimit::SafeLimit);
		}
	}

	*params->result = GetNumFoundPairsInThisThread();
}
