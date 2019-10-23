#include "stdafx.h"
#include "PrimeTables.h"
#include "RangeGen.h"

FORCEINLINE num64 IntegerCbrt(const num128 n)
{
	if (n < 8)
	{
		return 1;
	}

	unsigned long index;
	if (LowWord(n))
	{
		_BitScanReverse64(&index, LowWord(n));
	}
	else
	{
		_BitScanReverse64(&index, HighWord(n));
		index += 64;
	}

	index = (index * ((1048576 / 3) + 1)) >> 20;

	num64 result = num64(1) << index;
	for (num64 cur_bit = result >> 1; cur_bit > 0; cur_bit >>= 1)
	{
		const num64 k = result | cur_bit;
		if ((num128(k) * k) * k <= n)
		{
			result = k;
		}
	}
	return result;
}

num64 RunRanges(char* startFrom, char* stopAt, num64 largestPrimePower, num64 outputTaskSize = 0)
{
	num64 numbersProcessedTotal = 0;
	num64 prevNumbersProcessedTotal = 0;

	RangeData r;
	memset(&r, 0, sizeof(r));

	Factor stopAtFactors[MaxPrimeFactors + 1];
	memset(stopAtFactors, 0, sizeof(stopAtFactors));

	std::string range1 = startFrom;

	RangeGen::Init(startFrom, stopAt, &r, stopAtFactors, largestPrimePower);
	bool stop = false;
	for (;;)
	{
		stop |= RangeGen::HasReached(r, stopAtFactors);
		if (outputTaskSize && ((numbersProcessedTotal >= outputTaskSize) || stop))
		{
			std::cout << "/from " << range1 << " /to ";

			if (stop)
			{
				std::cout << stopAt;
			}
			else
			{
				std::stringstream s;
				const Factor(&mainFactors)[MaxPrimeFactors] = r.factors;
				for (num64 i = 0; i <= static_cast<num64>(r.last_factor_index); ++i)
				{
					if (i > 0)
					{
						s << '*';
					}
					s << mainFactors[i].p.Get();
					if (mainFactors[i].k > 1)
					{
						s << '^' << mainFactors[i].k;
					}
				}

				range1 = s.str();
				std::cout << range1;
			}

			if (largestPrimePower > 1)
			{
				std::cout << " /lpp " << largestPrimePower;
			}

			std::cout << " /task_size " << numbersProcessedTotal << std::endl;

			numbersProcessedTotal = 0;
		}

		if (stop)
			break;

		num128 prime_limit = 0;

		switch (largestPrimePower)
		{
		case 1:
			prime_limit = (SearchLimit::value - 1) / r.value;
			if (prime_limit > SearchLimit::MainPrimeTableBound)
			{
				prime_limit = SearchLimit::MainPrimeTableBound;
			}

			if (r.sum - r.value < r.value)
			{
				const num128 prime_limit2 = (r.sum - 1) / (r.value * 2 - r.sum);
				if (prime_limit > prime_limit2)
				{
					prime_limit = prime_limit2;
				}
			}
			break;

		case 2:
			prime_limit = static_cast<num64>(sqrt(Num128ToDouble(SearchLimit::value) / Num128ToDouble(r.value)));
			if (r.sum - r.value < r.value)
			{
				// r.sum * (p^2 + p + 1) = r.value * p^2 * 2
				// (r.sum - r.value * 2) * p^2 + r.sum * (p + 1) = 0
				// (r.value * 2 - r.sum) * p^2 - r.sum * (p + 1) = 0
				// (r.value * 2 - r.sum) * p^2 - r.sum * p - r.sum = 0
				// (r.value * 2 - r.sum) / r.sum * p^2 - p - 1 = 0
				// (r.value * 2 / r.sum - 1) * p^2 - p - 1 = 0
				const double A = Num128ToDouble(r.value * 2 - r.sum) / Num128ToDouble(r.sum);
				const num64 prime_limit2 = static_cast<num64>((1.0 + sqrt(1.0 + 4.0 * A)) / (2.0 * A));
				if (prime_limit > prime_limit2)
				{
					prime_limit = prime_limit2;
				}
			}
			break;

		case 3:
			prime_limit = IntegerCbrt(SearchLimit::value / r.value);
			if (r.sum - r.value < r.value)
			{
				// r.sum * (p^3 + p^2 + p + 1) = r.value * p^3 * 2
				// (r.value * 2 - r.sum) * p^3 = r.sum * (p^2 + p + 1)
				// A * p^3 = r.sum * (p^2 + p + 1)
				const num128 A = r.value * 2 - r.sum;

				// Lower bound
				// A * p^3 = r.sum * p^2 < r.sum * (p^2 + p + 1)
				// A * p^3 = r.sum * p^2
				// A * p = r.sum
				// p = r.sum / A
				double x1 = Num128ToDouble(r.sum) / Num128ToDouble(A);

				// Upper bound
				// A * p^3 = r.sum * p^2 * 2 > r.sum * (p^2 + p + 1)
				// A * p^3 = r.sum * p^2 * 2
				// A * p = r.sum * 2
				// p = r.sum * 2 / A
				double x2 = x1 * 2.0;
				// This upper bound holds when p^2 * 2 > p^2 + p + 1
				// p^2 * 2 > p^2 + p + 1
				// p^2 > p + 1
				// p^2 - p - 1 > 0 is true when p > (1 + sqrt(5)) / 2 = 1.6180339887498948482045868343656...

				double x;
				for (;;)
				{
					x = (x1 + x2) * 0.5;
					if ((static_cast<num64>(x1) == static_cast<num64>(x2)) || (x == x1) || (x == x2))
					{
						break;
					}
					if (Num128ToDouble(A) * x * x * x > Num128ToDouble(r.sum)* (x * x + x + 1))
					{
						x2 = x;
					}
					else
					{
						x1 = x;
					}
				}

				const num64 prime_limit2 = static_cast<num64>(x2);
				if (prime_limit > prime_limit2)
				{
					prime_limit = prime_limit2;
				}
			}
			break;
		}

		// Find the first prime "p" such that "p > prime_limit" is true
		num64 a = 0;
		num64 b = NumPrimes;
		const num64 t = LowWord(prime_limit);

		// Use prime counting function estimates to bound the initial search values
		if (t > 3)
		{
			unsigned long index;
			_BitScanReverse64(&index, t);
			if (index < 37)
			{
				// oeis.org sequence A185192
				static const num64 num_primes[37] = { 0, 2, 4, 6, 11, 18, 31, 54, 97, 172, 309, 564, 1028, 1900, 3512, 6542, 12251, 23000, 43390, 82025, 155611, 295947, 564163, 1077871, 2063689, 3957809, 7603553, 14630843, 28192750, 54400028, 105097565, 203280221, 393615806, 762939111, 1480206279, 2874398515, 5586502348 };
				const num64 b1 = num_primes[index];
				if (b > b1)
					b = b1;
			}

			double x = log(static_cast<double>(t)) - 1.0;
			for (int i = 0; i <= 1; ++i, x -= 0.1)
			{
				const num64 c = static_cast<num64>(t / x);
				if (a < c && c < b)
				{
					if (t < GetNthPrime(c))
						b = c;
					else
						a = c + 1;
				}
			}
		}

		// Run the search
		while (a < b)
		{
			const num64 c = (a + b) >> 1;
			if (t < GetNthPrime(c))
				b = c;
			else
				a = c + 1;
		}

		const num64 index_begin = r.index_start_prime;
		if (a > index_begin)
			numbersProcessedTotal += a - index_begin;

		stop |= !RangeGen::Iterate(r) || (RangeGen::cur_largest_prime_power != largestPrimePower);
	}

	return numbersProcessedTotal;
}

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		return 0;
	}

	SetPriorityClass(GetCurrentProcess(), IDLE_PRIORITY_CLASS);
	SetThreadPriorityBoost(GetCurrentThread(), true);
	SetThreadPriority(GetCurrentThread(), THREAD_MODE_BACKGROUND_BEGIN);
	SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_IDLE);

	PrimeTablesInit(0, 0, "");

	std::ifstream f(argv[1]);
	std::vector<char> data;
	std::vector<char*> params;
	data.reserve(256);
	params.reserve(8);

	std::string s;

	while (!f.eof())
	{
		std::getline(f, s);
		data.assign(s.begin(), s.end());

		char* c = data.data();

		params.clear();
		for (size_t i = 0, n = data.size(); i < n; ++i)
		{
			if (data[i] == ' ')
			{
				data[i] = '\0';
				params.emplace_back(c);
				c = data.data() + i + 1;
			}
		}
		params.emplace_back(c);

		num64 largestPrimePower = 1;
		if (!strcmp(params[4], "/lpp"))
			largestPrimePower = atoi(params[5]);

		const num64 numbersProcessedTotal = RunRanges(params[1], params[3], largestPrimePower);
		RunRanges(params[1], params[3], largestPrimePower, (numbersProcessedTotal / 10) + 1);
	}

	return 0;
}
