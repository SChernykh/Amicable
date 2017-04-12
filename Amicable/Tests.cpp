#include "stdafx.h"
#include "PrimeTables.h"
#include "Engine.h"
#include <fstream>

NOINLINE bool TestCheckPair()
{
	const char* files[] = { "c2_3.txt", "c2_4.txt", "c2_5.txt", "c2_6.txt", "c2_7.txt", "c2_8.txt", "c2_9.txt", "c2_10.txt", "c2_11.txt", "c2_12.txt", "c2_13.txt", "c2_14.txt", "c2_15.txt", "c2_16.txt", "c2_17.txt", "c2_18.txt", "c2_19.txt", "c2_20.txt" };
	for (const char* name : files)
	{
		std::ifstream f(name);
		while (!f.eof())
		{
			char buf[1024];
			buf[0] = '\0';
			f.getline(buf, sizeof(buf)); // type author year
			if (buf[0] == '\0')
				break;

			buf[0] = '\0';
			f.getline(buf, sizeof(buf)); // M and factorization

			const num128 m = atoi128(buf);
			if (m >= SearchLimit::value)
				break;

			buf[0] = '\0';
			f.getline(buf, sizeof(buf)); // N and factorization

			const num128 n = atoi128(buf);

			buf[0] = '\0';
			f.getline(buf, sizeof(buf)); // empty line

			const num64 oldNumPairs = GetNumFoundPairsInThisThread();
			const num128 sum = m + n;
			CheckPair128NoInline(m.lo, sum.lo, sum.hi);
			if (GetNumFoundPairsInThisThread() != oldNumPairs + 1)
			{
				std::cerr << "CheckPair128 didn't recognize " << m.lo << " as a valid amicable num64" << std::endl;
				return false;
			}
			if (!sum.hi)
			{
				CheckPairNoInline(m.lo, sum.lo);
				if (GetNumFoundPairsInThisThread() != oldNumPairs + 2)
				{
					std::cerr << "CheckPair didn't recognize " << m.lo << " as a valid amicable num64" << std::endl;
					return false;
				}
			}
		}
	}
	return true;
}

template<typename T>
NOINLINE bool TestGeneric(T test)
{
	const char* files[] = { "c2_3.txt", "c2_4.txt", "c2_5.txt", "c2_6.txt", "c2_7.txt", "c2_8.txt", "c2_9.txt", "c2_10.txt", "c2_11.txt", "c2_12.txt", "c2_13.txt", "c2_14.txt", "c2_15.txt", "c2_16.txt", "c2_17.txt", "c2_18.txt", "c2_19.txt", "c2_20.txt" };
	std::vector<std::pair<num64, num64>> factorization;
	for (const char* name : files)
	{
		std::ifstream f(name);
		while (!f.eof())
		{
			char buf[1024];
			buf[0] = '\0';
			f.getline(buf, sizeof(buf)); // type author year
			if (buf[0] == '\0')
				break;

			buf[0] = '\0';
			f.getline(buf, sizeof(buf)); // M and factorization

			num128 m = atoi128(buf);

			if (m < SearchLimit::value)
			{
				const size_t len = strlen(buf);
				const char* factor_begin = nullptr;
				factorization.clear();
				for (size_t i = 0; i <= len; ++i)
				{
					if (!isdigit(buf[i]))
					{
						if (factor_begin)
						{
							num64 p;
							if (sscanf_s(factor_begin, "%llu", &p) == 1)
							{
								num64 factor_power = 1;
								if (buf[i] == '^')
								{
									sscanf_s(buf + i + 1, "%llu", &factor_power);
								}
								factorization.emplace_back(std::pair<num64, num64>(p, factor_power));
							}
							factor_begin = nullptr;
						}
						if ((buf[i] == '=') || (buf[i] == '*'))
						{
							factor_begin = buf + i + 1;
						}
					}
				}
				std::sort(factorization.begin(), factorization.end(), [](const std::pair<num64, num64>& a, const std::pair<num64, num64>& b)
				{
					return a.first < b.first;
				});
				if (!test(m.lo, factorization))
				{
					return false;
				}
			}

			buf[0] = '\0';
			f.getline(buf, sizeof(buf)); // N and factorization

			buf[0] = '\0';
			f.getline(buf, sizeof(buf)); // empty line
		}
	}
	return true;
}

NOINLINE bool TestAmicableCandidates()
{
	return TestGeneric([](num64 m, const std::vector<std::pair<num64, num64>>& factorization)
	{
		const std::pair<num64, num64>& max_factor = factorization.back();
		if (max_factor.first >= SearchLimit::LinearLimit)
		{
			for (num64 i = 0; i < max_factor.second; ++i)
			{
				m /= max_factor.first;
			}
			auto it = std::lower_bound(CandidatesData.begin(), CandidatesData.end(), m, [](const AmicableCandidate& a, num64 b)
			{
				return a.value < b;
			});
			if ((it == CandidatesData.end()) || (it->value != m))
			{
				std::cerr << "LinearSearchData doesn't include " << m << " as a valid amicable candidate" << std::endl;
				return false;
			}

			const unsigned int mask = CandidatesDataMask[Mod385(max_factor.first + 1)];
			if (it->is_over_abundant_mask & mask)
			{
				std::cerr << "LinearSearchData has invalid is_over_abundant_mask for " << m << std::endl;
				return false;
			}
		}
		return true;
	});
}

NOINLINE bool TestMaximumSumOfDivisors3()
{
	double minEstimate = 1e10;
	int N = 1000;

	const byte* shift1 = NextPrimeShifts;
	num64 p1 = 2;
	for (int i = 0; i < N; ++i)
	{
		num64 a = p1 * p1 * p1;
		num64 sum = p1 * p1 * p1 + p1 * p1 + p1 + 1;
		num64 sumEstimate = MaximumSumOfDivisors3NoInline(a, p1, p1 * p1);
		// sumEstimate = a + (a_div_p0 + p0 * p0 + p0) * 2 + 1
		// sumEstimate = a + (a / p1 + p1 * p1 + p1) * 2 + 1
		// sumEstimate = a + (p1 * p1 + p1 * p1 + p1) * 2 + 1
		// sumEstimate = a + (p1 * p1 * 2 + p1) * 2 + 1
		// sumEstimate = a + p1 * p1 * 4 + p1 * 2 + 1
		// sumEstimate - sum = (a + p1^2 * 4 + p1 * 2 + 1) - (1 + p1 + p1^2 + p1^3)
		// sumEstimate - sum = (a + p1^2 * 4 + p1 * 2 + 1) - (1 + p1 + p1^2 + a)
		// sumEstimate - sum = (p1^2 * 4 + p1 * 2) - (p1 + p1^2)
		// sumEstimate - sum = p1^2 * 3 + p1
		// sumEstimate = sum + p1^2 * 3 + p1
		// therefore, sumEstimate > sum must always be true
		if (sumEstimate < sum)
		{
			std::cerr << "MaximumSumOfDivisors3 test failed" << std::endl;
			return false;
		}
		const double curEstimate = static_cast<double>(sumEstimate - sum) / sum;
		if (curEstimate < minEstimate)
			minEstimate = curEstimate;
		p1 += *shift1 ? *shift1 * ShiftMultiplier : 1ULL;
		shift1 += 2;
	}

	shift1 = NextPrimeShifts;
	p1 = 2;
	for (int i = 0; i < N; ++i)
	{
		const byte* shift2 = shift1;
		num64 p2 = p1 + (*shift2 ? *shift2 * ShiftMultiplier : 1ULL);
		shift2 += 2;
		for (int j = 0; j < N; ++j)
		{
			{
				num64 a = p1 * (p2 * p2);
				num64 sum = (p1 + 1) * (p2 * p2 + p2 + 1);
				num64 sumEstimate = MaximumSumOfDivisors3NoInline(a, p1, p2 * p2);
				// sumEstimate = a + (a_div_p0 + p0 * p0 + p0) * 2 + 1
				// sumEstimate = a + (p2 * p2 + p1 * p1 + p1) * 2 + 1
				// sum = (p1 + 1) * (1 + p2 + p2^2)
				// sum = p1 * (1 + p2 + p2^2) + (1 + p2 + p2^2)
				// sum = p1 + p1 * p2 + p1 * p2^2 + 1 + p2 + p2^2
				// sum = p1 + p1 * p2 + a + 1 + p2 + p2^2
				// sumEstimate - sum = (a + (p2 * p2 + p1 * p1 + p1) * 2 + 1) - (p1 + p1 * p2 + a + 1 + p2 + p2^2)
				// sumEstimate - sum = ((p2 * p2 + p1 * p1 + p1) * 2) - (p1 + p1 * p2 + p2 + p2^2)
				// sumEstimate - sum = (p2 * p2 * 2 + p1 * p1 * 2 + p1 * 2) - (p1 + p1 * p2 + p2 + p2^2)
				// sumEstimate - sum = (p2 * p2 * 2 + p1 * p1 * 2 + p1) - (p1 * p2 + p2 + p2^2)
				// sumEstimate - sum = (p2^2 * 2 + p1 * p1 * 2 + p1) - (p1 * p2 + p2 + p2^2)
				// sumEstimate - sum = (p2^2 + p1 * p1 * 2 + p1) - (p1 * p2 + p2)
				// sumEstimate - sum = p2^2 + p1 * p1 * 2 + p1 - p1 * p2 - p2
				// sumEstimate - sum = p2 * (p2 - p1) + p1 * p1 * 2 + p1 - p2
				// sumEstimate - sum = p2 * (p2 - p1) + p1 * p1 * 2 - (p2 - p1)
				// sumEstimate - sum = (p2 - p1) * (p2) + p1 * p1 * 2 - (p2 - p1)
				// sumEstimate - sum = (p2 - p1) * (p2 - 1) + p1 * p1 * 2
				// since p2 > p1 > 1 then all members are > 0
				// therefore, sumEstimate > sum must always be true
				if (sumEstimate < sum)
				{
					std::cerr << "MaximumSumOfDivisors3 test failed" << std::endl;
					return false;
				}
				const double curEstimate = static_cast<double>(sumEstimate - sum) / sum;
				if (curEstimate < minEstimate)
					minEstimate = curEstimate;
			}
			{
				num64 a = p1 * (p1 * p2);
				num64 sum = (p1 * p1 + p1 + 1) * (p2 + 1);
				num64 sumEstimate = MaximumSumOfDivisors3NoInline(a, p1, p1 * p2);
				// sumEstimate = a + (a_div_p0 + p0 * p0 + p0) * 2 + 1
				// sumEstimate = a + (p1 * p2 + p1 * p1 + p1) * 2 + 1
				// sum = (p1 * p1 + p1 + 1) * (p2 + 1)
				// sum = (p1 * p1 * p2) + (p1 * p2 + p2) + (p1 * p1 + p1 + 1)
				// sum = a + p1 * p2 + p2 + p1 * p1 + p1 + 1
				// sumEstimate - sum = (a + (p1 * p2 + p1 * p1 + p1) * 2 + 1) - (a + p1 * p2 + p2 + p1 * p1 + p1 + 1)
				// sumEstimate - sum = (p1 * p2 + p1 * p1 + p1) * 2 - (p1 * p2 + p2 + p1 * p1 + p1)
				// sumEstimate - sum = p1 * p2 * 2 + p1 * p1 * 2 + p1 * 2 - (p1 * p2 + p2 + p1 * p1 + p1)
				// sumEstimate - sum = p1 * p2 + p1 * p1 + p1 - p2
				// sumEstimate - sum = p2 * p1 + p1 * p1 + p1 - p2
				// sumEstimate - sum = p2 * (p1 - 1) + p1 * p1 + p1
				// since p1 > 1 then all members are > 0
				// therefore, sumEstimate > sum must always be true
				if (sumEstimate < sum)
				{
					std::cerr << "MaximumSumOfDivisors3 test failed" << std::endl;
					return false;
				}
				const double curEstimate = static_cast<double>(sumEstimate - sum) / sum;
				if (curEstimate < minEstimate)
					minEstimate = curEstimate;
			}
			p2 += *shift2 ? *shift2 * ShiftMultiplier : 1ULL;
			shift2 += 2;
		}
		p1 += *shift1 ? *shift1 * ShiftMultiplier : 1ULL;
		shift1 += 2;
	}

	shift1 = NextPrimeShifts;
	p1 = 2;
	for (int i = 0; i < N; ++i)
	{
		const byte* shift2 = shift1;
		num64 p2 = p1 + (*shift2 ? *shift2 * ShiftMultiplier : 1ULL);
		shift2 += 2;
		for (int j = 0; j < N; ++j)
		{
			const byte* shift3 = shift2;
			num64 p3 = p2 + (*shift3 ? *shift3 * ShiftMultiplier : 1ULL);
			shift3 += 2;
			for (int k = 0; k < N; ++k)
			{
				num64 a = p1 * p2 * p3;
				num64 sum = (p1 + 1) * (p2 + 1) * (p3 + 1);
				num64 sumEstimate = MaximumSumOfDivisors3NoInline(a, p1, p2 * p3);
				// sumEstimate - a - 1 = (p2 * p3 + p1 * p1 + p1) * 2
				// sum - a - 1 = p1 * p2 + p1 * p3 + p2 * p3 + p1 + p2 + p3
				//
				// sumEstimate - sum = (p2 * p3 + p1 * p1 + p1) * 2 - (p1 * p2 + p1 * p3 + p2 * p3 + p1 + p2 + p3)
				// sumEstimate - sum = p2 * p3 * 2 + p1 * p1 * 2 + p1 * 2 - (p1 * p2 + p1 * p3 + p2 * p3 + p1 + p2 + p3)
				// sumEstimate - sum = p2 * p3 + p1 * p1 * 2 + p1 - (p1 * p2 + p1 * p3 + p2 + p3)
				// sumEstimate - sum = p2 * p3 + p1 * p1 * 2 + p1 - p1 * p2 - p1 * p3 - p2 - p3
				// sumEstimate - sum = p2 * p3 - p1 * p3 + p1 * p1 * 2 + p1 - p1 * p2 - p2 - p3
				// sumEstimate - sum = (p2 - p1) * p3 + p1 * p1 * 2 + p1 - p1 * p2 - p2 - p3
				// sumEstimate - sum = (p2 - p1 - 1) * p3 + p1 * p1 * 2 + p1 - p1 * p2 - p2
				// sumEstimate - sum = ((p2 - p1 - 1) * p3 - p2) + p1 * (p1 * 2 + 1 - p2)
				//
				// let p2 = p1 + c1, then p2 - p1 = c1 >= 2
				// sumEstimate - sum = ((c1 - 1) * p3 - p1 - c1) + p1 * (p1 * 2 + 1 - p1 - c1)
				// sumEstimate - sum = ((c1 - 1) * p3 - p1 - c1) + p1 * (p1 + 1 - c1)
				// sumEstimate - sum = c1 * p3 - p3 - p1 - c1 + p1 * (p1 + 1 - c1)
				// sumEstimate - sum = c1 * (p3 - p1) - p3 - p1 - c1 + p1 * (p1 + 1)
				//
				// let p3 = p1 + c2, then p3 - p1 = c2 >= 4
				// sumEstimate - sum = c1 * c2 - p1 - c2 - p1 - c1 + p1 * p1 + p1
				// sumEstimate - sum = c1 * c2 - c2 - c1 + p1 * (p1 - 1)
				// sumEstimate - sum = p1 * (p1 - 1) + c1 * c2 - c2 - c1
				// sumEstimate - sum = p1 * (p1 - 1) + c2 * c1 - c2 - c1
				// sumEstimate - sum = p1 * (p1 - 1) + c2 * (c1 - 1) - c1
				//
				// c2 * (c1 - 1) - c1 > 0
				// c2 * (c1 - 1) > c1
				// c2 > c1 / (c1 - 1)
				// we know that c1 >= 2 and c2 >= 4 > 2
				// therefore, c1 / (c1 - 1) <= 2 and c2 > 2
				// c2 > c1 / (c1 - 1) is always true
				// therefore, sumEstimate > sum must always be true
				if (sumEstimate < sum)
				{
					std::cerr << "MaximumSumOfDivisors3 test failed" << std::endl;
					return false;
				}
				const double curEstimate = static_cast<double>(sumEstimate - sum) / sum;
				if (curEstimate < minEstimate)
					minEstimate = curEstimate;
				p3 += *shift3 ? *shift3 * ShiftMultiplier : 1ULL;
				shift3 += 2;
			}
			p2 += *shift2 ? *shift2 * ShiftMultiplier : 1ULL;
			shift2 += 2;
		}
		p1 += *shift1 ? *shift1 * ShiftMultiplier : 1ULL;
		shift1 += 2;
	}
	return true;
}

NOINLINE bool TestPrimeSieve()
{
	primesieve::PrimeSieve s;
	num64 primeCount = 0;
	num64 primeSum = 0;
	s.sieveTemplated(0, 10000000000, [&primeCount, &primeSum](num64 p){ ++primeCount; primeSum += p; });
	if (primeCount != 455052511)
	{
		std::cerr << "primesieve counted " << primeCount << " primes below 10^10 (correct value is 455052511)" << std::endl;
		return false;
	}
	if (primeSum != 2220822432581729238)
	{
		std::cerr << "primesieve counted that sum of primes below 10^10 is " << primeSum << " (correct value is 2220822432581729238)" << std::endl;
		return false;
	}
	return true;
}

NOINLINE bool TestNum128Division()
{
	num128 a, b, c;

	a.lo = num64(-1);
	a.hi = num64(-1);

	b = 1;
	c = a / b;
	if (c != a)
	{
		std::cerr << "(2^128 - 1) / 1 division check failed" << std::endl;
		return false;
	}

	b = num64(-1);
	c = a / b;
	if ((c.lo != 1) || (c.hi != 1))
	{
		std::cerr << "(2^128 - 1) / (2^64 - 1) division check failed" << std::endl;
		return false;
	}

	b.lo = 1;
	b.hi = 1;
	c = a / b;
	if ((c.lo != num64(-1)) || (c.hi != 0))
	{
		std::cerr << "(2^128 - 1) / (2^64 + 1) division check failed" << std::endl;
		return false;
	}

	b = 6700417;
	c = a / b;
	if ((c.lo != 2753074036095) || (c.hi != 2753074036095))
	{
		std::cerr << "(2^128 - 1) / 6700417 division check failed" << std::endl;
		return false;
	}

	b = a;
	c = a / b;
	if ((c.lo != 1) || (c.hi != 0))
	{
		std::cerr << "(2^128 - 1) / (2^128 - 1) division check failed" << std::endl;
		return false;
	}

	a.lo = num64(-2);
	a.hi = num64(-1);
	b.lo = 1;
	b.hi = 1;
	c = a / b;
	if ((c.lo != num64(-2)) || (c.hi != 0))
	{
		std::cerr << "(2^128 - 2) / (2^64 + 1) division check failed" << std::endl;
		return false;
	}

	a.lo = 0;
	a.hi = (num64(1) << 63) + 1;
	b.lo = 0;
	b.hi = 1;
	c = a / b;
	if ((c.lo != (num64(1) << 63) + 1) || (c.hi != 0))
	{
		std::cerr << "(2^127 + 1) / 2^64 division check failed" << std::endl;
		return false;
	}

	return true;
}
