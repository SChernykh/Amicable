#include "stdafx.h"
#include "PrimeTables.h"
#include "Engine.h"
#include "sprp64.h"
#include <fstream>

const char* FilesToCheck[] = { "c2_3.txt", "c2_4.txt", "c2_5.txt", "c2_6.txt", "c2_7.txt", "c2_8.txt", "c2_9.txt", "c2_10.txt", "c2_11.txt", "c2_12.txt", "c2_13.txt", "c2_14.txt", "c2_15.txt", "c2_16.txt", "c2_17.txt", "c2_18.txt", "c2_19.txt", "c2_20.txt", "c2_20_1.txt", "c2_20_2.txt", "c2_21.txt", "c2_21_1.txt", "c2_21_2.txt", "c2_21_3.txt", "c2_21_4.txt", "c2_22.txt" };

NOINLINE bool TestCheckPair()
{
	// (all pairs below the search limit; empty lines between the records are skipped, the count is printed)
	num64 numChecked = 0;
	for (const char* name : FilesToCheck)
	{
		std::ifstream f(name);
		while (!f.eof())
		{
			char buf[1024];
			buf[0] = '\0';
			f.getline(buf, sizeof(buf)); // type author year
			if (buf[0] == '\0')
			{
				if (f.eof() || f.fail())
					break;
				continue;
			}

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
			CheckPair128NoInline(m, sum);
			if (GetNumFoundPairsInThisThread() != oldNumPairs + 1)
			{
				std::cerr << "CheckPair128 didn't recognize " << m << " as a valid amicable number" << std::endl;
				return false;
			}
			++numChecked;
		}
	}
	std::cout << numChecked << " pairs...";
	return numChecked > 0;
}

template<typename T>
NOINLINE bool TestGeneric(T test)
{
	std::vector<std::pair<num64, num64>> factorization;
	for (const char* name : FilesToCheck)
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
				if (!test(m, factorization))
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
	return TestGeneric([](num128 m, const std::vector<std::pair<num64, num64>>& factorization)
	{
		const std::pair<num64, num64>& max_factor = factorization.back();
		if (max_factor.first >= SearchLimit::LinearLimit)
		{
			for (num64 i = 0; i < max_factor.second; ++i)
			{
				m /= max_factor.first;
			}
			auto it = std::lower_bound(CandidatesData.begin(), CandidatesData.end(), m, [](const AmicableCandidate& a, num128 b)
			{
				return b > a.GetValue();
			});
			if ((it == CandidatesData.end()) || (m != it->GetValue()))
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

	PrimeIterator it1;
	for (int i = 0; i < N; ++i, ++it1)
	{
		const num64 p1 = it1.Get();
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
	}

	it1 = PrimeIterator();
	for (int i = 0; i < N; ++i, ++it1)
	{
		const num64 p1 = it1.Get();
		PrimeIterator it2 = it1;
		++it2;
		for (int j = 0; j < N; ++j, ++it2)
		{
			const num64 p2 = it2.Get();
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
		}
	}

	it1 = PrimeIterator();
	for (int i = 0; i < N; ++i, ++it1)
	{
		const num64 p1 = it1.Get();
		PrimeIterator it2 = it1;
		++it2;
		for (int j = 0; j < N; ++j, ++it2)
		{
			const num64 p2 = it2.Get();
			PrimeIterator it3 = it2;
			++it3;
			for (int k = 0; k < N; ++k, ++it3)
			{
				const num64 p3 = it3.Get();
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
			}
		}
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

	a = NUM128_MAX;

	b = 1;
	c = a / b;
	if (c != a)
	{
		std::cerr << "(2^128 - 1) / 1 division check failed" << std::endl;
		return false;
	}

	b = num64(-1);
	c = a / b;
	if ((LowWord(c) != 1) || (HighWord(c) != 1))
	{
		std::cerr << "(2^128 - 1) / (2^64 - 1) division check failed" << std::endl;
		return false;
	}

	b = CombineNum128(1, 1);
	c = a / b;
	if ((LowWord(c) != num64(-1)) || (HighWord(c) != 0))
	{
		std::cerr << "(2^128 - 1) / (2^64 + 1) division check failed" << std::endl;
		return false;
	}

	b = 6700417;
	c = a / b;
	if ((LowWord(c) != 2753074036095) || (HighWord(c) != 2753074036095))
	{
		std::cerr << "(2^128 - 1) / 6700417 division check failed" << std::endl;
		return false;
	}

	b = a;
	c = a / b;
	if ((LowWord(c) != 1) || (HighWord(c) != 0))
	{
		std::cerr << "(2^128 - 1) / (2^128 - 1) division check failed" << std::endl;
		return false;
	}

	a = CombineNum128(num64(-2), num64(-1));
	b = CombineNum128(1, 1);
	c = a / b;
	if ((LowWord(c) != num64(-2)) || (HighWord(c) != 0))
	{
		std::cerr << "(2^128 - 2) / (2^64 + 1) division check failed" << std::endl;
		return false;
	}

	a = CombineNum128(0, (num64(1) << 63) + 1);
	b = CombineNum128(0, 1);
	c = a / b;
	if ((LowWord(c) != (num64(1) << 63) + 1) || (HighWord(c) != 0))
	{
		std::cerr << "(2^127 + 1) / 2^64 division check failed" << std::endl;
		return false;
	}

	return true;
}

// Regression test of InitialCheck128 (2026-11-08): the other member M = c * q^e of a pair whose largest prime power q^e has sigma(q^e) >= 2^64
// (q in the 128-bit trial division table, q <= cbrt(10^22), so e >= 3; c: primes below q) was rejected after its last trial division in
// the 128-bit loop (the cofactor 1 with the target sum 1: "n2 >= targetSum"). Synthetic pairs (n1 = sigma(M) - M, sigma(M)) through
// CheckPair128: n1 must be reported exactly once; with the target sum 3 sigma(M) (the cofactor 1 with the target sum 3) nothing. Controls:
// sigma(q^e) < 2^64 (decided by the 64-bit code).
NOINLINE bool TestCheckPairPowerCofactor()
{
	struct Case
	{
		num64 q;
		unsigned int e;
		num64 c;
	};
	const Case cases[] = {
		{ 3, 41, 2 }, { 569, 7, 2 }, { 1627, 6, 1 }, { 7151, 5, 1 }, { 65537, 4, 5 }, { 2642257, 3, 1 }, { 10000019, 3, 2 },
		// the largest prime of the 128-bit table (#1362209)
		{ 21544463, 3, 1 }, { 21544463, 3, 2 },
		// controls: sigma(q^3) < 2^64
		{ 1000003, 3, 1 }, { 1000003, 3, 10 },
	};

	auto sigmaSmall = [](num64 n)
	{
		num64 s = 1;
		for (num64 p = 2; p * p <= n; ++p)
		{
			if (n % p == 0)
			{
				num64 t = 1, pk = 1;
				while (n % p == 0)
				{
					n /= p;
					pk *= p;
					t += pk;
				}
				s *= t;
			}
		}
		if (n > 1)
		{
			s *= n + 1;
		}
		return s;
	};

	std::vector<num128> found;
	std::vector<num128>* oldLog = g_FoundNumbersLog;
	g_FoundNumbersLog = &found;
	bool ok = true;
	char buf[2][40];
	for (int pass = 0; pass < 2; ++pass)
	{
		for (const Case& c : cases)
		{
			num128 qe = 1;
			num128 sigmaQe = 1;
			for (unsigned int i = 0; i < c.e; ++i)
			{
				qe *= c.q;
				sigmaQe += qe;
			}
			const num128 M = qe * c.c;
			const num128 sigmaM = sigmaQe * sigmaSmall(c.c);
			// pass 0: the pair (n1, M); pass 1: the target sum 3 sigma(M) for the same M, not a match
			const num128 targetSum = (pass == 0) ? sigmaM : (sigmaM * 3);
			const num128 n1 = targetSum - M;
			found.clear();
			CheckPair128NoInline(n1, targetSum);
			const bool expected = (pass == 0);
			if ((found.size() != (expected ? 1U : 0U)) || (expected && (found[0] != n1)))
			{
				std::cerr << "CheckPair128: M = " << c.c << " * " << c.q << "^" << c.e << " = " << itoa128(M, buf[0], sizeof(buf[0])) << ", n1 = " << itoa128(n1, buf[1], sizeof(buf[1]))
					<< ": " << found.size() << " numbers reported, expected " << (expected ? 1 : 0) << std::endl;
				ok = false;
			}
		}
		// (the negative cases only with a correct check: before the fix they reached CheckPairInternal with n2 = 1)
		if (!ok)
		{
			break;
		}
	}
	g_FoundNumbersLog = oldLog;
	return ok;
}

// Regression test of OverAbundant (2026-11-08): the comparison (sum - value) (sigma(g) - g) >= value g is exact now (256-bit products); the
// 128-bit products wrapped around from value g >= 2^128 (m >= ~3e18 with a large g), and a wrapped "true" made the walk skip a range.
// Synthetic factorizations and sums (OverAbundant only uses the sum's divisibility by the prime powers of m); the expected results come from
// exact integers (a Python model of the function). The old code gave the wrong result in 9 of the 15 cases, a wrong "true" in 3 of them.
NOINLINE bool TestOverAbundantExact()
{
	struct Case
	{
		num64 p[3];
		unsigned int k[3];
		int n;
		const char* value;
		const char* sum;
		num64 sum_coeff;
		num64 max_factor;
		bool expected;
	};
	const Case cases[] = {
		{ { 2 }, { 70 }, 1, "1180591620717411303424", "2361183241434822606848", 1, 1, false }, // m = 2^70, sum = 2^71; old: true
		{ { 2 }, { 70 }, 1, "1180591620717411303424", "3541774862152233910272", 1, 1, true }, // m = 2^70, sum = 3 * 2^70; old: true
		{ { 2 }, { 64 }, 1, "18446744073709551616", "36893488147419103232", 1, 1, false }, // m = 2^64, sum = 2^65; old: true
		{ { 2 }, { 64 }, 1, "18446744073709551616", "36893488147419103232", 2, 2, false }, // m = 2^64, sum = 2^65, sum_coeff 2; old: true
		{ { 5 }, { 30 }, 1, "931322574615478515625", "1862645149230957031250", 1, 2, false }, // m = 5^30, sum = 2 * 5^30; old: false
		{ { 5 }, { 30 }, 1, "931322574615478515625", "2235174179077148437500", 10, 5, false }, // m = 5^30, sum = 12 * 5^29, sum_coeff 10; old: false
		{ { 2, 5 }, { 30, 17 }, 2, "819200000000000000000", "1638400000000000000000", 1, 2, true }, // m = 2^30 5^17, sum = 2^31 5^17; old: true
		{ { 2, 5, 7 }, { 3, 1, 1 }, 3, "280", "720", 2, 2, true }, // m = 280 (no overflow); old: true
		{ { 2, 7 }, { 5, 1 }, 2, "224", "504", 1, 1, true }, // m = 224 (no overflow); old: true
		{ { 2, 5 }, { 43, 12 }, 2, "2147483648000000000000", "10881568079872000000000", 1, 1, true }, // old: false
		{ { 2, 31 }, { 22, 10 }, 2, "3437770202596721557504", "1210552556542617326190592", 1, 2, true }, // old: false
		{ { 17, 19 }, { 13, 3 }, 2, "67935500727701821883", "41834137864112960300976336", 1, 1, true }, // old: false
		{ { 5, 43 }, { 2, 12 }, 2, "998990769931564410025", "250525906292667791181659475", 2, 2, true }, // old: false
		{ { 2, 19 }, { 47, 6 }, 2, "6621119129403646803968", "26164600472781249839104", 1, 1, true }, // old: false
		{ { 2, 7 }, { 31, 14 }, 2, "1456472958639540273152", "4048697585638607767470080", 1, 1, true }, // old: false
	};

	bool ok = true;
	for (const Case& c : cases)
	{
		Factor f[3];
		for (int i = 0; i < c.n; ++i)
		{
			f[i].p = PrimeIterator(c.p[i]);
			f[i].k = c.k[i];
			int index = 0;
			for (PrimeIterator it; it.Get() < c.p[i]; ++it)
			{
				++index;
			}
			f[i].index = index;
			f[i].p_inv = num64(0) - modular_inverse64(c.p[i]);
			f[i].p_inv128 = -modular_inverse128(c.p[i]);
		}
		const num128 value = atoi128(c.value);
		const num128 sum = atoi128(c.sum);
		byte result;
		if (c.max_factor == 1)
		{
			result = OverAbundant<1>(f, c.n - 1, value, sum, c.sum_coeff);
		}
		else if (c.max_factor == 2)
		{
			result = OverAbundant<2>(f, c.n - 1, value, sum, c.sum_coeff);
		}
		else
		{
			result = OverAbundant<5>(f, c.n - 1, value, sum, c.sum_coeff);
		}
		if ((result != 0) != c.expected)
		{
			std::cerr << "OverAbundant<" << c.max_factor << ">(m = " << c.value << ", sum = " << c.sum << ", sum_coeff = " << c.sum_coeff << ") = " << (result != 0) << ", expected " << c.expected << std::endl;
			ok = false;
		}
	}
	return ok;
}
