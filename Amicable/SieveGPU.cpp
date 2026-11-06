#include "stdafx.h"
#include "PrimeTables.h"
#include "RangeGen.h"
#include "SieveGPU.h"
#include "RangeReduction.h"
#include "sprp64.h"
#include <cmath>
#include <cstdint>
#include <cstring>

namespace
{
	const unsigned int W30[8] = { 1, 7, 11, 13, 17, 19, 23, 29 };

	int Class30(unsigned int x)
	{
		for (int c = 0; c < 8; ++c)
		{
			if (W30[c] == x % 30)
			{
				return c;
			}
		}
		return -1;
	}

	// a^-1 mod m for gcd(a, m) = 1, 1 < m < 2^32
	num64 InvMod(num64 a, num64 m)
	{
		int64_t t = 0;
		int64_t nt = 1;
		int64_t r = static_cast<int64_t>(m);
		int64_t nr = static_cast<int64_t>(a % m);
		while (nr)
		{
			const int64_t q = r / nr;
			const int64_t t1 = t - q * nt;
			t = nt;
			nt = t1;
			const int64_t r1 = r - q * nr;
			r = nr;
			nr = r1;
		}
		return static_cast<num64>((t < 0) ? (t + static_cast<int64_t>(m)) : t);
	}

	// m^-1 mod 2^32 for odd m
	unsigned int InvMod2_32(unsigned int m)
	{
		unsigned int x = m;
		for (int i = 0; i < 5; ++i)
		{
			x *= 2 - m * x;
		}
		return x;
	}

	long double LogSigmaRatio(long double q, int e)
	{
		if (e <= 0)
		{
			return 0.0L;
		}
		// sigma(q^e) / q^e = (q - q^-e) / (q - 1)
		return logl((q - powl(q, static_cast<long double>(-e))) / (q - 1.0L));
	}

	int FixedPoint(long double x)
	{
		return static_cast<int>(llroundl(x * static_cast<long double>(SIEVE_FX)));
	}

	unsigned int CountTrailingZeros128(const num128& x)
	{
		unsigned long index;
		if (LowWord(x))
		{
			_BitScanForward64(&index, LowWord(x));
			return index;
		}
		if (HighWord(x))
		{
			_BitScanForward64(&index, HighWord(x));
			return index + 64;
		}
		return 128;
	}

	unsigned int BitLength128(const num128& x)
	{
		unsigned long index;
		if (HighWord(x))
		{
			_BitScanReverse64(&index, HighWord(x));
			return index + 65;
		}
		if (LowWord(x))
		{
			_BitScanReverse64(&index, LowWord(x));
			return index + 1;
		}
		return 0;
	}

	// x mod q for 128-bit x and q < 2^32
	num64 Mod128(const num128& x, num64 q)
	{
		const num64 c64 = (num64(0) - q) % q;
		return ((HighWord(x) % q) * c64 + (LowWord(x) % q)) % q;
	}

	struct NumberAndSum
	{
		num128 N;
		num128 S;
	};

	// max sigma(x)/x over x < maxN made of consecutive primes after p (non-increasing exponents)
	void MaxSumRatio(PrimeIterator p, num64 maxPower, const num128& maxN, NumberAndSum cur, double& best)
	{
		++p;
		const num128 startSum = cur.S;
		for (num64 k = 1; k <= maxPower; ++k)
		{
			const num128 nextN = cur.N * p.Get();
			if (nextN >= maxN)
			{
				const double r = Num128ToDouble(cur.S) / Num128ToDouble(cur.N);
				if (r > best)
				{
					best = r;
				}
				return;
			}
			cur.N = nextN;
			cur.S = cur.S * p.Get() + startSum;
			MaxSumRatio(p, k, maxN, cur, best);
		}
	}

	double Li(double x)
	{
		if (x < 2.0)
		{
			return 0.0;
		}
		const double L = log(x);
		double sum = 0.0;
		double term = 1.0;
		for (int n = 1; n < 400; ++n)
		{
			term *= L / n;
			const double add = term / n;
			sum += add;
			if (add < 1e-15 * sum)
			{
				break;
			}
		}
		return 0.5772156649015329 + log(L) + sum;
	}

	double PrimePi(double x)
	{
		if (x < 2.0)
		{
			return 0.0;
		}
		if (x < 100.0)
		{
			static const int smallPrimes[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 };
			double n = 0;
			for (int p : smallPrimes)
			{
				if (p <= x)
				{
					n += 1;
				}
			}
			return n;
		}
		// first two terms of Riemann's R function
		return Li(x) - 0.5 * Li(sqrt(x));
	}

	struct PlannedTask
	{
		unsigned int x;
		unsigned int y;
		double hits;   // expected hits per segment
	};
}

SieveGPUConfig::SieveGPUConfig()
	: wgSize(1024)
	, segJ(1024)
	, chunkSegments(18)
	, hitsPerTask(16)
	, wheelThreshold(0)
		, numTiers(4)
		, highTierW(0.03)
		, killQMax(SieveKillQMax)
		, survLocal(SieveSurvLocal)
	{
		// Measured on typical work units (RTX 5070): a larger Q has fewer survivors but a costlier setup and more task visits per segment.
		// Q = 32749 only pays off for ranges with a narrow window (W < highTierW), see ChooseTier.
		tierQ[0] = 2039;
		tierQ[1] = 8191;
		tierQ[2] = 16381;
		tierQ[3] = 32749;
		tierMinPrimes[0] = 0;
		tierMinPrimes[1] = 2048;
		tierMinPrimes[2] = 16384;
		tierMinPrimes[3] = 16384;
	}

bool SieveGPUTables::Init(const SieveGPUConfig& config)
{
	cfg = config;
	killQMax = cfg.killQMax;
	// (segJ <= 4096: the kernel packs the row offset c << SIEVE_SEG_SHIFT of a row task above bit 16 of its geometry, below the wheel
	// flag in bit 31)
	if ((cfg.numTiers < 1) || (cfg.numTiers > SieveMaxTiers) || (cfg.segJ < 64) || (cfg.segJ > 4096) || (cfg.segJ & (cfg.segJ - 1)) || (cfg.wgSize < 32) || ((cfg.segJ * 8) % cfg.wgSize) || (cfg.chunkSegments < 1) ||
		(cfg.chunkSegments > SieveMaxChunkSegments) ||
		(cfg.segJ * 8 > cfg.wgSize * 32)) // the scan marks at most 32 positions per work-item
	{
		std::cerr << "Invalid sieve configuration" << std::endl;
		return false;
	}
	for (unsigned int t = 1; t < cfg.numTiers; ++t)
	{
		if ((cfg.tierQ[t] <= cfg.tierQ[t - 1]) || (cfg.tierMinPrimes[t] < cfg.tierMinPrimes[t - 1]))
		{
			std::cerr << "Sieve tiers must be increasing" << std::endl;
			return false;
		}
	}
	// wheel moduli must be >= 263 (their increments are stored in 18 bits)
	const unsigned int wheelThreshold = std::max(263U, cfg.wheelThreshold ? cfg.wheelThreshold : (cfg.segJ * 3 / 4));
	const unsigned int Qmax = cfg.tierQ[cfg.numTiers - 1];

	myPrimes.clear();
	for (PrimeIterator it(3); it.Get() <= Qmax; ++it)
	{
		myPrimes.push_back(static_cast<unsigned int>(it.Get()));
	}

	// T tables
	Ttab.assign(SieveT2Size + myPrimes.size() * SieveTStride, 0);
	for (int e = 0; e < SieveT2Size; ++e)
	{
		Ttab[e] = FixedPoint(LogSigmaRatio(2.0L, e));
	}
	for (size_t i = 0; i < myPrimes.size(); ++i)
	{
		for (int e = 0; e < SieveTStride; ++e)
		{
			Ttab[SieveT2Size + i * SieveTStride + e] = FixedPoint(LogSigmaRatio(static_cast<long double>(myPrimes[i]), e));
		}
	}

	// moduli
	primeRecs.clear();
	modRecs.clear();
	modRecs2.clear();
	myMaxPower.clear();
	myMaxTrunc = 0.0;
	unsigned int numRowModuli = 0;
	std::vector<unsigned int> tierNumPrimes(cfg.numTiers), tierNumModuli(cfg.numTiers), tierNumRowModuli(cfg.numTiers);
	std::vector<unsigned int> modulusPrime;
	std::vector<unsigned int> modulusPower;
	for (size_t i = 0; i < myPrimes.size(); ++i)
	{
		const unsigned int q = myPrimes[i];
		const unsigned int firstModulus = static_cast<unsigned int>(modRecs.size() / 4);
		unsigned int numModuli = 0;
		unsigned int maxPower = 1;

		// powers q^k <= SieveMaxModulus: k >= 1 for q >= 11, k >= 2 for q = 3, 5 (the first power is constant in a row) and q = 7 (constant
		// in a sub-sieve). The J0 tables of the kernel have the layout of the 30-wheel (stride q^k in a row for q = 7).
		num64 mk = q;
		unsigned int k = 1;
		for (;;)
		{
			if ((q > 7) || (k >= 2))
			{
				const unsigned int mm = static_cast<unsigned int>(mk);
				const unsigned int stride = (q > 5) ? mm : (mm / q);
				const unsigned int inv = static_cast<unsigned int>((q > 5) ? InvMod(30, mm) : InvMod(30 / q, stride));
				modRecs.push_back(mm);
				modRecs.push_back(0xFFFFFFFFU / mm);
				modRecs.push_back(stride);
				modRecs.push_back(inv);

				const bool wheel = (q > 7) && (mm >= wheelThreshold);
				// q; row moduli: INC index | the increment T(q, k) - T(q, k - 1) for q not dividing g << 8 (the GPU loads T only if q | g)
				const int inc0 = Ttab[SieveT2Size + i * SieveTStride + k] - Ttab[SieveT2Size + i * SieveTStride + k - 1];
				if (!wheel && ((numRowModuli >= 256) || (inc0 < 0) || (inc0 >= (1 << 24))))
				{
					std::cerr << "Row modulus " << mm << " does not fit into the modulus records" << std::endl;
					return false;
				}
				modRecs2.push_back(q);
				modRecs2.push_back(wheel ? 0xFFFFFFFFU : (numRowModuli | (static_cast<unsigned int>(inc0) << 8)));
				// type | class << 8 | prime index << 16 | k << 28 (the setup works per modulus: kernel_sieve.cl, SieveSetupModulus)
				if ((i >= 4096) || (k >= 16))
				{
					std::cerr << "Too many sieving primes or powers for the modulus records" << std::endl;
					return false;
				}
				modRecs2.push_back((wheel ? 1U : 0U) | (wheel ? (static_cast<unsigned int>(Class30(mm)) << 8) : 0U) | (static_cast<unsigned int>(i) << 16) | (k << 28));
				// q > 5: R(q), q = 3, 5: R(stride)
				modRecs2.push_back((q > 5) ? (0xFFFFFFFFU / q) : (0xFFFFFFFFU / stride));
				if (!wheel)
				{
					++numRowModuli;
				}
				modulusPrime.push_back(static_cast<unsigned int>(i));
				modulusPower.push_back(k);
				++numModuli;
				maxPower = k;
			}
			if (mk * q > SieveMaxModulus)
			{
				break;
			}
			mk *= q;
			++k;
		}
		myMaxPower.push_back(maxPower);

		// largest possible missing contribution if q^(maxPower + 1) divides A1 * p + B1
		const double trunc = -log1p(-pow(static_cast<double>(q), -static_cast<double>(maxPower + 1)));
		if (trunc > myMaxTrunc)
		{
			myMaxTrunc = trunc;
		}

		primeRecs.push_back(q);
		primeRecs.push_back(0xFFFFFFFFU / q);
		primeRecs.push_back(firstModulus | (numModuli << 16));
		primeRecs.push_back(static_cast<unsigned int>(SieveT2Size + i * SieveTStride));

		for (unsigned int t = 0; t < cfg.numTiers; ++t)
		{
			if (q <= cfg.tierQ[t])
			{
				tierNumPrimes[t] = static_cast<unsigned int>(i + 1);
				tierNumModuli[t] = static_cast<unsigned int>(modRecs.size() / 4);
				tierNumRowModuli[t] = numRowModuli;
			}
		}
	}
	maxModuli = static_cast<unsigned int>(modRecs.size() / 4);
	maxRowModuli = numRowModuli;
	if ((maxModuli > 4096) || (maxRowModuli > 4096))
	{
		std::cerr << "Too many sieve moduli" << std::endl;
		return false;
	}

	// Wheel tables (kernel_sieve.cl). The hits of a wheel modulus m in the sub-sieve s7 are the p == r (mod m), p == s7 (mod 7) that are
	// coprime to 30: they step by 7 m. For every sub-sieve, class u of m mod 30 and row c: the number dk of steps to the next hit, its row
	// cn and E = (dk * 7 W30[u] + w(c) - w(cn)) / 210, where w(c) = W30[c] + 30 tau(c) == s7 (mod 7) is the offset of row c: the next hit
	// is at j + dk * (m / 30) + E. Then for every class u and residue s mod 30 the number of steps by m to the first residue coprime to
	// 30, and its class (the kernel uses the class of 7 m).
	wheelTabs.assign(SieveSubSieves * 64 + 8 * 30, 0);
	for (unsigned int sub = 0; sub < SieveSubSieves; ++sub)
	{
		const unsigned int s7 = sub + 1;
		auto rowOffset = [s7](unsigned int row) { return static_cast<int>(W30[row] + 30 * (((s7 + 7 - W30[row] % 7) * 4) % 7)); };
		for (unsigned int u = 0; u < 8; ++u)
		{
			const unsigned int w7 = (7 * W30[u]) % 30;
			for (unsigned int row = 0; row < 8; ++row)
			{
				unsigned int dk = 1;
				while (Class30(W30[row] + dk * w7) < 0)
				{
					++dk;
				}
				const unsigned int cn = static_cast<unsigned int>(Class30(W30[row] + dk * w7));
				const int num = static_cast<int>(dk * 7 * W30[u]) + rowOffset(row) - rowOffset(cn);
				if ((num < 0) || (num % SieveWheel) || ((rowOffset(row) % 7) != static_cast<int>(s7)))
				{
					std::cerr << "Invalid wheel table" << std::endl;
					return false;
				}
				wheelTabs[sub * 64 + u * 8 + row] = dk | (cn << 8) | (static_cast<unsigned int>(num / SieveWheel) << 16);
			}
		}
	}
	for (unsigned int u = 0; u < 8; ++u)
	{
		const unsigned int w = W30[u];
		for (unsigned int s = 0; s < 30; ++s)
		{
			unsigned int t = 0;
			while (Class30(s + t * w) < 0)
			{
				++t;
			}
			wheelTabs[SieveSubSieves * 64 + u * 30 + s] = t | (static_cast<unsigned int>(Class30(s + t * w)) << 8);
		}
	}

	// tasks: row moduli are split into 8 rows (and slices of rows for small strides), wheel moduli are one task
	// per segment, in units of one iteration of the row hit loop (fitted to per-warp cycle counts on an RTX 5070, k9 work unit)
	const double C_task = 2.3, C_hit = 1.0, C_wheelTask = 4.4, C_wheelHit = 2.6;
	tasks.clear();
	tierInfo.clear();
	modTaskIdx.clear();
	modTaskList.clear();
	maxSlots = 0;
	maxSlotsMain = 1;
	for (unsigned int t = 0; t < cfg.numTiers; ++t)
	{
		struct TaskWithModulus
		{
			PlannedTask task;
			unsigned int mi;
		};
		std::vector<TaskWithModulus> rowTasks, wheelTasks;
		for (unsigned int mi = 0; mi < tierNumModuli[t]; ++mi)
		{
			const unsigned int m = modRecs[mi * 4];
			const unsigned int stride = modRecs[mi * 4 + 2];
			TaskWithModulus tm;
			tm.mi = mi;
			if (modRecs2[mi * 4 + 2] & 1)
			{
				// increment for q not dividing g (other cases are overridden on the GPU)
				const unsigned int pi = modulusPrime[mi];
				const unsigned int k = modulusPower[mi];
				const int inc = Ttab[SieveT2Size + pi * SieveTStride + k] - Ttab[SieveT2Size + pi * SieveTStride + k - 1];
				if ((inc <= 0) || (inc >= (1 << 18)) || (m / 30 >= (1U << 12)))
				{
					std::cerr << "Wheel modulus " << m << " is too small for the sieve configuration" << std::endl;
					return false;
				}
				tm.task.x = mi | (((modRecs2[mi * 4 + 2] >> 8) & 7) << 12) | 0x8000U;
				tm.task.y = (m / 30) | (static_cast<unsigned int>(inc) << 14);
				tm.task.hits = 8.0 * cfg.segJ / m;
				wheelTasks.push_back(tm);
			}
			else
			{
				// stride of the hits in a row of a sub-sieve: the stride of the row of the 30-wheel, but 7^(k-1) for q = 7 (all its hits are
				// in the sub-sieve p == r (mod 7))
				const unsigned int q = myPrimes[modulusPrime[mi]];
				const unsigned int rowStride = (q == 7) ? (stride / 7) : stride;
				unsigned int U = cfg.segJ / (rowStride * cfg.hitsPerTask);
				if (U < 1)
				{
					U = 1;
				}
				if (U > 255)
				{
					U = 255;
				}
				// the step of a task with q != 7 must be prime to 7: its hits in the row of the 30-wheel are distributed over the sub-sieves
				if ((q != 7) && ((U % 7) == 0))
				{
					--U;
				}
				if (U * rowStride >= 65536)
				{
					std::cerr << "Row modulus " << m << " is too large for the sieve configuration" << std::endl;
					return false;
				}
				// q = 3, 5: p mod q is the same in a row, so only the rows with p == r (mod q) have hits (4 rows for q = 3, 2 rows for
				// q = 5, the same rows for all powers of q): their tasks have the index k of the row among these rows (bit 24), the
				// kernel finds the row for the chunk. q = 7: bit 25.
				const unsigned int numRows = (q == 3) ? 4 : ((q == 5) ? 2 : 8);
				const unsigned int activeFlag = (q <= 5) ? (1U << 24) : ((q == 7) ? (1U << 25) : 0);
				for (unsigned int c = 0; c < numRows; ++c)
				{
					for (unsigned int u = 0; u < U; ++u)
					{
						tm.task.x = (modRecs2[mi * 4 + 1] & 0xFFU) | (c << 12) | (u << 16) | activeFlag;
						tm.task.y = U * rowStride;
						tm.task.hits = static_cast<double>(cfg.segJ) / (U * rowStride);
						rowTasks.push_back(tm);
					}
				}
			}
		}
		// A warp needs as many iterations of a hit loop as its busiest lane, the row and the wheel loop are different branches, and
		// all warps wait for the slowest one at the end of every segment. So tasks of one type with similar numbers of hits are
		// grouped into warp-slots (SieveWarpWidth consecutive work-items in one slot), and the warp-slots are distributed over the
		// warps, most expensive first, each to the warp with the smallest load so far (LPT scheduling).
		struct WarpSlot
		{
			double cost;
			size_t first;
			size_t count;
			int kind;        // 0 = row tasks, 1 = wheel tasks
			double maxHits;  // estimated hit loop iterations of the warp per segment
		};
		std::vector<TaskWithModulus> all;
		std::vector<WarpSlot> warpSlots;
		for (int kind = 0; kind < 2; ++kind)
		{
			std::vector<TaskWithModulus>& v = kind ? wheelTasks : rowTasks;
			std::stable_sort(v.begin(), v.end(), [](const TaskWithModulus& a, const TaskWithModulus& b) { return a.task.hits > b.task.hits; });
			if (kind == 0)
			{
				// Among tasks with the same number of loop iterations (ceil(hits)), interleave the moduli: some moduli have no hits in a
				// whole range (q | A1, or q = 3, 5 with p == 0 (mod q) as the root), and a warp-slot of one modulus would then be idle
				// while the others are not. Mixed warp-slots take the same time either way.
				std::vector<TaskWithModulus> mixed;
				mixed.reserve(v.size());
				for (size_t a = 0; a < v.size();)
				{
					size_t b = a;
					while ((b < v.size()) && (ceil(v[b].task.hits) == ceil(v[a].task.hits)))
					{
						++b;
					}
					std::vector<std::vector<TaskWithModulus>> byModulus;
					std::vector<unsigned int> moduli;
					for (size_t i = a; i < b; ++i)
					{
						size_t k = 0;
						while ((k < moduli.size()) && (moduli[k] != v[i].mi))
						{
							++k;
						}
						if (k == moduli.size())
						{
							moduli.push_back(v[i].mi);
							byModulus.emplace_back();
						}
						byModulus[k].push_back(v[i]);
					}
					for (size_t round = 0; mixed.size() < b; ++round)
					{
						for (const auto& list : byModulus)
						{
							if (round < list.size())
							{
								mixed.push_back(list[round]);
							}
						}
					}
					a = b;
				}
				v.swap(mixed);
			}
			for (size_t g = 0; g < v.size(); g += SieveWarpWidth)
			{
				const size_t n = std::min<size_t>(SieveWarpWidth, v.size() - g);
				// the first task of the group has the most hits
				double cost;
				double rowMaxHits;
				if (kind == 0)
				{
						cost = C_task + C_hit * ceil(v[g].task.hits);
						rowMaxHits = ceil(v[g].task.hits);
					}
				else
				{
					double maxHits;
					if (v[g].task.hits >= 1.0)
					{
						maxHits = ceil(v[g].task.hits) + 0.5;
					}
					else
					{
						// at most one hit per task: the warp iterates once if any of its tasks has a hit
						double none = 1.0;
						for (size_t k = 0; k < n; ++k)
						{
							none *= 1.0 - v[g + k].task.hits;
						}
						maxHits = 1.0 - none;
					}
						cost = C_wheelTask + C_wheelHit * maxHits;
						rowMaxHits = maxHits;
					}
					warpSlots.push_back({ cost, all.size(), n, kind, rowMaxHits });
				all.insert(all.end(), v.begin() + static_cast<std::ptrdiff_t>(g), v.begin() + static_cast<std::ptrdiff_t>(g + n));
			}
		}
		const unsigned int numWarps = cfg.wgSize / SieveWarpWidth;
		const unsigned int numSlots = std::max(1U, static_cast<unsigned int>((warpSlots.size() + numWarps - 1) / numWarps));
		std::vector<size_t> order(warpSlots.size());
		for (size_t i = 0; i < order.size(); ++i)
		{
			order[i] = i;
		}
		std::stable_sort(order.begin(), order.end(), [&warpSlots](size_t a, size_t b) { return warpSlots[a].cost > warpSlots[b].cost; });
		std::vector<double> warpLoad(numWarps, 0.0);
		std::vector<unsigned int> warpUsed(numWarps, 0);
		// development aid (AMICABLE_PRINT_WARP_FEATURES): per warp the numbers of row and wheel warp-slots and their hit loop iterations
		std::vector<double> warpFeatures(numWarps * 4, 0.0);

		const unsigned int offset = static_cast<unsigned int>(tasks.size() / 2);
		tasks.resize(tasks.size() + static_cast<size_t>(numSlots) * cfg.wgSize * 2, 0xFFFFFFFFU);
		std::vector<std::vector<unsigned int>> lists(tierNumModuli[t]);
		for (const size_t gi : order)
		{
			unsigned int w = numWarps;
			for (unsigned int k = 0; k < numWarps; ++k)
			{
				if ((warpUsed[k] < numSlots) && ((w == numWarps) || (warpLoad[k] < warpLoad[w])))
				{
					w = k;
				}
			}
			const unsigned int slot = warpUsed[w]++;
			warpLoad[w] += warpSlots[gi].cost;
			warpFeatures[w * 4 + warpSlots[gi].kind * 2] += 1.0;
			warpFeatures[w * 4 + warpSlots[gi].kind * 2 + 1] += warpSlots[gi].maxHits;
			for (size_t k = 0; k < warpSlots[gi].count; ++k)
			{
				const TaskWithModulus& tm = all[warpSlots[gi].first + k];
				const unsigned int lane = w * SieveWarpWidth + static_cast<unsigned int>(k);
				tasks[(offset + slot * cfg.wgSize + lane) * 2] = tm.task.x;
				tasks[(offset + slot * cfg.wgSize + lane) * 2 + 1] = tm.task.y;
				lists[tm.mi].push_back((slot << 16) | lane);
			}
		}
		if (getenv("AMICABLE_PRINT_WARP_FEATURES"))
		{
			std::cerr << "Warp features tier " << t << ":";
			for (unsigned int w = 0; w < numWarps; ++w)
			{
				std::cerr << " " << warpFeatures[w * 4] << "," << warpFeatures[w * 4 + 1] << "," << warpFeatures[w * 4 + 2] << "," << warpFeatures[w * 4 + 3];
			}
			std::cerr << std::endl;
		}
		if (t + 1 == cfg.numTiers)
		{
			myMaxWarpLoad = *std::max_element(warpLoad.begin(), warpLoad.end());
			myAvgWarpLoad = 0.0;
			for (const double l : warpLoad)
			{
				myAvgWarpLoad += l / numWarps;
			}
		}

		const unsigned int idxOffset = static_cast<unsigned int>(modTaskIdx.size() / 2);
		const unsigned int listOffset = static_cast<unsigned int>(modTaskList.size());
		for (unsigned int mi = 0; mi < tierNumModuli[t]; ++mi)
		{
			modTaskIdx.push_back(static_cast<unsigned int>(modTaskList.size()) - listOffset);
			modTaskIdx.push_back(static_cast<unsigned int>(lists[mi].size()));
			modTaskList.insert(modTaskList.end(), lists[mi].begin(), lists[mi].end());
		}

		tierInfo.push_back(tierNumPrimes[t]);
		tierInfo.push_back(tierNumModuli[t]);
		tierInfo.push_back(numSlots);
		tierInfo.push_back(offset);
		tierInfo.push_back(idxOffset);
		tierInfo.push_back(listOffset);
		tierInfo.push_back(tierNumRowModuli[t]);
		tierInfo.push_back(0);
		maxSlots = std::max(maxSlots, numSlots);
		if (t + 1 < cfg.numTiers)
		{
			maxSlotsMain = std::max(maxSlotsMain, numSlots);
		}
	}

	// kernel_sieve.cl: the initial task states share their memory with the sums of a segment and the two survivor lists
	// MEM4 (sums and survivor lists) holds J0TAB (9 words per row modulus) and STW (1 word per modulus) during the setup
	if (9 * std::max(1U, maxRowModuli) + maxModuli > 8 * cfg.segJ + 2 * cfg.survLocal)
	{
		std::cerr << "Too many sieve moduli for the segment size" << std::endl;
		return false;
	}
	// MEM4, INC, TBC, BASE, WHEEL, T0TAB, GF, OVR, SBN, SUB7, SUBW, UNDO + UNDON, WINS
	localMemSize = 4 * ((8 * cfg.segJ + 2 * cfg.survLocal) + std::max(1U, maxRowModuli) + 8 * 33 + 9 + SieveSubSieves * 64 + 240 + 16 + 1 + SieveMaxOverrides * 2 + 3 +
		2 + 16 + (cfg.killQMax ? SieveUndoCap : 1) + 1 + 2 * cfg.chunkSegments);

	// rho: largest sigma(x)/x for x < 2^L with all prime factors > Q
	for (unsigned int t = 0; t < cfg.numTiers; ++t)
	{
		for (unsigned int L = 0; L < 128; ++L)
		{
			double best = 1.0;
			if (L > 0)
			{
				NumberAndSum cur;
				cur.N = 1;
				cur.S = 1;
				MaxSumRatio(PrimeIterator(cfg.tierQ[t]), num64(-1), num128(1) << L, cur, best);
			}
			myLnRho[t][L] = log(best * (1.0 + 1e-12));
		}
	}

	return true;
}

void SieveGPUTables::LoExtraTable(std::vector<int>& out) const
{
	out.resize(static_cast<size_t>(cfg.numTiers) * 128);
	for (unsigned int t = 0; t < cfg.numTiers; ++t)
	{
		for (unsigned int L = 0; L < 128; ++L)
		{
			const double extra = myLnRho[t][L] + (L / 16) * myMaxTrunc;
			out[t * 128 + L] = static_cast<int>(ceil(extra * SIEVE_FX)) + 32 + SieveMargin;
		}
	}
}

unsigned int SieveGPUTables::ChooseTier(num64 numPrimes, double W) const
{
	const bool byW = (W >= 0.0) && (cfg.numTiers > 1);
	const unsigned int n = byW ? (cfg.numTiers - 1) : cfg.numTiers;
	unsigned int t = 0;
	while ((t + 1 < n) && (numPrimes >= cfg.tierMinPrimes[t + 1]))
	{
		++t;
	}
	if (byW && (W < cfg.highTierW) && (numPrimes >= cfg.tierMinPrimes[cfg.numTiers - 1]))
	{
		t = cfg.numTiers - 1;
	}
	return t;
}

namespace
{
	// Fast prime counting for EstimatePrimes (called for every range: the series in Li cost as much as the rest of the host work).
	// pi(x) for x < 10^4 exactly (table), above it Riemann's R function (2 terms) with Li(x) from the asymptotic expansion
	// x / Li(x) = L - 1 - 1/L - 3/L^2 - 13/L^3 - 71/L^4 - 461/L^5 - ...  (L = ln x; relative error ~5e-4 at 10^4, ~1e-6 at 10^9)
	const unsigned int SmallPiLimit = 10000;

	double SmallPi(unsigned int x)
	{
		static std::vector<unsigned short> table;
		if (table.empty())
		{
			std::vector<unsigned char> composite(SmallPiLimit, 0);
			table.resize(SmallPiLimit);
			unsigned short n = 0;
			for (unsigned int i = 0; i < SmallPiLimit; ++i)
			{
				if ((i >= 2) && !composite[i])
				{
					++n;
					for (unsigned int k = i * i; k < SmallPiLimit; k += i)
					{
						composite[k] = 1;
					}
				}
				table[i] = n;
			}
		}
		return table[(x < SmallPiLimit) ? x : (SmallPiLimit - 1)];
	}

	double FastLi(double x)
	{
		const double L = log(x);
		const double r = 1.0 / L;
		return x / (L - 1.0 - r * (1.0 + r * (3.0 + r * (13.0 + r * (71.0 + r * 461.0)))));
	}

	double FastPrimePi(double x)
	{
		if (x < SmallPiLimit)
		{
			return (x < 2.0) ? 0.0 : SmallPi(static_cast<unsigned int>(x));
		}
		return FastLi(x) - 0.5 * FastLi(sqrt(x));
	}
}

double SieveGPUTables::EstimatePrimes(num64 lo, num64 hi)
{
	if (hi < lo)
	{
		return 0.0;
	}
	const double a = static_cast<double>(lo);
	const double b = static_cast<double>(hi);
	if ((a >= SmallPiLimit) && (b - a < a * 0.25))
	{
		// narrow interval: integral of 1 / ln t (Simpson's rule; the difference of two estimates would lose the precision)
		return (b - a + 1.0) * (1.0 / log(a) + 4.0 / log((a + b) * 0.5) + 1.0 / log(b)) * (1.0 / 6.0);
	}
	const double n = FastPrimePi(b) - FastPrimePi(a - 1.0);
	return (n > 0.0) ? n : 0.0;
}

num64 SieveGPUTables::NumPositions(num64 p_lo, num64 p_hi) const
{
	if (p_hi < p_lo)
	{
		return 0;
	}
	const num64 P0 = (p_lo / SieveWheel) * SieveWheel;
	return ((p_hi - P0) / SieveWheel + 1) * SievePositionsPerJ;
}

num64 SieveGPUTables::NumChunks(num64 p_lo, num64 p_hi) const
{
	if (p_hi < p_lo)
	{
		return 0;
	}
	const num64 P0 = (p_lo / SieveWheel) * SieveWheel;
	const num64 numJ = (p_hi - P0) / SieveWheel + 1;
	const num64 numSegments = (numJ + cfg.segJ - 1) / cfg.segJ;
	return (numSegments + cfg.chunkSegments - 1) / cfg.chunkSegments;
}

void SieveGPUTables::ComputeReduction(num128 M0, num128 sumM0, const num64* primes, const unsigned int* exponents, unsigned int numFactors, num64 p_lo, num64 p_hi, num64 maxL, SieveReduction& r) const
{
	(void)p_lo;
	ComputeRangeReduction(M0, sumM0, primes, exponents, numFactors, p_hi, maxL, r);
}

bool SieveGPUTables::PrepareRange(const RangeData& r, num64 p_lo, num64 p_hi, unsigned int tier, bool adaptiveTier, num64* out) const
{
	num64 primes[MaxPrimeFactors];
	unsigned int exponents[MaxPrimeFactors];
	const unsigned int n = static_cast<unsigned int>(r.last_factor_index + 1);
	for (unsigned int i = 0; i < n; ++i)
	{
		primes[i] = r.factors[i].p.Get();
		exponents[i] = r.factors[i].k;
	}
	return PrepareRange(r.value, r.sum, primes, exponents, n, p_lo, p_hi, tier, adaptiveTier, out);
}

bool SieveGPUTables::PrepareRange(num128 M0, num128 sumM0, const num64* primes, const unsigned int* exponents, unsigned int numFactors, num64 p_lo, num64 p_hi, unsigned int tier, bool adaptiveTier, num64* out) const
{
	if ((p_hi >> (64 - SieveRangeIndexBits)) || (p_hi < p_lo) || (tier >= cfg.numTiers))
	{
		return false;
	}

	const num128 A = sumM0 - M0;
	const num128 B = sumM0;

	// g = gcd(A, B) = gcd(M0, sigma(M0)): the exponent of each prime of M0 in g is min(its exponent in M0, its exponent in B)
	num128 A1 = A;
	num128 B1 = B;
	unsigned int v2g = 0;
	unsigned int ngf = 0;
	unsigned int gf[SieveMaxGFactors] = {};
	for (unsigned int i = 0; i < numFactors; ++i)
	{
		const num64 q = primes[i];
		const unsigned int k = exponents[i];
		if (q == 2)
		{
			unsigned int v = CountTrailingZeros128(B1);
			if (v > k)
			{
				v = k;
			}
			v2g = v;
			A1 >>= v;
			B1 >>= v;
			continue;
		}

		// exact division test: x * q^-1 mod 2^128 <= x if and only if q divides x (because x < 2^128 / q)
		const num128 inv = -modular_inverse128(q);
		unsigned int e = 0;
		while (e < k)
		{
			const num128 y = B1 * inv;
			if (y > B1)
			{
				break;
			}
			B1 = y;
			A1 = A1 * inv;
			++e;
		}
		if ((e > 0) && (q <= cfg.tierQ[cfg.numTiers - 1]))
		{
			if (ngf >= SieveMaxGFactors)
			{
				return false;
			}
			gf[ngf++] = static_cast<unsigned int>(q) | (e << 24);
		}
	}

	const bool twoAdic = (LowWord(A1) & 1) && (LowWord(B1) & 1);
	const num64 P0 = (p_lo / SieveWheel) * SieveWheel;

	// window: ln(T/N) = ln(B1/A1) - ln(1 + (B1/A1 - 1) / (p + 1))
	const long double lnr = logl(static_cast<long double>(Num128ToDouble(B1)) / static_cast<long double>(Num128ToDouble(A1)));
	const int lnr_fx = static_cast<int>(llroundl(lnr * static_cast<long double>(SIEVE_FX)));
	const float rm1 = static_cast<float>(Num128ToDouble(M0) / Num128ToDouble(A));

	// The part of N without prime factors <= Q (rho) and the prime powers that are not sieved (slack) are divisors of N / F, where F is
	// the part with prime factors <= Q. F is a multiple of the Q-smooth part of g (and of 2^(v2g + 1) if N is always even), so the
	// bit length of N / F is at most L = bitlength(Nmax) - floor(log2(that part)).
	const num128 Nmax = A * p_hi + B;
	unsigned int L = BitLength128(Nmax);
	{
		double log2F = v2g + (twoAdic ? 1 : 0);
		for (unsigned int f = 0; f < ngf; ++f)
		{
			const unsigned int q = gf[f] & 0xFFFFFFU;
			if (q <= cfg.tierQ[tier])
			{
				log2F += (gf[f] >> 24) * log2(static_cast<double>(q));
			}
		}
		const unsigned int known = static_cast<unsigned int>(log2F);
		L = (L > known) ? (L - known) : 0;
	}
	if (L > 127)
	{
		L = 127;
	}

	// rho for the part of N without prime factors <= Q, slack for prime powers that are not sieved (each of them is > 2^16, so
	// at most L / 16 of them divide N), rounding of the table entries (0.5 per prime factor), margin
	int loExtra = 0;
	for (;;)
	{
		const double extra = myLnRho[tier][L] + (L / 16) * myMaxTrunc;
		loExtra = static_cast<int>(ceil(extra * SIEVE_FX)) + 32 + SieveMargin;
		if (!adaptiveTier || (tier + 1 >= cfg.numTiers))
		{
			break;
		}

		// smallest sum (all positions where A1 * p + B1 has no odd prime factors <= Q): g and the 2-adic part
		int base = Ttab[v2g + (twoAdic ? 1 : 0)];
		const unsigned int numPrimes = tierInfo[tier * 8];
		for (unsigned int f = 0; f < ngf; ++f)
		{
			const unsigned int q = gf[f] & 0xFFFFFFU;
			for (unsigned int i = 0; i < numPrimes; ++i)
			{
				if (myPrimes[i] == q)
				{
					base += Ttab[SieveT2Size + i * SieveTStride + (gf[f] >> 24)];
					break;
				}
			}
		}
		const double corr = log1p((Num128ToDouble(M0) / Num128ToDouble(A)) / (static_cast<double>(p_lo) + 1.0)) * SIEVE_FX;
		const int lo = lnr_fx - static_cast<int>(corr) - 1 - loExtra;
		if (base < lo)
		{
			break;
		}
		++tier;
	}

	unsigned int floatBits;
	memcpy(&floatBits, &rm1, sizeof(floatBits));

	out[0] = LowWord(A1);
	out[1] = HighWord(A1);
	out[2] = LowWord(B1);
	out[3] = HighWord(B1);
	out[4] = LowWord(M0);
	out[5] = HighWord(M0);
	out[6] = LowWord(sumM0);
	out[7] = HighWord(sumM0);
	out[8] = P0;
	out[9] = p_lo;
	out[10] = p_hi;
	out[11] = static_cast<num64>(static_cast<unsigned int>(lnr_fx)) | (static_cast<num64>(static_cast<unsigned int>(loExtra)) << 32);
	out[12] = static_cast<num64>(floatBits) | (static_cast<num64>(tier) << 32) | (static_cast<num64>(twoAdic ? 1 : 0) << 40) | (static_cast<num64>(v2g) << 48) | (static_cast<num64>(ngf) << 56);
	for (unsigned int i = 0; i < SieveMaxGFactors / 2; ++i)
	{
		out[13 + i] = static_cast<num64>(gf[i * 2]) | (static_cast<num64>(gf[i * 2 + 1]) << 32);
	}
	out[21] = 0;
	out[22] = 0;
	out[23] = 0;
	return true;
}

int SieveGPUTables::DirectSum(num128 N, unsigned int tier) const
{
	// exact value of the sum (natural log of sigma(F)/F for the part F of N with prime factors <= Q)
	int s = 0;
	const unsigned int v2 = CountTrailingZeros128(N);
	s += (v2 < SieveT2Size) ? Ttab[v2] : FixedPoint(LogSigmaRatio(2.0L, static_cast<int>(v2)));
	N >>= v2;
	const unsigned int numPrimes = tierInfo[tier * 8];
	for (unsigned int i = 0; i < numPrimes; ++i)
	{
		const num64 q = myPrimes[i];
		int e = 0;
		while (Mod128(N, q) == 0)
		{
			N /= q;
			++e;
		}
		s += (e < SieveTStride) ? Ttab[SieveT2Size + i * SieveTStride + e] : FixedPoint(LogSigmaRatio(static_cast<long double>(q), e));
	}
	return s;
}

int SieveGPUTables::ExpectedSum(const num64* rec, num64 p) const
{
	const num128 A1 = CombineNum128(rec[0], rec[1]);
	const num128 B1 = CombineNum128(rec[2], rec[3]);
	const unsigned int tier = static_cast<unsigned int>(rec[12] >> 32) & 0xFFU;
	const bool twoAdic = ((rec[12] >> 40) & 1) != 0;
	const unsigned int v2g = static_cast<unsigned int>(rec[12] >> 48) & 0xFFU;
	const unsigned int ngf = static_cast<unsigned int>(rec[12] >> 56);

	num128 x = A1 * p + B1;
	unsigned int v = 0;
	if (twoAdic)
	{
		v = CountTrailingZeros128(x);
		if (v > 32)
		{
			v = 32;
		}
	}
	int s = Ttab[v2g + v];

	const unsigned int numPrimes = tierInfo[tier * 8];
	for (unsigned int i = 0; i < numPrimes; ++i)
	{
		const num64 q = myPrimes[i];
		unsigned int eg = 0;
		for (unsigned int f = 0; f < ngf; ++f)
		{
			const unsigned int gf = static_cast<unsigned int>(rec[13 + f / 2] >> ((f & 1) * 32));
			if ((gf & 0xFFFFFFU) == q)
			{
				eg = gf >> 24;
			}
		}
		unsigned int e = 0;
		while ((e < myMaxPower[i]) && (Mod128(x, q) == 0))
		{
			x /= q;
			++e;
		}
		s += Ttab[SieveT2Size + i * SieveTStride + eg + e];

		// divisibility kill of the wheel primes (kernel_sieve.cl): exponent exactly 1 in A1 * p + B1 (q does not divide g), D(q) >= 7
		// not dividing p + 1
		const unsigned int mi = primeRecs[i * 4 + 2] & 0xFFFFU;
		if ((e == 1) && (eg == 0) && (q <= killQMax) && (tier + 1 == cfg.numTiers) && (modRecs2[mi * 4 + 2] & 1) && ((myMaxPower[i] > 1) || (Mod128(x, q) != 0)))
		{
			const num128 sumM0 = CombineNum128(rec[6], rec[7]);
			const num64 m = q + 1;
			const num64 D = m / GCD(m, Mod128(sumM0, m));
			// (kernel_sieve.cl, SieveKillD: if 7 | D, the sub-sieve p == 6 (mod 7) kills only if D >= 49)
			if ((D >= 7) && ((p + 1) % D) && (((D % 7) != 0) || ((p % 7) != 6) || (D >= 49)))
			{
				s += SieveKill;
			}
		}
	}
	return s;
}

bool SieveGPUTables::DivisibilityOK2(num128 M0, num128 sumM0, num128 g, num64 p, unsigned int maxQ, int mode) const
{
	num128 N = (sumM0 - M0) * p + sumM0;
	auto check = [&](num64 q, unsigned int v)
	{
		const num128 sq = RangeSigmaPrimePower(q, v);
		const num128 D = sq / RangeGcd128(sq, sumM0);
		return (D == 1) || ((num128(p + 1) % D) == 0);
	};
	const unsigned int v2 = CountTrailingZeros128(N);
	N >>= v2;
	if ((mode >= 2) && (v2 > 0) && !check(2, v2))
	{
		return false;
	}
	for (unsigned int i = 0; (i < myPrimes.size()) && (myPrimes[i] <= maxQ); ++i)
	{
		const num64 q = myPrimes[i];
		if (((g % num128(q)) == 0) || ((mode >= 3) && (q <= 255)))
		{
			continue;
		}
		unsigned int v = 0;
		while ((N % num128(q)) == 0)
		{
			N /= q;
			++v;
		}
		if ((v == 1) && !check(q, 1))
		{
			return false;
		}
	}
	return true;
}

bool SieveGPUTables::DivisibilityOK(num128 M0, num128 sumM0, num64 p, unsigned int Q) const
{
	num128 N = (sumM0 - M0) * p + sumM0;
	auto check = [&](num64 q, unsigned int v)
	{
		const num128 sq = RangeSigmaPrimePower(q, v);
		const num128 D = sq / RangeGcd128(sq, sumM0);
		return (D == 1) || ((num128(p + 1) % D) == 0);
	};
	const unsigned int v2 = CountTrailingZeros128(N);
	N >>= v2;
	if ((v2 > 0) && !check(2, v2))
	{
		return false;
	}
	for (unsigned int i = 0; (i < myPrimes.size()) && (myPrimes[i] <= Q); ++i)
	{
		const num64 q = myPrimes[i];
		unsigned int v = 0;
		while ((N % num128(q)) == 0)
		{
			N /= q;
			++v;
		}
		if ((v > 0) && !check(q, v))
		{
			return false;
		}
	}
	return true;
}

bool SieveGPUTables::InsideWindow(const num64* rec, num64 p, int S)
{
	const int lnr_fx = static_cast<int>(static_cast<unsigned int>(rec[11]));
	const int loExtra = static_cast<int>(static_cast<unsigned int>(rec[11] >> 32));
	float rm1;
	const unsigned int floatBits = static_cast<unsigned int>(rec[12]);
	memcpy(&rm1, &floatBits, sizeof(rm1));
	const double x = static_cast<double>(rm1) / (static_cast<double>(p) + 1.0);
	const double corr = log1p(x) * SIEVE_FX;
	const int lo = lnr_fx - static_cast<int>(corr) - 1 - loExtra;
	const int hi = lnr_fx - static_cast<int>(corr) + SieveMargin;
	return (lo <= S) && (S <= hi);
}

bool ParseKnownPair(const std::string& line2, const std::string& line3, KnownPair& pair, bool& bad, bool largePrimes)
{
	bad = false;
	const num128 m = atoi128(line2.c_str());
	const num128 n = atoi128(line3.c_str());
	num64 primes[64];
	unsigned int exps[64];
	unsigned int nf = 0;
	const char* s = strchr(line2.c_str(), '=');
	bool ok = (s != nullptr);
	if (ok)
	{
		++s;
		while ((*s >= '0') && (*s <= '9') && (nf < 64))
		{
			num64 q = 0;
			while ((*s >= '0') && (*s <= '9'))
			{
				q = q * 10 + static_cast<num64>(*s - '0');
				++s;
			}
			unsigned int e = 1;
			if (*s == '^')
			{
				++s;
				e = 0;
				while ((*s >= '0') && (*s <= '9'))
				{
					e = e * 10 + static_cast<unsigned int>(*s - '0');
					++s;
				}
			}
			primes[nf] = q;
			exps[nf] = e;
			++nf;
			if (*s == '*')
			{
				++s;
			}
			else
			{
				break;
			}
		}
	}
	// factorizations in the files are not always sorted
	for (unsigned int i = 1; i < nf; ++i)
	{
		for (unsigned int j = i; (j > 0) && (primes[j - 1] > primes[j]); --j)
		{
			std::swap(primes[j - 1], primes[j]);
			std::swap(exps[j - 1], exps[j]);
		}
	}
	if (!ok || (nf < 2) || (nf > 32) || (exps[nf - 1] != 1) || ((m % 6) == 0) ||
		(largePrimes ? (primes[nf - 1] <= SearchLimit::LinearLimit) : (primes[nf - 1] > SearchLimit::MainPrimeTableBound)))
	{
		return false;
	}

	pair.m = m;
	pair.n = n;
	pair.p = primes[nf - 1];
	pair.numFactors = nf - 1;
	pair.M0 = 1;
	pair.sumM0 = 1;
	for (unsigned int i = 0; i + 1 < nf; ++i)
	{
		pair.primes[i] = primes[i];
		pair.exponents[i] = exps[i];
		num128 pk = 1;
		num128 sp = 1;
		for (unsigned int k = 0; k < exps[i]; ++k)
		{
			pk *= primes[i];
			sp += pk;
		}
		pair.M0 *= pk;
		pair.sumM0 *= sp;
	}
	if ((pair.M0 * pair.p != m) || (pair.sumM0 * (pair.p + 1) != m + n))
	{
		bad = true;
		return false;
	}
	return true;
}

bool LoadKnownPairs(unsigned int stride, std::vector<KnownPair>& pairs, std::vector<num128>* allAmicable)
{
	const char* files[] = { "c2_3.txt", "c2_4.txt", "c2_5.txt", "c2_6.txt", "c2_7.txt", "c2_8.txt", "c2_9.txt", "c2_10.txt", "c2_11.txt", "c2_12.txt", "c2_13.txt", "c2_14.txt", "c2_15.txt", "c2_16.txt", "c2_17.txt", "c2_18.txt", "c2_19.txt", "c2_20.txt", "c2_20_1.txt", "c2_20_2.txt", "c2_21.txt", "c2_21_1.txt", "c2_21_2.txt", "c2_21_3.txt", "c2_21_4.txt", "c2_22.txt" };

	pairs.clear();
	if (allAmicable)
	{
		allAmicable->clear();
	}
	if (stride < 1)
	{
		stride = 1;
	}

	num64 numFirstStage = 0;
	for (const char* name : files)
	{
		std::ifstream f(name);
		std::string line1, line2, line3, line4;
		while (std::getline(f, line1) && std::getline(f, line2) && std::getline(f, line3))
		{
			std::getline(f, line4);

			const num128 m = atoi128(line2.c_str());
			if (m >= SearchLimit::value)
			{
				break;
			}
			if (allAmicable)
			{
				// both members: the search reports every amicable number m in its ranges, also the larger member of a pair
				allAmicable->push_back(m);
				allAmicable->push_back(atoi128(line3.c_str()));
			}

			KnownPair pair;
			bool bad = false;
			if (!ParseKnownPair(line2, line3, pair, bad))
			{
				if (bad)
				{
					std::cerr << "LoadKnownPairs: bad input for " << m << std::endl;
					return false;
				}
				continue;
			}
			if ((numFirstStage++ % stride) != 0)
			{
				continue;
			}
			pairs.push_back(pair);
		}
	}

	if (allAmicable)
	{
		std::sort(allAmicable->begin(), allAmicable->end());
	}
	std::cerr << "Loaded " << pairs.size() << " of " << numFirstStage << " first-stage pairs" << std::endl;
	return !pairs.empty();
}

bool TestSieveReduction(const SieveGPUTables& t, const std::vector<KnownPair>& pairs)
{
	// Every known pair must satisfy the necessary conditions of its range (with the range ending at its p, where T/N is largest)
	num64 numWithL = 0, numWithL235 = 0, numFixedL = 0;
	std::vector<num64> Lvalues;
	for (const KnownPair& pair : pairs)
	{
		SieveReduction red;
		t.ComputeReduction(pair.M0, pair.sumM0, pair.primes, pair.exponents, pair.numFactors, pair.p, pair.p, num64(1) << 62, red);
		bool ok = !red.impossible;
		if (ok && (red.Lfull > 1) && (((pair.p + 1) % red.Lfull) != 0))
		{
			ok = false;
		}
		if (ok && (red.L > 1) && (((pair.p + 1) % red.L) != 0))
		{
			ok = false;
		}
		if (ok && (red.L235 > 1) && (((pair.p + 1) % red.L235) != 0))
		{
			ok = false;
		}
		// the exponents of the primes of L in N must be the stated ones
		for (unsigned int i = 0; ok && (i < red.numLPrimes); ++i)
		{
			num128 x = pair.n;
			unsigned int e = 0;
			while ((x % num128(static_cast<num64>(red.lq[i]))) == 0)
			{
				x /= static_cast<num64>(red.lq[i]);
				++e;
			}
			if (e != red.le[i])
			{
				ok = false;
			}
		}
		if (!ok)
		{
			std::cerr << "TestSieveReduction: the conditions reject the known pair " << pair.m << " (p = " << pair.p << ", impossible = " << red.impossible << ", L = " << red.Lfull << ", L235 = " << red.L235 << ")" << std::endl;
			return false;
		}
		numWithL += (red.Lfull != 1) ? 1 : 0;
		numWithL235 += (red.L235 != 1) ? 1 : 0;
		numFixedL += red.numLPrimes ? 1 : 0;
	}
	std::cerr << "TestSieveReduction: " << pairs.size() << " pairs pass (" << numWithL << " with L > 1, " << numWithL235 << " with 2, 3, 5 conditions, " << numFixedL << " with primes of fixed exponent in L)" << std::endl;
	return !pairs.empty();
}

bool TestSieveGPUWindow(const SieveGPUTables& t, const std::vector<KnownPair>& pairs)
{
	num64 numChecks = 0;
	for (const KnownPair& pair : pairs)
	{
		for (unsigned int tier = 0; tier < t.cfg.numTiers; ++tier)
		{
			num64 rec[SieveRangeStride];
			if (!t.PrepareRange(pair.M0, pair.sumM0, pair.primes, pair.exponents, pair.numFactors, pair.p, pair.p, tier, false, rec))
			{
				std::cerr << "TestSieveGPUWindow: PrepareRange failed for " << pair.m << std::endl;
				return false;
			}
			++numChecks;
			// both the exact sum and the sum computed by the GPU (prime powers above the sieve limit are missing) must be inside
			if (!SieveGPUTables::InsideWindow(rec, pair.p, t.DirectSum(pair.n, tier)) || !SieveGPUTables::InsideWindow(rec, pair.p, t.ExpectedSum(rec, pair.p)))
			{
				std::cerr << "TestSieveGPUWindow: window rejects " << pair.m << " (tier " << tier << ")" << std::endl;
				return false;
			}
		}
	}

	std::cerr << "TestSieveGPUWindow: " << pairs.size() << " pairs, " << numChecks << " window checks passed" << std::endl;
	return !pairs.empty();
}
