#include "stdafx.h"
#include "PrimeTables.h"
#include "RangeGen.h"
#include "RangeReduction.h"
#include "sprp64.h"
#include <cmath>

// See RangeReduction.h. Shared by the GPU app (skips the ranges, sieves the progression) and the task generator (skips the
// ranges when it makes work units), so both must stay identical.

namespace
{
	long double LogSigmaRatio(long double q, int e)
	{
		if (e <= 0)
		{
			return 0.0L;
		}
		// sigma(q^e) / q^e = (q - q^-e) / (q - 1)
		return logl((q - powl(q, static_cast<long double>(-e))) / (q - 1.0L));
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
}

num128 RangeGcd128(num128 a, num128 b)
	{
		while (b != 0)
		{
			const num128 t = a % b;
			a = b;
			b = t;
		}
		return a;
	}

// 1 + q + ... + q^e (q^e divides M0, so this fits in 128 bits)
num128 RangeSigmaPrimePower(num64 q, unsigned int e)
	{
		num128 s = 1;
		num128 pk = 1;
		for (unsigned int k = 0; k < e; ++k)
		{
			pk *= q;
			s += pk;
		}
		return s;
	}

namespace
{
	num128 Gcd128(num128 a, num128 b)
	{
		return RangeGcd128(a, b);
	}

	num128 SigmaPrimePower(num64 q, unsigned int e)
	{
		return RangeSigmaPrimePower(q, e);
	}

	// sigma(q^e) / gcd(sigma(q^e), sigma(M0)), split into the part coprime to 30 and the 2-, 3-, 5-parts
	void DivisorCondition(const num128& sigmaQE, const num128& sumM0, num128& part30, num128& part235)
	{
		num128 d = sigmaQE / Gcd128(sigmaQE, sumM0);
		part235 = 1;
		const num64 wheelPrimes[3] = { 2, 3, 5 };
		for (num64 s : wheelPrimes)
		{
			while ((d % num128(s)) == 0)
			{
				d /= s;
				part235 *= s;
			}
		}
		part30 = d;
	}

	num128 Lcm128(const num128& a, const num128& b, bool& overflow)
	{
		if ((a == 0) || (b == 0))
		{
			overflow = true;
			return 0;
		}
		const num128 q = a / Gcd128(a, b);
		// q * b must fit in 64 bits for the kernel; larger values are reported as overflow
		if ((HighWord(q) != 0) || (HighWord(b) != 0) || ((LowWord(q) != 0) && (LowWord(b) > num64(-1) / LowWord(q))))
		{
			overflow = true;
			return 0;
		}
		return q * b;
	}
}

namespace
{
	// Prime powers q^a || L for primes q <= 65535 (L is coprime to 30); the rest (a product of larger primes) is returned in cofactor
	void FactorL(num64 L, num64* q, num64* qa, unsigned int* a, unsigned int& n, num64& cofactor)
	{
		n = 0;
		for (num64 d = 7; (d <= 65535) && (d * d <= L); d += 2)
		{
			if (L % d)
			{
				continue;
			}
			q[n] = d;
			qa[n] = 1;
			a[n] = 0;
			while ((L % d) == 0)
			{
				L /= d;
				qa[n] *= d;
				++a[n];
			}
			++n;
		}
		cofactor = 1;
		if (L > 1)
		{
			if (L <= 65535)
			{
				q[n] = L;
				qa[n] = L;
				a[n] = 1;
				++n;
			}
			else if (L < (num64(65537) * 65537))
			{
				// no prime factor below its square root: prime, but too large for the kernel
				cofactor = L;
			}
			else
			{
				cofactor = L;
			}
		}
	}
}

void ComputeRangeReduction(num128 M0, num128 sumM0, const num64* primes, const unsigned int* exponents, unsigned int numFactors, num64 p_hi, num64 maxL, RangeReduction& r)
{
	r.impossible = false;
	r.L = 1;
	r.Lfull = 1;
	r.L235 = 1;
	r.numForced = 0;
	r.numLPrimes = 0;
	r.W = 0;
	r.v2 = 0;
	r.fixed2 = true;

	const num128 A = sumM0 - M0;
	const num128 B = sumM0;

	// Known part K of N: the prime powers of g (all primes of M0, not only the sieved ones), later also primes with a fixed
	// exponent in the progression. ke = exponent of q in N (a lower bound, exact if kfixed)
	num128 A1 = A;
	num128 B1 = B;
	unsigned int v2g = 0;
	enum { MaxK = MaxPrimeFactors + 32 };
	num64 kq[MaxK];
	unsigned int ke[MaxK];
	unsigned int keg[MaxK]; // exponent of q in g
	bool kfixed[MaxK];
	unsigned int nk = 0;
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
		if (e > 0)
		{
			kq[nk] = q;
			ke[nk] = e;
			keg[nk] = e;
			kfixed[nk] = false;
			++nk;
		}
	}
	for (unsigned int i = 0; i < nk; ++i)
	{
		// q | A1 (then q does not divide B1): q never divides N1 = A1 * p + B1
		kfixed[i] = (A1 % num128(kq[i])) == 0;
	}
	const bool twoAdic = (LowWord(A1) & 1) && (LowWord(B1) & 1);
	unsigned int v2 = v2g + (twoAdic ? 1 : 0); // power of 2 in N (lower bound)
	bool fixed2 = !twoAdic;                    // N1 always odd: the power of 2 in N is exactly 2^v2g

	// largest T / N of the range (T / N = B (p + 1) / (A p + B) grows with p)
	const long double logHi = logl(static_cast<long double>(Num128ToDouble(B * (p_hi + 1)))) - logl(static_cast<long double>(Num128ToDouble(A * p_hi + B)));
	const long double eps = 1e-9L;
	const num128 c0 = B1 - A1; // N1 == B1 - A1 (mod q^a) for all p == -1 (mod q^a)

	num128 L = 1;
	bool overflow = false;
	for (int iteration = 0; iteration < 16; ++iteration)
	{
		long double logK = LogSigmaRatio(2.0L, static_cast<int>(v2));
		for (unsigned int i = 0; i < nk; ++i)
		{
			logK += LogSigmaRatio(static_cast<long double>(kq[i]), static_cast<int>(ke[i]));
		}
		const long double W = logHi - logK;
		r.W = static_cast<double>(W);
		if (W < -eps)
		{
			r.impossible = true;
			return;
		}

		// exponents that can not grow without exceeding the largest T / N
		for (unsigned int i = 0; i < nk; ++i)
		{
			if (!kfixed[i] && (LogSigmaRatio(static_cast<long double>(kq[i]), static_cast<int>(ke[i] + 1)) - LogSigmaRatio(static_cast<long double>(kq[i]), static_cast<int>(ke[i])) > W + eps))
			{
				kfixed[i] = true;
			}
		}
		if (!fixed2 && (LogSigmaRatio(2.0L, static_cast<int>(v2 + 1)) - LogSigmaRatio(2.0L, static_cast<int>(v2)) > W + eps))
		{
			fixed2 = true;
		}
		r.v2 = v2;
		r.fixed2 = fixed2;

		// every fixed q^e in N: sigma(q^e) | T = sigma(M0) * (p + 1)
		num128 newL = 1;
		num128 newL235 = 1;
		bool newOverflow = false;
		unsigned int numForced = 0;
		auto addCondition = [&](const num128& sigmaQE)
		{
			num128 part30, part235;
			DivisorCondition(sigmaQE, sumM0, part30, part235);
			bool o235 = false;
			const num128 l235 = Lcm128(newL235, part235, o235);
			if (!o235)
			{
				newL235 = l235;
			}
			if (!newOverflow)
			{
				newL = Lcm128(newL, part30, newOverflow);
			}
			++numForced;
		};
		for (unsigned int i = 0; i < nk; ++i)
		{
			if (kfixed[i] && (ke[i] > 0))
			{
				addCondition(SigmaPrimePower(kq[i], ke[i]));
			}
		}
		if (fixed2 && (v2 > 0))
		{
			addCondition(SigmaPrimePower(2, v2));
		}
		r.numForced = numForced;
		r.L235 = LowWord(newL235);
		L = newL;
		overflow = newOverflow;
		if (overflow)
		{
			break;
		}

		// In the progression p == -1 (mod q^a), q^a || L, N1 == B1 - A1 (mod q^a) is the same for all p: if its valuation v is
		// below a, the exponent of q in N is fixed (q^a with v >= a is dropped from L, the exponent of q is not fixed then)
		bool changed = false;
		num64 Lkept = 1;
		num64 pq[64], pqa[64];
		unsigned int pa[64];
		unsigned int np = 0;
		num64 cofactor = 1;
		FactorL(LowWord(L), pq, pqa, pa, np, cofactor);
		if (cofactor > 1)
		{
			// primes above 65535 (not factored): keep them if they don't divide N1 in the progression
			const num64 cm = LowWord(c0 % num128(cofactor));
			if (GCD(cm, cofactor) == 1)
			{
				Lkept *= cofactor;
			}
		}
		for (unsigned int f = 0; f < np; ++f)
		{
			const num64 q = pq[f];
			const num64 qa = pqa[f];
			const unsigned int a = pa[f];
			const num64 cm = LowWord(c0 % num128(qa));
			unsigned int v = 0;
			if (cm == 0)
			{
				v = a;
			}
			else
			{
				for (num64 x = cm; (x % q) == 0; x /= q)
				{
					++v;
				}
			}
			if (v >= a)
			{
				continue;
			}
			Lkept *= qa;
			unsigned int idx = nk;
			for (unsigned int i = 0; i < nk; ++i)
			{
				if (kq[i] == q)
				{
					idx = i;
				}
			}
			if (idx == nk)
			{
				if ((v == 0) || (nk >= MaxK))
				{
					continue; // q does not divide N
				}
				kq[nk] = q;
				ke[nk] = v;
				keg[nk] = 0;
				kfixed[nk] = true;
				++nk;
				changed = true;
			}
			else
			{
				const unsigned int total = keg[idx] + v;
				if (kfixed[idx] && (ke[idx] != total))
				{
					// the exponent was forced to a different value: no p satisfies both
					r.impossible = true;
					return;
				}
				if (!kfixed[idx] || (ke[idx] != total))
				{
					ke[idx] = total;
					kfixed[idx] = true;
					changed = true;
				}
			}
		}
		L = Lkept;
		if (!changed)
		{
			break;
		}
	}

	r.Lfull = overflow ? 0 : LowWord(L);
	if (overflow || (L == 1))
	{
		return;
	}

	// The kernel uses a divisor of L up to maxL (whole prime powers, so that the exponents of their primes in N stay fixed)
	num64 pq[64], pqa[64];
	unsigned int pa[64];
	unsigned int np = 0;
	num64 cofactor = 1;
	FactorL(LowWord(L), pq, pqa, pa, np, cofactor);
	num64 Luse = 1;
	for (unsigned int f = 0; f < np; ++f)
	{
		const num64 q = pq[f];
		const num64 qa = pqa[f];
		if ((qa <= maxL) && (Luse <= maxL / qa))
		{
			Luse *= qa;
			// the exponent of q in N (fixed in the progression; not in K means q does not divide N)
			if (r.numLPrimes < RangeReductionMaxLPrimes)
			{
				unsigned int e = 0;
				for (unsigned int i = 0; i < nk; ++i)
				{
					if (kq[i] == q)
					{
						e = ke[i];
					}
				}
				r.lq[r.numLPrimes] = static_cast<unsigned int>(q);
				r.le[r.numLPrimes] = e;
				++r.numLPrimes;
			}
		}
	}
	r.L = Luse;
}


bool IsRangeProvablyEmpty(const RangeData& r, num64 p_hi)
{
	num64 primes[MaxPrimeFactors];
	unsigned int exponents[MaxPrimeFactors];
	const unsigned int n = static_cast<unsigned int>(r.last_factor_index + 1);
	for (unsigned int i = 0; i < n; ++i)
	{
		primes[i] = r.factors[i].p.Get();
		exponents[i] = r.factors[i].k;
	}
	RangeReduction red;
	ComputeRangeReduction(r.value, r.sum, primes, exponents, n, p_hi, num64(1) << 62, red);
	return red.impossible;
}
