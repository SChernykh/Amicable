#pragma once

struct RangeData;

// Necessary conditions on p for every amicable number m = M0 * p of a range of the first search stage (p prime, larger than every
// prime factor of M0). With N = sigma(M0) * (p + 1) - M0 * p = A * p + B = g * (A1 * p + B1), g = gcd(A, B), and
// T = sigma(M0) * (p + 1), m is amicable if and only if sigma(N) = T.
// - sigma(N) / N >= sigma(K) / K for every divisor K of N, and T / N < B / A grows with p. K starts as g (times 2 if N1 is always
//   even): if its abundancy is above the largest T / N of the range, no p of the range gives an amicable number.
// - For a prime power q^e || K: if q | N1 would push the abundancy of K * q above the largest T / N (or q | A1, then q never divides
//   N1), the exponent of q in N is exactly e, so sigma(q^e) divides sigma(N) = T: D = sigma(q^e) / gcd(sigma(q^e), sigma(M0)) divides
//   p + 1. The same holds for the power of 2 in N when it is fixed (N1 always odd) or forced.
// - L = lcm of the parts of these D coprime to 30: every amicable p satisfies p == -1 (mod L). For p in this progression,
//   N1 == B1 - A1 (mod q^a) for q^a || L, so the exponent of q in N is known when the valuation of B1 - A1 is below a: K grows, and
//   the steps are repeated.
// Validated on all known first-stage pairs below 10^22 ("/test_reduction 1" in the GPU app).
enum
{
	RangeReductionMaxLPrimes = 8,
};

struct RangeReduction
{
	bool impossible;
	num64 L;          // progression modulus for the sieve (a divisor of Lfull up to maxL), 1 = no condition
	num64 Lfull;      // every amicable p satisfies p == -1 (mod Lfull); 0 = it does not fit in 64 bits
	num64 L235;       // lcm of the 2-, 3- and 5-parts of the conditions (not used yet)
	unsigned int numForced;
	unsigned int numLPrimes;                       // primes q | L: their exponent in N is the same for all p == -1 (mod L)
	unsigned int lq[RangeReductionMaxLPrimes];
	unsigned int le[RangeReductionMaxLPrimes];     // exponent of lq[i] in N (including its exponent in g)
	double W;         // log of the largest T / N minus the log of the abundancy of the known part K
	unsigned int v2;  // the power of 2 in N that W assumes (a lower bound)
	bool fixed2;      // true: no higher power of 2 fits in the window W (else the numbers with 2^(v2 + k) in N have a window
	                  // W - ln(sigma(2^(v2 + k)) / 2^(v2 + k)) + ln(sigma(2^v2) / 2^v2): OpenCL::PrepareSieveRange, band)
};

// M0 = product of primes[i]^exponents[i] (any order), sumM0 = sigma(M0), p <= p_hi
void ComputeRangeReduction(num128 M0, num128 sumM0, const num64* primes, const unsigned int* exponents, unsigned int numFactors, num64 p_hi, num64 maxL, RangeReduction& r);

// The range m = r.value * p, r.start_prime <= p <= p_hi has no amicable numbers (first stage only: the largest prime power is 1)
bool IsRangeProvablyEmpty(const RangeData& r, num64 p_hi);

num128 RangeGcd128(num128 a, num128 b);
num128 RangeSigmaPrimePower(num64 q, unsigned int e);
