#pragma once

#include <atomic>

struct RangeData;

// Arithmetic-progression sieve for the first search stage (SearchRange: m = M0 * p, p prime, p > largest prime factor of M0).
//
// For a fixed M0, N(p) = sigma(M0) * (p + 1) - M0 * p = A * p + B with A = sigma(M0) - M0 and B = sigma(M0) is linear in p,
// so every prime power q^k divides N(p) for exactly one residue class of p modulo q^k (or for none). Instead of trial-dividing
// every N(p), the sieve adds log(sigma(q^e) / q^e) for all prime powers q^e || N(p) with q <= Q, for a whole segment of p at once.
// sigma(N) = T = sigma(M0) * (p + 1) is only possible if the accumulated value S(p) satisfies
//
//     log(T / N) - log(rho_Q) <= S(p) <= log(T / N)
//
// where rho_Q is the largest possible sigma(x) / x for the part of N without prime factors <= Q. Only primes p that pass this
// test (about 0.1% of them) are checked exactly with CheckPair128.

extern bool g_UseRangeSieve;

// Testing/debugging knobs: force a sieve depth tier (0..2, -1 = automatic); with RANGE_SIEVE_VERIFY defined, every sieved sum
// is compared with trial division and mismatches are counted
extern int g_RangeSieveTierOverride;
// the layout of the range sieve: -1 = by the length of the range, 30 = the 30-wheel, 210 = the 210-wheel (6 sub-sieves by p mod 7)
extern int g_RangeSieveLayout;
// the vector code of the scan: -1 = the best one of the CPU (RangeSieveInit sets the level it uses), 0 = SSE2, 1 = AVX2, 2 = AVX-512
extern int g_RangeSieveScanLevel;
// the deepest tier (0..5: Q = 1021, 2039, 4093, 8191, 16381, 32749) for the ranges of each class of the window height W (below 1e-4,
// 1e-3, 1e-2, 0.03, 0.1, 0.3, 1, more)
extern int g_RangeSieveMaxTierByW[8];
extern std::atomic<num64> g_RangeSieveVerifyCount;
extern std::atomic<num64> g_RangeSieveVerifyErrors;

// SearchRange uses the sieve only for ranges with at least this many primes: below it, checking candidates one by one is faster.
// Default: 512 primes, the measured crossover point for the 10^22 search limit.
extern num64 g_RangeSieveMinPrimes;

// Must be called after PrimeTablesInit
void RangeSieveInit();

// Statistics of all calls of SearchRangeSieve: ranges, numbers in the ranges, sieve survivors (sums inside the window), survivors
// with a prime p (checked exactly), survivors whose sieve sum didn't fit the factors found (checked with CheckPair128, never seen),
// survivors where CheckSurvivor and CheckPair128 disagree (only counted in builds with RANGE_SIEVE_CHECK_COMPARE: every survivor is
// checked both ways); RangeSievePrintStats prints them to stderr
extern std::atomic<num64> g_RangeSieveStats[12];
// band sieve of the near-degenerate ranges (1 = on); the test band (TestRangeSieve: L0 and the top of the band for the next calls,
// hi = 0: the band of the range)
extern int g_RangeSieveBand;
extern int g_RangeSieveBandMinSpanLog2; // (/sieve_band_min k: the shortest range with a band where DecideAtMost2Factors can't decide all numbers, 2^k)
extern double g_RangeSieveTestBandL0;
extern double g_RangeSieveTestBandHi;
#ifdef RANGE_SIEVE_VERIFY_BAND
extern std::atomic<num64> g_bandVerifyCount;
extern std::atomic<num64> g_bandVerifyErrors;
#endif
void RangeSievePrintStats();

// Exact number of primes in [p_lo, p_hi], using the main prime table (p_hi must not exceed the table's bound)
num64 RangeSieveCountPrimes(num64 p_lo, num64 p_hi);

// Checks all primes p in [p_lo, p_hi] for the range r (same results as the per-candidate loop in SearchRange)
// (tierHint >= 0: the sieve tier 0..5 instead of the one chosen by the length and the window height of the range)
void SearchRangeSieve(const RangeData& r, num64 p_lo, num64 p_hi, unsigned int is_over_abundant_mask, num64 num_primes, int tierHint = -1);

// Same, for M0 given by its factorization (primes in increasing order)
void SearchRangeSieve(num128 M0, num128 sumM0, const num64* M0_primes, const unsigned int* M0_exponents, unsigned int M0_num_factors, num64 p_lo, num64 p_hi, unsigned int is_over_abundant_mask, num64 num_primes, int tierHint = -1);

// The tier of the last call of SearchRangeSieve in this thread, and the share of its time spent on the survivors (IsPrime,
// CheckSurvivor)
void RangeSieveLastCall(int& tier, double& checkShare);

// Self-tests (used by "/test")
bool TestRangeSieve();
