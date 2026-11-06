#pragma once

struct RangeData;

num64 SearchRange(const RangeData& r);
num64 SearchRangeSquared(const RangeData& r);
num64 SearchRangeCubed(const RangeData& r);

void SearchLargePrimes(volatile num64* SharedCounterForSearch, const num64 StartPrime, const num64 PrimeLimit, num64 &sharedCounterValue);

void CheckPair128NoInline(const num128 n1, const num128 targetSum);

// The rest of CheckPair128 for a partner n2 (odd, > 1 or 1) without prime factors below the prime with index startIndex of
// PrimeInverses128 (0: 3, startIndex >= CompileTimePrimesCount - 1), whose sigma must be targetSum: reports n1 if it is amicable
// Band sieve of a near-degenerate range (see BandSieve): the band primes are those of (lo, hi] (indexes [sub[0], sub[255]) of the
// tables of CheckCofactor, in 255 sub-bands [sub[s - 1], sub[s]) of equal numbers of primes); position k = (p - P0) / 30 * 8 +
// (class of p mod 30) has a byte: the sub-band (1..255) of the smallest band prime that divides A1 p + B1, 0 if none does. A number
// with L = N / D >= L0 whose smallest prime factor has to be <= hi (the band bound K L + K + 2) is not amicable if its sub-band is 0;
// otherwise that sub-band holds its smallest prime factor above lo. (15 sub-bands in a nibble: the scans of those sub-bands were 61%
// of t1b3.)
enum { BAND_SUBBANDS = 255 };
struct CofactorBand
{
	const unsigned char* subBands;
	num64 P0;
	double L0;
	double lo; // lo = max(L0, Q)
	double hi;
	unsigned int sub[BAND_SUBBANDS + 1];
};

// (band, p: the band sieve of the range of n1 = M0 * p, or nullptr; targetOrig: sigma(M0) * (p + 1), for the check of a report that
// used the band)
void CheckCofactor(const num128 n1, num128 targetSum, num128 n2, const unsigned int startIndex, const CofactorBand* band, num64 p, const num128& targetOrig);

// Fills the bytes of the band for the positions p in [P0, p_hi] (p coprime to 30, P0 a multiple of 30; subBands:
// ((p_hi - P0) / 30 + 1) * 8 bytes, zeroed by the caller) and band.sub; band.lo and band.hi must be set, hi below
// CofactorTableMaxPrime(). Returns the number of primes of the band.
unsigned int BandSieve(const num128& A1, const num128& B1, num64 P0, num64 p_hi, CofactorBand& band, unsigned char* subBands);
double CofactorTableMaxPrime();
// band.sub from band.lo and band.hi (as BandSieve), false if the band is outside the tables; the sub-band (1-255) of a prime of the
// band, 0 if it is not a prime of the band
bool CofactorBandSetup(CofactorBand& band);
unsigned int CofactorBandSubBandOf(const CofactorBand& band, num64 prime);
// the number of primes of the tables of CheckCofactor in (lo, hi]
double CofactorTablePrimeCount(double lo, double hi);
// The statistics of the band in this thread since the last call (numbers rejected by their sub-band 0, numbers with a sub-band,
// numbers whose band was not covered), then set to 0
void CofactorBandStats(num64* out3);
// (/test) the bounds L and E of CheckCofactor for n2, false if there are none; survivors with p == g_testBandP that pass the band check
bool CofactorBandBounds(const num128& n2, const num128& targetSum, double& L, double& E);
extern num64 g_testBandP;
extern num64 g_testBandPasses;
// The tables of CheckCofactor (after PrimeTablesInit; simdLevel: 0 = SSE2, 1 = AVX2, 2 = AVX-512)
void CheckCofactorInit(int simdLevel);
num64 MaximumSumOfDivisors3NoInline(const num64 a, const num64 p0, const num64 a_div_p0);
	
void SetNumFoundPairsInThisThread(unsigned int value);
unsigned int GetNumFoundPairsInThisThread();
// (tests, one thread) every number found is also appended here if set
extern std::vector<num128>* g_FoundNumbersLog;
// the over-abundance mask of a range (SearchRange)
unsigned int RangeOverAbundantMask(const RangeData& r);

extern bool g_PrintNumbers;
extern FILE* g_outputFile;
