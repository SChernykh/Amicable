# Amicable
This source code implements the leading-edge algorithm for running an exhaustive search for amicable pairs.

This is a special version for BOINC project https://sech.me/boinc/Amicable/

Performance when compiled with Visual Studio 2015 update 3 (profile-guided optimizations on) and running on Core i7-4770K:

- All amicable pairs < 10^9: 0.118 seconds (586 pairs)
- All amicable pairs < 10^11: 11.858 seconds (3340 pairs)
- All amicable pairs < 10^13: 21 minutes 23.53 seconds (17519 pairs)

## Search limit 10^22

The program is configured for the search limit 10^22 (`Definitions.cpp`, `PrimeTables.h`; `AmicableCandidate` stores 40-bit values
and sums).

The first stage (`/from ... /to ...`) uses an arithmetic-progression sieve (`RangeSieve.h`, `RangeSieve.cpp`) for the ranges with at
least 512 primes p instead of checking every number by trial division: for a fixed M0, N(p) = sigma(M0) * (p + 1) - M0 * p is linear
in p, so the prime factors q <= Q of N(p) are sieved for all p of a range at once, and only the primes p whose sieved sum can still
give sigma(N) = sigma(M0) * (p + 1) are checked exactly. Long ranges are shared by all threads, the sieve depth adapts to the measured
cost of the survivor checks, and ranges with a very narrow window also sieve the cofactors in a band (as in the OpenCL version).
On x86-64 the hot loops have SSE2, AVX2 and AVX-512 versions, chosen at run time; other architectures (ARM, RISC-V) compile the
portable code (`AMICABLE_PORTABLE` selects it on x86-64 too). The statistics of the sieve are printed to stderr at the end of a run.

## Command line options

- `/nosieve` - the original algorithm (trial division of every number)
- `/test` - self-tests (the `c2_*.txt` files must be in the current directory): all known pairs through `CheckPair`, and the range
  sieve around every known first-stage pair (the window of every tier, the survivor checks, the band, sampled full sieve runs)
- testing: `/sieve_tier t` (force the depth tier 0..5), `/sieve_wtiers a,b,...` (the deepest tier for each of the 8 classes of the
  window height), `/sieve_scan 0|1|2` (SSE2, AVX2, AVX-512: at most what the CPU supports; before `/test` on the command line),
  `/sieve_layout 30|210` (the wheel of the sieve), `/sieve_band 0|1` (the band sieve off/on), `/sieve_band_min k` (the band only in
  ranges of at least 2^k numbers, default 14)
