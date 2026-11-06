# Amicable
This source code implements the leading-edge algorithm for running an exhaustive search for amicable pairs.

This is an OpenCL version for BOINC project https://sech.me/boinc/Amicable/

## Search limit 10^22

The program is configured for the search limit 10^22 (`Definitions.cpp`, `PrimeTables.h`).

Both search stages use an arithmetic-progression sieve on the GPU (`kernel_sieve.cl`, `SieveGPU.h`, `SieveGPU.cpp`) instead of
checking every number by trial division:

- first stage (`/from ... /to ...`): for a fixed M0, N(p) = sigma(M0) * (p + 1) - M0 * p is linear in p, so the prime factors
  q <= Q of N(p) can be sieved for all p at once. Only the few primes p whose sieved sum of log(sigma(q^e)/q^e) can still give
  sigma(N) = sigma(M0) * (p + 1) are checked exactly (with the same kernels as before). Prime numbers p are not needed on the GPU,
  and the host only builds a small prime table (the full table up to 10^11 is only built for `/nosieve` and for `/lpp 2`, `/lpp 3`).
- large primes (`/lpr P1 P2`): every candidate c is sieved over p in [P1, P2] in the same way (M0 = c). Candidates are uploaded in
  batches, so GPU memory doesn't depend on the number of candidates.

## Command line options

- `/nosieve` - the original algorithm (trial division of every number)
- `/test [N]` - self-tests of the first stage sieve with every N-th known pair (`c2_*.txt` files must be in the current directory):
  all known pairs through phases 1-3, the sieve window, the sieve sums on the GPU against the host, and full sieve runs around known pairs
- `/test_lp [N]` - self-tests of the large primes sieve (known pairs with the largest prime factor above 10^11)
- `/test_reduction [N]` - only the host-side test of the range reduction (the necessary conditions that shrink the ranges)
- `/sieve_profile` - prints the GPU time of the sieve and of the survivor checks
- `/sieve_q Q1,Q2,...` - largest sieving prime of each tier, up to 4 tiers (default `2039,8191,16381,32749`)
- `/sieve_tier_min N1,N2,...` - a range uses tier t if it has at least Nt primes (default `0,2048,16384,16384`; the last tier only
  for ranges with a narrow window, see `/sieve_high_w`). A higher tier is used
  automatically when the constant part of the sum would let too many numbers through.
- `/sieve_high_w W` - the last tier only for ranges with a window narrower than W (default 0.03)
- `/sieve_wg N` - work-group size of the sieve (default 1024; halved while it is above the device limit)
- `/sieve_seg N` - positions per row in one segment, a power of 2 from 512 (smaller segments can't hold the task tables) to 4096 (default 1024; halved while it is above 32 positions
  per work-item or the local memory of the device is too small: 8 * N * 4 bytes + tables; 512 needs 29 KB)
- `/sieve_chunk N` - segments per chunk of a work-group (default 18, at most 128)
- `/sieve_hits N` - target number of hits per task (default 16)
- `/sieve_wheel N` - moduli >= N use wheel stepping (default 3/4 of the segment length)
- `/sieve_batch N` - 2^N positions per batch (default 35)
- `/sieve_dry_run` - only the host part (for profiling)
- `/nokill`, `/noreduce` - without the divisibility kill of the wheel primes, without the necessary conditions of the range reduction
  (for comparison)
- `/analyze R`, `/explore N` - development aids: survivors of every tier for positions sampled at rate R (host only); the
  factorization where a work unit of N numbers starting at `/from` ends

After the kernels are compiled, the app runs a quick self-test: known pairs through phases 1-3 (as before) and, with the sieve, full
sieve runs around 48 known first-stage pairs (17 of them through the band sieve) or 4 large primes pairs (`self_test_sieve.inl`).
A device or configuration that computes wrong sums fails it, and the app exits with "Self-test failed" instead of missing pairs.
