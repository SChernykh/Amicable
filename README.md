# Amicable
This source code implements the leading-edge algorithm for running an exhaustive search for amicable pairs.

Performance when compiled with Visual Studio 2015 update 3 (profile-guided optimizations on) and running on Core i7-4770K:

- All amicable pairs < 10^9: 0.13 seconds (586 pairs)
- All amicable pairs < 10^11: 12.78 seconds (3340 pairs)
- All amicable pairs < 10^13: 23 minutes 40.93 seconds (17519 pairs)
