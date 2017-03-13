# Amicable
This source code implements the leading-edge algorithm for running an exhaustive search for amicable pairs.

This is a special version for BOINC project https://sech.me/boinc/Amicable/

Performance when compiled with Visual Studio 2015 update 3 (profile-guided optimizations on) and running on Core i7-4770K:

- All amicable pairs < 10^9: 0.1226 seconds (586 pairs)
- All amicable pairs < 10^11: 12.318 seconds (3340 pairs)
- All amicable pairs < 10^13: 22 minutes 50.14 seconds (17519 pairs)
