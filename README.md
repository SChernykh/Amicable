# Amicable
This source code implements the leading-edge algorithm for running an exhaustive search for amicable pairs.

This is a special version for BOINC project https://sech.me/boinc/Amicable/

Performance when compiled with Visual Studio 2015 update 3 (profile-guided optimizations on) and running on Core i7-4770K:

- All amicable pairs < 10^9: 0.118 seconds (586 pairs)
- All amicable pairs < 10^11: 11.858 seconds (3340 pairs)
- All amicable pairs < 10^13: 21 minutes 23.53 seconds (17519 pairs)
