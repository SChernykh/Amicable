# Amicable
This source code implements the leading-edge algorithm for running an exhaustive search for amicable pairs.

Performance when compiled with Visual Studio 2015 update 3 (profile-guided optimizations on) and running on Core i7-4770K:

- All amicable pairs < 10^9: 0.118 seconds (586 pairs)
- All amicable pairs < 10^11: 11.858 seconds (3340 pairs)
- All amicable pairs < 10^13: 21 minutes 23.53 seconds (17519 pairs)

# Running Amicable

### No command line parameters
Outputs all amicable pairs up to the limit specified in Definitions.h

### /from factorization1 /to factorization2 [/lpp N]
Searches all amicable pairs with smallest number having factorization between factorization1 and factorization2 and having largest prime factor of the form p^N. N defaults to 1 if it's not specified.

For example, "/from 2^2\*11\*31 /to 2^2\*11\*37" will find this amicable pair among others:
```
32 TeRiele 1990
10020507332=2^2*11*31*61*83*1451
10306191676=2^2*11*31*1693*4463
```
Another example, "/from 5 /to 5^100" will find all amicable pairs coprime to 6, if search limit is large enough. For example, it will find:
```
X44 Walker&Einstein 2001
5480828320492525=5^2*7^2*11*13*19*31*17*23*103*1319
5786392931193875=5^3*7*11*13*19*31*37*43*61*809
```

### /lpr N1 N2
Searches all amicable pairs with largest prime factor in interval \[N1, N2\]. This is used to search amicable numbers which can't be found using "/from ... /to ..." command line because of very large prime factors in their factorization. For example, "/lpr 141000000 142000000" will find this amicable pair:
```
23 Poulet 1929
38453967088=2^4*17*141374879
40433215952=2^4*179*743*19001
```

### /threads N
Use N threads for searching. Use it to override default number of threads if the program doesn't use all CPU cores.

### /bench
Run a benchmark with current search limit specified in Definitions.h

### /test
Run a set of internal tests. You need to have *c2_3.txt, c2_4.txt, ..., c2_20.txt* files in current directory for the tests to run properly. These files can be downloaded https://sech.me/ap/ - go to "Download pairs from ... to ... digits", enter 3 and 20 and click download.
