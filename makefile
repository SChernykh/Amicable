all: Amicable/*.* primesieve/src/primesieve/*.cpp primesieve/include/*.* primesieve/include/primesieve/*.*
	mkdir -p release
	g++ -o release/amicable -static -std=c++11 -Wall -Wextra -Wno-unknown-pragmas -O3 -fprofile-dir=release -fprofile-generate -I primesieve/include -I Amicable -Wl,--whole-archive -lpthread -Wl,--no-whole-archive Amicable/*.cpp primesieve/src/primesieve/*.cpp
	release/amicable /instrument
	g++ -o release/amicable -static -std=c++11 -Wall -Wextra -Wno-unknown-pragmas -O3 -fprofile-dir=release -fprofile-use -fprofile-correction -I primesieve/include -I Amicable -Wl,--whole-archive -lpthread -Wl,--no-whole-archive Amicable/*.cpp primesieve/src/primesieve/*.cpp
	strip release/amicable

test:
	release/amicable /test

clean:
	rm -r release
