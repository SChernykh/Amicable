WARNINGS=-Wall -Wextra -Wpedantic -Wstrict-overflow=5 -Wshadow -Warray-bounds=2
INCLUDES=-I primesieve/include -I Amicable -I ../boinc/api -I ../boinc/lib
LIBS=-L ../boinc/api -L ../boinc/lib -Wl,--whole-archive -lboinc_api -lboinc -lpthread -Wl,--no-whole-archive

# Build profile-optimized, stripped and ready to use binary
all: Amicable/*.* primesieve/src/primesieve/*.cpp primesieve/include/*.* primesieve/include/primesieve/*.*
	mkdir -p release
	g++ -o release/amicable -static -std=c++11 $(WARNINGS) -O3 -fprofile-dir=release -fprofile-generate $(INCLUDES) $(LIBS) Amicable/*.cpp primesieve/src/primesieve/*.cpp
	release/amicable /instrument
	g++ -o release/amicable -static -std=c++11 $(WARNINGS) -O3 -fprofile-dir=release -fprofile-use -fprofile-correction $(INCLUDES) $(LIBS) Amicable/*.cpp primesieve/src/primesieve/*.cpp
	strip release/amicable

# Run clang's static analyzer + check for all clang's warnings
analyze: Amicable/*.* primesieve/src/primesieve/*.cpp primesieve/include/*.* primesieve/include/primesieve/*.*
	clang-3.9 -std=c++11 --analyze -Weverything -O3 -I primesieve/include $(INCLUDES) -pthread Amicable/*.cpp primesieve/src/primesieve/*.cpp

# Generate assembly output for the crucial file (Engine.cpp)
engine_asm:
	g++ -S -fverbose-asm -std=c++11 $(WARNINGS) -O3 -I primesieve/include -I Amicable -pthread Amicable/Engine.cpp

# Check that compiled binary doesn't miss any known amicable numbers, run some other internal tests
test:
	release/amicable /test

clean:
	rm -r release
