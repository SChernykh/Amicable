WARNINGS=-Wall -Wextra
INCLUDES=-I primesieve/include -I Amicable -I ../boinc/api -I ../boinc/lib
LIBS=-L ../boinc/api -L ../boinc/lib -Wl,--whole-archive -lboinc -lboinc_api -lboinc_opencl -lpthread -lrt -lOpenCL -Wl,--no-whole-archive

# Build profile-optimized, stripped and ready to use binary
all: Amicable/*.* primesieve/src/primesieve/*.cpp primesieve/include/*.* primesieve/include/primesieve/*.*
	mkdir -p release
	cd Amicable
	g++ -o release/amicable -static-libstdc++ -static-libgcc -std=c++17 $(WARNINGS) -O3 $(INCLUDES) Amicable/*.cpp primesieve/src/primesieve/*.cpp $(LIBS)
	strip release/amicable

# Check that compiled binary doesn't miss any known amicable numbers, run some other internal tests
test:
	release/amicable /test

clean:
	rm -r release
