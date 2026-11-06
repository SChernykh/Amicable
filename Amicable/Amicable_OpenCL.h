#pragma once

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

#include "SieveGPU.h"
#include <condition_variable>
#include <mutex>

extern FILE* g_outputFile;

// Command-line options for the sieve (see Amicable.cpp)
extern bool g_UseSieve;
extern bool g_SieveDryRun;
extern bool g_SieveProfile;
extern unsigned int g_SieveBatchShift;
extern bool g_RunTests;
extern double g_AnalyzeRate;
extern bool g_SieveReduce;
extern bool g_SieveKill;
extern num64 g_ExploreNumbers;
extern bool g_TestReductionOnly;
extern unsigned int g_TestStride;
extern SieveGPUConfig g_SieveConfig;

struct RangeData;
struct AmicableCandidate;

#pragma pack(push, 1)
struct RangeDataGPU
{
	RangeDataGPU(num128 _M0, num128 _sumM0, unsigned int _start_prime_index, unsigned int _size) : M0(_M0), sumM0(_sumM0), start_prime_index(_start_prime_index), size(_size) {}

	num128 M0;
	num128 sumM0;
	unsigned int start_prime_index;
	unsigned int size;
};
#pragma pack(pop)

static_assert(sizeof(RangeDataGPU) == sizeof(num64) * 5, "Invalid RangeDataGPU size");

struct LookupDataGPU
{
	LookupDataGPU() : range_index(0), range_offset(0) {}
	LookupDataGPU(unsigned int i, unsigned int offset) : range_index(i), range_offset(offset) {}

	unsigned int range_index;
	unsigned int range_offset;
};

static_assert(sizeof(LookupDataGPU) == sizeof(unsigned int) * 2, "Invalid LookupDataGPU size");

// Ranges and chunks of one sieve batch (first stage), see AddSieveRangeRecord and ProcessSieveBatch
struct SieveBatch
{
	std::vector<num64> ranges;
	std::vector<unsigned int> chunks;
	std::vector<num64> extra;
	num64 positions = 0;
};

// capacity of the buffer of found numbers (kernel.cl: NumberFound)
enum { AmicableFoundMax = 8192 };

class OpenCL
{
public:
	OpenCL(const char* preferences);
	~OpenCL();

	bool Run(int argc, char* argv[], char* startFrom, char* stopAt, unsigned int largestPrimePower, num64 startPrime, num64 primeLimit);

	bool Test();

	struct SPairToCheck
	{
		SPairToCheck() : M(0), targetSum(0) {}
		SPairToCheck(num128 m, num128 sum) : M(m), targetSum(sum) {}

		num128 M;
		num128 targetSum;
	};
	bool TestPairs(std::vector<SPairToCheck>& pairs);

	bool AddRange(const RangeData& r, num64 offset, bool& added);

	FORCEINLINE void operator()(num64 p)
	{
		myLargePrimes[myLargePrimesCount++] = p;
		if (myLargePrimesCount >= myLargePrimesMaxCount)
		{
			PassLargePrimesToThread();
		}
	}

private:
	int SetKernelSize(int size);

	bool GetPlatformID(cl_platform_id* clSelectedPlatformID);
	bool WaitForQueue(cl_command_queue queue, cl_event event);
	bool GetCounter(cl_command_queue queue, cl_mem buf, unsigned int &counter);
	bool ResetCounter(cl_command_queue queue, cl_mem buf);
	bool GetAndResetCounter(cl_command_queue queue, cl_mem buf, unsigned int &counter);
	bool GetAndResetCounter(cl_command_queue queue, cl_mem buf, num64 &counter);
	unsigned int GetMaxPhaseSize(num64 total_size, unsigned int max_size);

	bool RunRanges(char* startFrom, char* stopAt, num64& offset, unsigned int ChunkSizeShift);

	bool InitSieve(cl_device_id device, num64 localMemSize, num64 maxWorkGroupSize, std::string& buildOptions);
	bool CreateSieveBuffers();
	bool GetPrimeLimit(const RangeData& r, num64& p_lo, num64& p_hi) const;
	bool AddSieveRange(const RangeData& r);
	// W: height of the window above the known part (RangeReduction), -1 = unknown; it predicts the number of survivors
	bool AddSieveRangeRecord(const num64* record, num64 p_lo, num64 p_hi, double W = -1.0);
	bool ProcessSieveBatch(SieveBatch& batch);
	bool ProcessSieveBatch() { return WaitSieveBatch() && ProcessSieveBatch(mySieveBatch); }
	// The GPU work of a batch runs in a worker thread while the host prepares the next batch. SubmitSieveBatch waits for the previous
	// batch, saves a checkpoint at the first range of the submitted batch (if checkpoint is true and that range is known: all ranges
	// before it are done), then hands mySieveBatch to the worker. The worker is one persistent thread (SieveWorkerLoop, started with the
	// first batch, stopped in the destructor): a new thread for every batch meant thousands of short-lived threads per work unit.
	bool SubmitSieveBatch(bool checkpoint, num64 offset);
	bool WaitSieveBatch();
	void SieveWorkerLoop();
	// chunks [chunkBegin, killChunks) use mySieveRangesKill
	bool RunSieveLaunches(unsigned int chunkBegin, unsigned int chunkEnd, unsigned int killChunks, num64 positionsPerChunk, unsigned int numExtra, const SieveBatch* batch, num64& survivors, bool& overflow);
	bool FinishSieveBatch(bool& redo);
	bool RunSieveTests();
	bool TestSieveSums(const std::vector<struct KnownPair>& pairs);
	bool TestSieveKnownPairs(const std::vector<struct KnownPair>& pairs, const std::vector<num128>& allAmicable, bool verbose = true);
	bool TestSieveStartup();
	bool SetPhase2KernelArguments();
	bool TestSieveLargePrimesPair(const num128 m, num64 c, num64 sc, num64 p, unsigned int tier, const std::vector<num128>* allAmicable, num64& survivorsTotal);
	bool RunLargePrimes(num64 startPrime, num64 primeLimit);
	bool RunLargePrimesSieve(num64 startPrime, num64 primeLimit);
	bool UploadPrimeBitmap(num64 P1, num64 P2);
	bool SieveCandidates(const AmicableCandidate* candidates, unsigned int count, num64 P1, num64 P2, unsigned int tier, num64& survivorsTotal);
	bool TestSieveLargePrimes(unsigned int stride, const std::vector<num128>& allAmicable);

	void PassLargePrimesToThread();
	void ProcessLargePrimesThread();
	bool ProcessLargePrimes();

	bool ProcessNumbers();
	bool ProcessNumbersPhases2_3(unsigned int numbers_in_phase2);
	bool RunPhase3(unsigned int phase3_numbers_buf_count, int firstIteration);
	bool SortPhase2Numbers(unsigned int count);
	// (kernel.cl: PHASE2_BUCKETS) the numbers of a phase 2 run are sorted if there are at least Phase2SortMinCount
	enum { Phase2Buckets = 64, Phase2SortMinCount = 1 << 16 };
	enum { SieveBandEntry = 8, SieveMaxBandChunks = 64, SieveBandMaxOctaves = 4 };
	cl_kernel myPhase2BucketCount = nullptr;
	cl_kernel myPhase2BucketScatter = nullptr;
	cl_kernel myPhase2Prefilter = nullptr;
	cl_mem myPhase2BucketsBuf = nullptr;
	num64 mySieveStatPhase2Sorts = 0;
	std::vector<std::pair<int, cl_event>> myPhase2ProfileEvents;   // (/sieve_profile) kind, event
	double mySieveStatPhase2Times[4] = {};
	double mySieveStatPhase2Wall = 0.0;
	double mySieveStatCheckTimes[3] = {};   // (/sieve_profile) the survivor checks: first pass, SieveTestDeferred, prime p
	double mySieveStatPhase3Wall = 0.0;    // (/sieve_profile) host time of RunPhase3, and the numbers of its first iterations
	num64 mySieveStatPhase3Numbers = 0;
	double mySieveStatPhase3IterTime[20] = {};
	num64 mySieveStatPhase3IterNumbers[20] = {};
	num64 mySieveStatPhase3IterLaunches[20] = {};
	bool FlushPhase3Backlog();
	// numbers of phase 3 that are collected in the sieve batches (RunPhase3): when fewer than Phase3BacklogMinCount are left
	enum { Phase3BacklogMinCount = 1 << 16 };
	unsigned int Phase3BacklogCapacity() const { return myPhase2MaxNumbersCount / 4; }
	void CleanupRanges();

	bool SaveFoundNumbers();
	bool SaveProgressForRanges(const RangeData& r, num64 offset, num64 numbersProcessed = ~num64(0));
	// resets the counter of found numbers and sets the capacity of their buffer
	bool SetFoundCapacity(unsigned int capacity);
	bool SaveProgressForLargePrimes(num64 firstPrime, num64 lastPrime, num64 offset, bool queueEmpty);

	unsigned int myLargestPrimePower;

	cl_context myGPUContext;
	cl_program myProgram;
	cl_program myProgramKill = nullptr;        // first stage sieve with the divisibility kill (for the chunks of the last tier)
	cl_kernel mySieveRangesKill = nullptr;
	std::string mySieveOptions;                // sieve build options of myProgram
	std::string mySieveKillOptions;            // and of myProgramKill (empty: no second program)
	std::string myBuildOptionsKill;            // all build options of myProgramKill

	cl_command_queue myQueue;
	cl_kernel mySaveCounter;
	cl_kernel mySearchMultipleRanges;
	cl_kernel mySearchLargePrimes;
	cl_kernel myCheckPairs;
	cl_kernel myCheckPairPhase2;
	cl_kernel myCheckPairPhase3;
	cl_kernel myCheckPairPhase3Wide = nullptr;
	cl_mem myPhase3WideStatusBuf = nullptr;
	cl_mem mySpsp2Buf = nullptr;
	cl_kernel mySieveRanges;
	cl_kernel myCheckSieveSurvivors;
	cl_kernel mySieveLaunchDone;
	cl_kernel mySieveTestDeferred = nullptr;  // primality test of p for the numbers that passed phase 1 (first stage sieve)
	cl_kernel myCheckSieveSurvivorsPrime = nullptr; // second check of the survivors with a prime p
	cl_mem mySmallPrimesBuf;
	std::vector<cl_mem> myMainBuffers;
	cl_mem myPrimeInversesBuf;
	cl_mem myPrimeReciprocalsBuf;
	cl_mem myPQ_Buf;
	cl_mem myPowerOf2SumInverses_Buf;
	cl_mem myPowersOfP_128SumInverses_offsets_Buf;
	cl_mem myPowersOfP_128SumInverses_Buf;
	cl_mem myPrimeInverses_128_Buf;
	cl_mem mySumEstimates_128_Buf;
	cl_mem myLargePrimesBuf;
	cl_mem myRangesTable_Buf;
	cl_mem myRangesLookupTable_Buf;
	cl_mem myPhase1_offset_to_resume_buf;
	cl_mem myPhase2_numbers_count_buf;
	cl_mem myPhase2_numbers_buf;
	cl_mem myPhase3_numbers_count_buf;
	cl_mem myPhase3_numbers_buf;
	cl_mem myPhase3BacklogBuf = nullptr;
	unsigned int myPhase3BacklogCount = 0;
	int myPhase3BacklogIteration = 0;      // the first phase 3 iteration of the numbers in the backlog
	bool myPhase3BacklogEnabled = false;
	cl_mem myAmicable_numbers_count_buf;
	cl_mem myAmicable_numbers_data_buf;

	// Sieve (first stage and large primes)
	bool myUseSieve;
	bool myLargePrimesSieve;
	cl_uint mySieveChunkOffsetArg;
	unsigned int mySieveMaxCandidates;
	SieveGPUTables mySieve;
	cl_mem mySievePrimeRecsBuf;
	cl_mem mySieveModRecsBuf;
	cl_mem mySieveModRecs2Buf;
	cl_mem mySieveTtabBuf;
	cl_mem mySieveTasksBuf;
	cl_mem mySieveTierInfoBuf;
	cl_mem mySieveWheelTabsBuf;
	cl_mem mySieveModTaskIdxBuf;
	cl_mem mySieveModTaskListBuf;
	cl_mem mySieveRangesBuf;
	cl_mem mySieveChunksBuf;
	cl_mem mySurvivorCountBuf;
	cl_mem mySurvivorsBuf;
	cl_mem mySieveStatsBuf;
	cl_mem mySieveDeferredCountBuf = nullptr;
	cl_mem mySieveDeferredListBuf = nullptr;
	// band sieve (kernel_sieve.cl: SieveBand, kernel.cl: DecideByBand) for the ranges with W < SieveBandMaxW: record[21] = L0 | the index of
	// the smallest prime of g in the band << 32, record[22] = the indexes in smallPrimes [begin, end) of the band (the primes in (L0, L0 2^T]))
	cl_kernel mySieveBand = nullptr;
	cl_mem mySieveBandTableBuf = nullptr;    // SieveBandEntry uints per band chunk, for all launches of a round
	cl_mem mySieveBandBitmapBuf = nullptr;   // the index of the smallest band prime of every position of the band chunks of a launch
	std::vector<unsigned int> mySmallPrimesHost;
	num64 mySieveStatBandChunks = 0;
	double mySieveStatBandKernelTime = 0.0;
	void SetBandParameters(num64* record, double W, num64 p_hi, const num64* primes, unsigned int numPrimes) const;
	void SetTestBandParameters(num64* record, const struct KnownPair& pair) const;
	mutable num64 mySieveStatBandTestRanges = 0;
	bool TestBandBitmap(const SieveBatch* batch, const cl_uint* entries, cl_uint count, cl_uint bytes);    // survivor records (as many as the survivor buffer holds)
	cl_mem mySievePrimeCountBuf = nullptr;      // survivors with a prime p (SieveTestDeferred)
	// measured survivors / estimate of SieveSurvivorRateByW for the ranges with W < 1e-3 (almost every N with a rough cofactor
	// survives there) and for the others (launch sizes; the launches don't mix the two groups)
	double mySieveRateScale[2] = { 1.0, 1.0 };
	double mySievePhase2PerSurvivor[2] = { 0.3, 0.3 }; // numbers for phase 2 per survivor, per group of mySieveRateScale
	cl_mem myPrimeBitsBuf;
	std::vector<cl_uint> myPrimeBits;
	num64 myPrimeBitsP1;
	num64 myPrimeBitsP2;
	unsigned int mySieveMaxRanges;
	unsigned int mySieveMaxChunks;
	unsigned int mySurvivorCapacity;
	size_t mySieveCheckGlobalSize;
	num64 mySieveMaxPositions;
	SieveBatch mySieveBatch;            // the batch being built
	SieveBatch mySieveBatchRun;         // the batch of the worker thread
	std::thread mySieveWorker;
	std::mutex mySieveWorkerMutex;        // guards the three flags below and mySieveWorkerResult
	std::condition_variable mySieveWorkerCV;
	bool mySieveJobPending = false;       // a batch (mySieveBatchRun) is waiting for the worker
	bool mySieveJobDone = false;          // the worker has finished the last batch it was given
	bool mySieveWorkerStop = false;       // the worker ends (destructor)
	bool mySieveWorkerRunning = false;    // a batch was submitted and not waited for yet (main thread only)
	bool mySieveWorkerResult = true;
	RangeData* mySieveBatchFirstRange = nullptr; // first range of mySieveBatch (valid if mySieveBatchHasFirst)
	bool mySieveBatchHasFirst = false;
	num64 mySieveBatchStartNumbers = 0; // myNumbersProcessedTotal before mySieveBatchFirstRange
	double mySieveSurvivorRate;
	double mySievePhase2Rate;        // numbers for phase 2 per position (to limit the size of launches)
	num64 mySieveRoundPositions2;     // positions whose phase 2 numbers are still in the buffer
	num64 mySieveStatPositions;
	num64 mySieveStatSurvivors;
	num64 mySieveStatRanges;
	num64 mySieveStatChunks;
	num64 mySieveStatBatches;
	num64 mySieveStatPhase2;
	num64 mySieveStatOverflows;
	num64 mySieveStatLaunchedPositions = 0;  // including the launches that ran again
	num64 mySieveStatLaunches = 0;
	num64 mySieveStatRounds = 0;
	double myRedStat[4] = {};
	double myRedHist[64] = {};
	double myRedHistW[8] = {};
	num64 mySieveStatSkipped = 0;
	num64 myAnalyzeSamples = 0;
	num64 myAnalyzeCounts[SieveMaxTiers][6] = {};
	double mySieveStatTime;
	double mySieveStatSieveKernelTime;
	double mySieveStatCheckKernelTime;

	size_t myWorkGroupSize;

	unsigned int myMaxRangesCount;
	unsigned int myMaxLookupTableSize;
	
	unsigned int myPhase1MaxKernelSize;
	unsigned int myPhase2MinNumbersCount;
	unsigned int myPhase2MaxNumbersCount;

	num64 myOldTimerResolution;

	const unsigned char* myMainBufferData;
	num64 myMainBufferSize;
	num64 myMainBufferTotalSize;

	std::vector<RangeDataGPU> myRanges;
	unsigned int myTotalNumbersInRanges;
	num64 myNumbersProcessedTotal;
	num64 myAmicableNumbersFound;

	const char* myPreferences;
	std::string myBuildOptions;

	std::vector<LookupDataGPU> myLookupTable;

	std::streamsize myStdOldPrecision;

	struct SFoundPair
	{
		num64 M;
		num64 M_high;
		num64 N;
		num64 dummy;
	};

	std::vector<SFoundPair> myFoundPairs;

	num64* myLargePrimes;
	unsigned int myLargePrimesCount;
	unsigned int myLargePrimesMaxCount;
	num64 myLargePrimesStartOffset;

	Semaphore myLargePrimesReady;
	Semaphore myLargePrimesReceived;

	enum
	{
		ourZeroBufSize = 65536,

		// Large primes sieve: the interval of p is processed in parts of up to 30 * 2^LargePrimesBitmapShift numbers (bitmap of 32 MB)
		LargePrimesBitmapShift = 25,
	};
	static unsigned char* ourZeroBuf;
};

#if ((defined __cpp_variadic_templates) && (__cpp_variadic_templates >= 200704))

FORCEINLINE cl_int setKernelArgumentsWithOffset(cl_kernel, cl_uint) { return CL_SUCCESS; }

template<typename T0, typename... Arguments>
FORCEINLINE cl_int setKernelArgumentsWithOffset(cl_kernel kernel, cl_uint offset, const T0& arg0, const Arguments&... args)
{
	cl_int ciErrNum = clSetKernelArg(kernel, offset, sizeof(T0), &arg0);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	return setKernelArgumentsWithOffset(kernel, offset + 1, args...);
}

template<typename... Arguments>
FORCEINLINE cl_int setKernelArguments(cl_kernel kernel, const Arguments&... args)
{
	return setKernelArgumentsWithOffset(kernel, 0, args...);
}

#else
#include "OpenCL.inl"
#endif
