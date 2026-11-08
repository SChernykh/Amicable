#include <array>
#include "stdafx.h"
#include "Amicable_OpenCL.h"
#include "PrimeTables.h"
#include "RangeGen.h"
#include "sprp64.h"
#include "kernel.cl"
#include "kernel_sieve.cl"

PRAGMA_WARNING(push, 1)
PRAGMA_WARNING(disable : 4091 4917 4625 4626 5026 5027)
#include <boinc_api.h>
#include <boinc_opencl.h>
#include "Spsp2Table.h"
PRAGMA_WARNING(pop)

//#pragma optimize("", off)
//#undef FORCEINLINE
//#define FORCEINLINE NOINLINE

static const char* const CHECKPOINT_LOGICAL_NAME = "amicable_checkpoint";

FILE* g_outputFile = nullptr;

bool g_UseSieve = true;
bool g_SieveDryRun = false;
bool g_SieveProfile = false;
// Development aid: cycles spent in the host parts of the sieve path (printed with the sieve statistics): the preparation of the ranges
// (all threads; [1], [2] and [4] are parts of it), the reduction, PrepareRange, RangeGen::Iterate, EstimatePrimes, AddSieveRangeRecord, all
// of CommitSieveRange, and the main thread waiting for a block of ranges
static num64 g_HostCycles[8] = {};
unsigned int g_SieveBatchShift = 35;
// sieve batches on the GPU at once (SievePipe): 2 by default, 1 for one pipeline
int g_SievePipes = 2;
// helper threads that prepare the ranges of the first stage for the sieve (RunRangesSieve): 1 by default, 0 for none
int g_SieveThreads = 1;
bool g_RunTests = false;
double g_AnalyzeRate = 0;
static num128 Gcd128Host(num128 a, num128 b)
{
	while (b != 0)
	{
		const num128 t = a % b;
		a = b;
		b = t;
	}
	return a;
}
static const unsigned int SieveW30Host[8] = { 1, 7, 11, 13, 17, 19, 23, 29 };
bool g_SieveReduce = true;
bool g_SieveKill = true;
// band sieve for the ranges with W < 1e-4 (SetBandParameters); AMICABLE_NO_BAND=1 switches it off (for comparisons)
static const bool g_SieveBand = (getenv("AMICABLE_NO_BAND") == nullptr);
// (Development aid: AMICABLE_BAND_MAX_W=<W> sets the largest window of the band ranges. 1e-3, with the band for W >= 1e-4 only in
// dense ranges (SetBandParameters): 2^6*79*409*2287*3449.., 39% of its positions at W = 1e-4 - 1e-3, 11.7 -> 6.3 s, but the units of
// the 3^10 subtree 10% slower, e.g. 3^10*5*7*11*661*29131.. 30.1 -> 33.7 s: not the default.)
static const double SieveBandMaxW = getenv("AMICABLE_BAND_MAX_W") ? atof(getenv("AMICABLE_BAND_MAX_W")) : 1e-4;
num64 g_ExploreNumbers = 0;
bool g_TestReductionOnly = false;
unsigned int g_TestStride = 16;
SieveGPUConfig g_SieveConfig;

enum
{
	LogLevel = 0,

	// Limit phase 1 to 2^28 numbers
	Phase1Limit = 1 << 28,
};

#define DEVICE_TYPE_TO_USE CL_DEVICE_TYPE_GPU

#define LOG_ERROR(X) do { std::cerr << __FILE__ << ", line " << __LINE__ << ": " << X << std::endl << std::flush; } while (0)
#define LOG(LEVEL, X) do { IF_CONSTEXPR(LogLevel >= LEVEL) { std::cout << X; IF_CONSTEXPR(LEVEL > 0) { std::cout << std::endl; } } } while (0)
#define CL_CHECKED_CALL(FUNC, ...) do { const cl_int ciErrNum = FUNC(__VA_ARGS__); if (ciErrNum != CL_SUCCESS) { LOG_ERROR(#FUNC" returned error " << ciErrNum); return false; } } while (0)
#define CL_CHECKED_CALL_WITH_RESULT(result, FUNC, ...) do { cl_int ciErrNum; result = FUNC(__VA_ARGS__, &ciErrNum); if (ciErrNum != CL_SUCCESS) { LOG_ERROR(#FUNC" returned error " << ciErrNum); return false; } } while (0)

unsigned char* OpenCL::ourZeroBuf = nullptr;

OpenCL::OpenCL(const char* preferences)
	: myLargestPrimePower(1)
	, myGPUContext(nullptr)
	, myProgram(nullptr)
	, mySaveCounter(nullptr)
	, mySearchMultipleRanges(nullptr)
	, mySearchLargePrimes(nullptr)
	, myCheckPairs(nullptr)
	, mySmallPrimesBuf(nullptr)
	, myPrimeInversesBuf(nullptr)
	, myPrimeReciprocalsBuf(nullptr)
	, myPQ_Buf(nullptr)
	, myPowerOf2SumInverses_Buf(nullptr)
	, myPowersOfP_128SumInverses_offsets_Buf(nullptr)
	, myPowersOfP_128SumInverses_Buf(nullptr)
	, myPrimeInverses_128_Buf(nullptr)
	, mySumEstimates_128_Buf(nullptr)
	, myLargePrimesBuf(nullptr)
	, myRangesTable_Buf(nullptr)
	, myRangesLookupTable_Buf(nullptr)
	, myUseSieve(false)
	, myLargePrimesSieve(false)
	, mySieveChunkOffsetArg(2)
	, mySieveMaxCandidates(1U << 22)
	, mySievePrimeRecsBuf(nullptr)
	, mySieveModRecsBuf(nullptr)
	, mySieveModRecs2Buf(nullptr)
	, mySieveTtabBuf(nullptr)
	, mySieveTasksBuf(nullptr)
	, mySieveTierInfoBuf(nullptr)
	, mySieveWheelTabsBuf(nullptr)
	, mySieveModTaskIdxBuf(nullptr)
	, mySieveModTaskListBuf(nullptr)
	, myPrimeBitsBuf(nullptr)
	, myPrimeBitsP1(0)
	, myPrimeBitsP2(0)
	// (batches of up to 2^16 ranges. With bigger batches the rounds of RunSieveLaunches are sized by the phase 2 budget much more often
	// (a batch of 2^16 small ranges mostly fits in one or two rounds), the phase 2 buffer overflows more often and more positions are
	// sieved again: 3^11 unit 105 with 2^16 / 2^17 / 2^18 ranges 53.6 / 56.1 / 57.9 s (1.07 / 1.15 x the positions at 2^16 / 2^18),
	// unit 109 53.3 / 54.5 / 55.4 s; only the units of many tiny ranges gained 2% with 2^18, 3^11 units 122 and 252)
	, mySieveMaxRanges(65536)
	, mySieveMaxChunks(262144)
	, mySurvivorCapacity(1U << 25)
	, mySieveCheckGlobalSize(524288)
	, mySieveMaxPositions(0)
	, mySieveStatRanges(0)
	, myWorkGroupSize(128)
	, myMaxRangesCount(131072)
	, myMaxLookupTableSize((myMaxRangesCount + 1) * 4) // very conservative estimate to be safe
	, myMainBufferData(nullptr)
	, myMainBufferSize(0)
	, myMainBufferTotalSize(0)
	, myTotalNumbersInRanges(0)
	, myNumbersProcessedTotal(0)
	, myAmicableNumbersFound(0)
	, myPreferences(preferences)
	, myLargePrimes(nullptr)
	, myLargePrimesCount(0)
	, myLargePrimesMaxCount(1048576)
	, myLargePrimesStartOffset(0)
{
	SetKernelSize(21);

	myRanges.reserve(myMaxRangesCount);
	myLookupTable.reserve(myMaxLookupTableSize);
	myFoundPairs.reserve(8192);

	myStdOldPrecision = std::cout.precision();
	std::cout.precision(3);

	myOldTimerResolution = SetHighestTimerResolution();
}

OpenCL::~OpenCL()
{
	for (SievePipe& p : myPipes)
	{
		if (p.mySieveWorker.joinable())
		{
			{
				std::lock_guard<std::mutex> lock(p.mySieveWorkerMutex);
				p.mySieveWorkerStop = true;
			}
			p.mySieveWorkerCV.notify_all();
			p.mySieveWorker.join();
		}
	}
	SetTimerResoluion(myOldTimerResolution);

	std::cout.precision(myStdOldPrecision);

	// the queues, kernels and batch buffers of the pipelines
	for (int i = 0; i < MaxSievePipes; ++i)
	{
		PipeScope scope(i);
		ReleasePipe();
	}

	cl_mem sieveBuffers[] = { myPrimeBitsBuf, mySievePrimeRecsBuf, mySieveModRecsBuf, mySieveModRecs2Buf, mySieveTtabBuf, mySieveTasksBuf, mySieveTierInfoBuf, mySieveWheelTabsBuf, mySieveModTaskIdxBuf, mySieveModTaskListBuf };
	for (cl_mem buf : sieveBuffers)
	{
		if (buf) clReleaseMemObject(buf);
	}
	if (myProgramKill) clReleaseProgram(myProgramKill);

	clReleaseMemObject(mySmallPrimesBuf);
	for (cl_mem buf : myMainBuffers)
	{
		clReleaseMemObject(buf);
	}
	clReleaseMemObject(myPrimeInversesBuf);
	clReleaseMemObject(myPrimeReciprocalsBuf);
	clReleaseMemObject(myPQ_Buf);
	clReleaseMemObject(myPowerOf2SumInverses_Buf);
	clReleaseMemObject(myPowersOfP_128SumInverses_offsets_Buf);
	clReleaseMemObject(myPowersOfP_128SumInverses_Buf);
	clReleaseMemObject(myPrimeInverses_128_Buf);
	clReleaseMemObject(mySumEstimates_128_Buf);
	if (myLargePrimesBuf) clReleaseMemObject(myLargePrimesBuf);
	clReleaseMemObject(myRangesTable_Buf);
	clReleaseMemObject(myRangesLookupTable_Buf);
	clReleaseKernel(mySaveCounter);
	if (mySearchMultipleRanges) clReleaseKernel(mySearchMultipleRanges);
	if (mySearchLargePrimes) clReleaseKernel(mySearchLargePrimes);
	clReleaseKernel(myCheckPairs);
	if (mySpsp2Buf) clReleaseMemObject(mySpsp2Buf);
	clReleaseProgram(myProgram);
	clReleaseContext(myGPUContext);
}

static bool PlatformSupported(const char* platformName)
{
	return true;
}

template<size_t N>
static bool ParseIntegerFromXml(const char* xml, const char (&paramName)[N], int& value)
{
	const char* pos = strstr(xml, paramName);
	if (pos)
	{
		pos += N - 1;
		while (*pos && (*pos < '0'))
		{
			++pos;
		}
		if (('0' <= *pos) && (*pos <= '9'))
		{
			int k = 0;
			do
			{
				k = k * 10 + ((*pos) - '0');
				++pos;
			} while (('0' <= *pos) && (*pos <= '9'));
			value = k;
			return true;
		}
	}
	return false;
}

int OpenCL::SetKernelSize(int size)
{
	if (size < 16)
	{
		size = 16;
	}
	if (size > 23)
	{
		size = 23;
	}

	myPhase1MaxKernelSize = 1U << size;
	myPhase2MinNumbersCount = 1U << size;
	myPhase2MaxNumbersCount = myPhase2MinNumbersCount + myPhase1MaxKernelSize;

	return size;
}

static cl_mem clCreateBufferLogged(int line, cl_context context, cl_mem_flags flags, size_t size, void* host_ptr, cl_int* errcode_ret)
{
	static size_t total_allocated = 0;

	total_allocated += size;
	LOG(1, "Allocating " << size << " bytes on GPU (" __FILE__ ", line " << line << "), " << ((total_allocated + (1 << 20) - 1) >> 20) << " MB total");

	return clCreateBuffer(context, flags, size, host_ptr, errcode_ret);
}

bool OpenCL::Run(int argc, char* argv[], char* startFrom, char* stopAt, unsigned int largestPrimePower, num64 startPrime, num64 primeLimit)
{
	cl_platform_id platform = nullptr;
	cl_device_id device = nullptr;

	if (!boinc_is_standalone())
	{
		int retval = -1;
		for (int type = 0; retval && (type <= 3); ++type)
		{
			retval = boinc_get_opencl_ids(argc, argv, type, &device, &platform);
		}
		if (retval)
		{
			LOG_ERROR("Error: boinc_get_opencl_ids() failed with error " << retval);
			return false;
		}
	}
	else
	{
		if (!GetPlatformID(&platform))
		{
			return false;
		}

		if (clGetDeviceIDs(platform, DEVICE_TYPE_TO_USE, 1, &device, nullptr) != CL_SUCCESS)
		{
			// No GPU: use any OpenCL device (for example a CPU OpenCL runtime for testing)
			CL_CHECKED_CALL(clGetDeviceIDs, platform, CL_DEVICE_TYPE_ALL, 1, &device, nullptr);
		}
	}

	if (!(platform && device))
	{
		return false;
	}

	std::string vendor_specific_compiler_options;
	num64 MaxMemAllocSize = 0;
	num64 deviceLocalMemSize = 0;
	num64 deviceMaxWorkGroupSize = 0;
	{
		char platformName[1024] = {};
		CL_CHECKED_CALL(clGetPlatformInfo, platform, CL_PLATFORM_NAME, sizeof(platformName), &platformName, nullptr);
		if (!PlatformSupported(platformName))
		{
			LOG_ERROR("Platform '" << platformName << "' is not supported by this program");
			return false;
		}

		num64 global_memory_size = 0;
		num64 global_memory_cache_size = 0;
		num64 local_memory_size = 0;
		unsigned int constant_args = 0;
		num64 constant_buffer_size = 0;
		unsigned int max_compute_units = 0;
		unsigned int max_clock_frequency = 0;
		num64 max_work_group_size = 0;
		char deviceName[1024] = {};
		CL_CHECKED_CALL(clGetDeviceInfo, device, CL_DEVICE_NAME, sizeof(deviceName), deviceName, nullptr);
		CL_CHECKED_CALL(clGetDeviceInfo, device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(global_memory_size), &global_memory_size, nullptr);
		CL_CHECKED_CALL(clGetDeviceInfo, device, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, sizeof(global_memory_cache_size), &global_memory_cache_size, nullptr);
		CL_CHECKED_CALL(clGetDeviceInfo, device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(local_memory_size), &local_memory_size, nullptr);
		CL_CHECKED_CALL(clGetDeviceInfo, device, CL_DEVICE_MAX_CONSTANT_ARGS, sizeof(constant_args), &constant_args, nullptr);
		CL_CHECKED_CALL(clGetDeviceInfo, device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(constant_buffer_size), &constant_buffer_size, nullptr);
		CL_CHECKED_CALL(clGetDeviceInfo, device, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(MaxMemAllocSize), &MaxMemAllocSize, nullptr);
		CL_CHECKED_CALL(clGetDeviceInfo, device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(max_compute_units), &max_compute_units, nullptr);
		CL_CHECKED_CALL(clGetDeviceInfo, device, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(max_compute_units), &max_clock_frequency, nullptr);
		CL_CHECKED_CALL(clGetDeviceInfo, device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(max_work_group_size), &max_work_group_size, nullptr);

		size_t extensions_size = 0;
		CL_CHECKED_CALL(clGetDeviceInfo, device, CL_DEVICE_EXTENSIONS, 0, nullptr, &extensions_size);

		std::vector<char> extensions(extensions_size);
		CL_CHECKED_CALL(clGetDeviceInfo, device, CL_DEVICE_EXTENSIONS, extensions.size(), extensions.data(), nullptr);

		global_memory_size >>= 20;
		deviceLocalMemSize = local_memory_size;
		deviceMaxWorkGroupSize = max_work_group_size;

		LOG(2,
			"Device " << 0 << ": " << deviceName << std::endl <<
			"\t" << max_compute_units << " compute units" << std::endl <<
			"\t" << max_clock_frequency << " MHz" << std::endl <<
			"\t" << global_memory_size << " MB main memory" << std::endl <<
			"\t" << (global_memory_cache_size >> 10) << " KB cache" << std::endl <<
			"\t" << (local_memory_size >> 10) << " KB local memory" << std::endl <<
			"\t" << constant_args << " constant args max" << std::endl <<
			"\t" << (constant_buffer_size >> 10) << " KB constant memory" << std::endl <<
			"\t" << (MaxMemAllocSize >> 20) << " MB max allocation size" << std::endl <<
			"\t" << max_work_group_size << " max work group size" << std::endl <<
			"\t" << extensions.data() << std::endl
		);

		if (MaxMemAllocSize < (1 << 28))
		{
			LOG_ERROR("CL_DEVICE_MAX_MEM_ALLOC_SIZE is less than 256 MB, this GPU/driver is not fit to run this program");
			return false;
		}

		// NVIDIA: kernel size 22 by default, a phase 2 buffer of 2^23 numbers. Its overflows made the sieve run launches again (3^10 unit
		// 552 of the 3^14..3^10 search: 1.37 x its positions, 87.9 s with 21, 1.24 x and 80.9 s with 22). It takes ~300 MB more GPU memory
		// per sieve pipeline (RTX 5070 with one pipeline: 1.47 GB with 21, 1.78 GB with 22), so devices with less than 3 GB (3 GB cards
		// report 3072 MB) keep 21. <kernel_size_nvidia> in the preferences overrides it.
		if (strstr(platformName, "NVIDIA") && (global_memory_size >= 3072))
		{
			SetKernelSize(22);
		}

		if (myPreferences)
		{
			LOG_ERROR("Preferences:\n" << myPreferences << "\n");

			int kernel_size_amd;
			if ((strstr(platformName, "AMD") || strstr(platformName, "ATI")) && ParseIntegerFromXml(myPreferences, "<kernel_size_amd>", kernel_size_amd))
			{
				kernel_size_amd = SetKernelSize(kernel_size_amd);
				LOG_ERROR("Kernel size for AMD GPU has been set to " << kernel_size_amd);
			}

			int kernel_size_nvidia;
			if (strstr(platformName, "NVIDIA") && ParseIntegerFromXml(myPreferences, "<kernel_size_nvidia>", kernel_size_nvidia))
			{
				kernel_size_nvidia = SetKernelSize(kernel_size_nvidia);
				LOG_ERROR("Kernel size for NVIDIA GPU has been set to " << kernel_size_nvidia);
			}
		}

		while (myWorkGroupSize > max_work_group_size)
		{
			myWorkGroupSize >>= 1;
		}

		if (strstr(extensions.data(), "cl_nv_compiler_options"))
		{
			// NVIDIA
			vendor_specific_compiler_options += "-cl-nv-verbose ";
		}
	}

	myLargestPrimePower = largestPrimePower;

	cl_context_properties cps[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties) platform, 0 };
	CL_CHECKED_CALL_WITH_RESULT(myGPUContext, clCreateContext, cps, 1, &device, nullptr, nullptr);

	const char* kernel_sources[2] = { kernel_cl, kernel_sieve_cl };
	CL_CHECKED_CALL_WITH_RESULT(myProgram, clCreateProgramWithSource, myGPUContext, 2, kernel_sources, nullptr);

	num64 NumMainBufferChunks = 0;
	unsigned int ChunkSizeShift = 0;

	const bool IsLargePrimes = (startPrime && primeLimit);
	myMainBufferData = IsLargePrimes ? reinterpret_cast<const unsigned char*>(CandidatesData.data()) : reinterpret_cast<const unsigned char*>(PrimesCompact);
	myMainBufferTotalSize = IsLargePrimes ? static_cast<num64>((CandidatesData.size() + 1) * sizeof(AmicableCandidate)) : PrimesCompactAllocationSize;

	enum { MaxMainBufferSize = 512 << 20 };

	// The sieve is used for all ranges when the largest prime has exponent 1 (N(p) is linear in p only then), and for large primes
	myUseSieve = (IsLargePrimes || (largestPrimePower == 1)) && (g_UseSieve || g_RunTests);
	myLargePrimesSieve = myUseSieve && IsLargePrimes;
	mySieveChunkOffsetArg = myLargePrimesSieve ? 6 : 2;
	std::string sieve_compiler_options;
	if (myUseSieve)
	{
		if (!InitSieve(device, std::min<num64>(deviceLocalMemSize, 65536), deviceMaxWorkGroupSize, sieve_compiler_options))
		{
			return false;
		}
		if (myLargePrimesSieve)
		{
			sieve_compiler_options += " -D SIEVE_LARGE_PRIMES";
		}
		if (mySieve.cfg.wgSize > deviceMaxWorkGroupSize)
		{
			LOG_ERROR("Sieve work-group size " << mySieve.cfg.wgSize << " is not supported by this device (" << deviceMaxWorkGroupSize << " max)");
			return false;
		}

		// The sieve doesn't need prime numbers on the GPU
		myMainBufferData = nullptr;
		myMainBufferTotalSize = 0;
		myMainBufferSize = 0;
		NumMainBufferChunks = 1;
		ChunkSizeShift = 20;
	}
	else if (IsLargePrimes) {
		myMainBufferSize = myMainBufferTotalSize;
		constexpr num64 N = 1 << 20;
		myMainBufferSize = ((myMainBufferSize + N - 1) / N) * N;
	}
	else {
		myMainBufferSize = std::min<num64>(myMainBufferTotalSize, MaxMainBufferSize);
	}

	// Try to create a single continuous buffer for storing prime numbers
	// Even though its size is bigger than CL_DEVICE_MAX_MEM_ALLOC_SIZE for many devices,
	// they are still able to do big allocations sometimes
	if (!myUseSieve)
	{
		cl_int ciErrNum;
		cl_mem buf = clCreateBufferLogged(__LINE__, myGPUContext, CL_MEM_READ_ONLY, myMainBufferSize, 0, &ciErrNum);
		if (ciErrNum == CL_SUCCESS)
		{
			NumMainBufferChunks = 1;
			unsigned long index;
			_BitScanReverse64(&index, myMainBufferSize);
			ChunkSizeShift = index + 1;
			myMainBuffers.emplace_back(buf);
		}
		else
		{
			// If it failed, take into account maximum allocation size for this device
			unsigned long index;
			_BitScanReverse64(&index, MaxMemAllocSize);
			ChunkSizeShift = index;
			NumMainBufferChunks = (myMainBufferSize + (1ULL << ChunkSizeShift) - 1) >> ChunkSizeShift;

			myMainBuffers.reserve(NumMainBufferChunks);
			for (num64 totalSize = 0; totalSize < myMainBufferSize;)
			{
				const num64 BufSize = std::min<num64>(1ULL << ChunkSizeShift, myMainBufferSize - totalSize);
				cl_mem buf;
				CL_CHECKED_CALL_WITH_RESULT(buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, BufSize, 0);
				myMainBuffers.emplace_back(buf);
				totalSize += BufSize;
			}
		}
	}

	{
		char kernel_options[2048];
		char cBuildLog[10240];
		sprintf_s(kernel_options,
			"%s-cl-fast-relaxed-math -D PQ_STRIDE_SIZE=%u -D WORK_GROUP_SIZE=%u -D PHASE2_MAX_COUNT=%u -D LPP=%u -D NUM_DATA_CHUNKS=%u -D CHUNK_SIZE_SHIFT=%u -D SUM_ESTIMATES_128_SHIFT=%u -D MAIN_BUFFER_MASK=%u%s -Werror",
			vendor_specific_compiler_options.c_str(),
			SumEstimatesSize2_GPU,
			static_cast<unsigned int>(myWorkGroupSize),
			myPhase2MaxNumbersCount,
			myLargestPrimePower,
			static_cast<unsigned int>(NumMainBufferChunks),
			ChunkSizeShift,
			64 - SumEstimates128Shift,
			(MaxMainBufferSize / 2) - 1,
			sieve_compiler_options.c_str()
		);
		// Development aid: extra kernel build options (for example -D switches of experimental kernel variants)
		if (const char* extra = getenv("AMICABLE_CL_OPTIONS"))
		{
			strcat_s(kernel_options, " ");
			strcat_s(kernel_options, extra);
		}
		myBuildOptions = kernel_options;
		cl_int build_result = clBuildProgram(myProgram, 0, nullptr, kernel_options, nullptr, nullptr);
		if (build_result != CL_SUCCESS)
		{
			LOG_ERROR("clBuildProgram returned error " << build_result);
			LOG_ERROR("Build options: " << kernel_options);

			const cl_int get_info_result = clGetProgramBuildInfo(myProgram, device, CL_PROGRAM_BUILD_LOG, sizeof(cBuildLog), cBuildLog, nullptr);
			if (get_info_result == CL_SUCCESS)
			{
				LOG_ERROR("Build log:\n" << cBuildLog);
			}
			else
			{
				LOG_ERROR("Failed to get build log, clGetProgramBuildInfo returned error " << get_info_result);
			}

			LOG_ERROR("Trying to disable 'goto' and build again");

			strcat_s(kernel_options, " -D DISABLE_GOTO");
			myBuildOptions = kernel_options;
			build_result = clBuildProgram(myProgram, 0, nullptr, kernel_options, nullptr, nullptr);
		}

		const cl_int get_info_result2 = clGetProgramBuildInfo(myProgram, device, CL_PROGRAM_BUILD_LOG, sizeof(cBuildLog), cBuildLog, nullptr);
		if (get_info_result2 != CL_SUCCESS)
		{
			LOG_ERROR("Failed to get build log, clGetProgramBuildInfo returned error " << get_info_result2);
			cBuildLog[0] = '\0';
		}

		LOG(2, cBuildLog);

			if (build_result != CL_SUCCESS)
			{
				LOG_ERROR("clBuildProgram returned error " << build_result);
				LOG_ERROR("Build options: " << kernel_options);
				LOG_ERROR("Build log:\n" << cBuildLog);
				return false;
			}

			// the sieve with the divisibility kill (chunks of the last tier)
			if (!mySieveKillOptions.empty())
			{
				std::string killOptions = kernel_options;
				const size_t pos = killOptions.find(mySieveOptions);
				if (pos == std::string::npos)
				{
					LOG_ERROR("Sieve build options not found");
					return false;
				}
				killOptions.replace(pos, mySieveOptions.size(), mySieveKillOptions);
				CL_CHECKED_CALL_WITH_RESULT(myProgramKill, clCreateProgramWithSource, myGPUContext, 2, kernel_sources, nullptr);
				const cl_int killBuild = clBuildProgram(myProgramKill, 0, nullptr, killOptions.c_str(), nullptr, nullptr);
				if (killBuild != CL_SUCCESS)
				{
					clGetProgramBuildInfo(myProgramKill, device, CL_PROGRAM_BUILD_LOG, sizeof(cBuildLog), cBuildLog, nullptr);
					LOG_ERROR("clBuildProgram (sieve with the kill) returned error " << killBuild << "\nBuild options: " << killOptions << "\nBuild log:\n" << cBuildLog);
					return false;
				}
				myBuildOptionsKill = killOptions;
			}
		}

	CL_CHECKED_CALL_WITH_RESULT(mySaveCounter, clCreateKernel, myProgram, "SaveCounter");

	if (IsLargePrimes && !myUseSieve)
	{
		CL_CHECKED_CALL_WITH_RESULT(mySearchLargePrimes, clCreateKernel, myProgram, "SearchLargePrimes");
	}
	else if (!myUseSieve)
	{
		CL_CHECKED_CALL_WITH_RESULT(mySearchMultipleRanges, clCreateKernel, myProgram, "SearchMultipleRanges");
	}

	CL_CHECKED_CALL_WITH_RESULT(myCheckPairs, clCreateKernel, myProgram, "CheckPairs");
	// (the sieve kernels and the kernels of phases 2 and 3 belong to a pipeline: pipeline 0 here, the second one in CreateSecondPipe)
	if (!CreatePipeKernels())
	{
		return false;
	}

	CL_CHECKED_CALL_WITH_RESULT(P().myQueue, clCreateCommandQueue, myGPUContext, device, g_SieveProfile ? CL_QUEUE_PROFILING_ENABLE : 0);

	{
		std::vector<unsigned int> smallPrimes;
		smallPrimes.reserve(ReciprocalsTableSize128);
		for (unsigned int i = 0; i < ReciprocalsTableSize128; ++i)
		{
			smallPrimes.emplace_back(static_cast<unsigned int>(GetNthPrime(i)));
		}
		CL_CHECKED_CALL_WITH_RESULT(mySmallPrimesBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, ReciprocalsTableSize128 * sizeof(unsigned int), 0);
		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, mySmallPrimesBuf, CL_TRUE, 0, ReciprocalsTableSize128 * sizeof(unsigned int), smallPrimes.data(), 0, nullptr, nullptr);
		mySmallPrimesHost = smallPrimes;
	}

	const size_t primeInversesBufSize = sizeof(std::pair<num64, num64>) * ReciprocalsTableSize128;
	CL_CHECKED_CALL_WITH_RESULT(myPrimeInversesBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, primeInversesBufSize, 0);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, myPrimeInversesBuf, CL_TRUE, 0, primeInversesBufSize, PrimeInverses, 0, nullptr, nullptr);

	{
		const size_t primeReciprocalsBufSize = sizeof(std::pair<num64, num64>) * ReciprocalsTableSize128;
		CL_CHECKED_CALL_WITH_RESULT(myPrimeReciprocalsBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, primeReciprocalsBufSize, 0);

		std::vector<std::pair<num64, num64>> r;
		r.resize(ReciprocalsTableSize128);
		for (unsigned int i = 0; i < ReciprocalsTableSize128; ++i)
		{
			r[i].first = PrimeReciprocals[i].reciprocal;
			r[i].second = (GetNthPrime(i) << 32) + (num64(PrimeReciprocals[i].increment) << 8) + num64(PrimeReciprocals[i].shift);
		}
		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, myPrimeReciprocalsBuf, CL_TRUE, 0, primeReciprocalsBufSize, r.data(), 0, nullptr, nullptr);
	}

	{
		std::vector<std::pair<num64, num64>> pq;
		pq.resize(SumEstimatesSize + SumEstimatesSize * SumEstimatesSize2_GPU);

		for (unsigned int i = 0; i < SumEstimatesSize; ++i)
		{
			pq[i].first = SumEstimatesBeginP[i];
			pq[i].second = SumEstimatesBeginQ[i];
		}

		for (unsigned int i = 0; i < SumEstimatesSize; ++i)
		{
			for (unsigned int j = 0; j < SumEstimatesSize2_GPU; ++j)
			{
				pq[i * SumEstimatesSize2_GPU + j + SumEstimatesSize] = PQ[i][j * (SumEstimatesSize2 / SumEstimatesSize2_GPU)];
			}
		}

		// guard element which will stop the "while (N <= PQ[k].x) --k;" loop at k = 0
		pq[0].first = 0;

		const size_t PQ_BufSize = sizeof(std::pair<num64, num64>) * pq.size();
		CL_CHECKED_CALL_WITH_RESULT(myPQ_Buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, PQ_BufSize, 0);
		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, myPQ_Buf, CL_TRUE, 0, PQ_BufSize, pq.data(), 0, nullptr, nullptr);
	}

	{
		static num64 powerOf2Inverses[128][2];
		for (num64 i = 0; i < 64; ++i)
		{
			const num64 value = (i < 63) ? ((num64(1) << (i + 1)) - 1) : num64(-1);
			PRAGMA_WARNING(suppress : 4146)
			powerOf2Inverses[i][0] = -modular_inverse64(value);
			powerOf2Inverses[i][1] = num64(-1) / value;
		}
		for (num64 i = 0; i < 64; ++i)
		{
			powerOf2Inverses[i + 64][0] = LowWord(PowersOf2_128DivisibilityData[i]);
			powerOf2Inverses[i + 64][1] = HighWord(PowersOf2_128DivisibilityData[i]);
		}
		CL_CHECKED_CALL_WITH_RESULT(myPowerOf2SumInverses_Buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, sizeof(powerOf2Inverses), 0);
		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, myPowerOf2SumInverses_Buf, CL_TRUE, 0, sizeof(powerOf2Inverses), powerOf2Inverses, 0, nullptr, nullptr);
	}

	{
		std::vector<unsigned int> offsets;
		offsets.reserve(ReciprocalsTableSize128);
		for (unsigned int i = 0; i < ReciprocalsTableSize128; ++i)
		{
			offsets.emplace_back(static_cast<unsigned int>(PowersOfP_128DivisibilityData[i] - PowersOfP_128DivisibilityData_base));
		}

		CL_CHECKED_CALL_WITH_RESULT(myPowersOfP_128SumInverses_offsets_Buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, sizeof(unsigned int) * ReciprocalsTableSize128, 0);
		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, myPowersOfP_128SumInverses_offsets_Buf, CL_TRUE, 0, sizeof(unsigned int) * ReciprocalsTableSize128, offsets.data(), 0, nullptr, nullptr);

		std::vector<num64> inverses;
		inverses.resize(PowersOfP_128DivisibilityData_count * 4);
		for (unsigned int i = 0; i < PowersOfP_128DivisibilityData_count; ++i)
		{
			inverses[i * 4 + 0] = LowWord(PowersOfP_128DivisibilityData_base[i].inverse);
			inverses[i * 4 + 1] = HighWord(PowersOfP_128DivisibilityData_base[i].inverse);
			inverses[i * 4 + 2] = PowersOfP_128DivisibilityData_base[i].shift;
			inverses[i * 4 + 3] = PowersOfP_128DivisibilityData_base[i].shift_bits;
		}

		CL_CHECKED_CALL_WITH_RESULT(myPowersOfP_128SumInverses_Buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, sizeof(num64) * 4 * PowersOfP_128DivisibilityData_count, 0);
		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, myPowersOfP_128SumInverses_Buf, CL_TRUE, 0, sizeof(num64) * 4 * PowersOfP_128DivisibilityData_count, inverses.data(), 0, nullptr, nullptr);
	}

	CL_CHECKED_CALL_WITH_RESULT(myPrimeInverses_128_Buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, sizeof(num128) * ReciprocalsTableSize128, 0);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, myPrimeInverses_128_Buf, CL_TRUE, 0, sizeof(num128) * ReciprocalsTableSize128, PrimeInverses128, 0, nullptr, nullptr);

	CL_CHECKED_CALL_WITH_RESULT(mySumEstimates_128_Buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, sizeof(num64) * ReciprocalsTableSize128 / 16, 0);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, mySumEstimates_128_Buf, CL_TRUE, 0, sizeof(num64) * ReciprocalsTableSize128 / 16, SumEstimates128, 0, nullptr, nullptr);

	CL_CHECKED_CALL_WITH_RESULT(myRangesTable_Buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, sizeof(RangeDataGPU) * myMaxRangesCount, 0);
	CL_CHECKED_CALL_WITH_RESULT(myRangesLookupTable_Buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, sizeof(LookupDataGPU) * myMaxLookupTableSize, 0);

	ourZeroBuf = new unsigned char[ourZeroBufSize];
	memset(ourZeroBuf, 0, ourZeroBufSize);

	if (!CreatePipeBuffers())
	{
		return false;
	}

	if (myUseSieve && (!CreateSieveTables() || !CreateSievePipeBuffers()))
	{
		return false;
	}

	CL_CHECKED_CALL(clFinish, P().myQueue);

	if (Test() && (!myUseSieve || TestSieveStartup()))
	{
		LOG_ERROR("Self-test passed");
	}
	else
	{
		LOG_ERROR("Self-test failed");
		return false;
	}

	CL_CHECKED_CALL(setKernelArguments, mySaveCounter, P().myPhase2_numbers_count_buf);

	if (IsLargePrimes && !myUseSieve)
	{
		for (num64 totalSize = 0, i = 0; totalSize < myMainBufferSize; ++i)
		{
			const num64 BufSize = std::min<num64>(1ULL << ChunkSizeShift, myMainBufferSize - totalSize);
			CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, myMainBuffers[i], CL_TRUE, 0, BufSize, myMainBufferData + totalSize, 0, nullptr, nullptr);
			totalSize += BufSize;
		}

		CL_CHECKED_CALL_WITH_RESULT(myLargePrimesBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, sizeof(num64) * myLargePrimesMaxCount, 0);

		CL_CHECKED_CALL(setKernelArguments, mySearchLargePrimes,
			mySmallPrimesBuf,
			myPrimeInversesBuf,
			myPQ_Buf,
			myPowerOf2SumInverses_Buf,
			myPowersOfP_128SumInverses_offsets_Buf,
			myPowersOfP_128SumInverses_Buf,
			myPrimeInverses_128_Buf,
			mySumEstimates_128_Buf,
			myLargePrimesBuf,
			0U,
			0ULL,
			0U,
			0ULL,
			0ULL,
			P().myPhase1_offset_to_resume_buf,
			P().myPhase2_numbers_count_buf,
			P().myPhase2_numbers_buf,
			P().myAmicable_numbers_count_buf,
			P().myAmicable_numbers_data_buf
		);

		for (num64 i = 0; i < myMainBuffers.size(); ++i)
		{
			CL_CHECKED_CALL(clSetKernelArg, mySearchLargePrimes, static_cast<cl_uint>(i + 19), sizeof(cl_mem), &myMainBuffers[i]);
		}
	}
	else if (!myUseSieve)
	{
		CL_CHECKED_CALL(setKernelArguments, mySearchMultipleRanges,
			mySmallPrimesBuf,
			myPrimeInversesBuf,
			myPQ_Buf,
			myPowerOf2SumInverses_Buf,
			myPowersOfP_128SumInverses_offsets_Buf,
			myPowersOfP_128SumInverses_Buf,
			myPrimeInverses_128_Buf,
			mySumEstimates_128_Buf,
			myRangesTable_Buf,
			myRangesLookupTable_Buf,
			0U,
			0ULL,
			0U,
			P().myPhase1_offset_to_resume_buf,
			P().myPhase2_numbers_count_buf,
			P().myPhase2_numbers_buf,
			P().myAmicable_numbers_count_buf,
			P().myAmicable_numbers_data_buf
		);

		for (num64 i = 0; i < myMainBuffers.size(); ++i)
		{
			CL_CHECKED_CALL(clSetKernelArg, mySearchMultipleRanges, static_cast<cl_uint>(i + 18), sizeof(cl_mem), &myMainBuffers[i]);
		}
	}

	if (!SetPhase2KernelArguments())
	{
		return false;
	}

	if (g_RunTests)
	{
		return RunSieveTests();
	}

	Timer t;
	bool result;
	if (IsLargePrimes)
	{
		result = myUseSieve ? RunLargePrimesSieve(startPrime, primeLimit) : RunLargePrimes(startPrime, primeLimit);
	}
	else
	{
		result = true;
		if (myUseSieve)
		{
			// several batches on the GPU at once (SievePipe): pipeline i is created if the device has the memory for the buffers of i + 1
			// pipelines (~1.1 GB each with kernel size 21, ~1.4 GB with 22) and 1 GB for the rest
			cl_ulong memorySize = 0;
			clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(memorySize), &memorySize, nullptr);
			const cl_ulong pipeBytes = sizeof(num64) * 2 * cl_ulong(mySurvivorCapacity) + sizeof(num64) * 4 * cl_ulong(myPhase2MaxNumbersCount) * 2 +
				sizeof(cl_uint) * cl_ulong(myPhase2MaxNumbersCount) + sizeof(num64) * 4 * cl_ulong(Phase3BacklogCapacity()) +
				cl_ulong(SieveMaxBandChunks) * mySieve.cfg.chunkSegments * mySieve.cfg.segJ * 7 * 8 * sizeof(cl_uint) +
				sizeof(cl_uint) * SieveBandEntry * SieveMaxBandChunks * SieveMaxLaunchesPerRound + sizeof(num64) * SieveRangeStride * cl_ulong(mySieveMaxRanges) +
				sizeof(unsigned int) * 4 * cl_ulong(mySieveMaxChunks);
			for (int i = 1; (i < g_SievePipes) && (i < MaxSievePipes) && (memorySize >= pipeBytes * (i + 1) + (cl_ulong(1) << 30)); ++i)
			{
				if (!CreateExtraPipe(device, i))
				{
					// (what was created is released in the destructor)
					LOG_ERROR("Could not create sieve pipeline " << i << ", running with " << myNumPipes);
					break;
				}
			}
			num64 offset = 0;
			result = RunRanges(startFrom, stopAt, offset, ChunkSizeShift);
		}
		else
		{
			for (num64 offset = 0; offset < myMainBufferTotalSize; offset += myMainBufferSize)
			{
				if (!RunRanges(startFrom, stopAt, offset, ChunkSizeShift))
				{
					result = false;
					break;
				}
			}
		}
	}
	LOG(1, "\nWork unit finished in " << std::setprecision(6) << t.getElapsedTime() << " seconds, " << myNumbersProcessedTotal << " numbers processed\n");
	// (the statistics of all pipelines in pipeline 0)
	for (int i = 1; i < myNumPipes; ++i)
	{
		AddPipeStats(myPipes[i]);
	}
	if (myUseSieve)
	{
		std::cerr << "Sieve: " << mySieveStatRanges << " ranges, " << P().mySieveStatChunks << " chunks, " << P().mySieveStatPositions << " positions, " << P().mySieveStatSurvivors << " survivors (" << ((P().mySieveStatPositions > 0) ? (P().mySieveStatSurvivors * 100.0 / P().mySieveStatPositions) : 0.0) << "%), " << P().mySieveStatBatches << " batches, " << P().mySieveStatTime << " seconds in batches, " << P().mySieveStatPhase2 << " numbers in phase 2, " << P().mySieveStatOverflows << " overflows" << std::endl;
		std::cerr << "Host cycles (1e9): range preparation " << g_HostCycles[0] * 1e-9 << " (" << g_SieveThreads << " helper threads; reduction " << g_HostCycles[1] * 1e-9 << ", PrepareRange " << g_HostCycles[2] * 1e-9 << ", EstimatePrimes " << g_HostCycles[4] * 1e-9 << "), main thread: RangeGen::Iterate " << g_HostCycles[3] * 1e-9 << ", CommitSieveRange " << g_HostCycles[6] * 1e-9 << " (AddSieveRangeRecord " << g_HostCycles[5] * 1e-9 << "), waiting for the helpers " << g_HostCycles[7] * 1e-9 << std::endl;
		std::cerr << "Sieve launches: " << P().mySieveStatLaunches << " launches in " << P().mySieveStatRounds << " rounds, " << P().mySieveStatLaunchedPositions << " positions sieved (" << ((P().mySieveStatPositions > 0) ? (static_cast<double>(P().mySieveStatLaunchedPositions) / P().mySieveStatPositions) : 0.0) << " x)" << std::endl;
		if (g_SieveProfile)
		{
			std::cerr << "Sieve profile: " << P().mySieveStatSieveKernelTime << " seconds in SieveRanges, " << P().mySieveStatCheckKernelTime << " seconds in CheckSieveSurvivors, " << P().mySieveStatBandKernelTime << " seconds in SieveBand (" << P().mySieveStatBandChunks << " chunks), " << P().mySieveStatPhase2Sorts << " sorted phase 2 runs" << std::endl;
			std::cerr << "Phase 2 profile: " << P().mySieveStatPhase2Wall << " seconds (host), GPU: prefilter " << P().mySieveStatPhase2Times[0] << ", bucket count " << P().mySieveStatPhase2Times[1] << ", scatter and copy " << P().mySieveStatPhase2Times[2] << ", CheckPairPhase2 " << P().mySieveStatPhase2Times[3] << std::endl;
			std::cerr << "Survivor checks: first pass " << P().mySieveStatCheckTimes[0] << ", SieveTestDeferred " << P().mySieveStatCheckTimes[1] << ", prime p " << P().mySieveStatCheckTimes[2] << " seconds" << std::endl;
			std::cerr << "Phase 3 profile: " << P().mySieveStatPhase3Wall << " seconds (host), " << P().mySieveStatPhase3Numbers << " numbers; per iteration (seconds, numbers, runs):";
			for (int k = 1; k < 20; ++k)
			{
				if (P().mySieveStatPhase3IterLaunches[k])
				{
					std::cerr << " " << k << ": " << P().mySieveStatPhase3IterTime[k] << " " << P().mySieveStatPhase3IterNumbers[k] << " " << P().mySieveStatPhase3IterLaunches[k] << ";";
				}
			}
			std::cerr << std::endl;
		}
	}
	if (myUseSieve && !myLargePrimesSieve)
	{
		std::cerr << "Reduction: " << myRedStat[0] << " positions, " << myRedStat[1] << " in impossible ranges (" << (myRedStat[1] * 100.0 / myRedStat[0]) << "%), "
			<< myRedStat[2] << " after p == -1 (mod L) (" << (myRedStat[2] * 100.0 / myRedStat[0]) << "%), " << myRedStat[3] << " with the 2, 3, 5 parts (" << (myRedStat[3] * 100.0 / myRedStat[0]) << "%)" << std::endl;
		if (g_AnalyzeRate > 0)
		{
			std::cerr << "Analysis: " << myAnalyzeSamples << " positions sampled" << std::endl;
			for (unsigned int t = 0; t < mySieve.cfg.numTiers; ++t)
			{
				std::cerr << "  Q=" << mySieve.cfg.tierQ[t] << ": window " << (myAnalyzeCounts[t][0] * 100.0 / myAnalyzeSamples) << "%, full divisibility " << (myAnalyzeCounts[t][1] * 100.0 / myAnalyzeSamples)
					<< "%, e=1 255<q<=16381 " << (myAnalyzeCounts[t][2] * 100.0 / myAnalyzeSamples) << "%, +2-adic " << (myAnalyzeCounts[t][3] * 100.0 / myAnalyzeSamples)
					<< "%, e=1 q<=Q " << (myAnalyzeCounts[t][4] * 100.0 / myAnalyzeSamples) << "%, +2-adic " << (myAnalyzeCounts[t][5] * 100.0 / myAnalyzeSamples) << "%" << std::endl;
			}
		}
		std::cerr << "Positions by W (window height above the known part) <1e-4, 1e-3, 1e-2, 0.03, 0.1, 0.3, 1, more:";
		for (int b = 0; b < 8; ++b)
		{
			std::cerr << " " << (myRedHistW[b] * 100.0 / myRedStat[0]) << "%";
		}
		std::cerr << std::endl << mySieveStatSkipped << " ranges skipped (no amicable numbers possible)" << std::endl;
		std::cerr << "Positions by bit length of L:";
		for (int b = 0; b < 64; ++b)
		{
			if (myRedHist[b] > 0)
			{
				std::cerr << " " << b << ":" << (myRedHist[b] * 100.0 / myRedStat[0]) << "%";
			}
		}
		std::cerr << std::endl;
	}
	std::cerr << "Work unit finished in " << std::setprecision(6) << t.getElapsedTime() << " seconds" << std::endl;
	return result;
}

bool OpenCL::RunRanges(char* startFrom, char* stopAt, num64& offset, unsigned int ChunkSizeShift)
{
	RangeData r;
	memset(&r, 0, sizeof(r));

	Factor stopAtFactors[MaxPrimeFactors + 1];
	memset(stopAtFactors, 0, sizeof(stopAtFactors));

	// Few numbers, hard case for CPU
	// /from 2*5*7*11*23*41*2083 /to 2*5*7*11*23*41*2137 /task_size 851654026
	//char start[] = "2*5*7*11*23*41*2083";
	//char stop[]  = "2*5*7*11*23*41*2137";

	// A lot of numbers, easy case for CPU
	//char start[] = "2*5*7*11*43*31657";
	//char stop[]  = "2*5*7*11*43*95261";

	// Very easy case for CPU
	//char start[] = "2*5*7*11*89^2*840661";
	//char stop[]  = "2*5*7*11*97*107*317*3347";

	// Some random work unit, 3021.46 seconds on Intel Xeon E5-1650 v3 (6 cores, 12 threads)
	// 906509621271 numbers to check
	//char start[] = "2*5*7*210138413";
	//char stop[]  = "2*5*7*210493271";

	// Work unit with 15 amicable pairs
	// /from 2*5*7^2*37*41*53*4591 /to 2*5*7^2*37*41*1427*9787 /task_size 522154890671
	//char start[] = "2*5*7^2*37*41*53*4591";
	//char stop[]  = "2*5*7^2*37*41*1427*9787";

	// Work unit with 75% GPU load
	// /from 2*5*11*13*17*919*10861 /to 2*5*11*13*17*947*977*1973 /task_size 109093979525 

	char checkpoint_buf[256];
	std::string checkpoint_name;
	boinc_resolve_filename_s(CHECKPOINT_LOGICAL_NAME, checkpoint_name);
	FILE* checkpoint = boinc_fopen(checkpoint_name.c_str(), "r");
	if (checkpoint)
	{
		if (fgets(checkpoint_buf, sizeof(checkpoint_buf), checkpoint))
		{
			char* end_of_data = strchr(checkpoint_buf, '#');
			if (end_of_data)
			{
				*end_of_data = '\0';
				char* c = strchr(checkpoint_buf, ' ');
				if (c)
				{
					const num64 saved_NumbersProcessedTotal = StrToNumber(checkpoint_buf);
					do { ++c; } while (*c == ' ');
					const num64 saved_offset = StrToNumber(c);
					if (saved_offset >= offset)
					{
						myNumbersProcessedTotal = saved_NumbersProcessedTotal;
						const double f = myNumbersProcessedTotal / RangeGen::total_numbers_to_check;
						boinc_fraction_done((f < 1.0) ? f : 1.0);

						offset = saved_offset;
						do { ++c; } while (*c != ' ');
						do { ++c; } while (*c == ' ');
						startFrom = c;
					}
				}
			}
		}
		fclose(checkpoint);
	}

	const num64 n = std::min(myMainBufferSize, myMainBufferTotalSize - offset);
	for (num64 totalSize = 0, i = 0; totalSize < n; ++i)
	{
		const num64 BufSize = std::min<num64>(1ULL << ChunkSizeShift, n - totalSize);
		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, myMainBuffers[i], CL_TRUE, 0, BufSize, myMainBufferData + offset + totalSize, 0, nullptr, nullptr);
		totalSize += BufSize;
	}

	LOG(1, "\n---=== RUNNING RANGES WITH OFFSET " << offset << " ===---\n");

	RangeGen::Init(startFrom, stopAt, &r, stopAtFactors, myLargestPrimePower);
	if (!startFrom)
	{
		RangeGen::Iterate(r);
	}

	if (myUseSieve)
	{
		return RunRangesSieve(r, stopAtFactors, offset);
	}
	while (!RangeGen::HasReached(r, stopAtFactors))
	{
		bool range_added;
		if (!AddRange(r, offset, range_added))
		{
			return false;
		}
		if (!RangeGen::Iterate(r) || (RangeGen::cur_largest_prime_power != myLargestPrimePower))
		{
			break;
		}

		const double progress = static_cast<double>(myNumbersProcessedTotal) / RangeGen::total_numbers_to_check;
		boinc_fraction_done((progress < 1.0) ? progress : 1.0);

		if (range_added && (myTotalNumbersInRanges == 0))
		{
			unsigned int numbers_in_phase2 = 0;
			if (!GetCounter(P().myQueue, P().myPhase2_numbers_count_buf, numbers_in_phase2))
			{
				return false;
			}
			if (numbers_in_phase2 == 0)
			{
				if (!SaveProgressForRanges(r, offset))
				{
					return false;
				}
			}
		}
	}

	if ((myTotalNumbersInRanges > 0) && !ProcessNumbers())
	{
		return false;
	}

	unsigned int numbers_in_phase2 = 0;
	if (!GetCounter(P().myQueue, P().myPhase2_numbers_count_buf, numbers_in_phase2))
	{
		return false;
	}

	if ((numbers_in_phase2 > 0) && !ProcessNumbersPhases2_3(numbers_in_phase2))
	{
		return false;
	}
	CleanupRanges();

	return SaveFoundNumbers();
}

bool OpenCL::RunLargePrimes(num64 startPrime, num64 primeLimit)
{
	std::thread gpu_thread(&OpenCL::ProcessLargePrimesThread, this);

	if (!myLargePrimes)
	{
		myLargePrimes = reinterpret_cast<num64*>(AllocateSystemMemory(sizeof(num64) * myLargePrimesMaxCount, false));
	}

	primesieve::PrimeSieve sieve;

	char checkpoint_buf[256];
	std::string checkpoint_name;
	boinc_resolve_filename_s(CHECKPOINT_LOGICAL_NAME, checkpoint_name);
	FILE* checkpoint = boinc_fopen(checkpoint_name.c_str(), "r");
	if (checkpoint)
	{
		if (fgets(checkpoint_buf, sizeof(checkpoint_buf), checkpoint))
		{
			char* end_of_data = strchr(checkpoint_buf, '#');
			if (end_of_data)
			{
				*end_of_data = '\0';
				std::stringstream s;
				num64 firstPrime, lastPrime, offset;
				s << checkpoint_buf;
				s >> myNumbersProcessedTotal >> firstPrime >> lastPrime >> offset;
				if ((startPrime <= firstPrime) && (lastPrime <= primeLimit) && (firstPrime <= lastPrime) && IsPrime(firstPrime) && IsPrime(lastPrime))
				{
					myLargePrimesStartOffset = offset;
					sieve.sieveTemplated(firstPrime, lastPrime, *this);
					if (myLargePrimesCount > 0)
					{
						PassLargePrimesToThread();
					}
					startPrime = lastPrime + 1;
				}
			}
		}
		fclose(checkpoint);
	}

	sieve.sieveTemplated(startPrime, primeLimit, *this);
	if (myLargePrimesCount > 0)
	{
		PassLargePrimesToThread();
	}

	// Send 0 primes to GPU thread to make it exit
	PassLargePrimesToThread();
	gpu_thread.join();

	unsigned int numbers_in_phase2 = 0;
	if (!GetCounter(P().myQueue, P().myPhase2_numbers_count_buf, numbers_in_phase2))
	{
		return false;
	}

	if ((numbers_in_phase2 > 0) && !ProcessNumbersPhases2_3(numbers_in_phase2))
	{
		return false;
	}

	return SaveProgressForLargePrimes(0, 0, 0, true);
}

NOINLINE void OpenCL::PassLargePrimesToThread()
{
	if (!myLargePrimesReady.Signal())
	{
		const int error_code = errno;
		LOG_ERROR("Failed to signal semaphore, errno " << error_code);
		boinc_finish(-1);
	}

	if (!myLargePrimesReceived.Wait())
	{
		const int error_code = errno;
		LOG_ERROR("Failed to wait for semaphore, errno " << error_code);
		boinc_finish(-1);
	}

	myLargePrimesCount = 0;
	myLargePrimesStartOffset = 0;
}

void OpenCL::ProcessLargePrimesThread()
{
	if (!ProcessLargePrimes())
	{
		boinc_finish(-1);
	}
}

bool OpenCL::ProcessLargePrimes()
{
	while (myLargePrimesReady.Wait())
	{
		// If there are no primes to process it means we need to exit
		if (myLargePrimesCount == 0)
		{
			if (!myLargePrimesReceived.Signal())
			{
				const int error_code = errno;
				LOG_ERROR("Failed to signal semaphore, errno " << error_code);
				return false;
			}
			return true;
		}

		const unsigned int LargePrimesCount = myLargePrimesCount;
		const num64 FirstPrime = myLargePrimes[0];
		const num64 LastPrime = myLargePrimes[LargePrimesCount - 1];
		const num64 LargePrimesStartOffset = myLargePrimesStartOffset;
		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, myLargePrimesBuf, CL_TRUE, 0, sizeof(num64) * LargePrimesCount, myLargePrimes, 0, nullptr, nullptr);

		if (!myLargePrimesReceived.Signal())
		{
			const int error_code = errno;
			LOG_ERROR("Failed to signal semaphore, errno " << error_code);
			return false;
		}

		num64 LargePrimesCountReciprocal = 0;
		unsigned int LargePrimesCountIncrementAndShift = 0;

		if (LargePrimesCount & (LargePrimesCount - 1))
		{
			// Not a power of 2
			SReciprocal r;
			r.Init(LargePrimesCount);
			LargePrimesCountReciprocal = r.reciprocal;
			LargePrimesCountIncrementAndShift = static_cast<unsigned int>(r.increment + (r.shift << 1));
		}
		else
		{
			// Power of 2
			unsigned long index;
			_BitScanReverse64(&index, LargePrimesCount);
			LargePrimesCountIncrementAndShift = index;
		}

		const num64 LargestCandidate = LowWord(SearchLimit::value / FirstPrime);
		num64 CandidatesCount;
		{
			int a = 0;
			int b = static_cast<int>(CandidatesData.size());
			while (a < b)
			{
				const int c = (a + b) >> 1;

				if (CandidatesData[c].get_value() <= LargestCandidate)
				{
					a = c + 1;
				}
				else
				{
					b = c;
				}
			}
			CandidatesCount = static_cast<num64>(a);
		}

		const num64 TotalNumbersToCheck = CandidatesCount * LargePrimesCount;

		CL_CHECKED_CALL(clSetKernelArg, mySearchLargePrimes, 9, sizeof(LargePrimesCount), &LargePrimesCount);
		CL_CHECKED_CALL(clSetKernelArg, mySearchLargePrimes, 10, sizeof(LargePrimesCountReciprocal), &LargePrimesCountReciprocal);
		CL_CHECKED_CALL(clSetKernelArg, mySearchLargePrimes, 11, sizeof(LargePrimesCountIncrementAndShift), &LargePrimesCountIncrementAndShift);
		CL_CHECKED_CALL(clSetKernelArg, mySearchLargePrimes, 13, sizeof(TotalNumbersToCheck), &TotalNumbersToCheck);

		cl_event event = nullptr;
		unsigned int numbers_in_phase2_before = 0;
		unsigned int numbers_in_phase2_after = 0;
		num64 phase1_offset_to_resume = 0;
		num64 global_offset = LargePrimesStartOffset;

		while (global_offset < TotalNumbersToCheck)
		{
			if (!GetCounter(P().myQueue, P().myPhase2_numbers_count_buf, numbers_in_phase2_before))
			{
				return false;
			}

			//
			// Phase 1
			//
			{
				Timer t;

				CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().myPhase1_offset_to_resume_buf, CL_TRUE, 0, sizeof(num64), ourZeroBuf, 0, nullptr, nullptr);

				const num64 global_offset_before = global_offset;
				const unsigned int max_size_phase1 = GetMaxPhaseSize(TotalNumbersToCheck - global_offset, myPhase1MaxKernelSize);

				num64 N = global_offset + Phase1Limit;
				if (N > TotalNumbersToCheck)
				{
					N = TotalNumbersToCheck;
				}
				for (unsigned int k = 0; global_offset < N; ++k)
				{
					if (k)
					{
						CL_CHECKED_CALL(clFlush, P().myQueue);
					}
					size_t globalSizePhase1 = std::min<size_t>(N - global_offset, max_size_phase1);
					if (globalSizePhase1 & (myWorkGroupSize - 1))
					{
						globalSizePhase1 += myWorkGroupSize - (globalSizePhase1 & (myWorkGroupSize - 1));
					}
					CL_CHECKED_CALL(clSetKernelArg, mySearchLargePrimes, 12, sizeof(global_offset), &global_offset);
					CL_CHECKED_CALL(clEnqueueTask, P().myQueue, mySaveCounter, 0, nullptr, nullptr);
					global_offset += globalSizePhase1;
					CL_CHECKED_CALL(clEnqueueNDRangeKernel, P().myQueue, mySearchLargePrimes, 1, nullptr, &globalSizePhase1, &myWorkGroupSize, 0, nullptr, (global_offset >= N) ? &event : nullptr);
				}

				if (!WaitForQueue(P().myQueue, event) ||
					!GetCounter(P().myQueue, P().myPhase2_numbers_count_buf, numbers_in_phase2_after) ||
					!GetAndResetCounter(P().myQueue, P().myPhase1_offset_to_resume_buf, phase1_offset_to_resume))
				{
					return false;
				}

				const double phase1Time = t.getElapsedTime();
				const num64 new_numbers_in_phase2 = numbers_in_phase2_after - numbers_in_phase2_before;
				const num64 numbers_processed = phase1_offset_to_resume ? (phase1_offset_to_resume - global_offset_before) : (N - global_offset_before);

				myNumbersProcessedTotal += numbers_processed;
				const double f = myNumbersProcessedTotal / RangeGen::total_numbers_to_check;
				boinc_fraction_done((f < 1.0) ? f : 1.0);

				if (numbers_in_phase2_after > myPhase2MaxNumbersCount)
				{
					// Read counter value which was saved before executing the kernel that caused overflow
					unsigned int counters[2] = {};
					CL_CHECKED_CALL(clEnqueueReadBuffer, P().myQueue, P().myPhase2_numbers_count_buf, CL_TRUE, 0, sizeof(unsigned int) * 2, &counters, 0, nullptr, nullptr);
					numbers_in_phase2_after = counters[1];
				}

				LOG(2, "Phase 1: processed " << numbers_processed << " numbers in " << phase1Time << " seconds, " << new_numbers_in_phase2 << " numbers left (" << ((new_numbers_in_phase2 * 100.0) / numbers_processed) << "%), phase 2: " << numbers_in_phase2_before << " -> " << numbers_in_phase2_after << " numbers");

				if (phase1_offset_to_resume)
				{
					LOG(2, "[PERFORMANCE] Phase 2 buffer overflowed, phase 1 will be resumed starting with the chunk of numbers that caused overflow.\n[PERFORMANCE] Numbers will not be skipped. Consider increasing myPhase2MaxNumbersCount to improve performance.\n");
					global_offset = phase1_offset_to_resume;
				}
			}

			if ((numbers_in_phase2_after >= myPhase2MinNumbersCount) || phase1_offset_to_resume)
			{
				if ((numbers_in_phase2_after > 0) && !ProcessNumbersPhases2_3(static_cast<unsigned int>(numbers_in_phase2_after)))
				{
					return false;
				}
				numbers_in_phase2_after = 0;
			}

			if (!SaveProgressForLargePrimes(FirstPrime, LastPrime, global_offset, numbers_in_phase2_after == 0))
			{
				return false;
			}
		}
	}

	const int error_code = errno;
	LOG_ERROR("Failed to wait for semaphore, errno " << error_code);
	return false;
}

bool OpenCL::GetPlatformID(cl_platform_id* clSelectedPlatformID)
{
	char chBuffer[1024];
	cl_uint num_platforms;
	*clSelectedPlatformID = nullptr;

	CL_CHECKED_CALL(clGetPlatformIDs, 0, nullptr, &num_platforms);
	if (num_platforms == 0)
	{
		return false;
	}

	std::vector<cl_platform_id> clPlatformIDs(num_platforms);

	CL_CHECKED_CALL(clGetPlatformIDs, num_platforms, clPlatformIDs.data(), nullptr);
	for (cl_uint i = 0; i < num_platforms; ++i)
	{
		if (clGetPlatformInfo(clPlatformIDs[i], CL_PLATFORM_NAME, sizeof(chBuffer), &chBuffer, nullptr) != CL_SUCCESS)
		{
			continue;
		}
		LOG_ERROR("OpenCL platform available: " << chBuffer);

		if (!PlatformSupported(chBuffer))
		{
			continue;
		}

		cl_uint ciDeviceCount = 0;
		if (clGetDeviceIDs(clPlatformIDs[i], DEVICE_TYPE_TO_USE, 0, nullptr, &ciDeviceCount) != CL_SUCCESS)
		{
			continue;
		}

		if (ciDeviceCount > 0)
		{
			LOG(1, "Using '" << chBuffer << "' OpenCL platform\n");
			LOG_ERROR("Using '" << chBuffer << "' OpenCL platform\n");
			*clSelectedPlatformID = clPlatformIDs[i];
			break;
		}
	}

	if (*clSelectedPlatformID == nullptr)
	{
		*clSelectedPlatformID = clPlatformIDs[0];
	}

	return true;
}

bool OpenCL::WaitForQueue(cl_command_queue queue, cl_event event)
{
	CL_CHECKED_CALL(clFlush, queue);

	cl_int event_status;
	do
	{
		// It's highly unlikely that the queue will be empty right after it's been sent to GPU, so start iterations with 0.5 ms sleep, not with clGetEventInfo
		HiResSleep(0.5);
		CL_CHECKED_CALL(clGetEventInfo, event, CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(event_status), &event_status, nullptr);
	} while (event_status != CL_COMPLETE);

	return true;
}

bool OpenCL::GetCounter(cl_command_queue queue, cl_mem buf, unsigned int &counter)
{
	CL_CHECKED_CALL(clEnqueueReadBuffer, queue, buf, CL_TRUE, 0, sizeof(counter), &counter, 0, nullptr, nullptr);
	return true;
}

bool OpenCL::SetFoundCapacity(unsigned int capacity)
{
	const cl_uint data[2] = { 0, capacity };
	CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().myAmicable_numbers_count_buf, CL_TRUE, 0, sizeof(data), data, 0, nullptr, nullptr);
	return true;
}

bool OpenCL::ResetCounter(cl_command_queue queue, cl_mem buf)
{
	CL_CHECKED_CALL(clEnqueueWriteBuffer, queue, buf, CL_TRUE, 0, sizeof(num64), ourZeroBuf, 0, nullptr, nullptr);
	return true;
}

bool OpenCL::GetAndResetCounter(cl_command_queue queue, cl_mem buf, unsigned int &counter)
{
	CL_CHECKED_CALL(clEnqueueReadBuffer, queue, buf, CL_TRUE, 0, sizeof(counter), &counter, 0, nullptr, nullptr);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, queue, buf, CL_TRUE, 0, sizeof(counter), ourZeroBuf, 0, nullptr, nullptr);
	return true;
}

bool OpenCL::GetAndResetCounter(cl_command_queue queue, cl_mem buf, num64 &counter)
{
	CL_CHECKED_CALL(clEnqueueReadBuffer, queue, buf, CL_TRUE, 0, sizeof(counter), &counter, 0, nullptr, nullptr);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, queue, buf, CL_TRUE, 0, sizeof(counter), ourZeroBuf, 0, nullptr, nullptr);
	return true;
}

// Divide total_size into 2^k equal parts, each part's size is <= max_size and also a multiple of 1024
unsigned int OpenCL::GetMaxPhaseSize(num64 total_size, unsigned int max_size)
{
	for (unsigned int k = 1;; ++k)
	{
		const num64 result = (((total_size >> k) + 1023) / 1024) * 1024;
		if (result <= max_size)
		{
			return static_cast<unsigned int>(result ? result : 1024ULL);
		}
	}
}

bool OpenCL::Test()
{
	std::vector<SPairToCheck> pairs;

	IF_CONSTEXPR(LogLevel > 0)
	{
		const char* files[] = {
			"c2_3.txt", "c2_4.txt", "c2_5.txt", "c2_6.txt", "c2_7.txt", "c2_8.txt", "c2_9.txt", "c2_10.txt", "c2_11.txt", "c2_12.txt", "c2_13.txt", "c2_14.txt", "c2_15.txt", "c2_16.txt", "c2_17.txt", "c2_18.txt", "c2_19.txt", "c2_20.txt", "c2_21.txt"
		};
		for (const char* name : files)
		{
			std::ifstream f(name);
			while (!f.eof())
			{
				char buf[1024];
				buf[0] = '\0';
				f.getline(buf, sizeof(buf)); // type author year
				if (buf[0] == '\0')
					break;

				buf[0] = '\0';
				f.getline(buf, sizeof(buf)); // M and factorization

				const num128 m = atoi128(buf);
				if (m >= SearchLimit::value)
					break;

				buf[0] = '\0';
				f.getline(buf, sizeof(buf)); // N and factorization

				const num128 n = atoi128(buf);

				buf[0] = '\0';
				f.getline(buf, sizeof(buf)); // empty line

				pairs.emplace_back(m, m + n);
			}
			LOG(1, "Loaded " << name);
		}
	}
	else
	{
#include "self_test.inl"

		for (size_t i = 0; i < ARRAYSIZE(pair_list); ++i)
		{
			pairs.emplace_back(atoi128(pair_list[i][0]), atoi128(pair_list[i][1]));
		}
	}

	return TestPairs(pairs);
}

// The arguments of the phase 2 kernel for the search (TestPairs points them to its own buffers and restores them)
bool OpenCL::SetPhase2KernelArguments()
{
	CL_CHECKED_CALL(setKernelArguments, P().myCheckPairPhase2,
		mySmallPrimesBuf,
		P().myPhase2_numbers_buf,
		0,
		myPrimeInversesBuf,
		myPowersOfP_128SumInverses_offsets_Buf,
		myPowersOfP_128SumInverses_Buf,
		myPrimeInverses_128_Buf,
		mySumEstimates_128_Buf,
		myPQ_Buf,
		P().myPhase3_numbers_count_buf,
		P().myPhase3_numbers_buf,
		P().myAmicable_numbers_count_buf,
		P().myAmicable_numbers_data_buf
	);
	return true;
}

// Checks the given pairs (M, sigma(M)) with the same kernels as all numbers (phases 1-3), all of them must be found
bool OpenCL::TestPairs(std::vector<SPairToCheck>& pairs)
{
	const unsigned int total_numbers = static_cast<unsigned int>(pairs.size());
	LOG(1, "Total numbers: " << total_numbers);
	if (total_numbers == 0)
	{
		LOG(1, "Skipping test");
		return true;
	}

	Timer overall_timer;

	struct TmpBuf
	{
		FORCEINLINE ~TmpBuf() { clReleaseMemObject(mem); }
		FORCEINLINE operator cl_mem() { return mem; }

		cl_mem mem;
	} pairsToCheck, buf;

	pairs.resize(((pairs.size() + myWorkGroupSize - 1) / myWorkGroupSize) * myWorkGroupSize);
	CL_CHECKED_CALL_WITH_RESULT(pairsToCheck.mem, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(SPairToCheck) * pairs.size(), 0);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, pairsToCheck, CL_TRUE, 0, sizeof(SPairToCheck) * pairs.size(), pairs.data(), 0, nullptr, nullptr);

	CL_CHECKED_CALL_WITH_RESULT(buf.mem, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(num64) * 4 * pairs.size(), 0);
	// the found numbers go to buf (one entry per pair) during the test
	if (!SetFoundCapacity(static_cast<unsigned int>(pairs.size())))
	{
		return false;
	}
	struct RestoreCapacity
	{
		OpenCL* self;
		~RestoreCapacity() { self->SetFoundCapacity(AmicableFoundMax); }
	} restoreCapacity{ this };

	cl_event event = nullptr;

	unsigned int filtered_numbers_count = 0;
	unsigned int filtered_numbers2_count = 0;
	unsigned int amicable_numberscount = 0;

	//
	// Phase 1
	//
	{
		Timer t;

		const size_t globalSizePhase1 = pairs.size();
		CL_CHECKED_CALL(setKernelArguments, myCheckPairs, mySmallPrimesBuf, myPrimeInversesBuf, myPQ_Buf, myPowerOf2SumInverses_Buf, myPowersOfP_128SumInverses_offsets_Buf, myPowersOfP_128SumInverses_Buf, myPrimeInverses_128_Buf, mySumEstimates_128_Buf, pairsToCheck, P().myPhase2_numbers_count_buf, P().myPhase2_numbers_buf, P().myAmicable_numbers_count_buf, buf);
		CL_CHECKED_CALL(clEnqueueNDRangeKernel, P().myQueue, myCheckPairs, 1, nullptr, &globalSizePhase1, &myWorkGroupSize, 0, nullptr, &event);

		if (!WaitForQueue(P().myQueue, event) || !GetAndResetCounter(P().myQueue, P().myPhase2_numbers_count_buf, filtered_numbers_count) || !GetCounter(P().myQueue, P().myAmicable_numbers_count_buf, amicable_numberscount))
		{
			return false;
		}

		const double phase1Time = t.getElapsedTime();
		LOG(2, "Phase 1: " << phase1Time << " seconds, " << filtered_numbers_count << " numbers left (" << ((filtered_numbers_count * 100.0) / total_numbers) << "%), " << amicable_numberscount << " amicable numbers found");
	}

	//
	// Phase 2
	//
	{
		Timer t;

		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().myPhase2_numbers_buf, CL_TRUE, filtered_numbers_count * sizeof(num64) * 4, myWorkGroupSize * sizeof(num64) * 4, ourZeroBuf, 0, nullptr, nullptr);
		CL_CHECKED_CALL(setKernelArguments, P().myCheckPairPhase2,
			mySmallPrimesBuf,
			P().myPhase2_numbers_buf,
			0,
			myPrimeInversesBuf,
			myPowersOfP_128SumInverses_offsets_Buf,
			myPowersOfP_128SumInverses_Buf,
			myPrimeInverses_128_Buf,
			mySumEstimates_128_Buf,
			myPQ_Buf,
			P().myPhase3_numbers_count_buf,
			P().myPhase3_numbers_buf,
			P().myAmicable_numbers_count_buf,
			buf
		);

		size_t globalSizePhase2 = filtered_numbers_count;
		if (globalSizePhase2 & (myWorkGroupSize - 1))
		{
			globalSizePhase2 += myWorkGroupSize - (globalSizePhase2 & (myWorkGroupSize - 1));
		}
		CL_CHECKED_CALL(clEnqueueNDRangeKernel, P().myQueue, P().myCheckPairPhase2, 1, nullptr, &globalSizePhase2, &myWorkGroupSize, 0, nullptr, &event);

		if (!WaitForQueue(P().myQueue, event) || !GetAndResetCounter(P().myQueue, P().myPhase3_numbers_count_buf, filtered_numbers2_count) || !GetCounter(P().myQueue, P().myAmicable_numbers_count_buf, amicable_numberscount))
		{
			return false;
		}

		const double phase2Time = t.getElapsedTime();
		LOG(2, "Phase 2: " << phase2Time << " seconds, " << filtered_numbers2_count << " numbers left (" << ((filtered_numbers2_count * 100.0) / total_numbers) << "%), " << amicable_numberscount << " amicable numbers found");
	}

	//
	// Phase 3
	//
	cl_mem phase3_buffers[2] = {P().myPhase2_numbers_buf, P().myPhase3_numbers_buf};
	const size_t phase3_workGroupSize = myWorkGroupSize;

	unsigned int filtered_numbers3_count = filtered_numbers2_count;
	for (int phase3_iteration = 1; filtered_numbers3_count > 0; ++phase3_iteration)
	{
		Timer t1;

		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, phase3_buffers[phase3_iteration & 1], CL_TRUE, filtered_numbers3_count * sizeof(num64) * 4, phase3_workGroupSize * sizeof(num64) * 4, ourZeroBuf, 0, nullptr, nullptr);

		const unsigned int max_size_phase3 = GetMaxPhaseSize(filtered_numbers3_count, 1 << 18);
		size_t globalSizePhase3;
		for (unsigned int offset = 0; offset < filtered_numbers3_count; offset += static_cast<unsigned int>(globalSizePhase3))
		{
			globalSizePhase3 = std::min(filtered_numbers3_count - offset, max_size_phase3);
			if (globalSizePhase3 & (phase3_workGroupSize - 1))
			{
				globalSizePhase3 += phase3_workGroupSize - (globalSizePhase3 & (phase3_workGroupSize - 1));
			}
			CL_CHECKED_CALL(setKernelArguments, P().myCheckPairPhase3,
				phase3_buffers[phase3_iteration & 1],
				offset,
				myPrimeInversesBuf,
				myPrimeReciprocalsBuf,
				P().myPhase3_numbers_count_buf,
				phase3_buffers[(phase3_iteration - 1) & 1],
				P().myAmicable_numbers_count_buf,
				buf,
				phase3_iteration << 12,
				mySmallPrimesBuf,
				myPowersOfP_128SumInverses_offsets_Buf,
				myPowersOfP_128SumInverses_Buf,
				myPrimeInverses_128_Buf,
				P().myPhase3WideStatusBuf,
				0U
			);
			CL_CHECKED_CALL(clEnqueueNDRangeKernel, P().myQueue, P().myCheckPairPhase3, 1, nullptr, &globalSizePhase3, &phase3_workGroupSize, 0, nullptr, (offset + globalSizePhase3 >= filtered_numbers3_count) ? &event : nullptr);
		}

		if (!WaitForQueue(P().myQueue, event) || !GetAndResetCounter(P().myQueue, P().myPhase3_numbers_count_buf, filtered_numbers3_count) || !GetCounter(P().myQueue, P().myAmicable_numbers_count_buf, amicable_numberscount))
		{
			return false;
		}

		const double phase3Time = t1.getElapsedTime();
		LOG(2, "Phase 3, iteration " << phase3_iteration << ": " << phase3Time << " seconds, " << filtered_numbers3_count << " numbers left (" << ((filtered_numbers3_count * 100.0) / total_numbers) << "%), " << amicable_numberscount << " amicable numbers found");
	}

	if (!GetAndResetCounter(P().myQueue, P().myAmicable_numbers_count_buf, amicable_numberscount))
	{
		return false;
	}

	// Restore the arguments of the phase 2 kernel (they point to the temporary buffer now)
	if (!SetPhase2KernelArguments())
	{
		return false;
	}

	LOG(1, (static_cast<int>(total_numbers) - static_cast<int>(amicable_numberscount)) << " amicable numbers missed");

	const double dt1 = overall_timer.getElapsedTime();
	LOG(1, "OpenCL: " << dt1 << " seconds\n\n");

	pairs.resize(total_numbers);
	std::vector<SFoundPair> found_pairs(amicable_numberscount);
	if (amicable_numberscount > 0)
	{
		CL_CHECKED_CALL(clEnqueueReadBuffer, P().myQueue, buf, CL_TRUE, 0, sizeof(num64) * 4 * amicable_numberscount, found_pairs.data(), 0, nullptr, nullptr);
		for (const SFoundPair& pair : found_pairs)
		{
			const num128 M = CombineNum128(pair.M, pair.M_high);
			auto it = std::lower_bound(pairs.begin(), pairs.end(), M, [](const SPairToCheck& a, const num128 b) { return a.M < b; });
			if ((it != pairs.end()) && (it->M == M))
			{
				it->targetSum = NUM128_MAX;
			}
		}
	}

	num64 num_missed = 0;

	for (const SPairToCheck& pair : pairs)
	{
		if (pair.targetSum != NUM128_MAX)
		{
			if (num_missed == 0)
				LOG_ERROR("Fatal error: self-test returned false negative for these known amicable numbers:");

			++num_missed;

			if (num_missed < 10)
				LOG_ERROR(pair.M);
			else if (num_missed == 10)
				LOG_ERROR("... and more");
		}
	}

	return (num_missed == 0);
}

// Fast integer cube root. Returns num64 m such that m^3 <= n < (m+1)^3 for all n > 0
FORCEINLINE num64 IntegerCbrt(const num128 n)
{
	if (n < 8)
	{
		return 1;
	}

	unsigned long index;
	if (LowWord(n))
	{
		_BitScanReverse64(&index, LowWord(n));
	}
	else
	{
		_BitScanReverse64(&index, HighWord(n));
		index += 64;
	}

	index = (index * ((1048576 / 3) + 1)) >> 20;

	num64 result = num64(1) << index;
	for (num64 cur_bit = result >> 1; cur_bit > 0; cur_bit >>= 1)
	{
		const num64 k = result | cur_bit;
		if ((num128(k) * k) * k <= n)
		{
			result = k;
		}
	}
	return result;
}

bool OpenCL::AddRange(const RangeData& r, num64 offset, bool& added)
{
	added = false;

	num128 prime_limit = 0;

	switch (myLargestPrimePower)
	{
	case 1:
		prime_limit = (SearchLimit::value - 1) / r.value;
		if (prime_limit > SearchLimit::MainPrimeTableBound)
		{
			prime_limit = SearchLimit::MainPrimeTableBound;
		}

		if (r.sum - r.value < r.value)
		{
			const num128 prime_limit2 = (r.sum - 1) / (r.value * 2 - r.sum);
			if (prime_limit > prime_limit2)
			{
				prime_limit = prime_limit2;
			}
		}
		break;

	case 2:
		prime_limit = static_cast<num64>(sqrt(Num128ToDouble(SearchLimit::value) / Num128ToDouble(r.value)));
		if (r.sum - r.value < r.value)
		{
			// r.sum * (p^2 + p + 1) = r.value * p^2 * 2
			// (r.sum - r.value * 2) * p^2 + r.sum * (p + 1) = 0
			// (r.value * 2 - r.sum) * p^2 - r.sum * (p + 1) = 0
			// (r.value * 2 - r.sum) * p^2 - r.sum * p - r.sum = 0
			// (r.value * 2 - r.sum) / r.sum * p^2 - p - 1 = 0
			// (r.value * 2 / r.sum - 1) * p^2 - p - 1 = 0
			const double A = Num128ToDouble(r.value * 2 - r.sum) / Num128ToDouble(r.sum);
			const num64 prime_limit2 = static_cast<num64>((1.0 + sqrt(1.0 + 4.0 * A)) / (2.0 * A));
			if (prime_limit > prime_limit2)
			{
				prime_limit = prime_limit2;
			}
		}
		break;

	case 3:
		prime_limit = IntegerCbrt(SearchLimit::value / r.value);
		if (r.sum - r.value < r.value)
		{
			// r.sum * (p^3 + p^2 + p + 1) = r.value * p^3 * 2
			// (r.value * 2 - r.sum) * p^3 = r.sum * (p^2 + p + 1)
			// A * p^3 = r.sum * (p^2 + p + 1)
			const num128 A = r.value * 2 - r.sum;

			// Lower bound
			// A * p^3 = r.sum * p^2 < r.sum * (p^2 + p + 1)
			// A * p^3 = r.sum * p^2
			// A * p = r.sum
			// p = r.sum / A
			double x1 = Num128ToDouble(r.sum) / Num128ToDouble(A);

			// Upper bound
			// A * p^3 = r.sum * p^2 * 2 > r.sum * (p^2 + p + 1)
			// A * p^3 = r.sum * p^2 * 2
			// A * p = r.sum * 2
			// p = r.sum * 2 / A
			double x2 = x1 * 2.0;
			// This upper bound holds when p^2 * 2 > p^2 + p + 1
			// p^2 * 2 > p^2 + p + 1
			// p^2 > p + 1
			// p^2 - p - 1 > 0 is true when p > (1 + sqrt(5)) / 2 = 1.6180339887498948482045868343656...

			double x;
			for (;;)
			{
				x = (x1 + x2) * 0.5;
				if ((static_cast<num64>(x1) == static_cast<num64>(x2)) || (x == x1) || (x == x2))
				{
					break;
				}
				if (Num128ToDouble(A)*x*x*x > Num128ToDouble(r.sum)*(x*x + x + 1))
				{
					x2 = x;
				}
				else
				{
					x1 = x;
				}
			}

			const num64 prime_limit2 = static_cast<num64>(x2);
			if (prime_limit > prime_limit2)
			{
				prime_limit = prime_limit2;
			}
		}
		break;
	}

	// Find the first prime "p" such that "p > prime_limit" is true
	num64 a = 0;
	num64 b = NumPrimes;
	const num64 t = LowWord(prime_limit);

	// Use prime counting function estimates to bound the initial search values
	if (t > 3)
	{
		unsigned long index;
		_BitScanReverse64(&index, t);
		if (index < 37)
		{
			// oeis.org sequence A185192
			static const num64 num_primes[37] = { 0, 2, 4, 6, 11, 18, 31, 54, 97, 172, 309, 564, 1028, 1900, 3512, 6542, 12251, 23000, 43390, 82025, 155611, 295947, 564163, 1077871, 2063689, 3957809, 7603553, 14630843, 28192750, 54400028, 105097565, 203280221, 393615806, 762939111, 1480206279, 2874398515, 5586502348 };
			const num64 b1 = num_primes[index];
			if (b > b1)
				b = b1;
		}

		double x = log(static_cast<double>(t)) - 1.0;
		for (int i = 0; i <= 1; ++i, x -= 0.1)
		{
			const num64 c = static_cast<num64>(t / x);
			if (a < c && c < b)
			{
				if (t < GetNthPrime(c))
					b = c;
				else
					a = c + 1;
			}
		}
	}

	// Run the search
	while (a < b)
	{
		const num64 c = (a + b) >> 1;
		if (t < GetNthPrime(c))
			b = c;
		else
			a = c + 1;
	}

	const num64 index_begin = std::max<num64>(r.index_start_prime, offset / 2);
	if ((a > index_begin) && (index_begin < (offset + myMainBufferSize) / 2))
	{
		const unsigned int range_size = static_cast<unsigned int>(std::min(a - index_begin, (offset + myMainBufferSize) / 2 - index_begin));

		myRanges.emplace_back(r.value, r.sum, static_cast<unsigned int>(index_begin), range_size);
		myTotalNumbersInRanges += range_size;
		myNumbersProcessedTotal += range_size;
		added = true;

		if ((myTotalNumbersInRanges >= Phase1Limit) || (myRanges.size() >= myMaxRangesCount))
		{
			if (!ProcessNumbers())
			{
				return false;
			}
			CleanupRanges();
		}
	}

	return true;
}

bool OpenCL::ProcessNumbers()
{
	Timer total_time;

	LOG(2, "[PROCESS_NUMBERS] Processing " << myTotalNumbersInRanges << " numbers (" << myRanges.size() << " ranges)");

	const num64 averageRangeSize = myTotalNumbersInRanges / myRanges.size();
	unsigned long index;
	_BitScanReverse64(&index, averageRangeSize);
	const unsigned int lookupTableShift = index;

	myLookupTable.clear();
	for (unsigned int range_index = 0, range_offset = 0;;)
	{
		while ((range_index < myRanges.size()) && (range_offset >= myRanges[range_index].size))
		{
			range_offset -= static_cast<unsigned int>(myRanges[range_index].size);
			++range_index;
		}
		if (range_index >= myRanges.size())
		{
			break;
		}
		myLookupTable.emplace_back(range_index, range_offset);
		range_offset += 1 << lookupTableShift;
	}

	cl_event event = nullptr;
	unsigned int numbers_in_phase2_before = 0;
	unsigned int numbers_in_phase2_after = 0;
	num64 phase1_offset_to_resume = 0;
	num64 global_offset = 0;

	do
	{
		if (!GetCounter(P().myQueue, P().myPhase2_numbers_count_buf, numbers_in_phase2_before))
		{
			return false;
		}

		//
		// Phase 1
		//
		{
			Timer t;

			CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().myPhase1_offset_to_resume_buf, CL_TRUE, 0, sizeof(num64), ourZeroBuf, 0, nullptr, nullptr);

			CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, myRangesTable_Buf, CL_TRUE, 0, myRanges.size() * sizeof(RangeDataGPU), myRanges.data(), 0, nullptr, nullptr);
			CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, myRangesLookupTable_Buf, CL_TRUE, 0, myLookupTable.size() * sizeof(LookupDataGPU), myLookupTable.data(), 0, nullptr, nullptr);

			CL_CHECKED_CALL(clSetKernelArg, mySearchMultipleRanges, 10, sizeof(unsigned int), &lookupTableShift);
			CL_CHECKED_CALL(clSetKernelArg, mySearchMultipleRanges, 12, sizeof(unsigned int), &myTotalNumbersInRanges);

			const num64 global_offset_before = global_offset;
			const unsigned int max_size_phase1 = GetMaxPhaseSize(myTotalNumbersInRanges - global_offset, myPhase1MaxKernelSize);

			size_t globalSizePhase1;
			for (unsigned int k = 0; global_offset < myTotalNumbersInRanges; global_offset += globalSizePhase1, ++k)
			{
				if (k)
				{
					CL_CHECKED_CALL(clFlush, P().myQueue);
				}
				globalSizePhase1 = std::min<num64>(myTotalNumbersInRanges - global_offset, max_size_phase1);
				if (globalSizePhase1 & (myWorkGroupSize - 1))
				{
					globalSizePhase1 += myWorkGroupSize - (globalSizePhase1 & (myWorkGroupSize - 1));
				}
				CL_CHECKED_CALL(clSetKernelArg, mySearchMultipleRanges, 11, sizeof(global_offset), &global_offset);
				CL_CHECKED_CALL(clEnqueueTask, P().myQueue, mySaveCounter, 0, nullptr, nullptr);
				CL_CHECKED_CALL(clEnqueueNDRangeKernel, P().myQueue, mySearchMultipleRanges, 1, nullptr, &globalSizePhase1, &myWorkGroupSize, 0, nullptr, (global_offset + globalSizePhase1 >= myTotalNumbersInRanges) ? &event : nullptr);
			}

			if (!WaitForQueue(P().myQueue, event) ||
				!GetCounter(P().myQueue, P().myPhase2_numbers_count_buf, numbers_in_phase2_after) ||
				!GetAndResetCounter(P().myQueue, P().myPhase1_offset_to_resume_buf, phase1_offset_to_resume))
			{
				return false;
			}

			const double phase1Time = t.getElapsedTime();
			const num64 new_numbers_in_phase2 = numbers_in_phase2_after - numbers_in_phase2_before;
			const num64 numbers_processed = phase1_offset_to_resume ? (phase1_offset_to_resume - global_offset_before) : (myTotalNumbersInRanges - global_offset_before);

			if (numbers_in_phase2_after > myPhase2MaxNumbersCount)
			{
				// Read counter value which was saved before executing the kernel that caused overflow
				unsigned int counters[2] = {};
				CL_CHECKED_CALL(clEnqueueReadBuffer, P().myQueue, P().myPhase2_numbers_count_buf, CL_TRUE, 0, sizeof(unsigned int) * 2, &counters, 0, nullptr, nullptr);
				numbers_in_phase2_after = counters[1];
			}

			LOG(2, "Phase 1: processed " << numbers_processed << " numbers in " << phase1Time << " seconds, " << new_numbers_in_phase2 << " numbers left (" << ((new_numbers_in_phase2 * 100.0) / numbers_processed) << "%), phase 2: " << numbers_in_phase2_before << " -> " << numbers_in_phase2_after << " numbers");

			if (phase1_offset_to_resume)
			{
				LOG(2, "[PERFORMANCE] Phase 2 buffer overflowed, phase 1 will be resumed starting with the chunk of numbers that caused overflow.\n[PERFORMANCE] Numbers will not be skipped. Consider increasing myPhase2MaxNumbersCount to improve performance.\n");
				global_offset = phase1_offset_to_resume;
			}
		}

		if ((numbers_in_phase2_after >= myPhase2MinNumbersCount) || phase1_offset_to_resume)
		{
			if ((numbers_in_phase2_after > 0) && !ProcessNumbersPhases2_3(static_cast<unsigned int>(numbers_in_phase2_after)))
			{
				return false;
			}
			numbers_in_phase2_after = 0;
		}
	} while (phase1_offset_to_resume);

	LOG(2, "[PROCESS_NUMBERS] Finished " << myTotalNumbersInRanges << " numbers (" << myRanges.size() << " ranges) in " << total_time.getElapsedTime() << " seconds\n");
	return true;
}

bool OpenCL::ProcessNumbersPhases2_3(unsigned int numbers_in_phase2)
{
	cl_event event = nullptr;

	//
	// Phase 2
	//
	unsigned int phase3_numbers_buf_count;
	{
		Timer t;

		// (kernel.cl: Phase2Prefilter) the numbers of the form r^2 s and the others with D < 3 N^(2/3)
		{
			const size_t globalSize = (static_cast<size_t>(numbers_in_phase2) + myWorkGroupSize - 1) / myWorkGroupSize * myWorkGroupSize;
			CL_CHECKED_CALL(setKernelArguments, P().myPhase2Prefilter, P().myPhase2_numbers_buf, numbers_in_phase2, P().myAmicable_numbers_count_buf, P().myAmicable_numbers_data_buf);
			cl_event ep = nullptr;
			CL_CHECKED_CALL(clEnqueueNDRangeKernel, P().myQueue, P().myPhase2Prefilter, 1, nullptr, &globalSize, &myWorkGroupSize, 0, nullptr, g_SieveProfile ? &ep : nullptr);
			if (ep) P().myPhase2ProfileEvents.emplace_back(0, ep);
		}
		if ((numbers_in_phase2 >= Phase2SortMinCount) && !SortPhase2Numbers(numbers_in_phase2))
		{
			return false;
		}

		// debugging: AMICABLE_DUMP_PHASE2=<file> writes the numbers of the first phase 2 run with at least 2^20 numbers (32 bytes each,
		// CheckPairPhase2 input)
		static bool phase2Dumped = false;
		const char* dump2Path = getenv("AMICABLE_DUMP_PHASE2");
		if (dump2Path && !phase2Dumped && (numbers_in_phase2 >= (1U << 20)))
		{
			phase2Dumped = true;
			std::vector<num64> data(static_cast<size_t>(numbers_in_phase2) * 4);
			CL_CHECKED_CALL(clEnqueueReadBuffer, P().myQueue, P().myPhase2_numbers_buf, CL_TRUE, 0, data.size() * sizeof(num64), data.data(), 0, nullptr, nullptr);
			FILE* f = fopen(dump2Path, "wb");
			if (f)
			{
				fwrite(data.data(), sizeof(num64), data.size(), f);
				fclose(f);
			}
		}

		if (numbers_in_phase2 < myPhase2MaxNumbersCount)
		{
			CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().myPhase2_numbers_buf, CL_TRUE, numbers_in_phase2 * sizeof(num64) * 4, std::min(static_cast<unsigned int>(myWorkGroupSize), myPhase2MaxNumbersCount - numbers_in_phase2) * sizeof(num64) * 4, ourZeroBuf, 0, nullptr, nullptr);
		}

		const unsigned int max_size_phase2 = GetMaxPhaseSize(numbers_in_phase2, myPhase1MaxKernelSize / 4);
		LOG(2, "Phase 2: max kernel size = " << max_size_phase2 << " numbers");

		size_t globalSizePhase2;
		for (unsigned int k = 0, global_offset = 0; global_offset < numbers_in_phase2; global_offset += static_cast<unsigned int>(globalSizePhase2), ++k)
		{
			if (k)
			{
				CL_CHECKED_CALL(clFlush, P().myQueue);
			}
			globalSizePhase2 = std::min(numbers_in_phase2 - global_offset, max_size_phase2);
			if (globalSizePhase2 & (myWorkGroupSize - 1))
			{
				globalSizePhase2 += myWorkGroupSize - (globalSizePhase2 & (myWorkGroupSize - 1));
			}
			CL_CHECKED_CALL(clSetKernelArg, P().myCheckPairPhase2, 2, sizeof(global_offset), &global_offset);
			const bool last = (global_offset + static_cast<unsigned int>(globalSizePhase2) >= numbers_in_phase2);
			cl_event ek = nullptr;
			CL_CHECKED_CALL(clEnqueueNDRangeKernel, P().myQueue, P().myCheckPairPhase2, 1, nullptr, &globalSizePhase2, &myWorkGroupSize, 0, nullptr, last ? &event : (g_SieveProfile ? &ek : nullptr));
			if (ek) P().myPhase2ProfileEvents.emplace_back(3, ek);
			if (last && g_SieveProfile)
			{
				clRetainEvent(event);
				P().myPhase2ProfileEvents.emplace_back(3, event);
			}
		}

		if (!WaitForQueue(P().myQueue, event) || !GetAndResetCounter(P().myQueue, P().myPhase3_numbers_count_buf, phase3_numbers_buf_count))
		{
			return false;
		}
		// (/sieve_profile: GPU time of the parts of phase 2)
		for (const auto& pe : P().myPhase2ProfileEvents)
		{
			cl_ulong t0 = 0, t1 = 0;
			if ((clGetEventProfilingInfo(pe.second, CL_PROFILING_COMMAND_START, sizeof(t0), &t0, nullptr) == CL_SUCCESS) &&
				(clGetEventProfilingInfo(pe.second, CL_PROFILING_COMMAND_END, sizeof(t1), &t1, nullptr) == CL_SUCCESS))
			{
				P().mySieveStatPhase2Times[pe.first] += (t1 - t0) * 1e-9;
			}
			clReleaseEvent(pe.second);
		}
		P().myPhase2ProfileEvents.clear();
		P().mySieveStatPhase2Wall += t.getElapsedTime();

		const double phase2Time = t.getElapsedTime();
		LOG(2, "Phase 2: processed " << numbers_in_phase2 << " numbers in " << phase2Time << " seconds, " << phase3_numbers_buf_count << " numbers left (" << ((phase3_numbers_buf_count * 100.0) / numbers_in_phase2) << "%)");

		// debugging: AMICABLE_DUMP_PHASE3=<file> writes the first numbers that go to phase 3 (32 bytes each, CheckPairPhase3 input)
		// (AMICABLE_DUMP_PHASE3_MIN=<n>: the first run with at least n numbers)
		static bool phase3Dumped = false;
		const char* dumpPath = getenv("AMICABLE_DUMP_PHASE3");
		static const unsigned int dumpMin = getenv("AMICABLE_DUMP_PHASE3_MIN") ? static_cast<unsigned int>(atof(getenv("AMICABLE_DUMP_PHASE3_MIN"))) : 1;
		if (dumpPath && !phase3Dumped && (phase3_numbers_buf_count >= std::max(dumpMin, 1U)))
		{
			phase3Dumped = true;
			std::vector<num64> data(static_cast<size_t>(std::min(phase3_numbers_buf_count, 1U << 18)) * 4);
			CL_CHECKED_CALL(clEnqueueReadBuffer, P().myQueue, P().myPhase3_numbers_buf, CL_TRUE, 0, data.size() * sizeof(num64), data.data(), 0, nullptr, nullptr);
			FILE* f = fopen(dumpPath, "wb");
			if (f)
			{
				fwrite(data.data(), sizeof(num64), data.size(), f);
				fclose(f);
			}
		}
	}

	//
	// Phase 3
	//
	if ((phase3_numbers_buf_count > 0) && !RunPhase3(phase3_numbers_buf_count, 1))
	{
		return false;
	}

	if (!ResetCounter(P().myQueue, P().myPhase2_numbers_count_buf))
	{
		return false;
	}

	// (the backlog of phase 3 is processed when it is half full, and at the end of the sieve batch: FinishSieveBatch)
	if (P().myPhase3BacklogEnabled && (P().myPhase3BacklogCount >= Phase3BacklogCapacity() / 2) && !FlushPhase3Backlog())
	{
		return false;
	}

	return true;
}

// Phase 3 for "phase3_numbers_buf_count" numbers in the input buffer of iteration "firstIteration" (iteration k reads phase3_buffers[k & 1], the
// output of phase 2 is the input of iteration 1). In the sieve batches (P().myPhase3BacklogEnabled), the numbers that are left when there are
// few of them go to the backlog: they need many more primes (up to N^(1/3) or 2 * N / D), and their launches left most of the GPU idle
// (k4: iterations 4-7 took 40 s for 0.9% of the numbers of phase 3, a few thousand numbers per launch). The backlog runs them together.
bool OpenCL::RunPhase3(unsigned int phase3_numbers_buf_count, int firstIteration)
{
	struct Phase3Time
	{
		Timer t;
		double& total;
		~Phase3Time() { total += t.getElapsedTime(); }
	} phase3TotalTime{ Timer(), P().mySieveStatPhase3Wall };
	if (firstIteration == 1)
	{
		P().mySieveStatPhase3Numbers += phase3_numbers_buf_count;
	}
	cl_mem phase3_buffers[2] = { P().myPhase2_numbers_buf, P().myPhase3_numbers_buf };
	const size_t phase3_workGroupSize = myWorkGroupSize;

	unsigned int phase3_numbers_buf_count2 = phase3_numbers_buf_count;
	for (int phase3_iteration = firstIteration; phase3_numbers_buf_count2 > 0; ++phase3_iteration)
	{
		if (P().myPhase3BacklogEnabled && (phase3_iteration > firstIteration) && (phase3_numbers_buf_count2 < Phase3BacklogMinCount) &&
			(P().myPhase3BacklogCount + phase3_numbers_buf_count2 <= Phase3BacklogCapacity()))
		{
			CL_CHECKED_CALL(clEnqueueCopyBuffer, P().myQueue, phase3_buffers[phase3_iteration & 1], P().myPhase3BacklogBuf, 0, P().myPhase3BacklogCount * sizeof(num64) * 4, phase3_numbers_buf_count2 * sizeof(num64) * 4, 0, nullptr, nullptr);
			P().myPhase3BacklogIteration = (P().myPhase3BacklogCount == 0) ? phase3_iteration : std::min(P().myPhase3BacklogIteration, phase3_iteration);
			P().myPhase3BacklogCount += phase3_numbers_buf_count2;
			LOG(2, "Phase 3, iteration " << phase3_iteration << ": " << phase3_numbers_buf_count2 << " numbers moved to the backlog (" << P().myPhase3BacklogCount << " numbers)");
			return true;
		}

		Timer t1;

		if (phase3_numbers_buf_count2 < myPhase2MaxNumbersCount)
		{
			CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, phase3_buffers[phase3_iteration & 1], CL_FALSE, phase3_numbers_buf_count2 * sizeof(num64) * 4, std::min(static_cast<unsigned int>(phase3_workGroupSize), myPhase2MaxNumbersCount - phase3_numbers_buf_count2) * sizeof(num64) * 4, ourZeroBuf, 0, nullptr, nullptr);
		}

		const unsigned int max_size_phase3 = GetMaxPhaseSize(phase3_numbers_buf_count2, myPhase1MaxKernelSize / 16);
		// trial division up to prime index 4096 * 2^(iteration - 1) (it was 4096 * iteration: up to ~160 iterations, each a launch and a
		// blocking read, for the numbers whose cofactor has no small prime factors), never beyond the table of primes
		const int phase3MaxI = static_cast<int>(std::min<num64>(num64(4096) << std::min(phase3_iteration - 1, 16), ReciprocalsTableSize128 - 64));
		if ((phase3MaxI == static_cast<int>(ReciprocalsTableSize128 - 64)) && (phase3_iteration > 18))
		{
			LOG_ERROR("Phase 3: " << phase3_numbers_buf_count2 << " numbers need primes beyond the table");
			return false;
		}
		LOG(2, "Phase 3, iteration " << phase3_iteration << ": max kernel size = " << max_size_phase3 << " numbers");

		// the late iterations (few numbers with many primes each): CheckPairPhase3Wide (one warp per number) finds the first round of
		// 1024 primes where CheckPairPhase3 has something to do, CheckPairPhase3 continues there (kernel.cl)
		// (k4: phase 3 3.91 -> 1.08 s with the wide kernel from iteration 2, 1.73 s from iteration 4; AMICABLE_P3_WIDE_FROM=0 disables it)
		static const int wideFrom = getenv("AMICABLE_P3_WIDE_FROM") ? static_cast<int>(StrToNumber(getenv("AMICABLE_P3_WIDE_FROM"))) : 2;
		const cl_uint wideMode = ((wideFrom > 0) && (phase3_iteration >= wideFrom)) ? 1U : 0U;
		if (wideMode)
		{
			const cl_uint ones = 0xFFFFFFFFU;
			CL_CHECKED_CALL(clEnqueueFillBuffer, P().myQueue, P().myPhase3WideStatusBuf, &ones, sizeof(ones), 0, sizeof(cl_uint) * static_cast<size_t>(phase3_numbers_buf_count2), 0, nullptr, nullptr);
			const cl_uint wideCount = phase3_numbers_buf_count2;
			const cl_int tableLimit = static_cast<cl_int>(ReciprocalsTableSize128 - 64);
			CL_CHECKED_CALL(setKernelArguments, P().myCheckPairPhase3Wide, phase3_buffers[phase3_iteration & 1], wideCount, phase3MaxI, tableLimit, mySmallPrimesBuf, myPrimeInversesBuf, myPrimeInverses_128_Buf, P().myPhase3WideStatusBuf);
			const size_t wideGlobal = (static_cast<size_t>(wideCount) * 32 + myWorkGroupSize - 1) / myWorkGroupSize * myWorkGroupSize;
			CL_CHECKED_CALL(clEnqueueNDRangeKernel, P().myQueue, P().myCheckPairPhase3Wide, 1, nullptr, &wideGlobal, &phase3_workGroupSize, 0, nullptr, nullptr);
		}

		size_t globalSizePhase3;
		for (unsigned int k = 0, global_offset = 0; global_offset < phase3_numbers_buf_count2; global_offset += static_cast<unsigned int>(globalSizePhase3), ++k)
		{
			if (k)
			{
				CL_CHECKED_CALL(clFlush, P().myQueue);
			}
			globalSizePhase3 = std::min<num64>(phase3_numbers_buf_count2 - global_offset, max_size_phase3);
			if (globalSizePhase3 & (myWorkGroupSize - 1))
			{
				globalSizePhase3 += myWorkGroupSize - (globalSizePhase3 & (myWorkGroupSize - 1));
			}
			CL_CHECKED_CALL(setKernelArguments, P().myCheckPairPhase3,
				phase3_buffers[phase3_iteration & 1],
				global_offset,
				myPrimeInversesBuf,
				myPrimeReciprocalsBuf,
				P().myPhase3_numbers_count_buf,
				phase3_buffers[(phase3_iteration - 1) & 1],
				P().myAmicable_numbers_count_buf,
				P().myAmicable_numbers_data_buf,
				phase3MaxI,
				mySmallPrimesBuf,
				myPowersOfP_128SumInverses_offsets_Buf,
				myPowersOfP_128SumInverses_Buf,
				myPrimeInverses_128_Buf,
				P().myPhase3WideStatusBuf,
				wideMode
			);
			CL_CHECKED_CALL(clEnqueueNDRangeKernel, P().myQueue, P().myCheckPairPhase3, 1, nullptr, &globalSizePhase3, &phase3_workGroupSize, 0, nullptr, nullptr);
		}

		// The blocking read of the counter waits for the launches (in-order queue), the reset doesn't block (ourZeroBuf is static; the
		// padding above doesn't either). WaitForQueue polls with sleeps of 0.5 ms: with it and three blocking calls, every iteration took
		// ~1 ms more than its kernels (k4: ~1600 iterations, most of them with ~0.1 ms of GPU work).
		CL_CHECKED_CALL(clEnqueueReadBuffer, P().myQueue, P().myPhase3_numbers_count_buf, CL_TRUE, 0, sizeof(phase3_numbers_buf_count2), &phase3_numbers_buf_count2, 0, nullptr, nullptr);
		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().myPhase3_numbers_count_buf, CL_FALSE, 0, sizeof(phase3_numbers_buf_count2), ourZeroBuf, 0, nullptr, nullptr);

		const double phase3Time = t1.getElapsedTime();
		{
			const int k = std::min(phase3_iteration, 19);
			P().mySieveStatPhase3IterTime[k] += phase3Time;
			P().mySieveStatPhase3IterNumbers[k] += phase3_numbers_buf_count;
			++P().mySieveStatPhase3IterLaunches[k];
		}
		LOG(2, "Phase 3, iteration " << phase3_iteration << ": processed " << phase3_numbers_buf_count << " numbers in " << phase3Time << " seconds, " << phase3_numbers_buf_count2 << " numbers left (" << ((phase3_numbers_buf_count2 * 100.0) / phase3_numbers_buf_count) << "%)");
		phase3_numbers_buf_count = phase3_numbers_buf_count2;
	}

	LOG(2, "Phase 3 finished");
	return true;
}

// Sorts the numbers of phase 2 by the primes they need (kernel.cl: Phase2Bucket), in P().myPhase3_numbers_buf, and copies them back. Nothing is
// sorted if almost all numbers are in one bucket (most work units: the sort would only cost time).
bool OpenCL::SortPhase2Numbers(unsigned int count)
{
	cl_uint buckets[Phase2Buckets] = {};
	CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().myPhase2BucketsBuf, CL_TRUE, 0, sizeof(buckets), buckets, 0, nullptr, nullptr);
	const size_t globalSize = (static_cast<size_t>(count) + myWorkGroupSize - 1) / myWorkGroupSize * myWorkGroupSize;
	CL_CHECKED_CALL(setKernelArguments, P().myPhase2BucketCount, P().myPhase2_numbers_buf, count, P().myPhase2BucketsBuf);
	cl_event ec = nullptr;
	CL_CHECKED_CALL(clEnqueueNDRangeKernel, P().myQueue, P().myPhase2BucketCount, 1, nullptr, &globalSize, &myWorkGroupSize, 0, nullptr, g_SieveProfile ? &ec : nullptr);
	if (ec) P().myPhase2ProfileEvents.emplace_back(1, ec);
	CL_CHECKED_CALL(clEnqueueReadBuffer, P().myQueue, P().myPhase2BucketsBuf, CL_TRUE, 0, sizeof(buckets), buckets, 0, nullptr, nullptr);
	cl_uint largest = 0;
	cl_uint sum = 0;
	for (cl_uint& b : buckets)
	{
		largest = std::max(largest, b);
		const cl_uint n = b;
		b = sum;
		sum += n;
	}
	if (sum != count)
	{
		LOG_ERROR("Phase 2: " << sum << " numbers in the buckets, " << count << " expected");
		return false;
	}
	if (largest >= count - count / 16)
	{
		return true;
	}
	CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().myPhase2BucketsBuf, CL_TRUE, 0, sizeof(buckets), buckets, 0, nullptr, nullptr);
	CL_CHECKED_CALL(setKernelArguments, P().myPhase2BucketScatter, P().myPhase2_numbers_buf, count, P().myPhase2BucketsBuf, P().myPhase3_numbers_buf);
	cl_event es = nullptr;
	cl_event ey = nullptr;
	CL_CHECKED_CALL(clEnqueueNDRangeKernel, P().myQueue, P().myPhase2BucketScatter, 1, nullptr, &globalSize, &myWorkGroupSize, 0, nullptr, g_SieveProfile ? &es : nullptr);
	CL_CHECKED_CALL(clEnqueueCopyBuffer, P().myQueue, P().myPhase3_numbers_buf, P().myPhase2_numbers_buf, 0, 0, static_cast<size_t>(count) * sizeof(num64) * 4, 0, nullptr, g_SieveProfile ? &ey : nullptr);
	if (es) P().myPhase2ProfileEvents.emplace_back(2, es);
	if (ey) P().myPhase2ProfileEvents.emplace_back(2, ey);
	++P().mySieveStatPhase2Sorts;
	return true;
}

// Runs phase 3 for the numbers in the backlog (the buffers of phases 2 and 3 must be free: no numbers waiting for phase 2)
bool OpenCL::FlushPhase3Backlog()
{
	if (P().myPhase3BacklogCount == 0)
	{
		return true;
	}
	const unsigned int count = P().myPhase3BacklogCount;
	const int iteration = P().myPhase3BacklogIteration;
	P().myPhase3BacklogCount = 0;
	// (the numbers start at the first iteration of the backlog: the others continue at their own primes, kernel.cl: CheckPairPhase3)
	CL_CHECKED_CALL(clEnqueueCopyBuffer, P().myQueue, P().myPhase3BacklogBuf, (iteration & 1) ? P().myPhase3_numbers_buf : P().myPhase2_numbers_buf, 0, 0, count * sizeof(num64) * 4, 0, nullptr, nullptr);
	LOG(2, "Phase 3 backlog: " << count << " numbers from iteration " << iteration);
	const bool enabled = P().myPhase3BacklogEnabled;
	P().myPhase3BacklogEnabled = false;
	const bool result = RunPhase3(count, iteration);
	P().myPhase3BacklogEnabled = enabled;
	return result;
}

void OpenCL::CleanupRanges()
{
	myRanges.clear();
	myTotalNumbersInRanges = 0;
	myLookupTable.clear();
}

bool OpenCL::SaveFoundNumbers()
{
	unsigned int numbers_found = 0;
	if (!GetAndResetCounter(P().myQueue, P().myAmicable_numbers_count_buf, numbers_found))
	{
		return false;
	}

	if (numbers_found > AmicableFoundMax)
	{
		// (kernel.cl: NumberFound) numbers were lost: the work unit must not be reported as done
		LOG_ERROR("Too many numbers found (" << numbers_found << "), the buffer holds " << AmicableFoundMax);
		return false;
	}
	if (numbers_found > 0)
	{
		myFoundPairs.resize(numbers_found);
		CL_CHECKED_CALL(clEnqueueReadBuffer, P().myQueue, P().myAmicable_numbers_data_buf, CL_TRUE, 0, sizeof(num64) * 4 * numbers_found, myFoundPairs.data(), 0, nullptr, nullptr);
		for (const SFoundPair& k : myFoundPairs)
		{
			if ((k.N == 0) || IsPrime(k.N))
			{
				++myAmicableNumbersFound;
				char buf[40];
				fprintf(g_outputFile, "%s\n", itoa128(CombineNum128(k.M, k.M_high), buf, sizeof(buf)));
			}
		}
		fflush(g_outputFile);
	}

	const double progress = static_cast<double>(myNumbersProcessedTotal) / RangeGen::total_numbers_to_check;
	LOG(1, "Progress: " << (progress * 100.0) << "%, " << myAmicableNumbersFound << " amicable numbers found   \r");

	return true;
}

bool OpenCL::SaveProgressForRanges(const RangeData& r, num64 offset, num64 numbersProcessed)
{
	SieveRangeKey key;
	key.numFactors = static_cast<unsigned int>(r.last_factor_index + 1);
	for (unsigned int i = 0; i < key.numFactors; ++i)
	{
		key.primes[i] = r.factors[i].p.Get();
		key.exponents[i] = r.factors[i].k;
	}
	return SaveProgressForRanges(key, offset, numbersProcessed);
}

bool OpenCL::SaveProgressForRanges(const SieveRangeKey& r, num64 offset, num64 numbersProcessed)
{
	if (numbersProcessed == ~num64(0))
	{
		numbersProcessed = myNumbersProcessedTotal;
	}
	// This is only called when the queue is empty: there are no numbers in phases 1-3 on GPU and GPU is idle

	if (!SaveFoundNumbers())
	{
		return false;
	}

	if (boinc_time_to_checkpoint())
	{
		std::stringstream s;
		std::string checkpoint_name, data;

		boinc_resolve_filename_s(CHECKPOINT_LOGICAL_NAME, checkpoint_name);

		// First write how many numbers have been checked so far
		s << numbersProcessed << ' ' << offset << ' ';

		// Then current range
		for (unsigned int i = 0; i < r.numFactors; ++i)
		{
			if (i > 0)
			{
				s << '*';
			}
			s << r.primes[i];
			if (r.exponents[i] > 1)
			{
				s << '^' << r.exponents[i];
			}
		}

		// End of data marker
		s << '#';

		data = s.str();

		FILE* f = boinc_fopen(checkpoint_name.c_str(), "wb");
		if (f)
		{
			fwrite(data.c_str(), sizeof(char), data.length(), f);
			fclose(f);
		}

		boinc_checkpoint_completed();
	}

	return true;
}

bool OpenCL::SaveProgressForLargePrimes(num64 firstPrime, num64 lastPrime, num64 offset, bool queueEmpty)
{
	if (!SaveFoundNumbers())
	{
		return false;
	}

	if (queueEmpty && boinc_time_to_checkpoint())
	{
		std::stringstream s;
		std::string checkpoint_name, data;

		boinc_resolve_filename_s(CHECKPOINT_LOGICAL_NAME, checkpoint_name);

		s << myNumbersProcessedTotal << ' ' << firstPrime << ' ' << lastPrime << ' ' << offset << '#';

		data = s.str();

		FILE* f = boinc_fopen(checkpoint_name.c_str(), "wb");
		if (f)
		{
			fwrite(data.c_str(), sizeof(char), data.length(), f);
			fclose(f);
		}

		boinc_checkpoint_completed();
	}
		
	return true;
}

//-----------------------------------------------------------------------------------------------------------------------------------
// Arithmetic-progression sieve for the first stage (see SieveGPU.h)
//-----------------------------------------------------------------------------------------------------------------------------------

bool OpenCL::InitSieve(cl_device_id /*device*/, num64 localMemSize, num64 maxWorkGroupSize, std::string& options)
{
	SieveGPUConfig cfg = g_SieveConfig;
	// Devices with a smaller work-group limit (256 on many AMD and Intel GPUs): smaller work-groups, as for the other kernels, and
	// segments of at most 32 positions per work-item (the scan). Before, the sieve refused to run there.
	while ((cfg.wgSize > maxWorkGroupSize) && (cfg.wgSize > 32))
	{
		cfg.wgSize >>= 1;
	}
	while ((cfg.segJ * 8 > cfg.wgSize * 32) && (cfg.segJ > 64))
	{
		cfg.segJ >>= 1;
	}
	// The divisibility kill pays off only in the last tier (ranges with a narrow window, see ChooseTier): the ranges stage builds a
	// second program with the kill for their chunks (the kill code in the hit loops slows down the other chunks even when it is
	// not used), the large primes stage uses one program with the kill (it is only active in the last tier).
	const bool twoPrograms = g_SieveKill && !myLargePrimesSieve;
	cfg.killQMax = (g_SieveKill && !twoPrograms) ? SieveKillQMax : 0;
	// with the kill, the undo list needs local memory: smaller survivor lists (the kill also halves the survivors)
	cfg.survLocal = cfg.killQMax ? (SieveSurvLocal / 2) : SieveSurvLocal;
	for (;;)
	{
		if (!mySieve.Init(cfg))
		{
			return false;
		}

		// leave 1 KB for the compiler
		if (mySieve.localMemSize + 1024 <= localMemSize)
		{
			break;
		}
		if (cfg.segJ <= 128)
		{
			LOG_ERROR("Not enough local memory for the sieve: " << mySieve.localMemSize << " bytes needed, " << localMemSize << " bytes available");
			return false;
		}
		cfg.segJ >>= 1;
	}

	// ExpectedSum (self-tests) applies the kill in the last tier
	mySieve.killQMax = g_SieveKill ? SieveKillQMax : 0;

	unsigned long shift;
	_BitScanReverse64(&shift, cfg.segJ);

	auto makeOptions = [&](unsigned int survLocal, unsigned int killQMax, unsigned int maxSlots) -> std::string
	{
		char buf[512];
		sprintf_s(buf, " -D SIEVE_WG_SIZE=%u -D SIEVE_SEG_J=%u -D SIEVE_SEG_SHIFT=%u -D SIEVE_MAX_SLOTS=%u -D SIEVE_MAX_ROW_MODULI=%u -D SIEVE_RANGE_STRIDE=%u -D SIEVE_MARGIN=%u -D SIEVE_SURV_LOCAL=%u -D SIEVE_KILL_QMAX=%u -D SIEVE_UNDO_CAP=%u -D SIEVE_CHUNK_SEGMENTS=%u -D SIEVE_DEFERRED_MAX=%u -D SIEVE_KILL_TIER=%u",
			mySieve.cfg.wgSize,
			mySieve.cfg.segJ,
			static_cast<unsigned int>(shift),
			maxSlots,
			std::max(1U, mySieve.maxRowModuli),
			static_cast<unsigned int>(SieveRangeStride),
			static_cast<unsigned int>(SieveMargin),
			survLocal,
			killQMax,
			static_cast<unsigned int>(SieveUndoCap),
			mySieve.cfg.chunkSegments,
			mySurvivorCapacity,
			mySieve.cfg.numTiers - 1
		);
		return buf;
	};
	// With two programs, the main one never runs chunks of the last tier: its task slots are those of the other tiers (one slot less
	// means three registers less in the sieve loop of the kernel, which then has no spills there)
	options = makeOptions(mySieve.cfg.survLocal, mySieve.cfg.killQMax, twoPrograms ? mySieve.maxSlotsMain : mySieve.maxSlots);
	mySieveOptions = options;
	// (the kill variant has the undo list but smaller survivor lists: it needs less local memory than the main program)
	mySieveKillOptions = twoPrograms ? makeOptions(SieveSurvLocal / 2, SieveKillQMax, mySieve.maxSlots) : std::string();

	std::cerr << "Sieve: work-group size " << mySieve.cfg.wgSize << ", segment " << mySieve.cfg.segJ << " x 8 positions x 6 sub-sieves, " << mySieve.cfg.chunkSegments << " segments per chunk, " << mySieve.localMemSize << " bytes of local memory, tiers:";
	for (unsigned int t = 0; t < mySieve.cfg.numTiers; ++t)
	{
		std::cerr << " Q=" << mySieve.cfg.tierQ[t] << " (" << mySieve.tierInfo[t * 8 + 1] << " moduli, " << mySieve.tierInfo[t * 8 + 2] << " tasks per work-item, ranges >= " << mySieve.cfg.tierMinPrimes[t] << " primes)";
	}
	std::cerr << ", estimated load of the busiest warp " << std::setprecision(3) << mySieve.myMaxWarpLoad / mySieve.myAvgWarpLoad << " x average" << std::endl;

	// A batch is processed with one wait on the host, so it should be big (but the tables of ranges and chunks are limited too)
	mySieveMaxPositions = num64(1) << g_SieveBatchShift;
	return true;
}

thread_local int OpenCL::ourPipe = 0;

// The kernels of pipeline P(): the sieve kernels (first stage or large primes) and the kernels of phases 2 and 3
bool OpenCL::CreatePipeKernels()
{
	if (myUseSieve)
	{
		CL_CHECKED_CALL_WITH_RESULT(P().mySieveRanges, clCreateKernel, myProgram, "SieveRanges");
		if (myProgramKill)
		{
			CL_CHECKED_CALL_WITH_RESULT(P().mySieveRangesKill, clCreateKernel, myProgramKill, "SieveRanges");
		}
		CL_CHECKED_CALL_WITH_RESULT(P().myCheckSieveSurvivors, clCreateKernel, myProgram, "CheckSieveSurvivors");
		CL_CHECKED_CALL_WITH_RESULT(P().mySieveLaunchDone, clCreateKernel, myProgram, "SieveLaunchDone");
		if (!myLargePrimesSieve)
		{
			CL_CHECKED_CALL_WITH_RESULT(P().mySieveTestDeferred, clCreateKernel, myProgram, "SieveTestDeferred");
			CL_CHECKED_CALL_WITH_RESULT(P().mySieveBand, clCreateKernel, myProgram, "SieveBand");
			CL_CHECKED_CALL_WITH_RESULT(P().myCheckSieveSurvivorsPrime, clCreateKernel, myProgram, "CheckSieveSurvivors");
		}
	}
	CL_CHECKED_CALL_WITH_RESULT(P().myCheckPairPhase2, clCreateKernel, myProgram, "CheckPairPhase2");
	CL_CHECKED_CALL_WITH_RESULT(P().myCheckPairPhase3, clCreateKernel, myProgram, "CheckPairPhase3");
	CL_CHECKED_CALL_WITH_RESULT(P().myCheckPairPhase3Wide, clCreateKernel, myProgram, "CheckPairPhase3Wide");
	CL_CHECKED_CALL_WITH_RESULT(P().myPhase2BucketCount, clCreateKernel, myProgram, "Phase2BucketCount");
	CL_CHECKED_CALL_WITH_RESULT(P().myPhase2BucketScatter, clCreateKernel, myProgram, "Phase2BucketScatter");
	CL_CHECKED_CALL_WITH_RESULT(P().myPhase2Prefilter, clCreateKernel, myProgram, "Phase2Prefilter");
	CL_CHECKED_CALL_WITH_RESULT(P().myPhase2BucketsBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(cl_uint) * Phase2Buckets, 0);
	return true;
}

// The buffers of phases 1-3 and of the found numbers of pipeline P() (its queue must exist)
bool OpenCL::CreatePipeBuffers()
{
	CL_CHECKED_CALL_WITH_RESULT(P().myPhase1_offset_to_resume_buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(num64), 0);
	CL_CHECKED_CALL_WITH_RESULT(P().myPhase2_numbers_count_buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(unsigned int) * 2, 0);
	CL_CHECKED_CALL_WITH_RESULT(P().myPhase3_numbers_count_buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(num64), 0);
	CL_CHECKED_CALL_WITH_RESULT(P().myAmicable_numbers_count_buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(num64), 0);

	CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().myPhase1_offset_to_resume_buf, CL_TRUE, 0, sizeof(num64), ourZeroBuf, 0, nullptr, nullptr);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().myPhase2_numbers_count_buf, CL_TRUE, 0, sizeof(unsigned int) * 2, ourZeroBuf, 0, nullptr, nullptr);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().myPhase3_numbers_count_buf, CL_TRUE, 0, sizeof(num64), ourZeroBuf, 0, nullptr, nullptr);
	// counter, capacity of the buffer of found numbers (kernel.cl: NumberFound)
	if (!SetFoundCapacity(AmicableFoundMax))
	{
		return false;
	}

	CL_CHECKED_CALL_WITH_RESULT(P().myPhase2_numbers_buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(num64) * 4 * myPhase2MaxNumbersCount, 0);
	CL_CHECKED_CALL_WITH_RESULT(P().myPhase3_numbers_buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(num64) * 4 * myPhase2MaxNumbersCount, 0);
	CL_CHECKED_CALL_WITH_RESULT(P().myPhase3WideStatusBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(cl_uint) * myPhase2MaxNumbersCount, 0);
	if (myUseSieve && !myLargePrimesSieve)
	{
		CL_CHECKED_CALL_WITH_RESULT(P().myPhase3BacklogBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(num64) * 4 * Phase3BacklogCapacity(), 0);
	}
	CL_CHECKED_CALL_WITH_RESULT(P().myAmicable_numbers_data_buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(num64) * 4 * AmicableFoundMax, 0);
	return true;
}

// Pipeline i > 0 of RunRanges (first stage sieve): its own queue, kernels, batch buffers and kernel arguments, the context, the programs
// and the read-only tables of pipeline 0
bool OpenCL::CreateExtraPipe(cl_device_id device, int i)
{
	PipeScope scope(i);
	CL_CHECKED_CALL_WITH_RESULT(P().myQueue, clCreateCommandQueue, myGPUContext, device, g_SieveProfile ? CL_QUEUE_PROFILING_ENABLE : 0);
	if (!CreatePipeKernels() || !CreatePipeBuffers() || !CreateSievePipeBuffers() || !SetPhase2KernelArguments())
	{
		return false;
	}
	CL_CHECKED_CALL(clFinish, P().myQueue);
	P().mySieveBatchRun.ranges.reserve(static_cast<size_t>(mySieveMaxRanges) * SieveRangeStride);
	P().mySieveBatchRun.chunks.reserve(static_cast<size_t>(mySieveMaxChunks) * 4);
	myNumPipes = i + 1;
	return true;
}

// Releases the queue, the kernels and the buffers of pipeline P()
void OpenCL::ReleasePipe()
{
	SievePipe& p = P();
	const cl_kernel kernels[] = { p.mySieveRanges, p.mySieveRangesKill, p.myCheckSieveSurvivors, p.mySieveLaunchDone, p.mySieveTestDeferred,
		p.myCheckSieveSurvivorsPrime, p.mySieveBand, p.myCheckPairPhase2, p.myCheckPairPhase3, p.myCheckPairPhase3Wide, p.myPhase2BucketCount,
		p.myPhase2BucketScatter, p.myPhase2Prefilter };
	for (cl_kernel k : kernels)
	{
		if (k) clReleaseKernel(k);
	}
	const cl_mem buffers[] = { p.myPhase2BucketsBuf, p.myPhase1_offset_to_resume_buf, p.myPhase2_numbers_count_buf, p.myPhase2_numbers_buf,
		p.myPhase3_numbers_count_buf, p.myPhase3_numbers_buf, p.myPhase3WideStatusBuf, p.myPhase3BacklogBuf, p.myAmicable_numbers_count_buf,
		p.myAmicable_numbers_data_buf, p.mySieveRangesBuf, p.mySieveChunksBuf, p.mySurvivorCountBuf, p.mySurvivorsBuf, p.mySieveStatsBuf,
		p.mySieveDeferredCountBuf, p.mySieveDeferredListBuf, p.mySieveBandTableBuf, p.mySieveBandBitmapBuf, p.mySievePrimeCountBuf };
	for (cl_mem b : buffers)
	{
		if (b) clReleaseMemObject(b);
	}
	if (p.myQueue) clReleaseCommandQueue(p.myQueue);
}

// Adds the statistics of another pipeline to those of pipeline 0 (for the summary at the end)
void OpenCL::AddPipeStats(const SievePipe& o)
{
	SievePipe& p = myPipes[0];
	p.mySieveStatPositions += o.mySieveStatPositions;
	p.mySieveStatSurvivors += o.mySieveStatSurvivors;
	p.mySieveStatChunks += o.mySieveStatChunks;
	p.mySieveStatBatches += o.mySieveStatBatches;
	p.mySieveStatPhase2 += o.mySieveStatPhase2;
	p.mySieveStatOverflows += o.mySieveStatOverflows;
	p.mySieveStatLaunchedPositions += o.mySieveStatLaunchedPositions;
	p.mySieveStatLaunches += o.mySieveStatLaunches;
	p.mySieveStatRounds += o.mySieveStatRounds;
	p.mySieveStatBandChunks += o.mySieveStatBandChunks;
	p.mySieveStatPhase2Sorts += o.mySieveStatPhase2Sorts;
	p.mySieveStatTime += o.mySieveStatTime;
	p.mySieveStatSieveKernelTime += o.mySieveStatSieveKernelTime;
	p.mySieveStatCheckKernelTime += o.mySieveStatCheckKernelTime;
	p.mySieveStatBandKernelTime += o.mySieveStatBandKernelTime;
	for (int k = 0; k < 4; ++k) p.mySieveStatPhase2Times[k] += o.mySieveStatPhase2Times[k];
	p.mySieveStatPhase2Wall += o.mySieveStatPhase2Wall;
	for (int k = 0; k < 3; ++k) p.mySieveStatCheckTimes[k] += o.mySieveStatCheckTimes[k];
	p.mySieveStatPhase3Wall += o.mySieveStatPhase3Wall;
	p.mySieveStatPhase3Numbers += o.mySieveStatPhase3Numbers;
	for (int k = 0; k < 20; ++k)
	{
		p.mySieveStatPhase3IterTime[k] += o.mySieveStatPhase3IterTime[k];
		p.mySieveStatPhase3IterNumbers[k] += o.mySieveStatPhase3IterNumbers[k];
		p.mySieveStatPhase3IterLaunches[k] += o.mySieveStatPhase3IterLaunches[k];
	}
}

bool OpenCL::CreateSieveTables()
{
	struct TableBuf
	{
		cl_mem* buf;
		const void* data;
		size_t size;
	};
	const TableBuf tables[] = {
		{ &mySievePrimeRecsBuf, mySieve.primeRecs.data(), mySieve.primeRecs.size() * sizeof(unsigned int) },
		{ &mySieveModRecsBuf, mySieve.modRecs.data(), mySieve.modRecs.size() * sizeof(unsigned int) },
		{ &mySieveModRecs2Buf, mySieve.modRecs2.data(), mySieve.modRecs2.size() * sizeof(unsigned int) },
		{ &mySieveTtabBuf, mySieve.Ttab.data(), mySieve.Ttab.size() * sizeof(int) },
		{ &mySieveTasksBuf, mySieve.tasks.data(), mySieve.tasks.size() * sizeof(unsigned int) },
		{ &mySieveTierInfoBuf, mySieve.tierInfo.data(), mySieve.tierInfo.size() * sizeof(unsigned int) },
		{ &mySieveWheelTabsBuf, mySieve.wheelTabs.data(), mySieve.wheelTabs.size() * sizeof(unsigned int) },
		{ &mySieveModTaskIdxBuf, mySieve.modTaskIdx.data(), mySieve.modTaskIdx.size() * sizeof(unsigned int) },
		{ &mySieveModTaskListBuf, mySieve.modTaskList.data(), mySieve.modTaskList.size() * sizeof(unsigned int) },
	};
	for (const TableBuf& t : tables)
	{
		CL_CHECKED_CALL_WITH_RESULT(*t.buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, t.size, 0);
		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, *t.buf, CL_TRUE, 0, t.size, t.data, 0, nullptr, nullptr);
	}
	if (!myLargePrimesSieve)
	{
		// the base-2 strong pseudoprimes below 2^32 (Spsp2Table.h) for the primality test of p (kernel.cl: SieveIsPrimeList)
		CL_CHECKED_CALL_WITH_RESULT(mySpsp2Buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, sizeof(g_Spsp2Below2_32), 0);
		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, mySpsp2Buf, CL_TRUE, 0, sizeof(g_Spsp2Below2_32), g_Spsp2Below2_32, 0, nullptr, nullptr);
	}

	mySieveBatch.ranges.reserve(static_cast<size_t>(mySieveMaxRanges) * SieveRangeStride);
	mySieveBatch.chunks.reserve(static_cast<size_t>(mySieveMaxChunks) * 4);
	return true;
}

// The buffers of a sieve batch and the arguments of the sieve kernels of pipeline P() (after CreatePipeKernels and CreatePipeBuffers)
bool OpenCL::CreateSievePipeBuffers()
{
	CL_CHECKED_CALL_WITH_RESULT(P().mySurvivorCountBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(num64), 0);
	{
		const num64 zero64 = 0;
		CL_CHECKED_CALL_WITH_RESULT(P().mySieveDeferredCountBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(num64), 0);
		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().mySieveDeferredCountBuf, CL_TRUE, 0, sizeof(zero64), &zero64, 0, nullptr, nullptr);
		CL_CHECKED_CALL_WITH_RESULT(P().mySievePrimeCountBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(num64), 0);
		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().mySievePrimeCountBuf, CL_TRUE, 0, sizeof(zero64), &zero64, 0, nullptr, nullptr);
	}
	CL_CHECKED_CALL_WITH_RESULT(P().mySurvivorsBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(num64) * (myLargePrimesSieve ? 2 : 1) * mySurvivorCapacity, 0);
	const unsigned int zero = 0;

	if (myLargePrimesSieve)
	{
		// Candidates are uploaded in batches, every batch is one buffer (no more than CL_DEVICE_MAX_MEM_ALLOC_SIZE)
		CL_CHECKED_CALL_WITH_RESULT(P().mySieveRangesBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, sizeof(AmicableCandidate) * mySieveMaxCandidates + 16, 0);

		std::vector<int> loExtra;
		mySieve.LoExtraTable(loExtra);
		CL_CHECKED_CALL_WITH_RESULT(P().mySieveChunksBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, loExtra.size() * sizeof(int), 0);
		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().mySieveChunksBuf, CL_TRUE, 0, loExtra.size() * sizeof(int), loExtra.data(), 0, nullptr, nullptr);

		// bitmap of the primes of one interval of p (8 rows of 2^LargePrimesBitmapShift bits)
		CL_CHECKED_CALL_WITH_RESULT(myPrimeBitsBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, sizeof(cl_uint) * 8 * (static_cast<size_t>(1) << (LargePrimesBitmapShift - 5)), 0);

		const cl_ulong4 lpParams = {};
		const float limitF = static_cast<float>(Num128ToDouble(SearchLimit::value));
		const cl_uint wordsPerRow = 0;
		CL_CHECKED_CALL(setKernelArguments, P().mySieveRanges,
			P().mySieveRangesBuf,
			lpParams,
			limitF,
			P().mySieveChunksBuf,
			myPrimeBitsBuf,
			wordsPerRow,
			zero,
			mySievePrimeRecsBuf,
			mySieveModRecsBuf,
			mySieveModRecs2Buf,
			mySieveTtabBuf,
			mySieveTasksBuf,
			mySieveTierInfoBuf,
			mySieveWheelTabsBuf,
			mySieveModTaskIdxBuf,
			mySieveModTaskListBuf,
			P().mySurvivorCountBuf,
			P().mySurvivorsBuf,
			mySurvivorCapacity
		);
	}
	else
	{
	CL_CHECKED_CALL_WITH_RESULT(P().mySieveRangesBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, sizeof(num64) * SieveRangeStride * mySieveMaxRanges, 0);
	CL_CHECKED_CALL_WITH_RESULT(P().mySieveChunksBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, sizeof(unsigned int) * 4 * mySieveMaxChunks, 0);

	for (cl_kernel kernel : { P().mySieveRanges, P().mySieveRangesKill })
	{
		if (!kernel)
		{
			continue;
		}
		CL_CHECKED_CALL(setKernelArguments, kernel,
			P().mySieveRangesBuf,
			P().mySieveChunksBuf,
			zero,
			mySievePrimeRecsBuf,
			mySieveModRecsBuf,
			mySieveModRecs2Buf,
			mySieveTtabBuf,
			mySieveTasksBuf,
			mySieveTierInfoBuf,
			mySieveWheelTabsBuf,
			mySieveModTaskIdxBuf,
			mySieveModTaskListBuf,
			P().mySurvivorCountBuf,
			P().mySurvivorsBuf,
			mySurvivorCapacity
		);
	}
	}

	CL_CHECKED_CALL_WITH_RESULT(P().mySieveStatsBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(cl_uint) * (SieveStatsHeader + SieveMaxLaunchesPerRound), 0);
	CL_CHECKED_CALL(setKernelArguments, P().mySieveLaunchDone, P().mySurvivorCountBuf, P().mySieveStatsBuf, mySurvivorCapacity, 0U, P().myPhase2_numbers_count_buf, P().mySieveDeferredCountBuf, P().mySievePrimeCountBuf);

	if (myLargePrimesSieve)
	{
		// the large primes variant of CheckSieveSurvivors: p is prime, no deferred list
		CL_CHECKED_CALL(setKernelArguments, P().myCheckSieveSurvivors,
			P().mySurvivorsBuf,
			P().mySurvivorCountBuf,
			mySurvivorCapacity,
			P().mySieveRangesBuf,
			mySmallPrimesBuf,
			myPrimeInversesBuf,
			myPQ_Buf,
			myPowerOf2SumInverses_Buf,
			myPowersOfP_128SumInverses_offsets_Buf,
			myPowersOfP_128SumInverses_Buf,
			myPrimeInverses_128_Buf,
			mySumEstimates_128_Buf,
			P().myPhase1_offset_to_resume_buf,
			P().myPhase2_numbers_count_buf,
			P().myPhase2_numbers_buf,
			P().myAmicable_numbers_count_buf,
			P().myAmicable_numbers_data_buf
		);
	}
	else
	{
		CL_CHECKED_CALL_WITH_RESULT(P().mySieveDeferredListBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(num64) * mySurvivorCapacity, 0);
		CL_CHECKED_CALL_WITH_RESULT(P().mySieveBandTableBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, sizeof(cl_uint) * SieveBandEntry * SieveMaxBandChunks * SieveMaxLaunchesPerRound, 0);
		// (the band bitmap has the layout of the 30-wheel: 7 j of 30 for every j of 210, 8 positions each)
		CL_CHECKED_CALL_WITH_RESULT(P().mySieveBandBitmapBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, static_cast<size_t>(SieveMaxBandChunks) * mySieve.cfg.chunkSegments * mySieve.cfg.segJ * 7 * 8 * sizeof(cl_uint), 0);
		// first pass: all survivors (numbers that pass phase 1 are deferred), second pass: the deferred ones with a prime p
		for (int pass = 0; pass < 2; ++pass)
		{
			const cl_uint defer = (pass == 0) ? 1U : 0U;
			CL_CHECKED_CALL(setKernelArguments, (pass == 0) ? P().myCheckSieveSurvivors : P().myCheckSieveSurvivorsPrime,
				P().mySurvivorsBuf,
				(pass == 0) ? P().mySurvivorCountBuf : P().mySievePrimeCountBuf,
				mySurvivorCapacity,
				P().mySieveRangesBuf,
				mySmallPrimesBuf,
				myPrimeInversesBuf,
				myPQ_Buf,
				myPowerOf2SumInverses_Buf,
				myPowersOfP_128SumInverses_offsets_Buf,
				myPowersOfP_128SumInverses_Buf,
				myPrimeInverses_128_Buf,
				mySumEstimates_128_Buf,
				P().myPhase1_offset_to_resume_buf,
				P().myPhase2_numbers_count_buf,
				P().myPhase2_numbers_buf,
				P().myAmicable_numbers_count_buf,
				P().myAmicable_numbers_data_buf,
				P().mySieveDeferredCountBuf,
				P().mySieveDeferredListBuf,
				defer,
				P().mySieveBandTableBuf,
				0U,
				0U,
				P().mySieveBandBitmapBuf
			);
		}
		CL_CHECKED_CALL(setKernelArguments, P().mySieveBand, P().mySieveRangesBuf, P().mySieveBandTableBuf, 0U, 0U, 0U, 0U, mySmallPrimesBuf, P().mySieveBandBitmapBuf);
		const cl_uint spsp2Count = static_cast<cl_uint>(sizeof(g_Spsp2Below2_32) / sizeof(g_Spsp2Below2_32[0]));
		CL_CHECKED_CALL(setKernelArguments, P().mySieveTestDeferred,
			P().mySieveDeferredCountBuf,
			P().mySieveDeferredListBuf,
			P().mySurvivorsBuf,
			P().mySievePrimeCountBuf,
			mySpsp2Buf,
			spsp2Count
		);
	}
	return true;
}

// Range of p for M0 = r.value: the same limits as in AddRange (largest prime power 1)
bool OpenCL::GetPrimeLimit(num128 value, num128 sum, num64 start_prime, num64& p_lo, num64& p_hi)
{
	num128 prime_limit = (SearchLimit::value - 1) / value;
	if (prime_limit > SearchLimit::MainPrimeTableBound)
	{
		prime_limit = SearchLimit::MainPrimeTableBound;
	}

	if (sum - value < value)
	{
		const num128 prime_limit2 = (sum - 1) / (value * 2 - sum);
		if (prime_limit > prime_limit2)
		{
			prime_limit = prime_limit2;
		}
	}

	p_lo = start_prime;
	p_hi = LowWord(prime_limit);
	return (p_lo <= p_hi);
}

// The first stage with the sieve. RangeGen::Iterate is sequential: the ranges are generated here, in blocks of blockSize. Preparing a range
// for the sieve (PrepareSieveRange: prime limits, reduction, the record) is 3/4 of the host work per range, and in the units of many tiny
// ranges the main thread was the limit (3^11 unit 252 of the 3^14..3^10 search: 142.7 M ranges, the main thread busy all the time while the
// GPU waited for batches): g_SieveThreads helper threads prepare the blocks, and so does this thread while it waits for the oldest one. The
// blocks are committed in their order (CommitSieveRange), so the batches, the checkpoints and the results are the same as with one thread.
bool OpenCL::RunRangesSieve(RangeData& r, const Factor* stopAtFactors, num64 offset)
{
	const unsigned int numHelpers = static_cast<unsigned int>(std::max(g_SieveThreads, 0));
	const unsigned int blockSize = 4096;
	const size_t maxBlocks = 2 + 2 * static_cast<size_t>(numHelpers);   // generated and not committed yet

	std::vector<std::unique_ptr<SieveRangeBlock>> blocks;
	std::vector<SieveRangeBlock*> freeBlocks;
	std::deque<SieveRangeBlock*> order;   // generated and not committed yet, oldest first

	// The helper threads end with this function, also on an error: the pending blocks are dropped, a block being prepared is finished
	// first. (Declared after the blocks: the threads are joined before the blocks are freed.)
	struct Helpers
	{
		OpenCL& cl;
		std::vector<std::thread> threads;
		~Helpers()
		{
			{
				std::lock_guard<std::mutex> lock(cl.mySieveRangeMutex);
				cl.mySieveRangeStop = true;
				cl.mySieveRangePending.clear();
			}
			cl.mySieveRangeWorkCV.notify_all();
			for (std::thread& t : threads)
			{
				t.join();
			}
		}
	} helpers{ *this, {} };
	{
		std::lock_guard<std::mutex> lock(mySieveRangeMutex);
		mySieveRangeStop = false;
		mySieveRangePending.clear();
	}
	for (unsigned int i = 0; i < numHelpers; ++i)
	{
		helpers.threads.emplace_back(&OpenCL::SieveRangeHelperLoop, this);
	}

	bool generated = false;
	for (;;)
	{
		while (!generated && (order.size() < maxBlocks))
		{
			if (freeBlocks.empty())
			{
				blocks.emplace_back(new SieveRangeBlock);
				blocks.back()->ranges.reserve(blockSize);
				freeBlocks.push_back(blocks.back().get());
			}
			SieveRangeBlock* block = freeBlocks.back();
			freeBlocks.pop_back();
			generated = FillSieveRangeBlock(*block, r, stopAtFactors, blockSize);
			if (block->ranges.empty())
			{
				freeBlocks.push_back(block);
				break;
			}
			order.push_back(block);
			if (numHelpers > 0)
			{
				{
					std::lock_guard<std::mutex> lock(mySieveRangeMutex);
					mySieveRangePending.push_back(block);
				}
				mySieveRangeWorkCV.notify_one();
			}
		}
		if (order.empty())
		{
			break;
		}

		// the oldest block, prepared (by this thread if no other one has taken it: it prepares the pending blocks while it waits)
		SieveRangeBlock* block = order.front();
		order.pop_front();
		if (numHelpers == 0)
		{
			PrepareSieveRangeBlock(*block);
		}
		else
		{
			std::unique_lock<std::mutex> lock(mySieveRangeMutex);
			while (!block->prepared)
			{
				if (!mySieveRangePending.empty())
				{
					SieveRangeBlock* b = mySieveRangePending.front();
					mySieveRangePending.pop_front();
					lock.unlock();
					PrepareSieveRangeBlock(*b);
					lock.lock();
					b->prepared = true;
				}
				else
				{
					const num64 tw0 = ReadCycleCounter();
					mySieveRangeDoneCV.wait(lock);
					g_HostCycles[7] += ReadCycleCounter() - tw0;
				}
			}
		}

		const num64 tc0 = ReadCycleCounter();
		for (const SieveRangeBlock::Range& range : block->ranges)
		{
			bool explored = false;
			if (!CommitSieveRange(*block, range, offset, explored))
			{
				return false;
			}
			if (explored)
			{
				return WaitSieveBatch();
			}
		}
		g_HostCycles[6] += ReadCycleCounter() - tc0;
		for (int k = 0; k < 6; ++k)
		{
			g_HostCycles[k] += block->cycles[k];
		}
		freeBlocks.push_back(block);

		const double progress = static_cast<double>(myNumbersProcessedTotal) / RangeGen::total_numbers_to_check;
		boinc_fraction_done((progress < 1.0) ? progress : 1.0);
	}

	// (the last batch on pipeline 0 in this thread, after all batches on the GPU; then the found numbers of all pipelines)
	if (!WaitSieveBatch() || !ProcessSieveBatch(mySieveBatch))
	{
		return false;
	}
	for (int i = 0; i < myNumPipes; ++i)
	{
		PipeScope scope(i);
		if (!SaveFoundNumbers())
		{
			return false;
		}
	}
	return true;
}

// Generates up to count ranges from r on into the block (r is the next range afterwards). Returns true when there are no more ranges.
bool OpenCL::FillSieveRangeBlock(SieveRangeBlock& block, RangeData& r, const Factor* stopAtFactors, unsigned int count)
{
	block.ranges.clear();
	block.primes.clear();
	block.exponents.clear();
	std::fill(std::begin(block.cycles), std::end(block.cycles), num64(0));
	block.prepared = false;

	const num64 t0 = ReadCycleCounter();
	bool end = false;
	while (block.ranges.size() < count)
	{
		if (RangeGen::HasReached(r, stopAtFactors))
		{
			end = true;
			break;
		}
		block.ranges.emplace_back();
		SieveRangeBlock::Range& range = block.ranges.back();
		range.value = r.value;
		range.sum = r.sum;
		range.start_prime = r.start_prime;
		range.firstFactor = static_cast<unsigned int>(block.primes.size());
		range.numFactors = static_cast<unsigned int>(r.last_factor_index + 1);
		for (int i = 0; i <= r.last_factor_index; ++i)
		{
			block.primes.push_back(r.factors[i].p.Get());
			block.exponents.push_back(r.factors[i].k);
		}
		if (!RangeGen::Iterate(r) || (RangeGen::cur_largest_prime_power != myLargestPrimePower))
		{
			end = true;
			break;
		}
	}
	g_HostCycles[3] += ReadCycleCounter() - t0;
	return end;
}

void OpenCL::SieveRangeHelperLoop()
{
	// The rounding mode of the main thread (PrimeTables.cpp): the bounds of the windows and the reduction are computed with it (with the
	// default rounding of a new thread, 42 of 11.45e9 survivors of 3^11 unit 252 were different)
	ForceRoundUpFloatingPoint();

	std::unique_lock<std::mutex> lock(mySieveRangeMutex);
	for (;;)
	{
		mySieveRangeWorkCV.wait(lock, [this]() { return mySieveRangeStop || !mySieveRangePending.empty(); });
		if (mySieveRangeStop)
		{
			return;
		}
		SieveRangeBlock* block = mySieveRangePending.front();
		mySieveRangePending.pop_front();
		lock.unlock();
		PrepareSieveRangeBlock(*block);
		lock.lock();
		block->prepared = true;
		mySieveRangeDoneCV.notify_all();
	}
}

void OpenCL::PrepareSieveRangeBlock(SieveRangeBlock& block) const
{
	const num64 t0 = ReadCycleCounter();
	for (SieveRangeBlock::Range& range : block.ranges)
	{
		PrepareSieveRange(block, range);
	}
	block.cycles[0] += ReadCycleCounter() - t0;
}

// Everything about a range that doesn't depend on the batches. Helper threads run it: only read-only state (the tables of mySieve, the
// small primes; SieveGPUTables::EstimatePrimes builds its table thread-safely).
void OpenCL::PrepareSieveRange(SieveRangeBlock& block, SieveRangeBlock::Range& range) const
{
	range.status = SieveRangeBlock::NoPrimes;
	range.impossible = false;
	range.numPrimes = 0.0;
	range.W = -1.0;
	range.Wband = -1.0;
	range.Lfull = 0;
	range.L235 = 0;
	if (!GetPrimeLimit(range.value, range.sum, range.start_prime, range.p_lo, range.p_hi))
	{
		return;
	}

	const num64 te0 = ReadCycleCounter();
	range.numPrimes = SieveGPUTables::EstimatePrimes(range.p_lo, range.p_hi);
	block.cycles[4] += ReadCycleCounter() - te0;
	range.status = SieveRangeBlock::Counted;
	if (g_ExploreNumbers > 0)
	{
		return;
	}

	// necessary conditions (SieveReduction): ranges without possible amicable numbers are skipped and not counted (the task
	// generator doesn't count them in /task_size either)
	const num64* primes = block.primes.data() + range.firstFactor;
	const unsigned int* exponents = block.exponents.data() + range.firstFactor;
	SieveReduction red;
	const num64 tr0 = ReadCycleCounter();
	mySieve.ComputeReduction(range.value, range.sum, primes, exponents, range.numFactors, range.p_lo, range.p_hi, num64(1) << 62, red);
	block.cycles[1] += ReadCycleCounter() - tr0;
	range.impossible = red.impossible;
	if (red.impossible)
	{
		if (g_SieveReduce)
		{
			range.status = SieveRangeBlock::Skipped;
			return;
		}
	}
	else
	{
		range.W = red.W;
		range.Lfull = red.Lfull;
		range.L235 = red.L235;

		// The band is made for W. If W is too wide for it and the power of 2 in N is not fixed, the numbers with a higher power of 2
		// have a smaller window W - ln(sigma(2^v) / 2^v) + ln(sigma(2^v2) / 2^v2) (v2: the power that W assumes): the band is made for the
		// smallest of these windows above 0 if it is narrow enough (v73: 2^5*37*197*69557*p, W = 0.016 from 2^4 || N for p == 1 (mod 4),
		// W = 4e-6 for 2^5 || N, p == 7 (mod 8): 2.7e8 numbers of 1e10 went to phase 3, 15 s of 19 s; the numbers of the other powers
		// have L = N / D below the band and are not decided by it). Only for a window <= 1e-5 (L >= 1e5: such a number needs >= 2e4
		// primes of trial division without the band) of a power at most 3 above v2 (a fair share of the positions): the band costs the
		// same for all positions of the range, also if few of them have that power (k1 = 2^6*97*1777*7559..: W = 8.2e-4 for 2^7 || N, the
		// band made k1 3.7 -> 99 s and its phase 2 numbers 26.7M -> 24.8M)
		range.Wband = red.W;
		// (development aid: AMICABLE_BAND_CLASS=<mask> chooses the classes: 1 the powers of 2, 2 the odd primes; default both)
		static const unsigned int bandClasses = getenv("AMICABLE_BAND_CLASS") ? static_cast<unsigned int>(atoi(getenv("AMICABLE_BAND_CLASS"))) : 3U;
		if ((bandClasses & 1) && (red.W >= SieveBandMaxW) && !red.fixed2)
		{
			const long double base = logl(2.0L - ldexpl(1.0L, -static_cast<int>(red.v2)));
			double best = -1.0;
			unsigned int bestV = 0;
			for (unsigned int v = red.v2 + 1; v <= red.v2 + 16; ++v)
			{
				const double Wv = static_cast<double>(static_cast<long double>(red.W) - (logl(2.0L - ldexpl(1.0L, -static_cast<int>(v))) - base));
				if (!(Wv > 0.0))
				{
					break;
				}
				best = Wv;
				bestV = v;
			}
			if ((best > 0.0) && (best <= std::min(SieveBandMaxW, 1e-5)) && (bestV <= red.v2 + 3))
			{
				range.Wband = best;
			}
		}
		// The same for the odd primes q <= 13 outside the known part (not in g, not fixed by the conditions of the reduction, able to
		// divide N): the numbers with q^e || N (a progression of p) have the window W - ln(sigma(q^e) / q^e) (v73: 2*5^2*13*47*103*177823*p,
		// W = 0.18; 5 || N for p == 4 (mod 5) gives a window 2.5e-6, L = 4.0e5: 2.6e8 numbers of 1e10 went to phase 3, 25 s of 27 s).
		// (With q <= 97 and e <= 3, 72 windows per range, a band came up by chance in ~0.2% of the ranges of the typical units: k1
		// 3.63 -> 3.74 s, t3 2.52 -> 2.65 s; q <= 13 and e <= 2: 12 windows, the larger classes)
		if ((bandClasses & 2) && (red.W >= SieveBandMaxW) && !(range.Wband < red.W))
		{
			static const unsigned int oddPrimes[] = { 3, 5, 7, 11, 13 };
			const num128 Ag = range.sum - range.value;
			double best = -1.0;
			for (const unsigned int q : oddPrimes)
			{
				// q | g (q | M0 and q | sigma(M0)), or q | A and not B = sigma(M0) (then q never divides N = A p + B), or a fixed exponent
				const bool qA = (Ag % num128(static_cast<num64>(q))) == 0;
				const bool qB = (range.sum % num128(static_cast<num64>(q))) == 0;
				bool skip = qA;
				for (unsigned int k = 0; (k < range.numFactors) && !skip; ++k)
				{
					skip = (primes[k] == q) && qB;
				}
				for (unsigned int k = 0; (k < red.numLPrimes) && !skip; ++k)
				{
					skip = (red.lq[k] == q);
				}
				if (skip)
				{
					continue;
				}
				double qe = 1.0;
				double sq = 1.0;
				for (unsigned int e = 1; e <= 2; ++e)
				{
					qe *= q;
					sq += qe;
					const double Wq = red.W - log(sq / qe);
					if (!(Wq > 0.0))
					{
						break;
					}
					if ((Wq <= std::min(SieveBandMaxW, 1e-5)) && ((best < 0.0) || (Wq < best)))
					{
						best = Wq;
					}
				}
			}
			if (best > 0.0)
			{
				range.Wband = best;
			}
		}
	}
	if (g_AnalyzeRate > 0)
	{
		// (CommitSieveRange samples the range)
		return;
	}

	const num64 tp0 = ReadCycleCounter();
	const bool prepared = mySieve.PrepareRange(range.value, range.sum, primes, exponents, range.numFactors, range.p_lo, range.p_hi, mySieve.ChooseTier(static_cast<num64>(range.numPrimes), range.W), true, range.record);
	block.cycles[2] += ReadCycleCounter() - tp0;
	if (!prepared)
	{
		range.status = SieveRangeBlock::Failed;
		return;
	}
	SetBandParameters(range.record, range.Wband, range.p_lo, range.p_hi, primes, range.numFactors);
	range.status = SieveRangeBlock::Ready;
}

// The part of a range that depends on the batches, in the order of the ranges (main thread). The batch is submitted before a range that
// doesn't fit in it, so that ranges are not split between batches (unless a single range is bigger than a batch). Then all ranges before
// the first range of a batch are completely processed, which allows a checkpoint there.
bool OpenCL::CommitSieveRange(const SieveRangeBlock& block, const SieveRangeBlock::Range& range, num64 offset, bool& explored)
{
	num64 numSegments = 0;
	if (range.status != SieveRangeBlock::NoPrimes)
	{
		numSegments = ((range.p_hi - (range.p_lo / SieveWheel) * SieveWheel) / SieveWheel) / mySieve.cfg.segJ + 1;
	}
	const num64 numChunks = (numSegments + mySieve.cfg.chunkSegments - 1) / mySieve.cfg.chunkSegments;
	if (!mySieveBatch.chunks.empty() && (
		(mySieveBatch.ranges.size() / SieveRangeStride + 1 >= mySieveMaxRanges) ||
		(mySieveBatch.chunks.size() / 4 + numChunks + 1 >= mySieveMaxChunks) ||
		(mySieveBatch.positions + numSegments * mySieve.cfg.segJ * SievePositionsPerJ > mySieveMaxPositions)))
	{
		// the GPU processes the batch while the next one is prepared
		if (!SubmitSieveBatch(true, offset))
		{
			return false;
		}
	}
	const num64* primes = block.primes.data() + range.firstFactor;
	const unsigned int* exponents = block.exponents.data() + range.firstFactor;
	if (mySieveBatch.chunks.empty() && mySieveBatch.extra.empty())
	{
		mySieveBatchFirstRange.numFactors = range.numFactors;
		std::copy(primes, primes + range.numFactors, mySieveBatchFirstRange.primes);
		std::copy(exponents, exponents + range.numFactors, mySieveBatchFirstRange.exponents);
		mySieveBatchStartNumbers = myNumbersProcessedTotal;
		mySieveBatchHasFirst = true;
	}

	// Development aid: "/explore <numbers>" prints the factorization where a work unit of that many numbers starting at /from ends
	if (g_ExploreNumbers > 0)
	{
		if (range.status != SieveRangeBlock::NoPrimes)
		{
			myNumbersProcessedTotal += static_cast<num64>(range.numPrimes + 0.5);
		}
		if (myNumbersProcessedTotal >= g_ExploreNumbers)
		{
			std::stringstream f;
			for (unsigned int i = 0; i < range.numFactors; ++i)
			{
				if (i)
				{
					f << '*';
				}
				f << primes[i];
				if (exponents[i] > 1)
				{
					f << '^' << exponents[i];
				}
			}
			std::cerr << "Explore: /to " << f.str() << " /task_size " << myNumbersProcessedTotal << std::endl;
			explored = true;
		}
		return true;
	}
	if (range.status == SieveRangeBlock::NoPrimes)
	{
		return true;
	}

	const double positions = static_cast<double>(mySieve.NumPositions(range.p_lo, range.p_hi));
	myRedStat[0] += positions;
	if (range.impossible)
	{
		myRedStat[1] += positions;
		if (range.status == SieveRangeBlock::Skipped)
		{
			++mySieveStatSkipped;
			return true;
		}
		myNumbersProcessedTotal += static_cast<num64>(range.numPrimes + 0.5);
	}
	else
	{
		const double wb[8] = { 1e-4, 1e-3, 1e-2, 0.03, 0.1, 0.3, 1.0, 1e300 };
		int wi = 0;
		while (range.W >= wb[wi])
		{
			++wi;
		}
		myRedHistW[wi] += positions;

		// Development aid: AMICABLE_DUMP_RANGES=<file> writes every range (factorization of M0, p_lo, p_hi, W) for analysis
		static FILE* rangeDump = getenv("AMICABLE_DUMP_RANGES") ? fopen(getenv("AMICABLE_DUMP_RANGES"), "w") : nullptr;
		if (rangeDump)
		{
			for (unsigned int i = 0; i < range.numFactors; ++i)
			{
				fprintf(rangeDump, "%s%llu^%u", i ? "*" : "", primes[i], exponents[i]);
			}
			fprintf(rangeDump, " %llu %llu %.6e\n", range.p_lo, range.p_hi, range.W);
		}
		const double L = range.Lfull ? static_cast<double>(range.Lfull) : 1e30;
		myRedStat[2] += positions / L;
		myRedStat[3] += positions / (L * static_cast<double>(range.L235));
		unsigned long b = 0;
		if (range.Lfull > 1)
		{
			_BitScanReverse64(&b, range.Lfull);
			b += 1;
		}
		myRedHist[(range.Lfull == 0) ? 63 : b] += positions;
		myNumbersProcessedTotal += static_cast<num64>(range.numPrimes + 0.5);
	}

	// Development aid: "/analyze <rate>" samples positions at this rate and counts survivors of every tier, with and without the
	// divisibility conditions (host only)
	if (g_AnalyzeRate > 0)
	{
		const num64 p_lo = range.p_lo;
		const num64 p_hi = range.p_hi;
		static num64 rnd = 0x9E3779B97F4A7C15ULL;
		auto next = [&]() { rnd ^= rnd << 13; rnd ^= rnd >> 7; rnd ^= rnd << 17; return rnd; };
		const num64 P0 = (p_lo / 30) * 30;
		const num64 numJ = (p_hi - P0) / 30 + 1;
		const double expected = static_cast<double>(numJ * 8) * g_AnalyzeRate;
		num64 n = static_cast<num64>(expected);
		if (static_cast<double>(next() >> 11) * (1.0 / 9007199254740992.0) < expected - static_cast<double>(n))
		{
			++n;
		}
		if (n > 0)
		{
			num64 recs[SieveMaxTiers][SieveRangeStride];
			for (unsigned int t = 0; t < mySieve.cfg.numTiers; ++t)
			{
				if (!mySieve.PrepareRange(range.value, range.sum, primes, exponents, range.numFactors, p_lo, p_hi, t, false, recs[t]))
				{
					return false;
				}
			}
			for (num64 k = 0; k < n; ++k)
			{
				const num64 p = P0 + 30 * (next() % numJ) + SieveW30Host[next() % 8];
				if ((p < p_lo) || (p > p_hi))
				{
					continue;
				}
				++myAnalyzeSamples;
				for (unsigned int t = 0; t < mySieve.cfg.numTiers; ++t)
				{
					if (SieveGPUTables::InsideWindow(recs[t], p, mySieve.ExpectedSum(recs[t], p)))
					{
						++myAnalyzeCounts[t][0];
						if (mySieve.DivisibilityOK(range.value, range.sum, p, mySieve.cfg.tierQ[t]))
						{
							++myAnalyzeCounts[t][1];
						}
						// variants: exponent 1 only, q <= 255 / q <= Q, with and without the power of 2
						const num128 gg = Gcd128Host(range.value, range.sum);
						if (mySieve.DivisibilityOK2(range.value, range.sum, gg, p, std::min(16381U, mySieve.cfg.tierQ[t]), 3)) ++myAnalyzeCounts[t][2];
						if (mySieve.DivisibilityOK2(range.value, range.sum, gg, p, std::min(16381U, mySieve.cfg.tierQ[t]), 4)) ++myAnalyzeCounts[t][3];
						if (mySieve.DivisibilityOK2(range.value, range.sum, gg, p, mySieve.cfg.tierQ[t], 1)) ++myAnalyzeCounts[t][4];
						if (mySieve.DivisibilityOK2(range.value, range.sum, gg, p, mySieve.cfg.tierQ[t], 2)) ++myAnalyzeCounts[t][5];
					}
				}
			}
		}
		return true;
	}

	if (range.status != SieveRangeBlock::Ready)
	{
		LOG_ERROR("Failed to prepare a range for the sieve");
		return false;
	}

	++mySieveStatRanges;
	const num64 ta0 = ReadCycleCounter();
	const bool added = AddSieveRangeRecord(range.record, range.p_lo, range.p_hi, range.W);
	g_HostCycles[5] += ReadCycleCounter() - ta0;
	return added;
}

namespace
{
	// Expected survivors per position by the height W of the window above the known part of N (rough upper estimates from the
	// typical work units: small W means a window that contains every N whose cofactor has no small prime factors). The launches are
	// sized with them (times a factor measured in every round), so that a region with many survivors doesn't overflow the buffers.
	const unsigned int SieveNumWClasses = 7;
	const double SieveSurvivorRateByW[SieveNumWClasses] = { 0.5, 0.3, 0.01, 0.005, 0.003, 0.002, 0.001 };

	unsigned int SieveWClass(double W)
	{
		if (W < 0.0)
		{
			return SieveNumWClasses - 1;
		}
		const double bounds[SieveNumWClasses - 1] = { 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1 };
		unsigned int c = 0;
		while ((c < SieveNumWClasses - 1) && (W >= bounds[c]))
		{
			++c;
		}
		return c;
	}
}

// Development aid (AMICABLE_TEST_BAND): the bitmap of SieveBand for the band chunks of a launch, compared with the host (num128 arithmetic)
bool OpenCL::TestBandBitmap(const SieveBatch* batch, const cl_uint* entries, cl_uint count, cl_uint positions)
{
	std::vector<cl_uint> gpu(positions);
	CL_CHECKED_CALL(clEnqueueReadBuffer, P().myQueue, P().mySieveBandBitmapBuf, CL_TRUE, 0, static_cast<size_t>(positions) * sizeof(cl_uint), gpu.data(), 0, nullptr, nullptr);
	std::vector<cl_uint> host(positions, 0xFFFFFFFFU);
	auto inverse = [](num64 a, num64 m) -> num64
	{
		long long t = 0, nt = 1;
		long long r = static_cast<long long>(m), nr = static_cast<long long>(a % m);
		while (nr != 0)
		{
			const long long q = r / nr;
			const long long t2 = t - q * nt; t = nt; nt = t2;
			const long long r2 = r - q * nr; r = nr; nr = r2;
		}
		return static_cast<num64>((t < 0) ? (t + static_cast<long long>(m)) : t);
	};
	num64 hits = 0;
	for (cl_uint e = 0; e < count; ++e)
	{
		const cl_uint* be = entries + static_cast<size_t>(e) * SieveBandEntry;
		const num64* rec = batch->ranges.data() + static_cast<size_t>(be[0]) * SieveRangeStride;
		const num128 A1 = CombineNum128(rec[0], rec[1]);
		const num128 B1 = CombineNum128(rec[2], rec[3]);
		const num64 P0 = rec[8];
		for (cl_uint i = be[4]; i < be[5]; ++i)
		{
			const num64 r = mySmallPrimesHost[i];
			const num64 am = LowWord(A1 % num128(r));
			if (am == 0)
			{
				continue;
			}
			const num64 bmod = LowWord(B1 % num128(r));
			const num64 inv = inverse((30 * am) % r, r);
			for (unsigned int c = 0; c < 8; ++c)
			{
				const num64 w = (0x1D1713110D0B0701ULL >> (c * 8)) & 0xFF;
				const num64 t = (LowWord((num128(am) * ((P0 + w) % r)) % num128(r)) + bmod) % r;
				const num64 jc = LowWord((num128((r - t) % r) * inv) % num128(r));
				const num64 js = be[1] % r;
				for (num64 d = (jc + r - js) % r; d < be[2]; d += r)
				{
					// (SieveBand skips the positions whose p has a prime factor 11..127 in the entries of at least 2^20 j)
					if (be[2] >= (1U << 20))
					{
						const num64 p = P0 + 30 * (be[1] + d) + w;
						bool smallFactor = false;
						for (num64 q = 11; (q <= 127) && !smallFactor; q += 2)
						{
							smallFactor = ((q % 3) != 0) && ((q % 5) != 0) && ((q % 7) != 0) && ((p % q) == 0) && (p > 127);
						}
						if (smallFactor)
						{
							continue;
						}
					}
					cl_uint& v = host[be[3] + d * 8 + c];
					v = std::min<cl_uint>(v, i);
					++hits;
				}
			}
		}
	}
	num64 mismatches = 0;
	num64 found = 0;
	for (cl_uint k = 0; k < positions; ++k)
	{
		found += (gpu[k] != 0xFFFFFFFFU);
		if (gpu[k] != host[k])
		{
			if (mismatches < 5)
			{
				std::cerr << "Band test: position " << k << ": GPU " << gpu[k] << ", host " << host[k] << std::endl;
			}
			++mismatches;
		}
	}
	std::cerr << "Band test: " << count << " chunks, " << positions << " positions, " << hits << " hits, " << found << " positions with a band factor, " << mismatches << " mismatches" << std::endl;
	return mismatches == 0;
}

// The index in smallPrimes of the smallest prime factor q > 53 of g = gcd(A, B) = A / A1 (they are primes of M0) in the band [begin, end)
// (end if there is none): phase 1 doesn't always remove these primes from N, and the band sieve doesn't see them (it sieves A1 p + B1)
static unsigned int BandGMin(const num64* record, const num64* primes, unsigned int numPrimes, const std::vector<unsigned int>& smallPrimes, unsigned int begin, unsigned int end)
{
	const num128 A = CombineNum128(record[6], record[7]) - CombineNum128(record[4], record[5]);
	const num128 A1 = CombineNum128(record[0], record[1]);
	const num128 g = A / A1;
	unsigned int result = end;
	for (unsigned int k = 0; k < numPrimes; ++k)
	{
		const num64 q = primes[k];
		if ((q > 53) && (q < 0xFFFFFFFFULL) && ((g % num128(q)) == 0))
		{
			const unsigned int idx = static_cast<unsigned int>(std::lower_bound(smallPrimes.begin(), smallPrimes.end(), static_cast<unsigned int>(q)) - smallPrimes.begin());
			if ((idx >= begin) && (idx < result))
			{
				result = idx;
			}
		}
	}
	return result;
}

// TestSieveKnownPairs: band parameters with L0 just below N / D of the pair's own number (the cofactor of n after the primes <= 53 or
// <= 1619 of phase 1, D = sigma(N) - N - 1), so that the band sieve has to find the smallest prime factor of N for the pair to be found
// (kernel.cl: DecideByBand rejects the numbers without a factor in their band)
void OpenCL::SetTestBandParameters(num64* record, const KnownPair& pair) const
{
	record[21] = 0;
	record[22] = 0;
	if (!P().mySieveBand || mySmallPrimesHost.empty())
	{
		return;
	}
	double Lmin = 1e300;
	for (const unsigned int P : { 53U, 1619U })
	{
		num128 N = pair.n;
		num128 S = 1;
		for (size_t k = 0; (k < mySmallPrimesHost.size()) && (mySmallPrimesHost[k] <= P); ++k)
		{
			const num64 q = mySmallPrimesHost[k];
			num128 qe = 1;
			num128 s = 1;
			while ((N % num128(q)) == 0)
			{
				N /= num128(q);
				qe *= num128(q);
				s += qe;
			}
			S *= s;
		}
		const num128 T = (pair.m + pair.n) / S;
		if (T <= N + 1)
		{
			return;
		}
		Lmin = std::min(Lmin, Num128ToDouble(N) / Num128ToDouble(T - N - 1));
	}
	const double maxPrime = static_cast<double>(mySmallPrimesHost.back());
	if (!(Lmin >= 4.0) || (Lmin * 4.0 > maxPrime))
	{
		return;
	}
	const num64 L0 = static_cast<num64>(Lmin * 0.999);
	unsigned int T = 2;
	while ((T < 8) && (static_cast<double>(L0 << (T + 1)) <= maxPrime))
	{
		++T;
	}
	const unsigned int begin = static_cast<unsigned int>(std::upper_bound(mySmallPrimesHost.begin(), mySmallPrimesHost.end(), static_cast<unsigned int>(L0)) - mySmallPrimesHost.begin());
	const unsigned int end = static_cast<unsigned int>(std::upper_bound(mySmallPrimesHost.begin(), mySmallPrimesHost.end(), static_cast<unsigned int>(L0 << T)) - mySmallPrimesHost.begin());
	// (the band primes must be above the primes of phase 1: a band prime that divides A1 p + B1 then divides N)
	if ((L0 > 1619) && (end > begin) && (begin > 0) && (end < mySmallPrimesHost.size()))
	{
		record[21] = L0 | (static_cast<num64>(BandGMin(record, pair.primes, pair.numFactors, mySmallPrimesHost, begin, end)) << 32);
		record[22] = begin | (static_cast<num64>(end) << 32);
		++mySieveStatBandTestRanges;
	}
}

// Band sieve parameters of a range record (kernel_sieve.cl: SieveBand, kernel.cl: DecideByBand). Every number of the range has
// N / D >= 1 / (e^W - 1) (W is the largest ln(T/N) - ln(sigma(K)/K)), so the band starts at L0 = 0.999 / (e^W - 1) (>= 1e4 for
// W < 1e-4; it may contain primes of the sieve, then the band has their positions too). A number with L = N / D needs the band up to
// K L + K + 2 (K = floor(ln N / ln L) <= K0 = floor(ln Nmax / ln L0)). T / N grows with p, and so does the window W(p) of the numbers
// at p: they have L(p) ~ 1 / (e^W(p) - 1), close to L0 except near p_lo. The band ends at K0 Lx + K0 + 2 (Lx = L(px) + 1%, px: 1e-4 of
// the positions of the range are below it; the numbers there with a larger L continue from the end of the band in phases 1 - 3), at most
// at the octaves up to the bands K L + K + 2 <= 1.7 N^(1/3) + 8 that the numbers can need (DecideAtMost2Factors decides the numbers
// with D < 1.88 N^(2/3)), at most SieveBandMaxOctaves, and only within the table of primes. Before (v72) the band always had these
// octaves (16 L0 in the ranges with W ~ 0, where L0 - 4.2 L0 is enough: 1/2 of the hits of SieveBand, 4/5 of its primes). The choice
// only matters for the speed (DecideByBand checks that the band of each number is covered).
void OpenCL::SetBandParameters(num64* record, double W, num64 p_lo, num64 p_hi, const num64* primes, unsigned int numPrimes) const
{
	record[21] = 0;
	record[22] = 0;
	if (!g_SieveBand || !P().mySieveBand || !(W > 0.0) || (W >= SieveBandMaxW) || mySmallPrimesHost.empty())
	{
		return;
	}
	// (with AMICABLE_BAND_MAX_W above 1e-4) W >= 1e-4 only in dense ranges (record[12] bit 41, SieveGPUTables::PrepareRange: the positions
	// without small factors survive): the band costs the same for all positions, in the ranges with few survivors it costs more than the
	// trial division it saves (2^3*13*17*503*1167823.., W = 2.4e-4, 1% survivors: 1.6 -> 5.7 s without this condition)
	if ((W >= 1e-4) && !((record[12] >> 41) & 1))
	{
		return;
	}
	const double L0d = 0.999 / expm1(W);
	const double maxPrime = static_cast<double>(mySmallPrimesHost.back());
	if (L0d * 2.0 > maxPrime)
	{
		return;
	}
	const num64 L0 = static_cast<num64>(L0d);
	const double Nmax = Num128ToDouble(CombineNum128(record[0], record[1])) * static_cast<double>(p_hi) + Num128ToDouble(CombineNum128(record[2], record[3]));
	const double U = 1.7 * cbrt(Nmax) + 8.0;
	unsigned int T = 1;
	while ((T < SieveBandMaxOctaves) && (static_cast<double>(L0 << T) < U) && (static_cast<double>(L0 << (T + 1)) <= maxPrime))
	{
		++T;
	}
	num64 top = L0 << T;

	// the band that the numbers at p >= px need: W(px) = W - (ln(T/N)(p_hi) - ln(T/N)(px)), ln(T/N)(p) = -ln(a + (1 - a) / (p + 1)),
	// a = 1 - M0 / sigma(M0); N after the known power of 2 for K0 (BandEndForChunk: the same for the first p of each chunk)
	const unsigned int v2g = static_cast<unsigned int>((record[12] >> 48) & 0xFF);
	const double Nred = (Num128ToDouble(CombineNum128(record[6], record[7]) - CombineNum128(record[4], record[5])) * static_cast<double>(p_hi) +
		Num128ToDouble(CombineNum128(record[6], record[7]))) / ldexp(1.0, static_cast<int>(v2g));
	const double K0 = std::max(floor(log(Nred) / log(static_cast<double>(L0)) * 1.0001 + 0.0001), 2.0);
	if (p_hi > p_lo)
	{
		const double M0 = Num128ToDouble(CombineNum128(record[4], record[5]));
		const double S0 = Num128ToDouble(CombineNum128(record[6], record[7]));
		const double oneMinusA = M0 / S0;
		const double a = 1.0 - oneMinusA;
		const double px = static_cast<double>(p_lo) + static_cast<double>(p_hi - p_lo) * 1e-4;
		const double ph = static_cast<double>(p_hi);
		const double delta = log1p(oneMinusA * ((ph - px) / ((px + 1.0) * (ph + 1.0))) / (a + oneMinusA / (ph + 1.0)));
		const double Wx = W - delta;
		if (Wx > 0.0)
		{
			const double Lx = 1.01 / expm1(Wx);
			const double need = K0 * Lx + K0 + 2.0;
			if (need < static_cast<double>(top))
			{
				top = std::max(static_cast<num64>(need) + 1, L0 * 2);
			}
		}
	}

	const unsigned int begin = static_cast<unsigned int>(std::upper_bound(mySmallPrimesHost.begin(), mySmallPrimesHost.end(), static_cast<unsigned int>(L0)) - mySmallPrimesHost.begin());
	const unsigned int end = static_cast<unsigned int>(std::upper_bound(mySmallPrimesHost.begin(), mySmallPrimesHost.end(), static_cast<unsigned int>(top)) - mySmallPrimesHost.begin());
	if ((end <= begin) || (begin == 0) || (end >= mySmallPrimesHost.size()))
	{
		return;
	}
	record[21] = L0 | (static_cast<num64>(BandGMin(record, primes, numPrimes, mySmallPrimesHost, begin, end)) << 32);
	record[22] = begin | (static_cast<num64>(end) << 32);
	const float Wf = static_cast<float>(W);
	unsigned int Wbits;
	memcpy(&Wbits, &Wf, sizeof(Wbits));
	record[23] = Wbits | (static_cast<num64>(K0) << 32);
}

// The end (prime index) of the band for the chunk of a band range whose first position is P0 + 30 jFirst (record[23]: W and K0 of
// SetBandParameters): K0 L(pa) + K0 + 2 for the smallest p of the chunk (L(p) is largest there), within the band of the range. In a
// range with p_hi / p_lo ~ 1e3 (2^4*19*137*167*878699: L0 = 9.2e5, but L = 5.3e6 at the p with 1e-4 of the positions below it) the
// range needs a wide band, but only its first chunks do.
cl_uint OpenCL::BandEndForChunk(const num64* rec, num64 jFirst) const
{
	const cl_uint begin = static_cast<cl_uint>(rec[22]);
	const cl_uint end = static_cast<cl_uint>(rec[22] >> 32);
	// (only in ranges of 8 chunks or more: an entry of its own for the first chunks of a short range doubles the setup of SieveBand for
	// it, "t1" 2.6 chunks per range: 18.6 -> 20.3 s)
	if ((rec[23] == 0) || (mySieve.NumChunks(rec[9], rec[10]) < 8))
	{
		return end;
	}
	float Wf;
	const unsigned int Wbits = static_cast<unsigned int>(rec[23]);
	memcpy(&Wf, &Wbits, sizeof(Wf));
	const double K0 = static_cast<double>(rec[23] >> 32);
	const double M0 = Num128ToDouble(CombineNum128(rec[4], rec[5]));
	const double S0 = Num128ToDouble(CombineNum128(rec[6], rec[7]));
	const double oneMinusA = M0 / S0;
	const double a = 1.0 - oneMinusA;
	const double pa = std::max(static_cast<double>(rec[9]), static_cast<double>(rec[8]) + 30.0 * static_cast<double>(jFirst));
	const double ph = static_cast<double>(rec[10]);
	const double delta = (ph > pa) ? log1p(oneMinusA * ((ph - pa) / ((pa + 1.0) * (ph + 1.0))) / (a + oneMinusA / (ph + 1.0))) : 0.0;
	const double Wa = static_cast<double>(Wf) - delta;
	if (!(Wa > 0.0))
	{
		return end;
	}
	const double need = K0 * 1.01 / expm1(Wa) + K0 + 2.0;
	if (!(need < static_cast<double>(mySmallPrimesHost[end - 1])))
	{
		return end;
	}
	const cl_uint e = static_cast<cl_uint>(std::upper_bound(mySmallPrimesHost.begin(), mySmallPrimesHost.end(), static_cast<unsigned int>(need)) - mySmallPrimesHost.begin());
	return std::min(std::max(e, begin + 1), end);
}

bool OpenCL::AddSieveRangeRecord(const num64* record, num64 p_lo, num64 p_hi, double W)
{
	const num64 P0 = record[8];
	const num64 numJ = (p_hi - P0) / SieveWheel + 1;
	const num64 numSegments = (numJ + mySieve.cfg.segJ - 1) / mySieve.cfg.segJ;
	const num64 positionsPerSegment = static_cast<num64>(mySieve.cfg.segJ) * SievePositionsPerJ;

	num64 segment = 0;
	bool first = true;
	while (first || (segment < numSegments))
	{
		if ((mySieveBatch.ranges.size() / SieveRangeStride + 1 >= mySieveMaxRanges) ||
			(mySieveBatch.chunks.size() / 4 + 1 >= mySieveMaxChunks) ||
			(mySieveBatch.positions >= mySieveMaxPositions))
		{
			if (!SubmitSieveBatch(false, 0))
			{
				return false;
			}
		}

		const num64 rangeIndex = mySieveBatch.ranges.size() / SieveRangeStride;
		mySieveBatch.ranges.insert(mySieveBatch.ranges.end(), record, record + SieveRangeStride);

		if (first)
		{
			// p = 2, 3, 5, 7 are not in the sieve (only p coprime to 210 are), check them directly
			const num64 smallPrimes[4] = { 2, 3, 5, 7 };
			for (num64 p : smallPrimes)
			{
				if ((p_lo <= p) && (p <= p_hi))
				{
					mySieveBatch.extra.push_back((p << SieveRangeIndexBits) | rangeIndex);
				}
			}
			first = false;
		}

		while ((segment < numSegments) && (mySieveBatch.chunks.size() / 4 + 1 < mySieveMaxChunks) && (mySieveBatch.positions < mySieveMaxPositions))
		{
			const num64 n = std::min<num64>(mySieve.cfg.chunkSegments, numSegments - segment);
			mySieveBatch.chunks.push_back(static_cast<unsigned int>(rangeIndex));
			mySieveBatch.chunks.push_back(static_cast<unsigned int>(segment * mySieve.cfg.segJ));
			mySieveBatch.chunks.push_back(static_cast<unsigned int>(n));
			mySieveBatch.chunks.push_back(SieveWClass(W));
			mySieveBatch.positions += n * positionsPerSegment;
			segment += n;
		}
	}
	return true;
}

// Runs the sieve on chunks [chunkBegin, chunkEnd) of the current batch. Every launch of the sieve is followed by the survivor check
// and SieveLaunchDone on the device, so the host only waits once per round. If the survivors of a launch don't fit in the survivor
// buffer, only that launch is run again (with smaller launches). Chunks have positionsPerChunk positions each (large primes), or (0) the
// numbers of segments in mySieveBatch.chunks (ranges). numExtra survivors were written to the start of the survivor buffer by the host.
bool OpenCL::RunSieveLaunches(unsigned int chunkBegin, unsigned int chunkEnd, unsigned int killChunks, num64 positionsPerChunk, unsigned int numExtra, const SieveBatch* batch, num64& survivors, bool& overflow)
{
	survivors = 0;
	overflow = false;

	const size_t localSize = mySieve.cfg.wgSize;
	const size_t checkLocalSize = myWorkGroupSize;
	// Development aid: AMICABLE_CHECK_SIZE overrides the global size of CheckSieveSurvivors
	static const size_t checkSizeOverride = getenv("AMICABLE_CHECK_SIZE") ? static_cast<size_t>(StrToNumber(getenv("AMICABLE_CHECK_SIZE"))) : 0;
	const size_t checkGlobalSize = checkSizeOverride ? checkSizeOverride : mySieveCheckGlobalSize;
	std::vector<std::pair<unsigned int, unsigned int>> todo(1, std::make_pair(chunkBegin, chunkEnd));
	std::vector<std::pair<unsigned int, unsigned int>> launches;
	std::vector<cl_uint> stats;
	std::vector<cl_event> sieveEvents, checkEvents, bandEvents, deferEvents, check2Events;
	bool firstRound = true;
	// numbers in the phase 2 buffer (it is empty at the start of a batch: FinishSieveBatch processes it)
	cl_uint phase2Fill = 0;
	auto chunkPositions = [this, positionsPerChunk, batch](unsigned int c) -> num64
	{
		return positionsPerChunk ? positionsPerChunk : (static_cast<num64>(batch->chunks[c * 4 + 2]) * mySieve.cfg.segJ * SievePositionsPerJ);
	};
	// group of a chunk for the estimates (0: W < 1e-3, 1: the others and the large primes stage)
	auto chunkGroup = [batch](unsigned int c) -> unsigned int
	{
		return (batch && (batch->chunks[c * 4 + 3] < 2)) ? 0U : 1U;
	};
	// expected survivors of a chunk without the measured factor P().mySieveRateScale (ranges: by the W class of the range, large primes:
	// the measured rate)
	auto chunkSurvivors = [this, &chunkPositions, batch](unsigned int c) -> double
	{
		const double positions = static_cast<double>(chunkPositions(c));
		return batch ? (positions * SieveSurvivorRateByW[std::min(batch->chunks[c * 4 + 3], SieveNumWClasses - 1)]) : (positions * P().mySieveSurvivorRate);
	};
	auto chunkScale = [this, batch, &chunkGroup](unsigned int c) -> double
	{
		return batch ? P().mySieveRateScale[chunkGroup(c)] : 1.0;
	};
	// chunks of band ranges (SetBandParameters): at most SieveMaxBandChunks per launch (the bitmap)
	auto chunkBand = [batch](unsigned int c) -> unsigned int
	{
		return (batch && (batch->ranges[static_cast<size_t>(batch->chunks[c * 4]) * SieveRangeStride + 21] != 0)) ? 1U : 0U;
	};
	while (!todo.empty())
	{
		// Every launch has its own survivor buffer: it must not get more survivors than the buffer holds (1.5 times the expected number)
		// or more than 2^31 positions (the survivor counter must not wrap around). The launches of a round share the phase 2 buffer: the
		// round must not produce more numbers for phase 2 than the free part of it holds (1.5 times the expected number), else the
		// buffer is processed first.
		// (expected survivors with the measured factors per launch, expected phase 2 numbers per round)
		const double maxLaunchSurvivors = static_cast<double>(mySurvivorCapacity) / 1.5;
		if (phase2Fill > myPhase2MaxNumbersCount / 2)
		{
			P().mySieveStatPhase2 += phase2Fill;
			if (!ProcessNumbersPhases2_3(phase2Fill))
			{
				return false;
			}
			phase2Fill = 0;
			P().mySieveRoundPositions2 = 0;
		}
		const double roundBudget = (static_cast<double>(myPhase2MaxNumbersCount) - phase2Fill) / 1.5;
		auto chunkPhase2 = [this, batch, &chunkGroup](unsigned int c, double survivors) -> double
		{
			return survivors * (batch ? P().mySievePhase2PerSurvivor[chunkGroup(c)] : (P().mySievePhase2Rate / P().mySieveSurvivorRate));
		};
		launches.clear();
		std::vector<std::pair<unsigned int, unsigned int>> later;
		double roundPositions = 0.0;
		double roundPhase2 = 0.0;      // expected
		for (const auto& r : todo)
		{
			unsigned int c = r.first;
			while (c < r.second)
			{
				if ((launches.size() >= SieveMaxLaunchesPerRound) || (!launches.empty() && (roundPhase2 >= roundBudget)))
				{
					later.emplace_back(c, r.second);
					break;
				}
				// at most 2^31 positions per launch (the survivor counter must not wrap around), one kernel and one group per launch
				unsigned int n = 0;
				unsigned int nb = 0;
				num64 positions = 0;
				double expected = 0.0;
				while ((c + n < r.second) && ((n == 0) || ((positions < (num64(1) << 31)) && (expected + chunkSurvivors(c + n) * chunkScale(c + n) < maxLaunchSurvivors) &&
					(roundPhase2 + chunkPhase2(c, expected + chunkSurvivors(c + n) * chunkScale(c + n)) < roundBudget) && ((c + n < killChunks) == (c < killChunks)) &&
					(chunkGroup(c + n) == chunkGroup(c)) && (nb + chunkBand(c + n) <= SieveMaxBandChunks))))
				{
					positions += chunkPositions(c + n);
					expected += chunkSurvivors(c + n) * chunkScale(c + n);
					nb += chunkBand(c + n);
					++n;
				}
				launches.emplace_back(c, c + n);
				roundPositions += static_cast<double>(positions);
				roundPhase2 += chunkPhase2(c, expected);
				c += n;
			}
		}
		if (launches.empty() && firstRound && numExtra)
		{
			// only the extra survivors (p = 2, 3, 5)
			launches.emplace_back(chunkBegin, chunkBegin);
		}
		todo.clear();
		if (launches.empty())
		{
			todo.swap(later);
			continue;
		}

		// band table of the round: the band chunks of every launch (first entry, count, primes [begin, end), bytes of the bitmap)
		std::vector<cl_uint> bandTable;
		std::vector<std::array<cl_uint, 5>> bandLaunch(launches.size(), std::array<cl_uint, 5>{});
		if (P().mySieveBand && batch)
		{
			for (size_t i = 0; i < launches.size(); ++i)
			{
				std::array<cl_uint, 5>& bl = bandLaunch[i];
				bl[0] = static_cast<cl_uint>(bandTable.size() / SieveBandEntry);
				bl[2] = ~0U;
				for (unsigned int c = launches[i].first; c < launches[i].second; ++c)
				{
					const num64* rec = batch->ranges.data() + static_cast<size_t>(batch->chunks[c * 4]) * SieveRangeStride;
					if (rec[21] == 0)
					{
						continue;
					}
					// (the band table and bitmap have the layout of the 30-wheel: first j and number of j in units of 30)
					const cl_uint jCount = batch->chunks[c * 4 + 2] * mySieve.cfg.segJ * 7;
					// (the end of the band for this chunk: BandEndForChunk; the smallest prime of g stays, the kernel caps it at the end)
					const cl_uint entry[SieveBandEntry] = { batch->chunks[c * 4], batch->chunks[c * 4 + 1] * 7, jCount, bl[4], static_cast<cl_uint>(rec[22]), BandEndForChunk(rec, batch->chunks[c * 4 + 1] * 7), static_cast<cl_uint>(rec[21] >> 32), 0 };
					// The next chunk of the same range (its j follow those of the last entry, and so do its bits of the bitmap) extends that
					// entry: SieveBand computes A1, B1, their inverse and the first hit mod every band prime once per entry, and the chunks of
					// a range have the same A1, B1 and P0 (3^10.. units of "t1": 2.6 chunks per range, SieveBand 68% of the time)
					// (v73: only if the end of its band is within 1/32 of the band of that entry, which then covers it)
					cl_uint* last = (bl[1] > 0) ? (bandTable.data() + bandTable.size() - SieveBandEntry) : nullptr;
					if (last && (last[0] == entry[0]) && (last[1] + last[2] == entry[1]) && (entry[5] <= last[5]) && ((last[5] - entry[5]) * 32 <= last[5] - last[4]))
					{
						last[2] += jCount;
					}
					else
					{
						bandTable.insert(bandTable.end(), entry, entry + SieveBandEntry);
						++bl[1];
					}
					bl[4] += jCount * 8;
					bl[2] = std::min(bl[2], entry[4]);
					bl[3] = std::max(bl[3], entry[5]);
					++P().mySieveStatBandChunks;
				}
			}
			if (!bandTable.empty())
			{
				CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().mySieveBandTableBuf, CL_TRUE, 0, bandTable.size() * sizeof(cl_uint), bandTable.data(), 0, nullptr, nullptr);
			}
		}

		// statistics (survivors, overflow flags of the launches) and the survivor counter (which starts with the extra survivors)
		stats.assign(SieveStatsHeader + launches.size(), 0);
		const cl_uint startCount = firstRound ? numExtra : 0;
		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().mySieveStatsBuf, CL_FALSE, 0, stats.size() * sizeof(cl_uint), stats.data(), 0, nullptr, nullptr);
		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().mySurvivorCountBuf, CL_TRUE, 0, sizeof(startCount), &startCount, 0, nullptr, nullptr);

		for (size_t i = 0; i < launches.size(); ++i)
		{
			const unsigned int c = launches[i].first;
			const unsigned int n = launches[i].second - c;
			cl_event e1 = nullptr;
			cl_event e2 = nullptr;
			++P().mySieveStatLaunches;
			for (unsigned int k = c; k < c + n; ++k)
			{
				P().mySieveStatLaunchedPositions += positionsPerChunk ? positionsPerChunk : (static_cast<num64>(batch->chunks[k * 4 + 2]) * mySieve.cfg.segJ * SievePositionsPerJ);
			}
			if (n > 0)
			{
				const size_t globalSize = static_cast<size_t>(n) * localSize;
				const cl_kernel kernel = ((c < killChunks) && P().mySieveRangesKill) ? P().mySieveRangesKill : P().mySieveRanges;
				CL_CHECKED_CALL(clSetKernelArg, kernel, mySieveChunkOffsetArg, sizeof(c), &c);
				CL_CHECKED_CALL(clEnqueueNDRangeKernel, P().myQueue, kernel, 1, nullptr, &globalSize, &localSize, 0, nullptr, g_SieveProfile ? &e1 : nullptr);
			}
			if (P().mySieveBand)
			{
				const std::array<cl_uint, 5>& bl = bandLaunch[i];
				if (bl[1] > 0)
				{
					const cl_uint none = 0xFFFFFFFFU;
					CL_CHECKED_CALL(clEnqueueFillBuffer, P().myQueue, P().mySieveBandBitmapBuf, &none, sizeof(none), 0, static_cast<size_t>(bl[4]) * sizeof(cl_uint), 0, nullptr, nullptr);
					// SieveBand over slices of the entries of the launch, each with at most bandSliceBytes of the bitmap: its hits (atomic_min,
					// scattered) then mostly stay in the L2 cache. With all entries at once (up to 64 band chunks, ~180 MB of bitmap in "t1")
					// 3/4 of SieveBand was the hits (t1: SieveBand 36.8 s, 9.3 s without them). The prime range of a slice is that of its
					// entries. t1 (RTX 5070, 48 MB L2): 16 MB 48.9 s, 32 MB 46.3 s, 48 MB 47.4 s, 64 MB 49.5 s, unsliced 55.9 s.
					// (Development aid: AMICABLE_BAND_SLICE_MB=<n> sets the size of a slice.)
					static const size_t bandSliceBytes = static_cast<size_t>(getenv("AMICABLE_BAND_SLICE_MB") ? StrToNumber(getenv("AMICABLE_BAND_SLICE_MB")) : 32) << 20;
					const size_t bandLocal = 256;
					for (cl_uint first = 0; first < bl[1];)
					{
						const cl_uint* e0 = bandTable.data() + static_cast<size_t>(bl[0] + first) * SieveBandEntry;
						cl_uint count = 0;
						cl_uint pBegin = ~0U;
						cl_uint pEnd = 0;
						size_t bytes = 0;
						while (first + count < bl[1])
						{
							const cl_uint* e = e0 + static_cast<size_t>(count) * SieveBandEntry;
							const size_t eBytes = static_cast<size_t>(e[2]) * 8 * sizeof(cl_uint);
							if ((count > 0) && (bytes + eBytes > bandSliceBytes))
							{
								break;
							}
							bytes += eBytes;
							pBegin = std::min(pBegin, e[4]);
							pEnd = std::max(pEnd, e[5]);
							++count;
						}
						if (pEnd > pBegin)
						{
							CL_CHECKED_CALL(setKernelArgumentsWithOffset, P().mySieveBand, 2, bl[0] + first, count, pBegin, pEnd);
							const size_t bandGlobal = (static_cast<size_t>(pEnd - pBegin) + bandLocal - 1) / bandLocal * bandLocal;
							cl_event e5 = nullptr;
							CL_CHECKED_CALL(clEnqueueNDRangeKernel, P().myQueue, P().mySieveBand, 1, nullptr, &bandGlobal, &bandLocal, 0, nullptr, g_SieveProfile ? &e5 : nullptr);
							if (e5) bandEvents.push_back(e5);
						}
						first += count;
					}
					// Development aid: AMICABLE_TEST_BAND=<n> compares the bitmaps of the first n launches with band chunks with the host
					static int bandTests = getenv("AMICABLE_TEST_BAND") ? static_cast<int>(StrToNumber(getenv("AMICABLE_TEST_BAND"))) : 0;
					if (bandTests > 0)
					{
						--bandTests;
						if (!TestBandBitmap(batch, bandTable.data() + static_cast<size_t>(bl[0]) * SieveBandEntry, bl[1], bl[4]))
						{
							return false;
						}
					}
				}
				for (cl_kernel kc : { P().myCheckSieveSurvivors, P().myCheckSieveSurvivorsPrime })
				{
					CL_CHECKED_CALL(setKernelArgumentsWithOffset, kc, 21, bl[0], bl[1]);
				}
			}
			CL_CHECKED_CALL(clEnqueueNDRangeKernel, P().myQueue, P().myCheckSieveSurvivors, 1, nullptr, &checkGlobalSize, &checkLocalSize, 0, nullptr, g_SieveProfile ? &e2 : nullptr);
			if (P().mySieveTestDeferred)
			{
				cl_event e3 = nullptr;
				cl_event e4 = nullptr;
				CL_CHECKED_CALL(clEnqueueNDRangeKernel, P().myQueue, P().mySieveTestDeferred, 1, nullptr, &checkGlobalSize, &checkLocalSize, 0, nullptr, g_SieveProfile ? &e3 : nullptr);
				CL_CHECKED_CALL(clEnqueueNDRangeKernel, P().myQueue, P().myCheckSieveSurvivorsPrime, 1, nullptr, &checkGlobalSize, &checkLocalSize, 0, nullptr, g_SieveProfile ? &e4 : nullptr);
				if (e3) deferEvents.push_back(e3);
				if (e4) check2Events.push_back(e4);
			}
			const cl_uint launchIndex = static_cast<cl_uint>(i);
			CL_CHECKED_CALL(clSetKernelArg, P().mySieveLaunchDone, 3, sizeof(launchIndex), &launchIndex);
			CL_CHECKED_CALL(clEnqueueTask, P().myQueue, P().mySieveLaunchDone, 0, nullptr, nullptr);
			CL_CHECKED_CALL(clFlush, P().myQueue);
			if (e1) sieveEvents.push_back(e1);
			if (e2) checkEvents.push_back(e2);
		}

		++P().mySieveStatRounds;
		cl_event event = nullptr;
		CL_CHECKED_CALL(clEnqueueReadBuffer, P().myQueue, P().mySieveStatsBuf, CL_FALSE, 0, stats.size() * sizeof(cl_uint), stats.data(), 0, nullptr, &event);
		if (!WaitForQueue(P().myQueue, event))
		{
			return false;
		}
		clReleaseEvent(event);

			// survivors of the launches that didn't overflow; the overflowed launches are run again (their survivors that were checked
			// are checked again, found numbers are deduplicated at the end)
			const num64 roundMeasured = (static_cast<num64>(stats[1]) << 32) + stats[0];
			survivors += roundMeasured;
			unsigned int numOverflows[3] = {};
			num64 donePositions = 0;
			num64 minOverflowPositions = ~num64(0);
			// per group: survivors and their estimate (without the factor) of the launches that didn't overflow, and the largest factor
			// needed by an overflowed launch
			double groupMeasured[2] = {};
			double groupExpected[2] = {};
			double groupOverflowScale[2] = {};
			double groupChecked[2] = {};   // survivors checked by all launches (those that overflowed the survivor buffer: its capacity)
			// Numbers for phase 2 that the estimates predicted for the launches that didn't overflow (what sized the round). The estimates
			// are lowered after a round with enough survivors (or positions) measured, or with at least 1e5 predicted numbers: when they
			// are too high, the rounds they size are too small for the first condition (12574 rounds of 0.18 survivors in a batch of
			// 65535 tiny ranges after a dense range had raised them), so they would never come down again.
			double groupPredicted2[2] = {};
			unsigned int groupLaunches[2] = {};
			for (size_t i = 0; i < launches.size(); ++i)
			{
				const cl_uint flag = stats[SieveStatsHeader + i] & 3;
				const cl_uint launchSurvivors = stats[SieveStatsHeader + i] >> 2;
				num64 positions = 0;
				double expected = 0.0;
				for (unsigned int c = launches[i].first; c < launches[i].second; ++c)
				{
					positions += chunkPositions(c);
					expected += chunkSurvivors(c);
				}
				const unsigned int group = (launches[i].second > launches[i].first) ? chunkGroup(launches[i].first) : 1;
				++groupLaunches[group];
				groupChecked[group] += std::min<double>(launchSurvivors, mySurvivorCapacity);
				if (flag)
				{
					todo.push_back(launches[i]);
					++numOverflows[flag];
					if (flag == 1)
					{
						minOverflowPositions = std::min(minOverflowPositions, std::max<num64>(positions, 1));
						groupOverflowScale[group] = std::max(groupOverflowScale[group], 2.0 * static_cast<double>(launchSurvivors) / std::max(expected, 1.0));
					}
					else
					{
						// (phase 2 overflowed: all survivors of the launch were counted and checked, so they measure the survivor estimate.
						// Without them a dense range kept its estimate and the factor for phase 2 doubled instead, 13 times in one batch,
						// each time sieving the same launches again)
						groupMeasured[group] += launchSurvivors;
						groupExpected[group] += expected;
					}
				}
				else
				{
					donePositions += positions;
					groupMeasured[group] += launchSurvivors;
					groupExpected[group] += expected;
					if (launches[i].second > launches[i].first)
					{
						groupPredicted2[group] += chunkPhase2(launches[i].first, expected * chunkScale(launches[i].first));
					}
				}
			}

			// Measured survivors relative to the estimate: after an overflow at least twice what overflowed the buffer
			if (numOverflows[1])
			{
				LOG_ERROR("Sieve: survivor buffer overflow in " << numOverflows[1] << " of " << launches.size() << " launches, running them again");
					P().mySieveSurvivorRate = std::max(P().mySieveSurvivorRate * 2.0, 2.0 * static_cast<double>(mySurvivorCapacity) / static_cast<double>(minOverflowPositions));
					++P().mySieveStatOverflows;
				}
				else if ((donePositions >= (num64(1) << 24)) || (static_cast<double>(donePositions) * P().mySieveSurvivorRate >= 1e5))
				{
					const double measured = static_cast<double>(roundMeasured) / static_cast<double>(donePositions);
					P().mySieveSurvivorRate = std::max(std::max(measured, P().mySieveSurvivorRate * 0.5), 1e-6);
				}
				for (unsigned int group = 0; group < 2; ++group)
				{
					if (groupOverflowScale[group] > 0.0)
					{
						// (the counter of an overflowed launch holds all its survivors)
						P().mySieveRateScale[group] = std::max(P().mySieveRateScale[group] * 2.0, groupOverflowScale[group]);
					}
					else if ((groupExpected[group] > 1e5) || (groupPredicted2[group] >= 1e5))
					{
						P().mySieveRateScale[group] = std::max(std::max(groupMeasured[group] / std::max(groupExpected[group], 1e-30), P().mySieveRateScale[group] * 0.5), 1e-3);
					}
				}

			// Numbers for phase 2 (shared by the launches of a round): process them when the buffer is half full or overflowed (the
			// overflowed launches run again)
			cl_uint phase2Count = 0;
			CL_CHECKED_CALL(clEnqueueReadBuffer, P().myQueue, P().myPhase2_numbers_count_buf, CL_TRUE, 0, sizeof(phase2Count), &phase2Count, 0, nullptr, nullptr);
			const cl_uint phase2New = std::min<cl_uint>(phase2Count, myPhase2MaxNumbersCount) - std::min<cl_uint>(phase2Fill, myPhase2MaxNumbersCount);
			if ((roundPositions >= 16777216.0) || ((roundPositions > 0.0) && (roundPositions * P().mySievePhase2Rate >= 1e5)))
			{
				const double measured = static_cast<double>(phase2New) / roundPositions;
				P().mySievePhase2Rate = numOverflows[2] ? std::max(P().mySievePhase2Rate * 2.0, measured * 2.0) : std::max(std::max(measured, P().mySievePhase2Rate * 0.5), 1e-9);
			}
			else if (numOverflows[2])
			{
				P().mySievePhase2Rate *= 2.0;
			}
			// phase 2 numbers per checked survivor: measured for the group of the round (both groups if the round had both)
			{
				const double checked = groupChecked[0] + groupChecked[1];
				const double measured = (checked > 0.0) ? (static_cast<double>(phase2New) / checked) : 0.0;
				for (unsigned int group = 0; group < 2; ++group)
				{
					// (only the groups of the round's launches)
					if ((groupLaunches[group] == 0) || ((groupChecked[group] <= 0.0) && (groupPredicted2[group] < 1e5) && !numOverflows[2]))
					{
						continue;
					}
					double& ratio = P().mySievePhase2PerSurvivor[group];
					if (numOverflows[2])
					{
						ratio = std::max(ratio * 2.0, measured * 2.0);
					}
					else if ((checked >= 1e5) || (groupPredicted2[group] >= 1e5))
					{
						ratio = std::max(std::max(measured, ratio * 0.5), 1e-4);
					}
				}
			}
			if (numOverflows[2])
			{
				++P().mySieveStatOverflows;
			}
			phase2Fill = std::min<cl_uint>(phase2Count, myPhase2MaxNumbersCount);
			if (numOverflows[2] || (phase2Fill > myPhase2MaxNumbersCount / 2))
			{
				P().mySieveStatPhase2 += phase2Fill;
				if (!ProcessNumbersPhases2_3(phase2Fill))
				{
					return false;
				}
				phase2Fill = 0;
				P().mySieveRoundPositions2 = 0;
			}
			else
			{
				P().mySieveRoundPositions2 += donePositions;
			}
			todo.insert(todo.end(), later.begin(), later.end());
			firstRound = false;
		}

	// (other callers of the checks have no band table)
	if (P().mySieveBand)
	{
		for (cl_kernel kc : { P().myCheckSieveSurvivors, P().myCheckSieveSurvivorsPrime })
		{
			CL_CHECKED_CALL(setKernelArgumentsWithOffset, kc, 21, 0U, 0U);
		}
	}

	// /sieve_profile: GPU time of the sieve, of the survivor checks (the first pass, SieveTestDeferred and the pass with a prime p: all three
	// in P().mySieveStatCheckKernelTime) and of the band sieve
	const std::vector<cl_event>* eventLists[5] = { &sieveEvents, &checkEvents, &bandEvents, &deferEvents, &check2Events };
	double* eventTimes[5] = { &P().mySieveStatSieveKernelTime, &P().mySieveStatCheckKernelTime, &P().mySieveStatBandKernelTime, &P().mySieveStatCheckTimes[1], &P().mySieveStatCheckTimes[2] };
	for (int k = 0; k < 5; ++k)
	{
		for (cl_event e : *eventLists[k])
		{
			cl_ulong t0 = 0, t1 = 0;
			if ((clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_START, sizeof(t0), &t0, nullptr) == CL_SUCCESS) &&
				(clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_END, sizeof(t1), &t1, nullptr) == CL_SUCCESS))
			{
				*eventTimes[k] += (t1 - t0) * 1e-9;
				if (k == 1)
				{
					P().mySieveStatCheckTimes[0] += (t1 - t0) * 1e-9;
				}
				else if (k >= 3)
				{
					P().mySieveStatCheckKernelTime += (t1 - t0) * 1e-9;
				}
			}
			clReleaseEvent(e);
		}
	}
	return true;
}

// Finishes phases 2 and 3 for the numbers left by the survivor checks. Returns false in "redo" if phase 2 overflowed (some numbers were
// not saved, so the caller must run the batch again).
bool OpenCL::FinishSieveBatch(bool& redo)
{
	redo = false;
	P().mySieveRoundPositions2 = 0;
	unsigned int counters[2] = {};
	CL_CHECKED_CALL(clEnqueueReadBuffer, P().myQueue, P().myPhase2_numbers_count_buf, CL_TRUE, 0, sizeof(counters), counters, 0, nullptr, nullptr);
	unsigned int n = counters[0];
	P().mySieveStatPhase2 += n;
	if (n > myPhase2MaxNumbersCount)
	{
		redo = true;
		n = myPhase2MaxNumbersCount;
	}
	if (n > 0)
	{
		if (!ProcessNumbersPhases2_3(n))
		{
			return false;
		}
	}
	else if (!ResetCounter(P().myQueue, P().myPhase2_numbers_count_buf))
	{
		return false;
	}
	// (the progress of the batch is saved after this)
	return FlushPhase3Backlog();
}

bool OpenCL::ProcessSieveBatch(SieveBatch& batch)
{
	const unsigned int numChunks = static_cast<unsigned int>(batch.chunks.size() / 4);
	if ((numChunks == 0) && batch.extra.empty())
	{
		return true;
	}

	if (g_SieveDryRun)
	{
		// profiling of the host part: no GPU work
		P().mySieveStatPositions += batch.positions;
		P().mySieveStatChunks += numChunks;
		++P().mySieveStatBatches;
		batch.ranges.clear();
		batch.chunks.clear();
		batch.extra.clear();
		batch.positions = 0;
		return true;
	}

	Timer t;

	// Development aid: AMICABLE_DUMP_BATCH=<file> writes the sieve tables and batch number AMICABLE_DUMP_BATCH_INDEX (default 0) to a
	// file, for running the kernels outside the app (profilers without OpenCL support)
	if (const char* dumpName = getenv("AMICABLE_DUMP_BATCH"))
	{
		const char* idx = getenv("AMICABLE_DUMP_BATCH_INDEX");
		if (P().mySieveStatBatches == (idx ? StrToNumber(idx) : 0))
		{
			FILE* f = fopen(dumpName, "wb");
			auto put = [f](const void* data, size_t elementSize, size_t count)
			{
				const num64 n = count;
				fwrite(&n, sizeof(n), 1, f);
				fwrite(&elementSize, sizeof(num64), 1, f);
				fwrite(data, elementSize, count, f);
			};
			put(mySieve.primeRecs.data(), 4, mySieve.primeRecs.size());
			put(mySieve.modRecs.data(), 4, mySieve.modRecs.size());
			put(mySieve.modRecs2.data(), 4, mySieve.modRecs2.size());
			put(mySieve.Ttab.data(), 4, mySieve.Ttab.size());
			put(mySieve.tasks.data(), 4, mySieve.tasks.size());
			put(mySieve.tierInfo.data(), 4, mySieve.tierInfo.size());
			put(mySieve.wheelTabs.data(), 4, mySieve.wheelTabs.size());
			put(mySieve.modTaskIdx.data(), 4, mySieve.modTaskIdx.size());
			put(mySieve.modTaskList.data(), 4, mySieve.modTaskList.size());
			put(batch.ranges.data(), 8, batch.ranges.size());
			put(batch.chunks.data(), 4, batch.chunks.size());
			put(myBuildOptions.c_str(), 1, myBuildOptions.size());
			// tables of CheckSieveSurvivors (read back from the GPU)
			auto putBuf = [&](cl_mem b)
			{
				size_t size = 0;
				if (b)
				{
					clGetMemObjectInfo(b, CL_MEM_SIZE, sizeof(size), &size, nullptr);
				}
				std::vector<unsigned char> h(size);
				if (size)
				{
					clEnqueueReadBuffer(P().myQueue, b, CL_TRUE, 0, size, h.data(), 0, nullptr, nullptr);
				}
				put(h.data(), 1, size);
			};
			putBuf(mySmallPrimesBuf);
			putBuf(myPrimeInversesBuf);
			putBuf(myPQ_Buf);
			putBuf(myPowerOf2SumInverses_Buf);
			putBuf(myPowersOfP_128SumInverses_offsets_Buf);
			putBuf(myPowersOfP_128SumInverses_Buf);
			putBuf(myPrimeInverses_128_Buf);
			putBuf(mySumEstimates_128_Buf);
			putBuf(myPrimeReciprocalsBuf);
			fclose(f);
			std::cerr << "Dumped sieve batch " << P().mySieveStatBatches << " (" << numChunks << " chunks) to " << dumpName << std::endl;
		}
	}

	// Chunks of the last tier first (they use the sieve with the kill), and in both parts the chunks of the two groups of the launch
	// estimates (W class < 2 and the others, RunSieveLaunches: chunkGroup) together: a launch takes chunks of one kind only, and the
	// groups alternated (3^10 unit 20: 3766 launches of 19 chunks per batch, a few us between two launches). The order is kept within
	// each kind (the band entries of a launch must be sorted by range).
	unsigned int killChunks = 0;
	{
		std::vector<unsigned int> reordered;
		reordered.reserve(batch.chunks.size());
		for (int pass = 0; pass < 4; ++pass)
		{
			const bool killPass = (pass < 2);
			const bool group0Pass = ((pass & 1) == 0);
			if (!killPass || P().mySieveRangesKill)
			{
				for (size_t c = 0; c < batch.chunks.size(); c += 4)
				{
					const unsigned int tier = static_cast<unsigned int>(batch.ranges[static_cast<size_t>(batch.chunks[c]) * SieveRangeStride + 12] >> 32) & 0xFFU;
					const bool kill = P().mySieveRangesKill && (tier + 1 == mySieve.cfg.numTiers);
					const bool group0 = (batch.chunks[c + 3] < 2);
					if ((kill == killPass) && (group0 == group0Pass))
					{
						reordered.insert(reordered.end(), batch.chunks.begin() + static_cast<std::ptrdiff_t>(c), batch.chunks.begin() + static_cast<std::ptrdiff_t>(c + 4));
					}
				}
			}
			if (pass == 1)
			{
				killChunks = static_cast<unsigned int>(reordered.size() / 4);
			}
		}
		batch.chunks.swap(reordered);
	}

	if (!batch.ranges.empty())
	{
		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().mySieveRangesBuf, CL_FALSE, 0, batch.ranges.size() * sizeof(num64), batch.ranges.data(), 0, nullptr, nullptr);
	}
	if (numChunks > 0)
	{
		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().mySieveChunksBuf, CL_FALSE, 0, batch.chunks.size() * sizeof(unsigned int), batch.chunks.data(), 0, nullptr, nullptr);
	}

	num64 survivors = 0;
	// (phase 3: RunPhase3, FinishSieveBatch processes the backlog)
	P().myPhase3BacklogEnabled = (P().myPhase3BacklogBuf != nullptr);
	struct DisableBacklog
	{
		bool& enabled;
		~DisableBacklog() { enabled = false; }
	} disableBacklog{ P().myPhase3BacklogEnabled };
	for (;;)
	{
		// p = 2, 3, 5 are checked directly (they are not in the sieve)
		const unsigned int numExtra = static_cast<unsigned int>(batch.extra.size());
		if (numExtra > 0)
		{
			CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().mySurvivorsBuf, CL_FALSE, 0, numExtra * sizeof(num64), batch.extra.data(), 0, nullptr, nullptr);
		}

		bool overflow = false;
		if (!RunSieveLaunches(0, numChunks, killChunks, 0, numExtra, &batch, survivors, overflow))
		{
			return false;
		}

		bool redo = false;
		if (!FinishSieveBatch(redo))
		{
			return false;
		}

		if (!overflow && !redo)
		{
			break;
		}

		// Too many survivors in one launch (or too many numbers for phase 2): run the whole batch again with shorter launches.
		// Numbers found twice are removed at the end.
		LOG_ERROR("Sieve batch overflow, running it again" << (redo ? " (phase 2)" : ""));
		P().mySieveSurvivorRate *= 4.0;
	}

	P().mySieveStatPositions += batch.positions;
	P().mySieveStatSurvivors += survivors;
	P().mySieveStatChunks += numChunks;
	++P().mySieveStatBatches;

	batch.ranges.clear();
	batch.chunks.clear();
	batch.extra.clear();
	batch.positions = 0;

	P().mySieveStatTime += t.getElapsedTime();
	return true;
}

// Waits for the batch of pipeline i (main thread)
bool OpenCL::WaitPipe(int i)
{
	SievePipe& p = myPipes[i];
	if (!p.mySieveWorkerRunning)
	{
		return true;
	}
	{
		std::unique_lock<std::mutex> lock(p.mySieveWorkerMutex);
		p.mySieveWorkerCV.wait(lock, [&p]() { return p.mySieveJobDone; });
	}
	p.mySieveWorkerRunning = false;
	return p.mySieveWorkerResult;
}

// Waits for all batches on the GPU
bool OpenCL::WaitSieveBatch()
{
	bool result = true;
	while (!myPipeOrder.empty())
	{
		const int i = myPipeOrder.front();
		myPipeOrder.pop_front();
		result = WaitPipe(i) && result;
	}
	return result;
}

void OpenCL::SieveWorkerLoop(int pipe)
{
	ourPipe = pipe;
	SievePipe& p = myPipes[pipe];
	std::unique_lock<std::mutex> lock(p.mySieveWorkerMutex);
	for (;;)
	{
		p.mySieveWorkerCV.wait(lock, [&p]() { return p.mySieveJobPending || p.mySieveWorkerStop; });
		if (!p.mySieveJobPending)
		{
			// (stop: only between batches, a submitted batch is always finished first)
			return;
		}
		p.mySieveJobPending = false;
		lock.unlock();
		const bool result = ProcessSieveBatch(p.mySieveBatchRun);
		lock.lock();
		p.mySieveWorkerResult = result;
		p.mySieveJobDone = true;
		p.mySieveWorkerCV.notify_all();
	}
}

bool OpenCL::SubmitSieveBatch(bool checkpoint, num64 offset)
{
	// a free pipeline, else the one of the oldest batch when that batch is done
	int pipe = -1;
	for (int i = 0; i < myNumPipes; ++i)
	{
		if (!myPipes[i].mySieveWorkerRunning)
		{
			pipe = i;
			break;
		}
	}
	if (pipe < 0)
	{
		pipe = myPipeOrder.front();
		myPipeOrder.pop_front();
		if (!WaitPipe(pipe))
		{
			return false;
		}
		// All ranges before the first range of the oldest batch still on the GPU (or of the submitted batch) are done: the found numbers
		// of the finished batch are saved, with a checkpoint at that range. (The batches are finished in the order they were submitted,
		// so the found numbers in output.txt are those of the ranges before the checkpoint, as with one pipeline.)
		if (checkpoint)
		{
			const SievePipe* next = myPipeOrder.empty() ? nullptr : &myPipes[myPipeOrder.front()];
			const SieveRangeKey* r = next ? (next->myRunHasFirst ? &next->myRunFirstRange : nullptr) : (mySieveBatchHasFirst ? &mySieveBatchFirstRange : nullptr);
			PipeScope scope(pipe);
			if (r && !SaveProgressForRanges(*r, offset, next ? next->myRunStartNumbers : mySieveBatchStartNumbers))
			{
				return false;
			}
		}
	}

	SievePipe& p = myPipes[pipe];
	p.myRunHasFirst = mySieveBatchHasFirst;
	if (p.myRunHasFirst)
	{
		p.myRunFirstRange = mySieveBatchFirstRange;
		p.myRunStartNumbers = mySieveBatchStartNumbers;
	}
	mySieveBatchHasFirst = false;
	std::swap(p.mySieveBatchRun, mySieveBatch);
	if (!p.mySieveWorker.joinable())
	{
		p.mySieveWorker = std::thread(&OpenCL::SieveWorkerLoop, this, pipe);
	}
	{
		std::lock_guard<std::mutex> lock(p.mySieveWorkerMutex);
		p.mySieveJobPending = true;
		p.mySieveJobDone = false;
	}
	p.mySieveWorkerCV.notify_all();
	p.mySieveWorkerRunning = true;
	myPipeOrder.push_back(pipe);
	return true;
}

//-----------------------------------------------------------------------------------------------------------------------------------
// Sieve self-tests ("/test")
//-----------------------------------------------------------------------------------------------------------------------------------

bool OpenCL::RunSieveTests()
{
	if (!myUseSieve)
	{
		LOG_ERROR("Sieve tests need the sieve to be enabled");
		return false;
	}

	Timer t;
	std::vector<KnownPair> pairs;
	std::vector<num128> allAmicable;
	const unsigned int stride = g_TestStride;
	if (!LoadKnownPairs(stride, pairs, &allAmicable))
	{
		LOG_ERROR("No known pairs loaded (c2_*.txt files are needed)");
		return false;
	}

	// Necessary conditions of the ranges (host only); "/test_reduction" runs only this test
	bool ok = TestSieveReduction(mySieve, pairs);
	std::cerr << "Reduction test: " << t.getElapsedTime() << " seconds" << std::endl;
	if (g_TestReductionOnly || !ok)
	{
		std::cerr << (ok ? "All sieve tests passed" : "Sieve tests FAILED") << " (" << t.getElapsedTime() << " seconds)" << std::endl;
		return ok;
	}

	// All known pairs below the search limit (not only the first stage) through the kernels for phases 1-3
	{
		Timer t1;
		std::vector<SPairToCheck> batch;
		num64 numChecked = 0;
		std::vector<KnownPair> all;
		std::vector<num128> allM;
		const char* files[] = { "c2_3.txt", "c2_4.txt", "c2_5.txt", "c2_6.txt", "c2_7.txt", "c2_8.txt", "c2_9.txt", "c2_10.txt", "c2_11.txt", "c2_12.txt", "c2_13.txt", "c2_14.txt", "c2_15.txt", "c2_16.txt", "c2_17.txt", "c2_18.txt", "c2_19.txt", "c2_20.txt", "c2_20_1.txt", "c2_20_2.txt", "c2_21.txt", "c2_21_1.txt", "c2_21_2.txt", "c2_21_3.txt", "c2_21_4.txt", "c2_22.txt" };
		num64 k = 0;
		for (const char* name : files)
		{
			std::ifstream f(name);
			std::string line1, line2, line3, line4;
			while (ok && std::getline(f, line1) && std::getline(f, line2) && std::getline(f, line3))
			{
				std::getline(f, line4);
				const num128 m = atoi128(line2.c_str());
				if (m >= SearchLimit::value)
				{
					break;
				}
				if ((k++ % stride) != 0)
				{
					continue;
				}
				batch.emplace_back(m, m + atoi128(line3.c_str()));
				if (batch.size() >= (1U << 20))
				{
					std::sort(batch.begin(), batch.end(), [](const SPairToCheck& a, const SPairToCheck& b) { return a.M < b.M; });
					numChecked += batch.size();
					ok = TestPairs(batch);
					batch.clear();
				}
			}
		}
		if (ok && !batch.empty())
		{
			std::sort(batch.begin(), batch.end(), [](const SPairToCheck& a, const SPairToCheck& b) { return a.M < b.M; });
			numChecked += batch.size();
			ok = TestPairs(batch);
		}
		std::cerr << "Phases 1-3 test: " << numChecked << " pairs " << (ok ? "found" : "NOT found") << " (" << t1.getElapsedTime() << " seconds)" << std::endl;
	}

	if (myLargePrimesSieve)
	{
		Timer t1;
		ok = ok && TestSieveLargePrimes(stride, allAmicable);
		std::cerr << "Large primes test: " << t1.getElapsedTime() << " seconds" << std::endl;
		std::cerr << (ok ? "All sieve tests passed" : "Sieve tests FAILED") << " (" << t.getElapsedTime() << " seconds)" << std::endl;
		return ok;
	}

	ok = ok && TestSieveGPUWindow(mySieve, pairs);
	std::cerr << "Window test: " << t.getElapsedTime() << " seconds" << std::endl;

	if (ok)
	{
		Timer t1;
		ok = TestSieveSums(pairs);
		std::cerr << "Sum test: " << t1.getElapsedTime() << " seconds" << std::endl;
	}

	if (ok)
	{
		Timer t1;
		ok = TestSieveKnownPairs(pairs, allAmicable);
		std::cerr << "Known pairs test: " << t1.getElapsedTime() << " seconds" << std::endl;
	}

	std::cerr << (ok ? "All sieve tests passed" : "Sieve tests FAILED") << " (" << t.getElapsedTime() << " seconds)" << std::endl;
	return ok;
}

// Runs the sieve with "-D SIEVE_DUMP" on short ranges around known pairs and compares the sum of every position with the host
bool OpenCL::TestSieveSums(const std::vector<KnownPair>& pairs)
{
	struct TmpProgram
	{
		~TmpProgram() { if (k) clReleaseKernel(k); if (p) clReleaseProgram(p); if (dump) clReleaseMemObject(dump); }
		cl_program p = nullptr;
		cl_kernel k = nullptr;
		cl_mem dump = nullptr;
	} tmp;

	const char* sources[2] = { kernel_cl, kernel_sieve_cl };
	CL_CHECKED_CALL_WITH_RESULT(tmp.p, clCreateProgramWithSource, myGPUContext, 2, sources, nullptr);
	const std::string options = myBuildOptions + " -D SIEVE_DUMP";
	if (clBuildProgram(tmp.p, 0, nullptr, options.c_str(), nullptr, nullptr) != CL_SUCCESS)
	{
		LOG_ERROR("Failed to build the sieve test program");
		return false;
	}
	CL_CHECKED_CALL_WITH_RESULT(tmp.k, clCreateKernel, tmp.p, "SieveRanges");
	// the sieve with the kill (ranges of the last tier)
	struct TmpKernel
	{
		~TmpKernel() { if (k) clReleaseKernel(k); if (p) clReleaseProgram(p); }
		cl_program p = nullptr;
		cl_kernel k = nullptr;
	} tmpKill;
	if (!myBuildOptionsKill.empty())
	{
		CL_CHECKED_CALL_WITH_RESULT(tmpKill.p, clCreateProgramWithSource, myGPUContext, 2, sources, nullptr);
		const std::string killOptions = myBuildOptionsKill + " -D SIEVE_DUMP";
		if (clBuildProgram(tmpKill.p, 0, nullptr, killOptions.c_str(), nullptr, nullptr) != CL_SUCCESS)
		{
			LOG_ERROR("Failed to build the sieve test program (kill)");
			return false;
		}
		CL_CHECKED_CALL_WITH_RESULT(tmpKill.k, clCreateKernel, tmpKill.p, "SieveRanges");
	}

	// one segment per test range: the sums of all 6 sub-sieves (kernel_sieve.cl, SIEVE_DUMP)
	const unsigned int segPos = mySieve.cfg.segJ * SievePositionsPerJ;
	const unsigned int maxTestRanges = 1024;
	CL_CHECKED_CALL_WITH_RESULT(tmp.dump, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, (static_cast<size_t>(maxTestRanges) * segPos + maxTestRanges * 2) * sizeof(int), 0);

	const unsigned int zero = 0;
	const unsigned int capacity = 0;
	for (cl_kernel kernel : { tmp.k, tmpKill.k })
	{
		if (kernel)
		{
			CL_CHECKED_CALL(setKernelArguments, kernel, P().mySieveRangesBuf, P().mySieveChunksBuf, zero, mySievePrimeRecsBuf, mySieveModRecsBuf, mySieveModRecs2Buf, mySieveTtabBuf, mySieveTasksBuf, mySieveTierInfoBuf, mySieveWheelTabsBuf, mySieveModTaskIdxBuf, mySieveModTaskListBuf, P().mySurvivorCountBuf, tmp.dump, capacity);
		}
	}

	// one segment per test range: p in [P0, P0 + 210 * segJ) around the known p. Every 16th range fills the segment, the others have the
	// width of a segment of the 30-wheel (the last segment of a chunk ends at p_hi: SieveRanges sieves only its first j there; the check
	// of every position on the host is the slow part of the test)
	const num64 halfWidthFull = static_cast<num64>(mySieve.cfg.segJ) * (SieveWheel / 2) - SieveWheel;
	const num64 halfWidthNarrow = static_cast<num64>(mySieve.cfg.segJ) * 15 - 30;

	// a sample of the pairs, plus up to 1024 pairs where a prime factor of g has wheel moduli (their increments are overridden)
	std::vector<size_t> selected;
	{
		const size_t step = std::max<size_t>(1, pairs.size() / 4096);
		for (size_t i = 0; i < pairs.size(); i += step)
		{
			selected.push_back(i);
		}
		size_t numOverrides = 0;
		for (size_t i = 0; (i < pairs.size()) && (numOverrides < 1024); ++i)
		{
			const KnownPair& pair = pairs[i];
			num64 rec[SieveRangeStride];
			if (!mySieve.PrepareRange(pair.M0, pair.sumM0, pair.primes, pair.exponents, pair.numFactors, pair.p, pair.p, mySieve.cfg.numTiers - 1, false, rec))
			{
				continue;
			}
			const unsigned int ngf = static_cast<unsigned int>(rec[12] >> 56);
			for (unsigned int f = 0; f < ngf; ++f)
			{
				const unsigned int q = static_cast<unsigned int>(rec[13 + f / 2] >> ((f & 1) * 32)) & 0xFFFFFFU;
				if (q * 4 >= mySieve.cfg.segJ * 3)
				{
					selected.push_back(i);
					++numOverrides;
					break;
				}
			}
		}
		std::cerr << "TestSieveSums: " << numOverrides << " ranges with overridden increments" << std::endl;
	}

	num64 numPositions = 0;
	num64 numRanges = 0;
	std::vector<int> dump;
	std::vector<num64> records;
	std::vector<unsigned int> chunks;
	std::vector<std::pair<num64, num64>> limits;
	for (size_t first = 0; first < selected.size(); first += maxTestRanges)
	{
		records.clear();
		chunks.clear();
		limits.clear();
		for (size_t k = first; (k < selected.size()) && (k < first + maxTestRanges); ++k)
		{
			const size_t i = selected[k];
			const KnownPair& pair = pairs[i];
			const num64 halfWidth = ((k % 16) == 0) ? halfWidthFull : halfWidthNarrow;
			num64 p_lo = (pair.p > halfWidth) ? (pair.p - halfWidth) : 2;
			if (p_lo <= pair.primes[pair.numFactors - 1])
			{
				p_lo = pair.primes[pair.numFactors - 1] + 1;
			}
			const num64 p_hi = p_lo + halfWidth * 2;
			num64 rec[SieveRangeStride];
			if (!mySieve.PrepareRange(pair.M0, pair.sumM0, pair.primes, pair.exponents, pair.numFactors, p_lo, p_hi, static_cast<unsigned int>(i % mySieve.cfg.numTiers), false, rec))
			{
				LOG_ERROR("PrepareRange failed for " << pair.m);
				return false;
			}
			const unsigned int ri = static_cast<unsigned int>(records.size() / SieveRangeStride);
			records.insert(records.end(), rec, rec + SieveRangeStride);
			chunks.push_back(ri);
			chunks.push_back(0);
			chunks.push_back(1);
			chunks.push_back(0);
			limits.emplace_back(p_lo, p_hi);
		}

		const unsigned int n = static_cast<unsigned int>(limits.size());
		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().mySieveRangesBuf, CL_TRUE, 0, records.size() * sizeof(num64), records.data(), 0, nullptr, nullptr);
		CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().mySieveChunksBuf, CL_TRUE, 0, chunks.size() * sizeof(unsigned int), chunks.data(), 0, nullptr, nullptr);
		const size_t localSize = mySieve.cfg.wgSize;
		const size_t globalSize = static_cast<size_t>(n) * localSize;
		CL_CHECKED_CALL(clEnqueueNDRangeKernel, P().myQueue, tmp.k, 1, nullptr, &globalSize, &localSize, 0, nullptr, nullptr);
		CL_CHECKED_CALL(clFinish, P().myQueue);
		dump.resize(static_cast<size_t>(n) * segPos + n * 2);
		CL_CHECKED_CALL(clEnqueueReadBuffer, P().myQueue, tmp.dump, CL_TRUE, 0, dump.size() * sizeof(int), dump.data(), 0, nullptr, nullptr);
		// ranges of the last tier: the sums of the sieve with the kill
		if (tmpKill.k)
		{
			std::vector<int> dumpKill(dump.size());
			CL_CHECKED_CALL(clEnqueueNDRangeKernel, P().myQueue, tmpKill.k, 1, nullptr, &globalSize, &localSize, 0, nullptr, nullptr);
			CL_CHECKED_CALL(clFinish, P().myQueue);
			CL_CHECKED_CALL(clEnqueueReadBuffer, P().myQueue, tmp.dump, CL_TRUE, 0, dumpKill.size() * sizeof(int), dumpKill.data(), 0, nullptr, nullptr);
			for (unsigned int r = 0; r < n; ++r)
			{
				const unsigned int tier = static_cast<unsigned int>(records[static_cast<size_t>(r) * SieveRangeStride + 12] >> 32) & 0xFFU;
				if (tier + 1 == mySieve.cfg.numTiers)
				{
					std::copy(dumpKill.begin() + static_cast<std::ptrdiff_t>(r) * segPos, dumpKill.begin() + static_cast<std::ptrdiff_t>(r + 1) * segPos, dump.begin() + static_cast<std::ptrdiff_t>(r) * segPos);
				}
			}
		}

		for (unsigned int r = 0; r < n; ++r)
		{
			const num64* rec = records.data() + static_cast<size_t>(r) * SieveRangeStride;
			const num64 P0 = rec[8];
			for (unsigned int idx = 0; idx < segPos; idx += 1)
			{
				// sub-sieve s7 = p mod 7, row c: p = P0 + 210 j + w(c) with w(c) == W30[c] (mod 30), w(c) == s7 (mod 7)
				static const unsigned int W30[8] = { 1, 7, 11, 13, 17, 19, 23, 29 };
				const unsigned int s7 = idx / (mySieve.cfg.segJ * 8) + 1;
				const unsigned int c = (idx / mySieve.cfg.segJ) % 8;
				const unsigned int j = idx % mySieve.cfg.segJ;
				const num64 p = P0 + SieveWheel * j + W30[c] + 30 * (((s7 + 7 - W30[c] % 7) * 4) % 7);
				if ((p < limits[r].first) || (p > limits[r].second))
				{
					continue;
				}
				++numPositions;
				const int expected = mySieve.ExpectedSum(rec, p);
				const int actual = dump[static_cast<size_t>(r) * segPos + idx];
				if (actual != expected)
				{
					LOG_ERROR("Sieve sum mismatch: M0 = " << CombineNum128(rec[4], rec[5]) << ", p = " << p << ", tier " << ((rec[12] >> 32) & 0xFF) << ", GPU " << actual << ", expected " << expected);
					return false;
				}
			}

			// window of the segment: must contain the window of every p in it
			const int lo = dump[static_cast<size_t>(n) * segPos + r * 2];
			const int hi = dump[static_cast<size_t>(n) * segPos + r * 2 + 1];
			for (num64 p : { limits[r].first, limits[r].second })
			{
				const int lnr_fx = static_cast<int>(static_cast<unsigned int>(rec[11]));
				float rm1;
				const unsigned int bits = static_cast<unsigned int>(rec[12]);
				memcpy(&rm1, &bits, sizeof(rm1));
				const double corr = log1p(static_cast<double>(rm1) / (static_cast<double>(p) + 1.0)) * SIEVE_FX;
				const int exact = lnr_fx - static_cast<int>(corr);
				if ((exact - static_cast<int>(static_cast<unsigned int>(rec[11] >> 32)) + 64 < lo) || (exact - 64 > hi))
				{
					LOG_ERROR("Sieve window mismatch: M0 = " << CombineNum128(rec[4], rec[5]) << ", p = " << p << ", window [" << lo << ", " << hi << "], exact log(T/N) " << exact);
					return false;
				}
			}
		}
		numRanges += n;
	}

	std::cerr << "TestSieveSums: " << numRanges << " ranges, " << numPositions << " positions, all sums are correct" << std::endl;
	return numPositions > 0;
}

// Runs the complete sieve pipeline on short ranges around known pairs, all of them must be found
bool OpenCL::TestSieveKnownPairs(const std::vector<KnownPair>& pairs, const std::vector<num128>& allAmicable, bool verbose)
{
	const num64 halfWidth = 50000;
	const size_t batchSize = 2048;
	num64 numFound = 0;
	std::vector<SFoundPair> found;
	std::vector<num128> foundM;
	for (size_t first = 0; first < pairs.size(); first += batchSize)
	{
		const size_t last = std::min(pairs.size(), first + batchSize);
		for (size_t i = first; i < last; ++i)
		{
			const KnownPair& pair = pairs[i];
			num64 p_lo = (pair.p > halfWidth) ? (pair.p - halfWidth) : 2;
			if (p_lo <= pair.primes[pair.numFactors - 1])
			{
				p_lo = pair.primes[pair.numFactors - 1] + 1;
			}
			num128 limit = (SearchLimit::value - 1) / pair.M0;
			if (limit > SearchLimit::MainPrimeTableBound)
			{
				limit = SearchLimit::MainPrimeTableBound;
			}
			num64 p_hi = pair.p + halfWidth;
			if (p_hi > LowWord(limit))
			{
				p_hi = LowWord(limit);
			}
			num64 rec[SieveRangeStride];
			const bool prepared = mySieve.PrepareRange(pair.M0, pair.sumM0, pair.primes, pair.exponents, pair.numFactors, p_lo, p_hi, static_cast<unsigned int>(i % mySieve.cfg.numTiers), false, rec);
			if (prepared)
			{
				SetTestBandParameters(rec, pair);
			}
			if (!prepared || !AddSieveRangeRecord(rec, p_lo, p_hi))
			{
				LOG_ERROR("Failed to add a test range for " << pair.m);
				return false;
			}
		}
		if (!ProcessSieveBatch())
		{
			return false;
		}

		unsigned int numbers_found = 0;
		if (!GetAndResetCounter(P().myQueue, P().myAmicable_numbers_count_buf, numbers_found))
		{
			return false;
		}
		if (numbers_found > 8192)
		{
			LOG_ERROR("Too many numbers found");
			return false;
		}
		found.resize(numbers_found);
		if (numbers_found > 0)
		{
			CL_CHECKED_CALL(clEnqueueReadBuffer, P().myQueue, P().myAmicable_numbers_data_buf, CL_TRUE, 0, sizeof(num64) * 4 * numbers_found, found.data(), 0, nullptr, nullptr);
		}
		foundM.clear();
		for (const SFoundPair& k : found)
		{
			if ((k.N == 0) || IsPrime(k.N))
			{
				const num128 M = CombineNum128(k.M, k.M_high);
				if (verbose && !std::binary_search(allAmicable.begin(), allAmicable.end(), M))
				{
					// not in the c2_*.txt files: report it (the lists above 10^21 are not complete)
					std::cerr << "Sieve test found " << M << " which is not in the known pairs, check it" << std::endl;
				}
				foundM.push_back(M);
			}
		}
		std::sort(foundM.begin(), foundM.end());
		for (size_t i = first; i < last; ++i)
		{
			if (!std::binary_search(foundM.begin(), foundM.end(), pairs[i].m))
			{
				LOG_ERROR("Sieve test missed amicable number " << pairs[i].m << " (tier " << (i % mySieve.cfg.numTiers) << ")");
				return false;
			}
			++numFound;
		}
	}

	if (verbose)
	{
		std::cerr << "TestSieveKnownPairs: all " << numFound << " pairs found (" << P().mySieveStatPositions << " positions, " << P().mySieveStatSurvivors << " survivors, " << mySieveStatBandTestRanges << " ranges with the band sieve)" << std::endl;
	}
	return numFound > 0;
}

//-----------------------------------------------------------------------------------------------------------------------------------
// Large primes stage with the sieve: every candidate c (m = c * p) is a range like in the first stage, with M0 = c and p in [P1, P2].
// The GPU derives the range parameters from the candidate itself (there are billions of candidates).
//-----------------------------------------------------------------------------------------------------------------------------------

bool OpenCL::SieveCandidates(const AmicableCandidate* candidates, unsigned int count, num64 P1, num64 P2, unsigned int tier, num64& survivorsTotal)
{
	CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, P().mySieveRangesBuf, CL_FALSE, 0, sizeof(AmicableCandidate) * count, candidates, 0, nullptr, nullptr);

	const num64 P0 = (P1 / SieveWheel) * SieveWheel;
	const num64 numSegments = ((P2 - P0) / SieveWheel) / mySieve.cfg.segJ + 1;
	const num64 chunksPerCandidate = (numSegments + mySieve.cfg.chunkSegments - 1) / mySieve.cfg.chunkSegments;
	if (chunksPerCandidate * count >= (num64(1) << 31))
	{
		LOG_ERROR("Too many chunks per batch");
		return false;
	}

	cl_ulong4 lpParams;
	lpParams.s[0] = P1;
	lpParams.s[1] = P2;
	lpParams.s[2] = P0;
	lpParams.s[3] = chunksPerCandidate | (static_cast<num64>(tier) << 32) | (static_cast<num64>(mySieve.cfg.chunkSegments) << 40);
	CL_CHECKED_CALL(clSetKernelArg, P().mySieveRanges, 1, sizeof(lpParams), &lpParams);
	if ((P1 != myPrimeBitsP1) || (P2 != myPrimeBitsP2))
	{
		LOG_ERROR("The bitmap of primes doesn't match the interval");
		return false;
	}

	const num64 positionsPerChunk = std::min<num64>(static_cast<num64>(mySieve.cfg.chunkSegments), numSegments) * mySieve.cfg.segJ * SievePositionsPerJ;
	for (;;)
	{
		num64 survivors = 0;
		bool overflow = false;
		if (!RunSieveLaunches(0, static_cast<unsigned int>(count * chunksPerCandidate), 0, positionsPerChunk, 0, nullptr, survivors, overflow))
		{
			return false;
		}

		bool redo = false;
		if (!FinishSieveBatch(redo))
		{
			return false;
		}
		if (!overflow && !redo)
		{
			survivorsTotal += survivors;
			break;
		}
		LOG_ERROR("Sieve batch overflow, running it again" << (redo ? " (phase 2)" : ""));
		P().mySieveSurvivorRate *= 4.0;
	}
	return true;
}

// Bitmap of the primes in [P1, P2] (P2 - P1 < 30 * (2^LargePrimesBitmapShift - 8)): row c, bit j <=> P0 + 30 * j + W30[c] is prime,
// P0 = P1 rounded down to a multiple of 210 (as in SieveCandidates and the kernel; the bitmap keeps the layout of the 30-wheel)
bool OpenCL::UploadPrimeBitmap(num64 P1, num64 P2)
{
	const num64 P0 = (P1 / SieveWheel) * SieveWheel;
	const num64 numJ = (P2 - P0) / 30 + 1;
	if ((P1 < 7) || (P2 < P1) || (numJ > (num64(1) << LargePrimesBitmapShift)))
	{
		LOG_ERROR("Invalid interval for the bitmap of primes");
		return false;
	}
	const unsigned int wordsPerRow = static_cast<unsigned int>((numJ + 31) / 32);
	myPrimeBits.assign(static_cast<size_t>(wordsPerRow) * 8, 0);

	struct SetBits
	{
		cl_uint* bits;
		num64 P0;
		unsigned int wordsPerRow;
		FORCEINLINE void operator()(num64 p)
		{
			static const signed char cls[30] = { -1, 0, -1, -1, -1, -1, -1, 1, -1, -1, -1, 2, -1, 3, -1, -1, -1, 4, -1, 5, -1, -1, -1, 6, -1, -1, -1, -1, -1, 7 };
			const num64 d = p - P0;
			const num64 j = d / 30;
			const int c = cls[d - j * 30];
			if (c >= 0)
			{
				bits[static_cast<size_t>(c) * wordsPerRow + (j >> 5)] |= 1U << (j & 31);
			}
		}
	} setBits = { myPrimeBits.data(), P0, wordsPerRow };

	primesieve::PrimeSieve sieve;
	sieve.sieveTemplated(P1, P2, setBits);

	CL_CHECKED_CALL(clEnqueueWriteBuffer, P().myQueue, myPrimeBitsBuf, CL_TRUE, 0, myPrimeBits.size() * sizeof(cl_uint), myPrimeBits.data(), 0, nullptr, nullptr);
	CL_CHECKED_CALL(clSetKernelArg, P().mySieveRanges, 5, sizeof(wordsPerRow), &wordsPerRow);
	myPrimeBitsP1 = P1;
	myPrimeBitsP2 = P2;
	return true;
}

bool OpenCL::RunLargePrimesSieve(num64 startPrime, num64 primeLimit)
{
	if ((primeLimit < startPrime) || (startPrime < 7))
	{
		return SaveFoundNumbers();
	}
	if (startPrime < 11)
	{
		// the sieve has no positions for p = 7 (only p coprime to 210)
		LOG_ERROR("The large primes sieve needs p >= 11");
		return false;
	}

	// The interval of p is split into parts with a bitmap of primes on the GPU for each of them
	// (P0 is up to 209 below P1: the bitmap has 2^LargePrimesBitmapShift j of 30)
	const num64 partSize = ((num64(1) << LargePrimesBitmapShift) - 8) * 30;

	// Resume from a checkpoint: "numbers startPrime primeLimit part_start candidate_offset S#"
	num64 partStart = startPrime;
	num64 start = 0;
	{
		std::string checkpoint_name;
		boinc_resolve_filename_s(CHECKPOINT_LOGICAL_NAME, checkpoint_name);
		FILE* checkpoint = boinc_fopen(checkpoint_name.c_str(), "r");
		if (checkpoint)
		{
			char buf[256];
			if (fgets(buf, sizeof(buf), checkpoint) && strstr(buf, "S#"))
			{
				std::stringstream s;
				s << buf;
				num64 numbers = 0, firstPrime = 0, lastPrime = 0, part = 0, offset = 0;
				s >> numbers >> firstPrime >> lastPrime >> part >> offset;
				if ((firstPrime == startPrime) && (lastPrime == primeLimit) && (part >= startPrime) && (part <= primeLimit))
				{
					myNumbersProcessedTotal = numbers;
					partStart = part;
					start = offset;
				}
			}
			fclose(checkpoint);
		}
	}

	Timer t;
	num64 survivorsTotal = 0;
	for (; partStart <= primeLimit; partStart = (primeLimit - partStart > partSize) ? (partStart + partSize + 1) : (primeLimit + 1), start = 0)
	{
		const num64 P1 = partStart;
		const num64 P2 = (primeLimit - partStart > partSize) ? (partStart + partSize) : primeLimit;

		// candidates c with c * P1 < SearchLimit::value
		const num64 maxC = LowWord((SearchLimit::value - 1) / P1);
		size_t numCandidates;
		{
			size_t a = 0;
			size_t b = CandidatesData.size();
			while (a < b)
			{
				const size_t c = (a + b) / 2;
				if (CandidatesData[c].get_value() <= maxC)
				{
					a = c + 1;
				}
				else
				{
					b = c;
				}
			}
			numCandidates = a;
		}

		const double primesPerCandidate = SieveGPUTables::EstimatePrimes(P1, P2);
		const unsigned int tier = mySieve.ChooseTier(static_cast<num64>(primesPerCandidate));
		const num64 positionsPerCandidate = mySieve.NumPositions(P1, P2);
		std::cerr << "Large primes sieve: p in [" << P1 << ", " << P2 << "], " << numCandidates << " candidates, tier " << tier << " (Q = " << mySieve.cfg.tierQ[tier] << ")" << std::endl;

		if (!UploadPrimeBitmap(P1, P2))
		{
			return false;
		}

		while (start < numCandidates)
		{
			num64 count = std::min<num64>(numCandidates - start, mySieveMaxCandidates);
			const num64 byPositions = std::max<num64>(1, mySieveMaxPositions / std::max<num64>(1, positionsPerCandidate));
			count = std::min(count, byPositions);

			const num64 survivorsBefore = survivorsTotal;
			if (!SieveCandidates(CandidatesData.data() + start, static_cast<unsigned int>(count), P1, P2, tier, survivorsTotal))
			{
				return false;
			}

			// Adapt the size of the launches to the observed number of survivors
			const num64 positions = count * positionsPerCandidate;
			if (positions >= (num64(1) << 24))
			{
				const double rate = static_cast<double>(survivorsTotal - survivorsBefore) / positions;
				P().mySieveSurvivorRate = std::max(std::max(rate, P().mySieveSurvivorRate * 0.9), 1e-5);
			}
			P().mySieveStatPositions += positions;
			P().mySieveStatChunks += count;
			++P().mySieveStatBatches;

			start += count;
			myNumbersProcessedTotal += static_cast<num64>(count * primesPerCandidate + 0.5);
			const double progress = static_cast<double>(myNumbersProcessedTotal) / RangeGen::total_numbers_to_check;
			boinc_fraction_done((progress < 1.0) ? progress : 1.0);

			if (!SaveFoundNumbers())
			{
				return false;
			}
			if (boinc_time_to_checkpoint())
			{
				std::stringstream s;
				std::string checkpoint_name;
				boinc_resolve_filename_s(CHECKPOINT_LOGICAL_NAME, checkpoint_name);
				s << myNumbersProcessedTotal << ' ' << startPrime << ' ' << primeLimit << ' ' << P1 << ' ' << start << " S#";
				const std::string data = s.str();
				FILE* f = boinc_fopen(checkpoint_name.c_str(), "wb");
				if (f)
				{
					fwrite(data.c_str(), sizeof(char), data.length(), f);
					fclose(f);
				}
				boinc_checkpoint_completed();
			}
		}
	}

	P().mySieveStatSurvivors += survivorsTotal;
	P().mySieveStatTime += t.getElapsedTime();
	return SaveFoundNumbers();
}

// Large primes stage: known pairs m = c * p with p > SearchLimit::LinearLimit must be found by the sieve on a short interval around p
// One pair m = c p of the large primes stage: the candidate c alone, p in [p - W, p + W] (W covers several chunks with one segment
// each); allAmicable (if not null): numbers found that are not in it are reported
bool OpenCL::TestSieveLargePrimesPair(const num128 m, num64 c, num64 sc, num64 p, unsigned int tier, const std::vector<num128>* allAmicable, num64& survivorsTotal)
{
	const AmicableCandidate cand(c, sc);
	const num64 W = 200000;
	if (!UploadPrimeBitmap(p - W, p + W) || !SieveCandidates(&cand, 1, p - W, p + W, tier, survivorsTotal))
	{
		return false;
	}

	unsigned int numbers_found = 0;
	if (!GetAndResetCounter(P().myQueue, P().myAmicable_numbers_count_buf, numbers_found))
	{
		return false;
	}
	std::vector<SFoundPair> found(numbers_found);
	if (numbers_found > 0)
	{
		CL_CHECKED_CALL(clEnqueueReadBuffer, P().myQueue, P().myAmicable_numbers_data_buf, CL_TRUE, 0, sizeof(num64) * 4 * numbers_found, found.data(), 0, nullptr, nullptr);
	}
	bool ok = false;
	for (const SFoundPair& k : found)
	{
		if ((k.N == 0) || IsPrime(k.N))
		{
			const num128 M = CombineNum128(k.M, k.M_high);
			if (allAmicable && !std::binary_search(allAmicable->begin(), allAmicable->end(), M))
			{
				// not in the c2_*.txt files: report it (the lists above 10^21 are not complete)
				std::cerr << "Large primes sieve test found " << M << " which is not in the known pairs, check it" << std::endl;
			}
			ok = ok || (M == m);
		}
	}
	if (!ok)
	{
		LOG_ERROR("Large primes sieve test missed " << m << " (tier " << tier << ")");
	}
	return ok;
}

// The sieve part of the startup self-test: the sieve kernels as this device's compiler built them (a compiler or driver that breaks
// them would lose amicable numbers without any error). Sieve runs around known pairs that are part of the app (self_test_sieve.inl,
// no c2 files needed): the first-stage pairs through every tier (the last one with the divisibility kill), some of them only through
// the band sieve (TestSieveKnownPairs puts the band just below their own N/D), numbers of 64 and 128 bits, p above and below 2^32;
// with /lpr the pairs of the large primes stage. All of them must be found. The statistics are reset afterwards.
bool OpenCL::TestSieveStartup()
{
#include "self_test_sieve.inl"
	// (the batches run phase 2: its arguments are set after the self-tests, or by TestPairs if it had pairs)
	if (!SetPhase2KernelArguments())
	{
		return false;
	}
	bool ok = true;
	if (myLargePrimesSieve)
	{
		num64 survivors = 0;
		for (size_t i = 0; ok && (i < ARRAYSIZE(sieve_test_pairs_lp)); ++i)
		{
			KnownPair pair;
			bool bad = false;
			if (!ParseKnownPair(sieve_test_pairs_lp[i][0], sieve_test_pairs_lp[i][1], pair, bad, true) || HighWord(pair.M0) || HighWord(pair.sumM0))
			{
				LOG_ERROR("Bad self-test pair " << sieve_test_pairs_lp[i][0]);
				return false;
			}
			ok = TestSieveLargePrimesPair(pair.m, LowWord(pair.M0), LowWord(pair.sumM0), pair.p, static_cast<unsigned int>(i % mySieve.cfg.numTiers), nullptr, survivors);
		}
	}
	else
	{
		std::vector<KnownPair> pairs;
		for (size_t i = 0; i < ARRAYSIZE(sieve_test_pairs); ++i)
		{
			KnownPair pair;
			bool bad = false;
			if (!ParseKnownPair(sieve_test_pairs[i][0], sieve_test_pairs[i][1], pair, bad))
			{
				LOG_ERROR("Bad self-test pair " << sieve_test_pairs[i][0]);
				return false;
			}
			pairs.push_back(pair);
		}
		const std::vector<num128> none;
		ok = TestSieveKnownPairs(pairs, none, false);
		LOG(1, "Sieve self-test: " << pairs.size() << " pairs, " << mySieveStatBandTestRanges << " of them through the band");
	}

	// (the statistics of the work unit)
	P().mySieveStatPositions = 0;
	P().mySieveStatSurvivors = 0;
	mySieveStatRanges = 0;
	P().mySieveStatChunks = 0;
	P().mySieveStatBatches = 0;
	P().mySieveStatOverflows = 0;
	P().mySieveStatTime = 0.0;
	P().mySieveStatSieveKernelTime = 0.0;
	P().mySieveStatCheckKernelTime = 0.0;
	mySieveStatBandTestRanges = 0;
	return ok;
}

bool OpenCL::TestSieveLargePrimes(unsigned int stride, const std::vector<num128>& allAmicable)
{
	const char* files[] = { "c2_3.txt", "c2_4.txt", "c2_5.txt", "c2_6.txt", "c2_7.txt", "c2_8.txt", "c2_9.txt", "c2_10.txt", "c2_11.txt", "c2_12.txt", "c2_13.txt", "c2_14.txt", "c2_15.txt", "c2_16.txt", "c2_17.txt", "c2_18.txt", "c2_19.txt", "c2_20.txt", "c2_20_1.txt", "c2_20_2.txt", "c2_21.txt", "c2_21_1.txt", "c2_21_2.txt", "c2_21_3.txt", "c2_21_4.txt", "c2_22.txt" };
	num64 numPairs = 0;
	num64 numTested = 0;
	num64 survivorsTotal = 0;
	std::vector<SFoundPair> found;
	for (const char* name : files)
	{
		std::ifstream f(name);
		std::string line1, line2, line3, line4;
		while (std::getline(f, line1) && std::getline(f, line2) && std::getline(f, line3))
		{
			std::getline(f, line4);
			const num128 m = atoi128(line2.c_str());
			if (m >= SearchLimit::value)
			{
				break;
			}

			// largest prime factor of m with exponent 1
			const char* s = strchr(line2.c_str(), '=');
			if (!s)
			{
				continue;
			}
			num64 p = 0;
			bool single = false;
			num128 c = 1, sc = 1;
			std::vector<std::pair<num64, unsigned int>> factors;
			for (++s; (*s >= '0') && (*s <= '9');)
			{
				num64 q = 0;
				while ((*s >= '0') && (*s <= '9'))
				{
					q = q * 10 + static_cast<num64>(*s - '0');
					++s;
				}
				unsigned int e = 1;
				if (*s == '^')
				{
					e = static_cast<unsigned int>(StrToNumber(++s));
					while ((*s >= '0') && (*s <= '9'))
					{
						++s;
					}
				}
				factors.emplace_back(q, e);
				if (*s == '*')
				{
					++s;
				}
				else
				{
					break;
				}
			}
			for (const auto& fe : factors)
			{
				if (fe.first > p)
				{
					p = fe.first;
					single = (fe.second == 1);
				}
			}
			if (!single || (p <= SearchLimit::LinearLimit) || ((m % 6) == 0))
			{
				continue;
			}
			if ((numPairs++ % stride) != 0)
			{
				continue;
			}
			for (const auto& fe : factors)
			{
				if (fe.first == p)
				{
					continue;
				}
				num128 pk = 1, sp = 1;
				for (unsigned int k = 0; k < fe.second; ++k)
				{
					pk *= fe.first;
					sp += pk;
				}
				c *= pk;
				sc *= sp;
			}
			if ((c * p != m) || HighWord(c) || HighWord(sc))
			{
				continue;
			}

			const unsigned int tier = static_cast<unsigned int>(numTested % mySieve.cfg.numTiers);
			if (!TestSieveLargePrimesPair(m, LowWord(c), LowWord(sc), p, tier, &allAmicable, survivorsTotal))
			{
				return false;
			}
			++numTested;
		}
	}

	std::cerr << "TestSieveLargePrimes: all " << numTested << " of " << numPairs << " large primes pairs found (" << survivorsTotal << " survivors)" << std::endl;
	return numTested > 0;
}
