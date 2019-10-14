#include "stdafx.h"
#include "OpenCL.h"
#include "PrimeTables.h"
#include "RangeGen.h"
#include "sprp64.h"
#include "kernel.cl.inl"

PRAGMA_WARNING(push, 1)
PRAGMA_WARNING(disable : 4091 4917 4625 4626 5026 5027)
#include <boinc_api.h>
#include <boinc_opencl.h>
PRAGMA_WARNING(pop)

//#pragma optimize("", off)
//#undef FORCEINLINE
//#define FORCEINLINE NOINLINE

static const char* const CHECKPOINT_LOGICAL_NAME = "amicable_checkpoint";

FILE* g_outputFile = nullptr;

enum
{
	LogLevel = 0,

	// Limit phase 1 to 2^28 numbers
	Phase1Limit = 1 << 28,
};

#define DEVICE_TYPE_TO_USE CL_DEVICE_TYPE_GPU

#define LOG_ERROR(X) std::cerr << __FILE__ << ", line " << __LINE__ << ": " << X << std::endl << std::flush;
#define LOG(LEVEL, X) IF_CONSTEXPR(LogLevel >= LEVEL) { std::cout << X; IF_CONSTEXPR(LEVEL > 0) { std::cout << std::endl; } }
#define CL_CHECKED_CALL(FUNC, ...) { const cl_int ciErrNum = FUNC(__VA_ARGS__); if (ciErrNum != CL_SUCCESS) { LOG_ERROR(#FUNC" returned error " << ciErrNum); return false; } }
#define CL_CHECKED_CALL_WITH_RESULT(result, FUNC, ...) { cl_int ciErrNum; result = FUNC(__VA_ARGS__, &ciErrNum); if (ciErrNum != CL_SUCCESS) { LOG_ERROR(#FUNC" returned error " << ciErrNum); return false; } }

unsigned char* OpenCL::ourZeroBuf = nullptr;

OpenCL::OpenCL(const char* preferences)
	: myLargestPrimePower(1)
	, myGPUContext(nullptr)
	, myProgram(nullptr)
	, myQueue(nullptr)
	, mySaveCounter(nullptr)
	, mySearchMultipleRanges(nullptr)
	, mySearchLargePrimes(nullptr)
	, myCheckPairs(nullptr)
	, myCheckPairPhase2(nullptr)
	, myCheckPairPhase3(nullptr)
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
	, myPhase1_offset_to_resume_buf(nullptr)
	, myPhase2_numbers_count_buf(nullptr)
	, myPhase2_numbers_buf(nullptr)
	, myPhase3_numbers_count_buf(nullptr)
	, myPhase3_numbers_buf(nullptr)
	, myAmicable_numbers_count_buf(nullptr)
	, myAmicable_numbers_data_buf(nullptr)
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
	SetTimerResoluion(myOldTimerResolution);

	std::cout.precision(myStdOldPrecision);

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
	clReleaseMemObject(myPhase1_offset_to_resume_buf);
	clReleaseMemObject(myPhase2_numbers_count_buf);
	clReleaseMemObject(myPhase2_numbers_buf);
	clReleaseMemObject(myPhase3_numbers_count_buf);
	clReleaseMemObject(myPhase3_numbers_buf);
	clReleaseMemObject(myAmicable_numbers_count_buf);
	clReleaseMemObject(myAmicable_numbers_data_buf);
	clReleaseKernel(mySaveCounter);
	if (mySearchMultipleRanges) clReleaseKernel(mySearchMultipleRanges);
	if (mySearchLargePrimes) clReleaseKernel(mySearchLargePrimes);
	clReleaseKernel(myCheckPairs);
	clReleaseKernel(myCheckPairPhase2);
	clReleaseKernel(myCheckPairPhase3);
	clReleaseCommandQueue(myQueue);
	clReleaseProgram(myProgram);
	clReleaseContext(myGPUContext);
}

static unsigned int crc32(const char *message)
{
	unsigned int crc = 0xFFFFFFFF;
	for (int i = 0; message[i]; ++i)
	{
		crc ^= static_cast<unsigned int>(message[i]);
		for (int j = 7; j >= 0; --j)
		{
			PRAGMA_WARNING(suppress : 4146)
			const unsigned int mask = -(crc & 1);
			crc = (crc >> 1) ^ (0xEDB88320 & mask);
		}
	}
	return ~crc;
}

static bool PlatformSupported(const char* platformName)
{
	// Sorry Intel, but your OpenCL drivers are buggy
	if (strstr(platformName, "Intel") || strstr(platformName, "INTEL") || strstr(platformName, "intel"))
	{
		return false;
	}

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

		CL_CHECKED_CALL(clGetDeviceIDs, platform, DEVICE_TYPE_TO_USE, 1, &device, nullptr);
	}

	if (!(platform && device))
	{
		return false;
	}

	std::string vendor_specific_compiler_options;
	num64 MaxMemAllocSize = 0;
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
		char extensions[1024] = {};
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
		CL_CHECKED_CALL(clGetDeviceInfo, device, CL_DEVICE_EXTENSIONS, sizeof(extensions), extensions, nullptr);

		global_memory_size >>= 20;

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
			"\t" << extensions << std::endl
		);

		if (MaxMemAllocSize < (1 << 28))
		{
			LOG_ERROR("CL_DEVICE_MAX_MEM_ALLOC_SIZE is less than 256 MB, this GPU/driver is not fit to run this program");
			return false;
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

		if (strstr(extensions, "cl_nv_compiler_options"))
		{
			// NVIDIA
			vendor_specific_compiler_options += "-cl-nv-verbose ";
		}
	}

	myLargestPrimePower = largestPrimePower;

	cl_context_properties cps[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties) platform, 0 };
	CL_CHECKED_CALL_WITH_RESULT(myGPUContext, clCreateContext, cps, 1, &device, nullptr, nullptr);

	if (crc32(kernel_cl) != kernel_cl_crc32)
	{
		LOG_ERROR("OpenCL kernel checksum failed");
		return false;
	}

	CL_CHECKED_CALL_WITH_RESULT(myProgram, clCreateProgramWithSource, myGPUContext, 1, &kernel_cl, nullptr);

	num64 NumMainBufferChunks = 0;
	unsigned int ChunkSizeShift = 0;

	const bool IsLargePrimes = (startPrime && primeLimit);
	myMainBufferData = IsLargePrimes ? reinterpret_cast<const unsigned char*>(CandidatesData.data()) : reinterpret_cast<const unsigned char*>(PrimesCompact);
	myMainBufferTotalSize = IsLargePrimes ? static_cast<unsigned int>(CandidatesData.capacity() * AmicableCandidate::PackedSize) : PrimesCompactAllocationSize;
	myMainBufferSize = std::min<num64>(myMainBufferTotalSize, 512 << 20);

	// Try to create a single continuous buffer for storing prime numbers
	// Even though its size is bigger than CL_DEVICE_MAX_MEM_ALLOC_SIZE for many devices,
	// they are still able to do big allocations sometimes
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
		}
	}

	{
		char kernel_options[1024];
		char cBuildLog[10240];
		sprintf_s(kernel_options,
			"%s-cl-fast-relaxed-math -D PQ_STRIDE_SIZE=%u -D WORK_GROUP_SIZE=%u -D PHASE2_MAX_COUNT=%u -D LPP=%u -D NUM_DATA_CHUNKS=%u -D CHUNK_SIZE_SHIFT=%u -D SUM_ESTIMATES_128_SHIFT=%u -Werror",
			vendor_specific_compiler_options.c_str(),
			SumEstimatesSize2_GPU,
			static_cast<unsigned int>(myWorkGroupSize),
			myPhase2MaxNumbersCount,
			myLargestPrimePower,
			static_cast<unsigned int>(NumMainBufferChunks),
			ChunkSizeShift,
			64 - SumEstimates128Shift
		);
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
	}

	CL_CHECKED_CALL_WITH_RESULT(mySaveCounter, clCreateKernel, myProgram, "SaveCounter");

	if (IsLargePrimes)
	{
		CL_CHECKED_CALL_WITH_RESULT(mySearchLargePrimes, clCreateKernel, myProgram, "SearchLargePrimes");
	}
	else
	{
		CL_CHECKED_CALL_WITH_RESULT(mySearchMultipleRanges, clCreateKernel, myProgram, "SearchMultipleRanges");
	}

	CL_CHECKED_CALL_WITH_RESULT(myCheckPairs, clCreateKernel, myProgram, "CheckPairs");
	CL_CHECKED_CALL_WITH_RESULT(myCheckPairPhase2, clCreateKernel, myProgram, "CheckPairPhase2");
	CL_CHECKED_CALL_WITH_RESULT(myCheckPairPhase3, clCreateKernel, myProgram, "CheckPairPhase3");

	CL_CHECKED_CALL_WITH_RESULT(myQueue, clCreateCommandQueue, myGPUContext, device, 0);

	if (!myMainBuffers.empty())
	{
		CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myMainBuffers.front(), CL_TRUE, 0, myMainBufferSize, myMainBufferData, 0, nullptr, nullptr);
	}
	else
	{
		myMainBuffers.reserve(NumMainBufferChunks);
		for (num64 totalSize = 0; totalSize < myMainBufferSize;)
		{
			const num64 BufSize = std::min<num64>(1ULL << ChunkSizeShift, myMainBufferSize - totalSize);
			cl_mem buf;
			CL_CHECKED_CALL_WITH_RESULT(buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, BufSize, 0);
			CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, buf, CL_TRUE, 0, BufSize, myMainBufferData + totalSize, 0, nullptr, nullptr);
			myMainBuffers.emplace_back(buf);
			totalSize += BufSize;
		}
	}

	{
		std::vector<unsigned int> smallPrimes;
		smallPrimes.reserve(ReciprocalsTableSize128);
		for (unsigned int i = 0; i < ReciprocalsTableSize128; ++i)
		{
			smallPrimes.emplace_back(static_cast<unsigned int>(GetNthPrime(i)));
		}
		CL_CHECKED_CALL_WITH_RESULT(mySmallPrimesBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, ReciprocalsTableSize128 * sizeof(unsigned int), 0);
		CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, mySmallPrimesBuf, CL_TRUE, 0, ReciprocalsTableSize128 * sizeof(unsigned int), smallPrimes.data(), 0, nullptr, nullptr);
	}

	const size_t primeInversesBufSize = sizeof(std::pair<num64, num64>) * ReciprocalsTableSize128;
	CL_CHECKED_CALL_WITH_RESULT(myPrimeInversesBuf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, primeInversesBufSize, 0);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPrimeInversesBuf, CL_TRUE, 0, primeInversesBufSize, PrimeInverses, 0, nullptr, nullptr);

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
		CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPrimeReciprocalsBuf, CL_TRUE, 0, primeReciprocalsBufSize, r.data(), 0, nullptr, nullptr);
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
		CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPQ_Buf, CL_TRUE, 0, PQ_BufSize, pq.data(), 0, nullptr, nullptr);
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
		CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPowerOf2SumInverses_Buf, CL_TRUE, 0, sizeof(powerOf2Inverses), powerOf2Inverses, 0, nullptr, nullptr);
	}

	{
		std::vector<unsigned int> offsets;
		offsets.reserve(ReciprocalsTableSize128);
		for (unsigned int i = 0; i < ReciprocalsTableSize128; ++i)
		{
			offsets.emplace_back(static_cast<unsigned int>(PowersOfP_128DivisibilityData[i] - PowersOfP_128DivisibilityData_base));
		}

		CL_CHECKED_CALL_WITH_RESULT(myPowersOfP_128SumInverses_offsets_Buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, sizeof(unsigned int) * ReciprocalsTableSize128, 0);
		CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPowersOfP_128SumInverses_offsets_Buf, CL_TRUE, 0, sizeof(unsigned int) * ReciprocalsTableSize128, offsets.data(), 0, nullptr, nullptr);

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
		CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPowersOfP_128SumInverses_Buf, CL_TRUE, 0, sizeof(num64) * 4 * PowersOfP_128DivisibilityData_count, inverses.data(), 0, nullptr, nullptr);
	}

	CL_CHECKED_CALL_WITH_RESULT(myPrimeInverses_128_Buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, sizeof(num128) * ReciprocalsTableSize128, 0);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPrimeInverses_128_Buf, CL_TRUE, 0, sizeof(num128) * ReciprocalsTableSize128, PrimeInverses128, 0, nullptr, nullptr);

	CL_CHECKED_CALL_WITH_RESULT(mySumEstimates_128_Buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, sizeof(num64) * ReciprocalsTableSize128 / 16, 0);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, mySumEstimates_128_Buf, CL_TRUE, 0, sizeof(num64) * ReciprocalsTableSize128 / 16, SumEstimates128, 0, nullptr, nullptr);

	CL_CHECKED_CALL_WITH_RESULT(myRangesTable_Buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, sizeof(RangeDataGPU) * myMaxRangesCount, 0);
	CL_CHECKED_CALL_WITH_RESULT(myRangesLookupTable_Buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_ONLY, sizeof(LookupDataGPU) * myMaxLookupTableSize, 0);

	CL_CHECKED_CALL_WITH_RESULT(myPhase1_offset_to_resume_buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(num64), 0);
	CL_CHECKED_CALL_WITH_RESULT(myPhase2_numbers_count_buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(unsigned int) * 2, 0);
	CL_CHECKED_CALL_WITH_RESULT(myPhase3_numbers_count_buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(num64), 0);
	CL_CHECKED_CALL_WITH_RESULT(myAmicable_numbers_count_buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(num64), 0);

	ourZeroBuf = new unsigned char[ourZeroBufSize];
	memset(ourZeroBuf, 0, ourZeroBufSize);

	CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPhase1_offset_to_resume_buf, CL_TRUE, 0, sizeof(num64), ourZeroBuf, 0, nullptr, nullptr);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPhase2_numbers_count_buf, CL_TRUE, 0, sizeof(unsigned int) * 2, ourZeroBuf, 0, nullptr, nullptr);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPhase3_numbers_count_buf, CL_TRUE, 0, sizeof(num64), ourZeroBuf, 0, nullptr, nullptr);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myAmicable_numbers_count_buf, CL_TRUE, 0, sizeof(num64), ourZeroBuf, 0, nullptr, nullptr);

	CL_CHECKED_CALL_WITH_RESULT(myPhase2_numbers_buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(num64) * 4 * myPhase2MaxNumbersCount, 0);
	CL_CHECKED_CALL_WITH_RESULT(myPhase3_numbers_buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(num64) * 4 * myPhase2MaxNumbersCount, 0);
	CL_CHECKED_CALL_WITH_RESULT(myAmicable_numbers_data_buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(num64) * 4 * 8192, 0);

	CL_CHECKED_CALL(clFinish, myQueue);

	CL_CHECKED_CALL(setKernelArguments, mySaveCounter, myPhase2_numbers_count_buf);

	if (IsLargePrimes)
	{
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
			CandidatesDataHighBitOffsets,
			myPhase1_offset_to_resume_buf,
			myPhase2_numbers_count_buf,
			myPhase2_numbers_buf,
			myAmicable_numbers_count_buf,
			myAmicable_numbers_data_buf
		);

		for (num64 i = 0; i < myMainBuffers.size(); ++i)
		{
			CL_CHECKED_CALL(clSetKernelArg, mySearchLargePrimes, static_cast<cl_uint>(i + 20), sizeof(cl_mem), &myMainBuffers[i]);
		}
	}
	else
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
			myPhase1_offset_to_resume_buf,
			myPhase2_numbers_count_buf,
			myPhase2_numbers_buf,
			myAmicable_numbers_count_buf,
			myAmicable_numbers_data_buf
		);

		for (num64 i = 0; i < myMainBuffers.size(); ++i)
		{
			CL_CHECKED_CALL(clSetKernelArg, mySearchMultipleRanges, static_cast<cl_uint>(i + 18), sizeof(cl_mem), &myMainBuffers[i]);
		}
	}

	CL_CHECKED_CALL(setKernelArguments, myCheckPairPhase2,
		mySmallPrimesBuf,
		myPhase2_numbers_buf,
		0,
		myPrimeInversesBuf,
		myPowersOfP_128SumInverses_offsets_Buf,
		myPowersOfP_128SumInverses_Buf,
		myPrimeInverses_128_Buf,
		mySumEstimates_128_Buf,
		myPQ_Buf,
		myPhase3_numbers_count_buf,
		myPhase3_numbers_buf,
		myAmicable_numbers_count_buf,
		myAmicable_numbers_data_buf
	);

	IF_CONSTEXPR(LogLevel > 0)
	{
		if (!Test())
		{
			return false;
		}
	}

	Timer t;
	const bool result = IsLargePrimes ? RunLargePrimes(startPrime, primeLimit) : RunRanges(startFrom, stopAt);
	LOG(1, "\nWork unit finished in " << std::setprecision(6) << t.getElapsedTime() << " seconds, " << myNumbersProcessedTotal << " numbers processed\n");
	return result;
}

bool OpenCL::RunRanges(char* startFrom, char* stopAt)
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
				char* first_token_end = strchr(checkpoint_buf, ' ');
				if (first_token_end)
				{
					myNumbersProcessedTotal = StrToNumber(checkpoint_buf);
					const double f = myNumbersProcessedTotal / RangeGen::total_numbers_to_check;
					boinc_fraction_done((f < 1.0) ? f : 1.0);
					do
					{
						++first_token_end;
					} while (*first_token_end == ' ');
					startFrom = first_token_end;
				}
			}
		}
		fclose(checkpoint);
	}

	RangeGen::Init(startFrom, stopAt, &r, stopAtFactors, myLargestPrimePower);
	if (!startFrom)
	{
		RangeGen::Iterate(r);
	}
	while (!RangeGen::HasReached(r, stopAtFactors))
	{
		if (!AddRange(r))
		{
			return false;
		}
		if (!RangeGen::Iterate(r) || (RangeGen::cur_largest_prime_power != myLargestPrimePower))
		{
			break;
		}

		const double progress = static_cast<double>(myNumbersProcessedTotal) / RangeGen::total_numbers_to_check;
		boinc_fraction_done((progress < 1.0) ? progress : 1.0);

		if (myTotalNumbersInRanges == 0)
		{
			unsigned int numbers_in_phase2 = 0;
			if (!GetCounter(myQueue, myPhase2_numbers_count_buf, numbers_in_phase2))
			{
				return false;
			}
			if (numbers_in_phase2 == 0)
			{
				if (!SaveProgressForRanges(r))
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
	if (!GetCounter(myQueue, myPhase2_numbers_count_buf, numbers_in_phase2))
	{
		return false;
	}

	if ((numbers_in_phase2 > 0) && !ProcessNumbersPhases2_3(numbers_in_phase2))
	{
		return false;
	}
	CleanupRanges();

	return SaveProgressForRanges(r);
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
	if (!GetCounter(myQueue, myPhase2_numbers_count_buf, numbers_in_phase2))
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
		CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myLargePrimesBuf, CL_TRUE, 0, sizeof(num64) * LargePrimesCount, myLargePrimes, 0, nullptr, nullptr);

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
			const std::pair<unsigned int, unsigned int>* packedCandidates = reinterpret_cast<const std::pair<unsigned int, unsigned int>*>(CandidatesData.data());
			int a = 0;
			int b = static_cast<int>(CandidatesData.size());
			while (a < b)
			{
				const int c = (a + b) >> 1;

				num64 value = packedCandidates[c].first;
				if (c >= CandidatesDataHighBitOffsets.first)
				{
					value |= 0x100000000ULL;
				}

				if (value <= LargestCandidate)
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
			if (!GetCounter(myQueue, myPhase2_numbers_count_buf, numbers_in_phase2_before))
			{
				return false;
			}

			//
			// Phase 1
			//
			{
				Timer t;

				CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPhase1_offset_to_resume_buf, CL_TRUE, 0, sizeof(num64), ourZeroBuf, 0, nullptr, nullptr);

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
						CL_CHECKED_CALL(clFlush, myQueue);
					}
					size_t globalSizePhase1 = std::min<size_t>(N - global_offset, max_size_phase1);
					if (globalSizePhase1 & (myWorkGroupSize - 1))
					{
						globalSizePhase1 += myWorkGroupSize - (globalSizePhase1 & (myWorkGroupSize - 1));
					}
					CL_CHECKED_CALL(clSetKernelArg, mySearchLargePrimes, 12, sizeof(global_offset), &global_offset);
					CL_CHECKED_CALL(clEnqueueTask, myQueue, mySaveCounter, 0, nullptr, nullptr);
					global_offset += globalSizePhase1;
					CL_CHECKED_CALL(clEnqueueNDRangeKernel, myQueue, mySearchLargePrimes, 1, nullptr, &globalSizePhase1, &myWorkGroupSize, 0, nullptr, (global_offset >= N) ? &event : nullptr);
				}

				if (!WaitForQueue(myQueue, event) ||
					!GetCounter(myQueue, myPhase2_numbers_count_buf, numbers_in_phase2_after) ||
					!GetAndResetCounter(myQueue, myPhase1_offset_to_resume_buf, phase1_offset_to_resume))
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
					CL_CHECKED_CALL(clEnqueueReadBuffer, myQueue, myPhase2_numbers_count_buf, CL_TRUE, 0, sizeof(unsigned int) * 2, &counters, 0, nullptr, nullptr);
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
	struct SPairToCheck
	{
		SPairToCheck() : M(0), targetSum(0) {}
		SPairToCheck(num128 m, num128 sum) : M(m), targetSum(sum) {}

		num128 M;
		num128 targetSum;
	};
	std::vector<SPairToCheck> pairs;

	const char* files[] = {
		//"c2_1.txt"
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

	const unsigned int total_numbers = static_cast<unsigned int>(pairs.size());
	LOG(1, "Total numbers: " << total_numbers);
	if (total_numbers == 0)
	{
		LOG(1, "Skipping test");
		return true;
	}

	Timer overall_timer;

	cl_mem pairsToCheck_Buf;
	pairs.resize(((pairs.size() + myWorkGroupSize - 1) / myWorkGroupSize) * myWorkGroupSize);
	CL_CHECKED_CALL_WITH_RESULT(pairsToCheck_Buf, clCreateBufferLogged, __LINE__, myGPUContext, CL_MEM_READ_WRITE, sizeof(SPairToCheck) * pairs.size(), 0);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, pairsToCheck_Buf, CL_TRUE, 0, sizeof(SPairToCheck) * pairs.size(), pairs.data(), 0, nullptr, nullptr);

	struct TmpBuf
	{
		explicit TmpBuf(cl_context gpuContext) { mem = clCreateBuffer(gpuContext, CL_MEM_READ_WRITE, sizeof(num64) * 4 * 5000000, 0, nullptr); }
		~TmpBuf() { clReleaseMemObject(mem); }

		TmpBuf(const TmpBuf&) = delete;
		TmpBuf(TmpBuf&&) = delete;
		TmpBuf& operator=(const TmpBuf&) = delete;
		TmpBuf& operator=(TmpBuf&&) = delete;

		operator cl_mem() { return mem; }

		cl_mem mem;
	} buf(myGPUContext);

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
		CL_CHECKED_CALL(setKernelArguments, myCheckPairs, mySmallPrimesBuf, myPrimeInversesBuf, myPQ_Buf, myPowerOf2SumInverses_Buf, myPowersOfP_128SumInverses_offsets_Buf, myPowersOfP_128SumInverses_Buf, myPrimeInverses_128_Buf, mySumEstimates_128_Buf, pairsToCheck_Buf, myPhase2_numbers_count_buf, myPhase2_numbers_buf, myAmicable_numbers_count_buf, buf);
		CL_CHECKED_CALL(clEnqueueNDRangeKernel, myQueue, myCheckPairs, 1, nullptr, &globalSizePhase1, &myWorkGroupSize, 0, nullptr, &event);

		if (!WaitForQueue(myQueue, event) || !GetAndResetCounter(myQueue, myPhase2_numbers_count_buf, filtered_numbers_count) || !GetCounter(myQueue, myAmicable_numbers_count_buf, amicable_numberscount))
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

		CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPhase2_numbers_buf, CL_TRUE, filtered_numbers_count * sizeof(num64) * 4, myWorkGroupSize * sizeof(num64) * 4, ourZeroBuf, 0, nullptr, nullptr);
		CL_CHECKED_CALL(setKernelArguments, myCheckPairPhase2,
			mySmallPrimesBuf,
			myPhase2_numbers_buf,
			0,
			myPrimeInversesBuf,
			myPowersOfP_128SumInverses_offsets_Buf,
			myPowersOfP_128SumInverses_Buf,
			myPrimeInverses_128_Buf,
			mySumEstimates_128_Buf,
			myPQ_Buf,
			myPhase3_numbers_count_buf,
			myPhase3_numbers_buf,
			myAmicable_numbers_count_buf,
			buf
		);

		size_t globalSizePhase2 = filtered_numbers_count;
		if (globalSizePhase2 & (myWorkGroupSize - 1))
		{
			globalSizePhase2 += myWorkGroupSize - (globalSizePhase2 & (myWorkGroupSize - 1));
		}
		CL_CHECKED_CALL(clEnqueueNDRangeKernel, myQueue, myCheckPairPhase2, 1, nullptr, &globalSizePhase2, &myWorkGroupSize, 0, nullptr, &event);

		if (!WaitForQueue(myQueue, event) || !GetAndResetCounter(myQueue, myPhase3_numbers_count_buf, filtered_numbers2_count) || !GetCounter(myQueue, myAmicable_numbers_count_buf, amicable_numberscount))
		{
			return false;
		}

		const double phase2Time = t.getElapsedTime();
		LOG(2, "Phase 2: " << phase2Time << " seconds, " << filtered_numbers2_count << " numbers left (" << ((filtered_numbers2_count * 100.0) / total_numbers) << "%), " << amicable_numberscount << " amicable numbers found");
	}

	//
	// Phase 3
	//
	cl_mem phase3_buffers[2] = {myPhase2_numbers_buf, myPhase3_numbers_buf};
	const size_t phase3_workGroupSize = myWorkGroupSize;

	unsigned int filtered_numbers3_count = filtered_numbers2_count;
	for (int phase3_iteration = 1; filtered_numbers3_count > 0; ++phase3_iteration)
	{
		Timer t1;

		CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, phase3_buffers[phase3_iteration & 1], CL_TRUE, filtered_numbers3_count * sizeof(num64) * 4, phase3_workGroupSize * sizeof(num64) * 4, ourZeroBuf, 0, nullptr, nullptr);

		const unsigned int max_size_phase3 = GetMaxPhaseSize(filtered_numbers3_count, 1 << 18);
		size_t globalSizePhase3;
		for (unsigned int offset = 0; offset < filtered_numbers3_count; offset += static_cast<unsigned int>(globalSizePhase3))
		{
			globalSizePhase3 = std::min(filtered_numbers3_count - offset, max_size_phase3);
			if (globalSizePhase3 & (phase3_workGroupSize - 1))
			{
				globalSizePhase3 += phase3_workGroupSize - (globalSizePhase3 & (phase3_workGroupSize - 1));
			}
			CL_CHECKED_CALL(setKernelArguments, myCheckPairPhase3,
				phase3_buffers[phase3_iteration & 1],
				offset,
				myPrimeInversesBuf,
				myPrimeReciprocalsBuf,
				myPhase3_numbers_count_buf,
				phase3_buffers[(phase3_iteration - 1) & 1],
				myAmicable_numbers_count_buf,
				buf,
				phase3_iteration << 12,
				mySmallPrimesBuf,
				myPowersOfP_128SumInverses_offsets_Buf,
				myPowersOfP_128SumInverses_Buf,
				myPrimeInverses_128_Buf
			);
			CL_CHECKED_CALL(clEnqueueNDRangeKernel, myQueue, myCheckPairPhase3, 1, nullptr, &globalSizePhase3, &phase3_workGroupSize, 0, nullptr, (offset + globalSizePhase3 >= filtered_numbers3_count) ? &event : nullptr);
		}

		if (!WaitForQueue(myQueue, event) || !GetAndResetCounter(myQueue, myPhase3_numbers_count_buf, filtered_numbers3_count) || !GetCounter(myQueue, myAmicable_numbers_count_buf, amicable_numberscount))
		{
			return false;
		}

		const double phase3Time = t1.getElapsedTime();
		LOG(2, "Phase 3, iteration " << phase3_iteration << ": " << phase3Time << " seconds, " << filtered_numbers3_count << " numbers left (" << ((filtered_numbers3_count * 100.0) / total_numbers) << "%), " << amicable_numberscount << " amicable numbers found");
	}

	if (!GetAndResetCounter(myQueue, myAmicable_numbers_count_buf, amicable_numberscount))
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
		CL_CHECKED_CALL(clEnqueueReadBuffer, myQueue, buf, CL_TRUE, 0, sizeof(num64) * 4 * amicable_numberscount, found_pairs.data(), 0, nullptr, nullptr);
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

	for (const SPairToCheck& pair : pairs)
	{
		if (pair.targetSum != NUM128_MAX)
		{
			LOG(1, "First missed number: " << pair.M);
			return false;
		}
	}

	return true;
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

bool OpenCL::AddRange(const RangeData& r)
{
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

	unsigned int index_begin = r.index_start_prime;

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

	if ((a > index_begin) && (index_begin < myMainBufferSize / 2))
	{
		const unsigned int range_size = static_cast<unsigned int>(std::min(a - index_begin, myMainBufferSize / 2 - index_begin));

		myRanges.emplace_back(r.value, r.sum, index_begin, range_size);
		myTotalNumbersInRanges += range_size;
		myNumbersProcessedTotal += range_size;

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
		if (!GetCounter(myQueue, myPhase2_numbers_count_buf, numbers_in_phase2_before))
		{
			return false;
		}

		//
		// Phase 1
		//
		{
			Timer t;

			CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPhase1_offset_to_resume_buf, CL_TRUE, 0, sizeof(num64), ourZeroBuf, 0, nullptr, nullptr);

			CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myRangesTable_Buf, CL_TRUE, 0, myRanges.size() * sizeof(RangeDataGPU), myRanges.data(), 0, nullptr, nullptr);
			CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myRangesLookupTable_Buf, CL_TRUE, 0, myLookupTable.size() * sizeof(LookupDataGPU), myLookupTable.data(), 0, nullptr, nullptr);

			CL_CHECKED_CALL(clSetKernelArg, mySearchMultipleRanges, 10, sizeof(unsigned int), &lookupTableShift);
			CL_CHECKED_CALL(clSetKernelArg, mySearchMultipleRanges, 12, sizeof(unsigned int), &myTotalNumbersInRanges);

			const num64 global_offset_before = global_offset;
			const unsigned int max_size_phase1 = GetMaxPhaseSize(myTotalNumbersInRanges - global_offset, myPhase1MaxKernelSize);

			size_t globalSizePhase1;
			for (unsigned int k = 0; global_offset < myTotalNumbersInRanges; global_offset += globalSizePhase1, ++k)
			{
				if (k)
				{
					CL_CHECKED_CALL(clFlush, myQueue);
				}
				globalSizePhase1 = std::min<num64>(myTotalNumbersInRanges - global_offset, max_size_phase1);
				if (globalSizePhase1 & (myWorkGroupSize - 1))
				{
					globalSizePhase1 += myWorkGroupSize - (globalSizePhase1 & (myWorkGroupSize - 1));
				}
				CL_CHECKED_CALL(clSetKernelArg, mySearchMultipleRanges, 11, sizeof(global_offset), &global_offset);
				CL_CHECKED_CALL(clEnqueueTask, myQueue, mySaveCounter, 0, nullptr, nullptr);
				CL_CHECKED_CALL(clEnqueueNDRangeKernel, myQueue, mySearchMultipleRanges, 1, nullptr, &globalSizePhase1, &myWorkGroupSize, 0, nullptr, (global_offset + globalSizePhase1 >= myTotalNumbersInRanges) ? &event : nullptr);
			}

			if (!WaitForQueue(myQueue, event) ||
				!GetCounter(myQueue, myPhase2_numbers_count_buf, numbers_in_phase2_after) ||
				!GetAndResetCounter(myQueue, myPhase1_offset_to_resume_buf, phase1_offset_to_resume))
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
				CL_CHECKED_CALL(clEnqueueReadBuffer, myQueue, myPhase2_numbers_count_buf, CL_TRUE, 0, sizeof(unsigned int) * 2, &counters, 0, nullptr, nullptr);
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

		if (numbers_in_phase2 < myPhase2MaxNumbersCount)
		{
			CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPhase2_numbers_buf, CL_TRUE, numbers_in_phase2 * sizeof(num64) * 4, std::min(static_cast<unsigned int>(myWorkGroupSize), myPhase2MaxNumbersCount - numbers_in_phase2) * sizeof(num64) * 4, ourZeroBuf, 0, nullptr, nullptr);
		}

		const unsigned int max_size_phase2 = GetMaxPhaseSize(numbers_in_phase2, myPhase1MaxKernelSize / 4);
		LOG(2, "Phase 2: max kernel size = " << max_size_phase2 << " numbers");

		size_t globalSizePhase2;
		for (unsigned int k = 0, global_offset = 0; global_offset < numbers_in_phase2; global_offset += static_cast<unsigned int>(globalSizePhase2), ++k)
		{
			if (k)
			{
				CL_CHECKED_CALL(clFlush, myQueue);
			}
			globalSizePhase2 = std::min(numbers_in_phase2 - global_offset, max_size_phase2);
			if (globalSizePhase2 & (myWorkGroupSize - 1))
			{
				globalSizePhase2 += myWorkGroupSize - (globalSizePhase2 & (myWorkGroupSize - 1));
			}
			CL_CHECKED_CALL(clSetKernelArg, myCheckPairPhase2, 2, sizeof(global_offset), &global_offset);
			CL_CHECKED_CALL(clEnqueueNDRangeKernel, myQueue, myCheckPairPhase2, 1, nullptr, &globalSizePhase2, &myWorkGroupSize, 0, nullptr, (global_offset + static_cast<unsigned int>(globalSizePhase2) >= numbers_in_phase2) ? &event : nullptr);
		}

		if (!WaitForQueue(myQueue, event) || !GetAndResetCounter(myQueue, myPhase3_numbers_count_buf, phase3_numbers_buf_count))
		{
			return false;
		}

		const double phase2Time = t.getElapsedTime();
		LOG(2, "Phase 2: processed " << numbers_in_phase2 << " numbers in " << phase2Time << " seconds, " << phase3_numbers_buf_count << " numbers left (" << ((phase3_numbers_buf_count * 100.0) / numbers_in_phase2) << "%)");
	}

	//
	// Phase 3
	//
	if (phase3_numbers_buf_count > 0)
	{
		cl_mem phase3_buffers[2] = { myPhase2_numbers_buf, myPhase3_numbers_buf };
		const size_t phase3_workGroupSize = myWorkGroupSize;

		unsigned int phase3_numbers_buf_count2 = phase3_numbers_buf_count;
		for (int phase3_iteration = 1; phase3_numbers_buf_count2 > 0; ++phase3_iteration)
		{
			Timer t1;

			if (phase3_numbers_buf_count2 < myPhase2MaxNumbersCount)
			{
				CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, phase3_buffers[phase3_iteration & 1], CL_TRUE, phase3_numbers_buf_count2 * sizeof(num64) * 4, std::min(static_cast<unsigned int>(phase3_workGroupSize), myPhase2MaxNumbersCount - phase3_numbers_buf_count2) * sizeof(num64) * 4, ourZeroBuf, 0, nullptr, nullptr);
			}

			const unsigned int max_size_phase3 = GetMaxPhaseSize(phase3_numbers_buf_count2, myPhase1MaxKernelSize / 16);
			LOG(2, "Phase 3, iteration " << phase3_iteration << ": max kernel size = " << max_size_phase3 << " numbers");

			size_t globalSizePhase3;
			for (unsigned int k = 0, global_offset = 0; global_offset < phase3_numbers_buf_count2; global_offset += static_cast<unsigned int>(globalSizePhase3), ++k)
			{
				if (k)
				{
					CL_CHECKED_CALL(clFlush, myQueue);
				}
				globalSizePhase3 = std::min<num64>(phase3_numbers_buf_count2 - global_offset, max_size_phase3);
				if (globalSizePhase3 & (myWorkGroupSize - 1))
				{
					globalSizePhase3 += myWorkGroupSize - (globalSizePhase3 & (myWorkGroupSize - 1));
				}
				CL_CHECKED_CALL(setKernelArguments, myCheckPairPhase3,
					phase3_buffers[phase3_iteration & 1],
					global_offset,
					myPrimeInversesBuf,
					myPrimeReciprocalsBuf,
					myPhase3_numbers_count_buf,
					phase3_buffers[(phase3_iteration - 1) & 1],
					myAmicable_numbers_count_buf,
					myAmicable_numbers_data_buf,
					phase3_iteration << 12,
					mySmallPrimesBuf,
					myPowersOfP_128SumInverses_offsets_Buf,
					myPowersOfP_128SumInverses_Buf,
					myPrimeInverses_128_Buf
				);
				CL_CHECKED_CALL(clEnqueueNDRangeKernel, myQueue, myCheckPairPhase3, 1, nullptr, &globalSizePhase3, &phase3_workGroupSize, 0, nullptr, (global_offset + static_cast<unsigned int>(globalSizePhase3) >= phase3_numbers_buf_count2) ? &event : nullptr);
			}

			if (!WaitForQueue(myQueue, event) || !GetAndResetCounter(myQueue, myPhase3_numbers_count_buf, phase3_numbers_buf_count2))
			{
				return false;
			}

			const double phase3Time = t1.getElapsedTime();
			LOG(2, "Phase 3, iteration " << phase3_iteration << ": processed " << phase3_numbers_buf_count << " numbers in " << phase3Time << " seconds, " << phase3_numbers_buf_count2 << " numbers left (" << ((phase3_numbers_buf_count2 * 100.0) / phase3_numbers_buf_count) << "%)");
			phase3_numbers_buf_count = phase3_numbers_buf_count2;
		}

		LOG(2, "Phase 3 finished");
	}

	if (!ResetCounter(myQueue, myPhase2_numbers_count_buf))
	{
		return false;
	}

	return true;
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
	if (!GetAndResetCounter(myQueue, myAmicable_numbers_count_buf, numbers_found))
	{
		return false;
	}

	if (numbers_found > 0)
	{
		myFoundPairs.resize(numbers_found);
		CL_CHECKED_CALL(clEnqueueReadBuffer, myQueue, myAmicable_numbers_data_buf, CL_TRUE, 0, sizeof(num64) * 4 * numbers_found, myFoundPairs.data(), 0, nullptr, nullptr);
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

bool OpenCL::SaveProgressForRanges(const RangeData& r)
{
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
		s << myNumbersProcessedTotal << ' ';

		// Then current range
		const Factor(&mainFactors)[MaxPrimeFactors] = r.factors;
		for (num64 i = 0; i <= static_cast<num64>(r.last_factor_index); ++i)
		{
			if (i > 0)
			{
				s << '*';
			}
			s << mainFactors[i].p.Get();
			if (mainFactors[i].k > 1)
			{
				s << '^' << mainFactors[i].k;
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
