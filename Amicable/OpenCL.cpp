#include "stdafx.h"
#include "OpenCL.h"
#include "PrimeTables.h"
#include "RangeGen.h"
#include "sprp64.h"
#include "inverses128.h"
#include "kernel.cl.inl"
#include <iomanip>
#include <sstream>

PRAGMA_WARNING(push, 1)
PRAGMA_WARNING(disable : 4091 4917)
#include <boinc_api.h>
#include <boinc_opencl.h>
PRAGMA_WARNING(pop)

static const char* const CHECKPOINT_LOGICAL_NAME = "amicable_checkpoint";

FILE* g_outputFile = nullptr;

enum
{
	LogLevel = 0,
};

#define DEVICE_TYPE_TO_USE CL_DEVICE_TYPE_GPU

#define LOG_ERROR(X) std::cerr << __FILE__ << ", line " << __LINE__ << ": " << X << std::endl << std::flush;
#define LOG(LEVEL, X) IF_CONSTEXPR(LogLevel >= LEVEL) { std::cout << X; IF_CONSTEXPR(LEVEL > 0) { std::cout << std::endl; } }
#define CL_CHECKED_CALL(FUNC, ...) { const cl_int ciErrNum = FUNC(__VA_ARGS__); if (ciErrNum != CL_SUCCESS) { LOG_ERROR(#FUNC" returned error " << ciErrNum); return false; } }
#define CL_CHECKED_CALL_WITH_RESULT(result, FUNC, ...) { cl_int ciErrNum; result = FUNC(__VA_ARGS__, &ciErrNum); if (ciErrNum != CL_SUCCESS) { LOG_ERROR(#FUNC" returned error " << ciErrNum); return false; } }

const unsigned char OpenCL::ourZeroBuf[16384] = {};

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
	, myPowersOfP_128SumInverses_Buf(nullptr)
	, myPrimeInverses_128_Buf(nullptr)
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
	, myMaxRangesCount(8192)
	, myMaxLookupTableSize((myMaxRangesCount + 1) * 4) // very conservative estimate to be safe
	, myTotalNumbersInRanges(0)
	, myNumbersProcessedTotal(0)
	, myAmicableNumbersFound(0)
	, myPreferences(preferences)
	, myLargePrimes(nullptr)
	, myLargePrimesCount(0)
	, myLargePrimesMaxCount(1048576)
{
	SetKernelSize(20);

	myRanges.reserve(myMaxRangesCount);
	myLookupTable.reserve(myMaxLookupTableSize);
	myBufUlong2.reserve(8192);

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
	clReleaseMemObject(myPowersOfP_128SumInverses_Buf);
	clReleaseMemObject(myPrimeInverses_128_Buf);
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
	if (size > 21)
	{
		size = 21;
	}

	myPhase1MaxKernelSize = 1U << size;
	myPhase2MinNumbersCount = 1U << size;
	myPhase2MaxNumbersCount = myPhase2MinNumbersCount + myPhase1MaxKernelSize;

	return size;
}

bool OpenCL::Run(int argc, char* argv[], char* startFrom, char* stopAt, unsigned int largestPrimePower, number startPrime, number primeLimit)
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
	number MaxMemAllocSize = 0;
	{
		char platformName[1024] = {};
		CL_CHECKED_CALL(clGetPlatformInfo, platform, CL_PLATFORM_NAME, sizeof(platformName), &platformName, nullptr);
		if (!PlatformSupported(platformName))
		{
			LOG_ERROR("Platform '" << platformName << "' is not supported by this program");
			return false;
		}

		number global_memory_size = 0;
		number global_memory_cache_size = 0;
		number local_memory_size = 0;
		unsigned int constant_args = 0;
		number constant_buffer_size = 0;
		unsigned int max_compute_units = 0;
		unsigned int max_clock_frequency = 0;
		number max_work_group_size = 0;
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

		// If GPU has 2 GB of memory or more, increase internal buffers
		if (global_memory_size >= 2048)
		{
			myMaxRangesCount *= 8;
			myPhase1MaxKernelSize *= 2;
			myPhase2MinNumbersCount *= 2;

			myMaxLookupTableSize = (myMaxRangesCount + 1) * 4;
			myPhase2MaxNumbersCount = myPhase2MinNumbersCount + myPhase1MaxKernelSize;
			myRanges.reserve(myMaxRangesCount);
			myLookupTable.reserve(myMaxLookupTableSize);
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

	unsigned int NumMainBufferChunks = 0;
	unsigned int ChunkSizeShift = 0;

	const bool IsLargePrimes = (startPrime && primeLimit);
	const unsigned char* MainBufferData = IsLargePrimes ? reinterpret_cast<const unsigned char*>(CandidatesData.data()) : reinterpret_cast<const unsigned char*>(PrimesCompact);
	const unsigned int MainBufferSize = IsLargePrimes ? static_cast<unsigned int>(CandidatesData.capacity() * sizeof(AmicableCandidate)) : PrimesCompactAllocationSize;

	// Try to create a single continuous buffer for storing prime numbers
	// Even though its size is bigger than CL_DEVICE_MAX_MEM_ALLOC_SIZE for many devices,
	// they are still able to do big allocations sometimes
	{
		cl_int ciErrNum;
		cl_mem buf = clCreateBuffer(myGPUContext, CL_MEM_READ_ONLY, MainBufferSize, 0, &ciErrNum);
		if (ciErrNum == CL_SUCCESS)
		{
			NumMainBufferChunks = 1;
			unsigned long index;
			_BitScanReverse64(&index, MainBufferSize);
			ChunkSizeShift = index + 1;
			myMainBuffers.emplace_back(buf);
		}
		else
		{
			// If it failed, take into account maximum allocation size for this device
			unsigned long index;
			_BitScanReverse64(&index, MaxMemAllocSize);
			ChunkSizeShift = index;
			NumMainBufferChunks = (MainBufferSize + (1 << ChunkSizeShift) - 1) >> ChunkSizeShift;
		}
	}

	{
		char kernel_options[1024];
		char cBuildLog[10240];
		sprintf_s(kernel_options,
			"%s-cl-fast-relaxed-math -D PQ_STRIDE_SIZE=%u -D WORK_GROUP_SIZE=%u -D PHASE2_MAX_COUNT=%u -D LPP=%u -D NUM_DATA_CHUNKS=%u -D CHUNK_SIZE_SHIFT=%u -Werror",
			vendor_specific_compiler_options.c_str(),
			SumEstimatesSize2_GPU,
			static_cast<unsigned int>(myWorkGroupSize),
			myPhase2MaxNumbersCount,
			myLargestPrimePower,
			NumMainBufferChunks,
			ChunkSizeShift
		);
		cl_int build_result = clBuildProgram(myProgram, 0, nullptr, kernel_options, nullptr, nullptr);

		CL_CHECKED_CALL(clGetProgramBuildInfo, myProgram, device, CL_PROGRAM_BUILD_LOG, sizeof(cBuildLog), cBuildLog, nullptr);
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
		CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myMainBuffers.front(), CL_TRUE, 0, MainBufferSize, MainBufferData, 0, nullptr, nullptr);
	}
	else
	{
		myMainBuffers.reserve(NumMainBufferChunks);
		for (unsigned int totalSize = 0; totalSize < MainBufferSize;)
		{
			const unsigned int BufSize = std::min<unsigned int>(1U << ChunkSizeShift, MainBufferSize - totalSize);
			cl_mem buf;
			CL_CHECKED_CALL_WITH_RESULT(buf, clCreateBuffer, myGPUContext, CL_MEM_READ_ONLY, BufSize, 0);
			CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, buf, CL_TRUE, 0, BufSize, MainBufferData + totalSize, 0, nullptr, nullptr);
			myMainBuffers.emplace_back(buf);
			totalSize += BufSize;
		}
	}

	{
		std::vector<unsigned int> smallPrimes;
		smallPrimes.reserve(ReciprocalsTableSize);
		for (int i = 0; i < ReciprocalsTableSize; ++i)
		{
			smallPrimes.emplace_back(static_cast<unsigned int>(GetNthPrime(i)));
		}
		CL_CHECKED_CALL_WITH_RESULT(mySmallPrimesBuf, clCreateBuffer, myGPUContext, CL_MEM_READ_ONLY, ReciprocalsTableSize * sizeof(unsigned int), 0);
		CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, mySmallPrimesBuf, CL_TRUE, 0, ReciprocalsTableSize * sizeof(unsigned int), smallPrimes.data(), 0, nullptr, nullptr);
	}

	const size_t primeInversesBufSize = sizeof(std::pair<number, number>) * ReciprocalsTableSize;
	CL_CHECKED_CALL_WITH_RESULT(myPrimeInversesBuf, clCreateBuffer, myGPUContext, CL_MEM_READ_ONLY, primeInversesBufSize, 0);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPrimeInversesBuf, CL_TRUE, 0, primeInversesBufSize, PrimeInverses, 0, nullptr, nullptr);

	{
		const size_t primeReciprocalsBufSize = sizeof(std::pair<number, number>) * ReciprocalsTableSize;
		CL_CHECKED_CALL_WITH_RESULT(myPrimeReciprocalsBuf, clCreateBuffer, myGPUContext, CL_MEM_READ_ONLY, primeReciprocalsBufSize, 0);

		std::vector<std::pair<number, number>> r;
		r.resize(ReciprocalsTableSize);
		for (unsigned int i = 0; i < ReciprocalsTableSize; ++i)
		{
			r[i].first = PrimeReciprocals[i].reciprocal;
			r[i].second = (GetNthPrime(i) << 32) + (number(PrimeReciprocals[i].increment) << 8) + number(PrimeReciprocals[i].shift);
		}
		CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPrimeReciprocalsBuf, CL_TRUE, 0, primeReciprocalsBufSize, r.data(), 0, nullptr, nullptr);
	}

	{
		std::vector<std::pair<number, number>> pq;
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

		const size_t PQ_BufSize = sizeof(std::pair<number, number>) * pq.size();
		CL_CHECKED_CALL_WITH_RESULT(myPQ_Buf, clCreateBuffer, myGPUContext, CL_MEM_READ_ONLY, PQ_BufSize, 0);
		CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPQ_Buf, CL_TRUE, 0, PQ_BufSize, pq.data(), 0, nullptr, nullptr);
	}

	{
		static number powerOf2Inverses[128][2];
		for (number i = 0; i < 64; ++i)
		{
			const number value = (i < 63) ? ((number(1) << (i + 1)) - 1) : number(-1);
			PRAGMA_WARNING(suppress : 4146)
			powerOf2Inverses[i][0] = -modular_inverse64(value);
			powerOf2Inverses[i][1] = number(-1) / value;
		}
		for (number i = 0; i < 64; ++i)
		{
			powerOf2Inverses[i + 64][0] = PowersOf2_128DivisibilityData[i][0][0];
			powerOf2Inverses[i + 64][1] = PowersOf2_128DivisibilityData[i][0][1];
		}
		CL_CHECKED_CALL_WITH_RESULT(myPowerOf2SumInverses_Buf, clCreateBuffer, myGPUContext, CL_MEM_READ_ONLY, sizeof(powerOf2Inverses), 0);
		CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPowerOf2SumInverses_Buf, CL_TRUE, 0, sizeof(powerOf2Inverses), powerOf2Inverses, 0, nullptr, nullptr);
	}

	{
		static number buf[3][64][4];
		for (number i = 0; i < 3; ++i)
		{
			for (number j = 0; j < 64; ++j)
			{
				buf[i][j][0] = PowersOfP_128DivisibilityData[i][j].inverse[0];
				buf[i][j][1] = PowersOfP_128DivisibilityData[i][j].inverse[1];
				buf[i][j][2] = PowersOfP_128DivisibilityData[i][j].shift;
				buf[i][j][3] = 0;
			}
		}
		CL_CHECKED_CALL_WITH_RESULT(myPowersOfP_128SumInverses_Buf, clCreateBuffer, myGPUContext, CL_MEM_READ_ONLY, sizeof(buf), 0);
		CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPowersOfP_128SumInverses_Buf, CL_TRUE, 0, sizeof(buf), buf, 0, nullptr, nullptr);
	}

	{
		static number buf[3][2];
		for (number i = 0; i < 3; ++i)
		{
			buf[i][0] = PrimeInverses_128[i][0][0];
			buf[i][1] = PrimeInverses_128[i][0][1];
		}
		CL_CHECKED_CALL_WITH_RESULT(myPrimeInverses_128_Buf, clCreateBuffer, myGPUContext, CL_MEM_READ_ONLY, sizeof(buf), 0);
		CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPrimeInverses_128_Buf, CL_TRUE, 0, sizeof(buf), buf, 0, nullptr, nullptr);
	}

	CL_CHECKED_CALL_WITH_RESULT(myRangesTable_Buf, clCreateBuffer, myGPUContext, CL_MEM_READ_ONLY, sizeof(RangeDataGPU) * myMaxRangesCount, 0);
	CL_CHECKED_CALL_WITH_RESULT(myRangesLookupTable_Buf, clCreateBuffer, myGPUContext, CL_MEM_READ_ONLY, sizeof(LookupDataGPU) * myMaxLookupTableSize, 0);

	CL_CHECKED_CALL_WITH_RESULT(myPhase1_offset_to_resume_buf, clCreateBuffer, myGPUContext, CL_MEM_READ_WRITE, sizeof(number), 0);
	CL_CHECKED_CALL_WITH_RESULT(myPhase2_numbers_count_buf, clCreateBuffer, myGPUContext, CL_MEM_READ_WRITE, sizeof(unsigned int) * 2, 0);
	CL_CHECKED_CALL_WITH_RESULT(myPhase3_numbers_count_buf, clCreateBuffer, myGPUContext, CL_MEM_READ_WRITE, sizeof(number), 0);
	CL_CHECKED_CALL_WITH_RESULT(myAmicable_numbers_count_buf, clCreateBuffer, myGPUContext, CL_MEM_READ_WRITE, sizeof(number), 0);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPhase1_offset_to_resume_buf, CL_TRUE, 0, sizeof(number), &ourZeroBuf, 0, nullptr, nullptr);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPhase2_numbers_count_buf, CL_TRUE, 0, sizeof(unsigned int) * 2, &ourZeroBuf, 0, nullptr, nullptr);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPhase3_numbers_count_buf, CL_TRUE, 0, sizeof(number), &ourZeroBuf, 0, nullptr, nullptr);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myAmicable_numbers_count_buf, CL_TRUE, 0, sizeof(number), &ourZeroBuf, 0, nullptr, nullptr);

	CL_CHECKED_CALL_WITH_RESULT(myPhase2_numbers_buf, clCreateBuffer, myGPUContext, CL_MEM_READ_WRITE, sizeof(number) * 4 * myPhase2MaxNumbersCount, 0);
	CL_CHECKED_CALL_WITH_RESULT(myPhase3_numbers_buf, clCreateBuffer, myGPUContext, CL_MEM_READ_WRITE, sizeof(number) * 4 * myPhase2MaxNumbersCount, 0);
	CL_CHECKED_CALL_WITH_RESULT(myAmicable_numbers_data_buf, clCreateBuffer, myGPUContext, CL_MEM_READ_WRITE, sizeof(std::pair<number, number>) * 8192, 0);

	CL_CHECKED_CALL(clFinish, myQueue);

	CL_CHECKED_CALL(setKernelArguments, mySaveCounter, myPhase2_numbers_count_buf);

	if (IsLargePrimes)
	{
		CL_CHECKED_CALL_WITH_RESULT(myLargePrimesBuf, clCreateBuffer, myGPUContext, CL_MEM_READ_ONLY, sizeof(number) * myLargePrimesMaxCount, 0);

		CL_CHECKED_CALL(setKernelArguments, mySearchLargePrimes,
			mySmallPrimesBuf,
			myPrimeInversesBuf,
			myPQ_Buf,
			myPowerOf2SumInverses_Buf,
			myPowersOfP_128SumInverses_Buf,
			myPrimeInverses_128_Buf,
			myLargePrimesBuf,
			0U,
			0ULL,
			0U,
			0ULL,
			0ULL,
			myPhase1_offset_to_resume_buf,
			myPhase2_numbers_count_buf,
			myPhase2_numbers_buf,
			myAmicable_numbers_count_buf,
			myAmicable_numbers_data_buf
		);

		for (number i = 0; i < myMainBuffers.size(); ++i)
		{
			CL_CHECKED_CALL(clSetKernelArg, mySearchLargePrimes, static_cast<cl_uint>(i + 17), sizeof(cl_mem), &myMainBuffers[i]);
		}
	}
	else
	{
		CL_CHECKED_CALL(setKernelArguments, mySearchMultipleRanges,
			mySmallPrimesBuf,
			myPrimeInversesBuf,
			myPQ_Buf,
			myPowerOf2SumInverses_Buf,
			myPowersOfP_128SumInverses_Buf,
			myPrimeInverses_128_Buf,
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

		for (number i = 0; i < myMainBuffers.size(); ++i)
		{
			CL_CHECKED_CALL(clSetKernelArg, mySearchMultipleRanges, static_cast<cl_uint>(i + 16), sizeof(cl_mem), &myMainBuffers[i]);
		}
	}

	CL_CHECKED_CALL(setKernelArguments, myCheckPairPhase2, mySmallPrimesBuf, myPhase2_numbers_buf, 0, myPrimeInversesBuf, myPQ_Buf, myPhase3_numbers_count_buf, myPhase3_numbers_buf, myAmicable_numbers_count_buf, myAmicable_numbers_data_buf);

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
	Factor stopAtFactors[MaxPrimeFactors + 1] = {};

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

	if (!SaveProgressForRanges(r))
	{
		return false;
	}
	return true;
}

bool OpenCL::RunLargePrimes(number startPrime, number primeLimit)
{
	std::thread gpu_thread(&OpenCL::ProcessLargePrimesThread, this);

	if (!myLargePrimes)
	{
		myLargePrimes = reinterpret_cast<number*>(AllocateSystemMemory(sizeof(number) * myLargePrimesMaxCount, false));
	}

	primesieve::PrimeSieve sieve;
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

	return true;
}

NOINLINE void OpenCL::PassLargePrimesToThread()
{
	if (!myLargePrimesReady.Signal())
	{
		LOG_ERROR("Failed to signal semaphore");
		boinc_finish(-1);
	}

	if (!myLargePrimesReceived.Wait())
	{
		LOG_ERROR("Failed to wait for semaphore");
		boinc_finish(-1);
	}

	myLargePrimesCount = 0;
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
				LOG_ERROR("Failed to signal semaphore");
				return false;
			}
			return true;
		}

		const unsigned int LargePrimesCount = myLargePrimesCount;
		CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myLargePrimesBuf, CL_TRUE, 0, sizeof(number) * LargePrimesCount, myLargePrimes, 0, nullptr, nullptr);

		if (!myLargePrimesReceived.Signal())
		{
			LOG_ERROR("Failed to signal semaphore");
			return false;
		}

		number LargePrimesCountReciprocal = 0;
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

		const number LargestCandidate = SearchLimit::value / myLargePrimes[0];
		auto it = std::lower_bound(CandidatesData.begin(), CandidatesData.end(), LargestCandidate, [](const AmicableCandidate& a, number b) { return a.value <= b; });
		const number CandidatesCount = static_cast<number>(it - CandidatesData.begin());

		const number TotalNumbersToCheck = CandidatesCount * LargePrimesCount;

		CL_CHECKED_CALL(clSetKernelArg, mySearchLargePrimes, 7, sizeof(LargePrimesCount), &LargePrimesCount);
		CL_CHECKED_CALL(clSetKernelArg, mySearchLargePrimes, 8, sizeof(LargePrimesCountReciprocal), &LargePrimesCountReciprocal);
		CL_CHECKED_CALL(clSetKernelArg, mySearchLargePrimes, 9, sizeof(LargePrimesCountIncrementAndShift), &LargePrimesCountIncrementAndShift);
		CL_CHECKED_CALL(clSetKernelArg, mySearchLargePrimes, 11, sizeof(TotalNumbersToCheck), &TotalNumbersToCheck);

		cl_event event = nullptr;
		unsigned int numbers_in_phase2_before = 0;
		unsigned int numbers_in_phase2_after = 0;
		number phase1_offset_to_resume = 0;
		number global_offset = 0;

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

				CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPhase1_offset_to_resume_buf, CL_TRUE, 0, sizeof(number), &ourZeroBuf, 0, nullptr, nullptr);

				const number global_offset_before = global_offset;
				const unsigned int max_size_phase1 = GetMaxPhaseSize(TotalNumbersToCheck - global_offset, myPhase1MaxKernelSize);

				size_t globalSizePhase1;
				for (unsigned int k = 0; global_offset < TotalNumbersToCheck; ++k)
				{
					if (k)
					{
						CL_CHECKED_CALL(clFlush, myQueue);
					}
					globalSizePhase1 = std::min<number>(TotalNumbersToCheck - global_offset, max_size_phase1);
					if (globalSizePhase1 & (myWorkGroupSize - 1))
					{
						globalSizePhase1 += myWorkGroupSize - (globalSizePhase1 & (myWorkGroupSize - 1));
					}
					CL_CHECKED_CALL(clSetKernelArg, mySearchLargePrimes, 10, sizeof(global_offset), &global_offset);
					CL_CHECKED_CALL(clEnqueueTask, myQueue, mySaveCounter, 0, nullptr, nullptr);

					const number k1 = (global_offset >> 27);
					global_offset += globalSizePhase1;
					const number k2 = (global_offset >> 27);

					CL_CHECKED_CALL(clEnqueueNDRangeKernel, myQueue, mySearchLargePrimes, 1, nullptr, &globalSizePhase1, &myWorkGroupSize, 0, nullptr, ((global_offset >= TotalNumbersToCheck) || (k1 != k2)) ? &event : nullptr);

					// Wait for queue every 2^27 numbers
					if ((global_offset < TotalNumbersToCheck) && (k1 != k2))
					{
						LOG(2, "Phase 1: k = " << k << ", offset = " << global_offset << ", waiting for queue");
						if (!WaitForQueue(myQueue, event) || !GetCounter(myQueue, myPhase2_numbers_count_buf, numbers_in_phase2_after))
						{
							return false;
						}
						LOG(2, "Phase 1: k = " << k << ", offset = " << global_offset << ", finished waiting for queue, phase 2: " << numbers_in_phase2_after << " numbers");

						// Exit the loop if phase 2 is already full
						if (numbers_in_phase2_after > myPhase2MaxNumbersCount)
						{
							break;
						}
					}
				}

				if (!WaitForQueue(myQueue, event) ||
					!GetCounter(myQueue, myPhase2_numbers_count_buf, numbers_in_phase2_after) ||
					!GetAndResetCounter(myQueue, myPhase1_offset_to_resume_buf, phase1_offset_to_resume))
				{
					return false;
				}

				const double phase1Time = t.getElapsedTime();
				const number new_numbers_in_phase2 = numbers_in_phase2_after - numbers_in_phase2_before;
				const number numbers_processed = phase1_offset_to_resume ? (phase1_offset_to_resume - global_offset_before) : (TotalNumbersToCheck - global_offset_before);

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

			if (numbers_in_phase2_after >= myPhase2MinNumbersCount)
			{
				if (!ProcessNumbersPhases2_3(static_cast<unsigned int>(numbers_in_phase2_after)))
				{
					return false;
				}
			}
		} while (phase1_offset_to_resume);

		myNumbersProcessedTotal += TotalNumbersToCheck;
	}

	LOG_ERROR("Failed to wait for semaphore");
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

bool OpenCL::GetCounter(cl_command_queue queue, cl_mem buf, number &counter)
{
	CL_CHECKED_CALL(clEnqueueReadBuffer, queue, buf, CL_TRUE, 0, sizeof(counter), &counter, 0, nullptr, nullptr);
	return true;
}

bool OpenCL::ResetCounter(cl_command_queue queue, cl_mem buf)
{
	CL_CHECKED_CALL(clEnqueueWriteBuffer, queue, buf, CL_TRUE, 0, sizeof(number), &ourZeroBuf, 0, nullptr, nullptr);
	return true;
}

bool OpenCL::GetAndResetCounter(cl_command_queue queue, cl_mem buf, unsigned int &counter)
{
	CL_CHECKED_CALL(clEnqueueReadBuffer, queue, buf, CL_TRUE, 0, sizeof(counter), &counter, 0, nullptr, nullptr);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, queue, buf, CL_TRUE, 0, sizeof(counter), &ourZeroBuf, 0, nullptr, nullptr);
	return true;
}

bool OpenCL::GetAndResetCounter(cl_command_queue queue, cl_mem buf, number &counter)
{
	CL_CHECKED_CALL(clEnqueueReadBuffer, queue, buf, CL_TRUE, 0, sizeof(counter), &counter, 0, nullptr, nullptr);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, queue, buf, CL_TRUE, 0, sizeof(counter), &ourZeroBuf, 0, nullptr, nullptr);
	return true;
}

// Divide total_size into 2^k equal parts, each part's size is <= max_size and also a multiple of 1024
unsigned int OpenCL::GetMaxPhaseSize(number total_size, unsigned int max_size)
{
	for (unsigned int k = 1;; ++k)
	{
		const number result = (((total_size >> k) + 1023) / 1024) * 1024;
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
		SPairToCheck() : M(0), targetSumLow(0), targetSumHigh(0), found(0) {}
		SPairToCheck(number m, number sumLow, number sumHigh) : M(m), targetSumLow(sumLow), targetSumHigh(sumHigh), found(0) {}

		number M;
		number targetSumLow;
		number targetSumHigh;
		number found;
	};
	std::vector<SPairToCheck> pairs;

	const char* files[] = {
		//"c2_1.txt"
		"c2_3.txt", "c2_4.txt", "c2_5.txt", "c2_6.txt", "c2_7.txt", "c2_8.txt", "c2_9.txt", "c2_10.txt", "c2_11.txt", "c2_12.txt", "c2_13.txt", "c2_14.txt", "c2_15.txt", "c2_16.txt", "c2_17.txt", "c2_18.txt", "c2_19.txt", "c2_20.txt"
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

			number m[2];
			atoi128(buf, m[0], m[1]);

			if (m[1] || (m[0] >= SearchLimit::value))
				break;

			buf[0] = '\0';
			f.getline(buf, sizeof(buf)); // N and factorization

			number n[2];
			atoi128(buf, n[0], n[1]);

			buf[0] = '\0';
			f.getline(buf, sizeof(buf)); // empty line

			number sum[2];
			add128(m[0], m[1], n[0], n[1], sum, sum + 1);
			pairs.emplace_back(m[0], sum[0], sum[1]);
		}
		LOG(1, "Loaded " << name);
	}

	const unsigned int total_numbers = static_cast<unsigned int>(pairs.size());
	LOG(1, "Total numbers: " << total_numbers);

	Timer overall_timer;

	cl_mem pairsToCheck_Buf;
	pairs.resize(((pairs.size() + myWorkGroupSize - 1) / myWorkGroupSize) * myWorkGroupSize);
	CL_CHECKED_CALL_WITH_RESULT(pairsToCheck_Buf, clCreateBuffer, myGPUContext, CL_MEM_READ_WRITE, sizeof(SPairToCheck) * pairs.size(), 0);
	CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, pairsToCheck_Buf, CL_TRUE, 0, sizeof(SPairToCheck) * pairs.size(), pairs.data(), 0, nullptr, nullptr);

	struct TmpBuf
	{
		explicit TmpBuf(cl_context gpuContext) { mem = clCreateBuffer(gpuContext, CL_MEM_READ_WRITE, sizeof(std::pair<number, number>) * 2500000, 0, nullptr); }
		~TmpBuf() { clReleaseMemObject(mem); }

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
		CL_CHECKED_CALL(setKernelArguments, myCheckPairs, mySmallPrimesBuf, myPrimeInversesBuf, myPQ_Buf, myPowerOf2SumInverses_Buf, myPowersOfP_128SumInverses_Buf, myPrimeInverses_128_Buf, pairsToCheck_Buf, myPhase2_numbers_count_buf, myPhase2_numbers_buf, myAmicable_numbers_count_buf, buf);
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

		CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPhase2_numbers_buf, CL_TRUE, filtered_numbers_count * sizeof(number) * 4, myWorkGroupSize * sizeof(number) * 4, &ourZeroBuf, 0, nullptr, nullptr);
		CL_CHECKED_CALL(setKernelArguments, myCheckPairPhase2, mySmallPrimesBuf, myPhase2_numbers_buf, 0, myPrimeInversesBuf, myPQ_Buf, myPhase3_numbers_count_buf, myPhase3_numbers_buf, myAmicable_numbers_count_buf, buf);

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

		CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, phase3_buffers[phase3_iteration & 1], CL_TRUE, filtered_numbers3_count * sizeof(number) * 4, phase3_workGroupSize * sizeof(number) * 4, &ourZeroBuf, 0, nullptr, nullptr);

		const unsigned int max_size_phase3 = GetMaxPhaseSize(filtered_numbers3_count, 1 << 18);
		size_t globalSizePhase3;
		for (unsigned int offset = 0; offset < filtered_numbers3_count; offset += static_cast<unsigned int>(globalSizePhase3))
		{
			globalSizePhase3 = std::min(filtered_numbers3_count - offset, max_size_phase3);
			if (globalSizePhase3 & (phase3_workGroupSize - 1))
			{
				globalSizePhase3 += phase3_workGroupSize - (globalSizePhase3 & (phase3_workGroupSize - 1));
			}
			CL_CHECKED_CALL(setKernelArguments, myCheckPairPhase3, phase3_buffers[phase3_iteration & 1], offset, myPrimeInversesBuf, myPrimeReciprocalsBuf, myPhase3_numbers_count_buf, phase3_buffers[(phase3_iteration - 1) & 1], myAmicable_numbers_count_buf, buf, phase3_iteration << 12);
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
	if (amicable_numberscount > 0)
	{
		std::vector<std::pair<number, number>> found_pairs(amicable_numberscount);
		CL_CHECKED_CALL(clEnqueueReadBuffer, myQueue, buf, CL_TRUE, 0, sizeof(std::pair<number, number>) * amicable_numberscount, found_pairs.data(), 0, nullptr, nullptr);
		for (const std::pair<number, number>& pair : found_pairs)
		{
			auto it = std::lower_bound(pairs.begin(), pairs.end(), pair.first, [](const SPairToCheck& a, const number b) { return a.M < b; });
			if ((it != pairs.end()) && (it->M == pair.first))
			{
				it->found = 1;
			}
		}
	}

	for (const SPairToCheck& pair : pairs)
	{
		if (!pair.found)
		{
			LOG(1, "First missed number: " << pair.M);
			return false;
		}
	}

	return true;
}

// Fast integer cube root. Returns number m such that m^3 <= n < (m+1)^3 for all n > 0
FORCEINLINE number IntegerCbrt(const number n)
{
	if (n < 8)
	{
		return 1;
	}

	unsigned long index;
	_BitScanReverse64(&index, n);

	index = (index * ((65536 / 3) + 1)) >> 16;

	number result = number(1) << index;
	for (number cur_bit = result >> 1; cur_bit > 0; cur_bit >>= 1)
	{
		const number k = result | cur_bit;
		if ((k <= 2642245) && (k * k * k <= n))
		{
			result = k;
		}
	}
	return result;
}

bool OpenCL::AddRange(const RangeData& r)
{
	number prime_limit = 0;

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
			const number prime_limit2 = (r.sum - 1) / (r.value * 2 - r.sum);
			if (prime_limit > prime_limit2)
			{
				prime_limit = prime_limit2;
			}
		}
		break;

	case 2:
		prime_limit = static_cast<number>(sqrt(static_cast<double>(SearchLimit::value) / r.value));
		if (r.sum - r.value < r.value)
		{
			// r.sum * (p^2 + p + 1) = r.value * p^2 * 2
			// (r.sum - r.value * 2) * p^2 + r.sum * (p + 1) = 0
			// (r.value * 2 - r.sum) * p^2 - r.sum * (p + 1) = 0
			// (r.value * 2 - r.sum) * p^2 - r.sum * p - r.sum = 0
			// (r.value * 2 - r.sum) / r.sum * p^2 - p - 1 = 0
			// (r.value * 2 / r.sum - 1) * p^2 - p - 1 = 0
			const double A = static_cast<double>(r.value * 2 - r.sum) / r.sum;
			const number prime_limit2 = static_cast<number>((1.0 + sqrt(1.0 + 4.0 * A)) / (2.0 * A));
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
			const number A = r.value * 2 - r.sum;

			// Lower bound
			// A * p^3 = r.sum * p^2 < r.sum * (p^2 + p + 1)
			// A * p^3 = r.sum * p^2
			// A * p = r.sum
			// p = r.sum / A
			double x1 = static_cast<double>(r.sum) / A;

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
				if ((static_cast<number>(x1) == static_cast<number>(x2)) || (x == x1) || (x == x2))
				{
					break;
				}
				if (A*x*x*x > r.sum*(x*x + x + 1))
				{
					x2 = x;
				}
				else
				{
					x1 = x;
				}
			}

			const number prime_limit2 = static_cast<number>(x2);
			if (prime_limit > prime_limit2)
			{
				prime_limit = prime_limit2;
			}
		}
		break;
	}

	unsigned int index_begin = r.index_start_prime;

	// Find the first prime "p" such that "p > prime_limit" is true
	unsigned int a = 0;
	unsigned int b = NumPrimes;
	do
	{
		const unsigned int c = (a + b) >> 1;
		if (GetNthPrime(c) > prime_limit)
		{
			b = c;
		}
		else
		{
			a = c + 1;
		}
	} while (a < b);

	if (a > index_begin)
	{
		const unsigned int range_size = a - index_begin;

		myRanges.emplace_back(r.value, r.sum, index_begin, range_size);
		myTotalNumbersInRanges += range_size;
		myNumbersProcessedTotal += range_size;

		if ((myTotalNumbersInRanges >= (1 << 27)) || (myRanges.size() >= myMaxRangesCount))
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

	const number averageRangeSize = myTotalNumbersInRanges / myRanges.size();
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
	number numbers_in_phase2_before = 0;
	number numbers_in_phase2_after = 0;
	number phase1_offset_to_resume = 0;
	number global_offset = 0;

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

			CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPhase1_offset_to_resume_buf, CL_TRUE, 0, sizeof(number), &ourZeroBuf, 0, nullptr, nullptr);

			CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myRangesTable_Buf, CL_TRUE, 0, myRanges.size() * sizeof(RangeDataGPU), myRanges.data(), 0, nullptr, nullptr);
			CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myRangesLookupTable_Buf, CL_TRUE, 0, myLookupTable.size() * sizeof(LookupDataGPU), myLookupTable.data(), 0, nullptr, nullptr);

			CL_CHECKED_CALL(clSetKernelArg, mySearchMultipleRanges, 8, sizeof(unsigned int), &lookupTableShift);
			CL_CHECKED_CALL(clSetKernelArg, mySearchMultipleRanges, 10, sizeof(unsigned int), &myTotalNumbersInRanges);

			const number global_offset_before = global_offset;
			const unsigned int max_size_phase1 = GetMaxPhaseSize(myTotalNumbersInRanges - global_offset, myPhase1MaxKernelSize);

			size_t globalSizePhase1;
			for (unsigned int k = 0; global_offset < myTotalNumbersInRanges; global_offset += globalSizePhase1, ++k)
			{
				if (k)
				{
					CL_CHECKED_CALL(clFlush, myQueue);
				}
				globalSizePhase1 = std::min<number>(myTotalNumbersInRanges - global_offset, max_size_phase1);
				if (globalSizePhase1 & (myWorkGroupSize - 1))
				{
					globalSizePhase1 += myWorkGroupSize - (globalSizePhase1 & (myWorkGroupSize - 1));
				}
				CL_CHECKED_CALL(clSetKernelArg, mySearchMultipleRanges, 9, sizeof(global_offset), &global_offset);
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
			const number new_numbers_in_phase2 = numbers_in_phase2_after - numbers_in_phase2_before;
			const number numbers_processed = phase1_offset_to_resume ? (phase1_offset_to_resume - global_offset_before) : (myTotalNumbersInRanges - global_offset_before);

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

		if (numbers_in_phase2_after >= myPhase2MinNumbersCount)
		{
			if (!ProcessNumbersPhases2_3(static_cast<unsigned int>(numbers_in_phase2_after)))
			{
				return false;
			}
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
			CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, myPhase2_numbers_buf, CL_TRUE, numbers_in_phase2 * sizeof(number) * 4, std::min(static_cast<unsigned int>(myWorkGroupSize), myPhase2MaxNumbersCount - numbers_in_phase2) * sizeof(number) * 4, &ourZeroBuf, 0, nullptr, nullptr);
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
				CL_CHECKED_CALL(clEnqueueWriteBuffer, myQueue, phase3_buffers[phase3_iteration & 1], CL_TRUE, phase3_numbers_buf_count2 * sizeof(number) * 4, std::min(static_cast<unsigned int>(phase3_workGroupSize), myPhase2MaxNumbersCount - phase3_numbers_buf_count2) * sizeof(number) * 4, &ourZeroBuf, 0, nullptr, nullptr);
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
				globalSizePhase3 = std::min<number>(phase3_numbers_buf_count2 - global_offset, max_size_phase3);
				if (globalSizePhase3 & (myWorkGroupSize - 1))
				{
					globalSizePhase3 += myWorkGroupSize - (globalSizePhase3 & (myWorkGroupSize - 1));
				}
				CL_CHECKED_CALL(setKernelArguments, myCheckPairPhase3, phase3_buffers[phase3_iteration & 1], global_offset, myPrimeInversesBuf, myPrimeReciprocalsBuf, myPhase3_numbers_count_buf, phase3_buffers[(phase3_iteration - 1) & 1], myAmicable_numbers_count_buf, myAmicable_numbers_data_buf, phase3_iteration << 12);
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
		myBufUlong2.resize(numbers_found);
		CL_CHECKED_CALL(clEnqueueReadBuffer, myQueue, myAmicable_numbers_data_buf, CL_TRUE, 0, sizeof(std::pair<number, number>) * numbers_found, myBufUlong2.data(), 0, nullptr, nullptr);
		for (const auto& k : myBufUlong2)
		{
			if ((k.second == 0) || IsPrime(k.second))
			{
				++myAmicableNumbersFound;
				fprintf(g_outputFile, "%llu\n", k.first);
			}
		}
		fflush(g_outputFile);
	}

	return true;
}

bool OpenCL::SaveProgressForRanges(const RangeData& r)
{
	// This is only called when the queue is empty: there are no numbers in phases 1-3 on GPU and GPU is idle

	if (!SaveFoundNumbers())
	{
		return false;
	}

	const double progress = static_cast<double>(myNumbersProcessedTotal) / RangeGen::total_numbers_to_check;
	LOG(1, "Progress: " << (progress * 100.0) << "%, " << myAmicableNumbersFound << " amicable numbers found   \r");

	if (boinc_time_to_checkpoint())
	{
		std::stringstream s;
		std::string checkpoint_name, data;

		boinc_resolve_filename_s(CHECKPOINT_LOGICAL_NAME, checkpoint_name);

		// First write how many numbers have been checked so far
		s << myNumbersProcessedTotal << ' ';

		// Then current range
		const Factor(&mainFactors)[MaxPrimeFactors] = r.factors;
		for (number i = 0; i <= static_cast<number>(r.last_factor_index); ++i)
		{
			if (i > 0)
			{
				s << '*';
			}
			s << mainFactors[i].p;
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
