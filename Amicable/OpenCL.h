#pragma once

#include <CL/opencl.h>

extern FILE* g_outputFile;

struct RangeData;

struct RangeDataGPU
{
	RangeDataGPU(number _M0, number _sumM0, unsigned int _start_prime_index, unsigned int _size) : M0(_M0), sumM0(_sumM0), start_prime_index(_start_prime_index), size(_size) {}

	number M0;
	number sumM0;
	number start_prime_index;
	number size;
};

static_assert(sizeof(RangeDataGPU) == sizeof(number) * 4, "Invalid RangeDataGPU size");

struct LookupDataGPU
{
	LookupDataGPU() : range_index(0), range_offset(0) {}
	LookupDataGPU(unsigned int i, unsigned int offset) : range_index(i), range_offset(offset) {}

	unsigned int range_index;
	unsigned int range_offset;
};

static_assert(sizeof(LookupDataGPU) == sizeof(unsigned int) * 2, "Invalid LookupDataGPU size");

class OpenCL
{
public:
	OpenCL(const char* preferences);
	~OpenCL();

	bool Run(int argc, char* argv[], char* startFrom, char* stopAt, unsigned int largestPrimePower);

	bool Test();

	bool AddRange(const RangeData& r);

private:
	int SetKernelSize(int size);

	bool GetPlatformID(cl_platform_id* clSelectedPlatformID);
	bool WaitForQueue(cl_command_queue queue, cl_event event);
	bool GetCounter(cl_command_queue queue, cl_mem buf, unsigned int &counter);
	bool ResetCounter(cl_command_queue queue, cl_mem buf);
	bool GetAndResetCounter(cl_command_queue queue, cl_mem buf, unsigned int &counter);
	unsigned int GetMaxPhaseSize(unsigned int total_size, unsigned int max_size);

	bool ProcessNumbers();
	bool ProcessNumbersPhases2_3(unsigned int numbers_in_phase2);
	void CleanupRanges();

	bool SaveProgress(const RangeData& r);

	unsigned int myLargestPrimePower;

	cl_context myGPUContext;
	cl_program myProgram;

	cl_command_queue myQueue;
	cl_kernel mySaveCounter;
	cl_kernel mySearchMultipleRanges;
	cl_kernel myCheckPairs;
	cl_kernel myCheckPairPhase2;
	cl_kernel myCheckPairPhase3;
	cl_mem mySmallPrimesBuf;
	std::vector<cl_mem> myPrimesBuffers;
	cl_mem myPrimeInversesBuf;
	cl_mem myPrimeReciprocalsBuf;
	cl_mem myPQ_Buf;
	cl_mem myPowerOf2SumInverses_Buf;
	cl_mem myPowersOfP_128SumInverses_Buf;
	cl_mem myPrimeInverses_128_Buf;
	cl_mem myRangesTable_Buf;
	cl_mem myRangesLookupTable_Buf;
	cl_mem myPhase1_offset_to_resume_buf;
	cl_mem myPhase2_numbers_count_buf;
	cl_mem myPhase2_numbers_buf;
	cl_mem myPhase3_numbers_count_buf;
	cl_mem myPhase3_numbers_buf;
	cl_mem myAmicable_numbers_count_buf;
	cl_mem myAmicable_numbers_data_buf;

	size_t myWorkGroupSize;

	unsigned int myMaxRangesCount;
	unsigned int myMaxLookupTableSize;
	
	unsigned int myPhase1MaxKernelSize;
	unsigned int myPhase2MinNumbersCount;
	unsigned int myPhase2MaxNumbersCount;

	number myOldTimerResolution;

	std::vector<RangeDataGPU> myRanges;
	number myTotalNumbersInRanges;
	number myNumbersProcessedTotal;
	number myAmicableNumbersFound;

	const char* myPreferences;

	std::vector<LookupDataGPU> myLookupTable;

	std::streamsize myStdOldPrecision;

	std::vector<std::pair<number, number>> myBufUlong2;

	static const unsigned char ourZeroBuf[16384];
};

#if ((defined __cpp_variadic_templates) && (__cpp_variadic_templates >= 200704))

FORCEINLINE cl_int setKernelArgumentsWithOffset(cl_kernel, cl_uint) { return CL_SUCCESS; }

template<typename T0, typename... Arguments>
FORCEINLINE cl_int setKernelArgumentsWithOffset(cl_kernel kernel, cl_uint offset, const T0& arg0, Arguments... args)
{
	cl_int ciErrNum = clSetKernelArg(kernel, offset, sizeof(T0), &arg0);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	return setKernelArgumentsWithOffset(kernel, offset + 1, args...);
}

template<typename... Arguments>
FORCEINLINE cl_int setKernelArguments(cl_kernel kernel, Arguments... args)
{
	return setKernelArgumentsWithOffset(kernel, 0, args...);
}

#else
#include "OpenCL.inl"
#endif
