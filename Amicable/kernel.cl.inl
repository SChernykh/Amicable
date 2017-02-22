#pragma once

static const char* kernel_cl = "#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable\n"\
"\n"\
"//-----------------------------------------------------------------------------------------------------------------------------------------------------------\n"\
"// Check amicability for a range of numbers of the form \"M0*p\" where p is prime\n"\
"// p runs through all values in primes[p_offset]...primes[p_offset + get_global_size(0) - 1] (inclusive)\n"\
"// primes stores prime numbers in compact form: 4 primes per each 8 bytes\n"\
"// primeInverses is an array of multiplicative inverses of first 192725 prime numbers, this is used for checking divisibility by doing just one multiplication\n"\
"// primeReciprocals is an array of reciprocals of first 192725 prime numbers, this is used for calculating N / p (not just checking divisibility)\n"\
"// PQ corresponds to tables P and Q in lemma 2.1 from \"Computation of All the Amicable Pairs Below 10^10 By H.J.J.te Riele\": http://www.ams.org/journals/mcom/1986-47-175/S0025-5718-1986-0842142-3/S0025-5718-1986-0842142-3.pdf\n"\
"// The only difference is that we calculate exact (hence better) upper bounds for S(m)/m instead of inexact estimates\n"\
"//-----------------------------------------------------------------------------------------------------------------------------------------------------------\n"\
"\n"\
"// A: smaller member of an amicable pair\n"\
"// B: if it's not zero, it must be prime in order for A to be amicable number\n"\
"void NumberFound(__global uint* amicable_numbers_count, __global ulong2* amicable_numbers_data, ulong A, ulong B)\n"\
"{\n"\
"	const uint index = atom_inc(amicable_numbers_count);\n"\
"	amicable_numbers_data[index].x = A; \n"\
"	amicable_numbers_data[index].y = B;\n"\
"}\n"\
"\n"\
"// 32-bit floats don't have enough precision\n"\
"// 64-bit floats are not supported by all OpenCL devices\n"\
"// So we have to use integer arithmetic to calculate square roots\n"\
"\n"\
"// Returns number x such that x^2 <= n < (x+1)^2\n"\
"uint IntegerSquareRoot(const ulong n)\n"\
"{\n"\
"	uint result = 1;\n"\
"	result <<= ((63 - clz(n)) >> 1);\n"\
"	for (uint cur_bit = result >> 1; cur_bit > 0; cur_bit >>= 1)\n"\
"	{\n"\
"		const uint result1 = result | cur_bit;\n"\
"		const ulong k = result1;\n"\
"		if (k * k <= n)\n"\
"		{\n"\
"			result = result1;\n"\
"		}\n"\
"	}\n"\
"	return result;\n"\
"}\n"\
"\n"\
"// Returns number x such that x^2 <= n < (x+1)^2\n"\
"ulong IntegerSquareRoot128(const ulong2 n)\n"\
"{\n"\
"	ulong result = 1;\n"\
"	const ulong highest_bit_index = n.y ? (127 - clz(n.y)) : (63 - clz(n.x));\n"\
"	result <<= (highest_bit_index >> 1);\n"\
"	for (ulong cur_bit = result >> 1; cur_bit > 0; cur_bit >>= 1)\n"\
"	{\n"\
"		const ulong k = result | cur_bit;\n"\
"		ulong2 k2;\n"\
"		k2.x = k * k;\n"\
"		k2.y = mul_hi(k, k);\n"\
"		if ((k2.y < n.y) || ((k2.y == n.y) && (k2.x <= n.x)))\n"\
"		{\n"\
"			result = k;\n"\
"		}\n"\
"	}\n"\
"	return result;\n"\
"}\n"\
"\n"\
"void CheckPairPhase1(\n"\
"	__global const uint* smallPrimes,\n"\
"	__global const ulong2* primeInverses,\n"\
"	__global const ulong2* PQ,\n"\
"	__constant const ulong2* PowerOf2SumInverses,\n"\
"	__constant const ulong4* PowersOfP_128SumInverses,\n"\
"	__constant const ulong2* PrimeInverses_128,\n"\
"	const ulong M,\n"\
"	ulong targetSum,\n"\
"	ulong targetSumHigh,\n"\
"	const uint global_offset,\n"\
"	__global uint* phase1_offset_to_resume_after_overflow,\n"\
"	__global uint* phase2_numbers_count,\n"\
"	__global ulong4* phase2_numbers_data,\n"\
"	__global uint* amicable_numbers_count,\n"\
"	__global ulong2* amicable_numbers_data\n"\
")\n"\
"{\n"\
"	ulong N = targetSum - M;\n"\
"	ulong N_high = targetSumHigh - ((targetSum < M) ? 1 : 0);\n"\
"\n"\
"	// Divide out power of 2 using fast bitwise operations\n"\
"	// All numbers in the same range have the same parity, so threads won't diverge here\n"\
"	if ((N & 1) == 0)\n"\
"	{\n"\
"		if (N == 0)\n"\
"		{\n"\
"			// N = 2^64*k, k >= 0. It's not amicable for k = 0, 1, 2 and too large for k > 2.\n"\
"			return;\n"\
"		}\n"\
"\n"\
"		const uint powerOf2 = 63 - clz(N ^ (N - 1));\n"\
"\n"\
"		// M and targetSum (128-bit) come in ranges where they both grow monotonically,\n"\
"		// so there will be exactly one divergence in each range and one divergence on a border between ranges\n"\
"		// Threads don't diverge 99.9% of the time here\n"\
"		if (targetSumHigh > 0)\n"\
"		{\n"\
"			const ulong2 inv = PowerOf2SumInverses[powerOf2 + 64];\n"\
"			if (mul_hi(targetSum, inv.x) + targetSum * inv.y + targetSumHigh * inv.x != 0)\n"\
"			{\n"\
"				return;\n"\
"			}\n"\
"\n"\
"			targetSum = targetSum * inv.x;\n"\
"			targetSumHigh = 0;\n"\
"\n"\
"			N = (N >> powerOf2) | (N_high << (64 - powerOf2));\n"\
"		}\n"\
"		else\n"\
"		{\n"\
"			targetSum *= PowerOf2SumInverses[powerOf2].x;\n"\
"			if (targetSum > PowerOf2SumInverses[powerOf2].y)\n"\
"			{\n"\
"				return;\n"\
"			}\n"\
"			N >>= powerOf2;\n"\
"		}\n"\
"	}\n"\
"\n"\
"	int i = 1;\n"\
"\n"\
"	// M and targetSum (128-bit) come in ranges where they both grow monotonically,\n"\
"	// so there will be exactly one divergence in each range and one divergence on a border between ranges\n"\
"	// Threads don't diverge 99.9% of the time here\n"\
"	if (targetSumHigh > 0)\n"\
"	{\n"\
"		do\n"\
"		{\n"\
"			const ulong2 inv = PrimeInverses_128[i - 1];\n"\
"			if (mul_hi(N, inv.x) + N * inv.y + N_high * inv.x == 0)\n"\
"			{\n"\
"				ulong q = N * inv.x;\n"\
"				int powerOfP = 0;\n"\
"				do\n"\
"				{\n"\
"					++powerOfP;\n"\
"					N = q;\n"\
"					q *= primeInverses[i].x;\n"\
"				} while (q <= primeInverses[i].y);\n"\
"\n"\
"				const ulong4 inv2 = PowersOfP_128SumInverses[(i - 1) * 64 + powerOfP];\n"\
"				if (inv2.z != 0)\n"\
"				{\n"\
"					if (targetSum & ((1 << inv2.z) - 1))\n"\
"					{\n"\
"						return;\n"\
"					}\n"\
"					targetSum = (targetSum >> inv2.z) | (targetSumHigh << (64 - inv2.z));\n"\
"					targetSumHigh >>= inv2.z;\n"\
"				}\n"\
"\n"\
"				if (mul_hi(targetSum, inv2.x) + targetSum * inv2.y + targetSumHigh * inv2.x != 0)\n"\
"				{\n"\
"					return;\n"\
"				}\n"\
"\n"\
"				targetSum *= inv2.x;\n"\
"				targetSumHigh = 0;\n"\
"\n"\
"				const ulong p2 = (0x31190900 >> (i << 3)) & 0xFF;\n"\
"				// i = 1: (0x31190900 >>  8) & 0xFF =  9 = 3^2\n"\
"				// i = 2: (0x31190900 >> 16) & 0xFF = 25 = 5^2\n"\
"				// i = 3: (0x31190900 >> 24) & 0xFF = 49 = 7^2\n"\
"				if (p2 > N)\n"\
"				{\n"\
"					const ulong sumN = ((N > 1) ? N : 0) + 1;\n"\
"					if (sumN == targetSum)\n"\
"					{\n"\
"						NumberFound(amicable_numbers_count, amicable_numbers_data, M, 0);\n"\
"					}\n"\
"					return;\n"\
"				}\n"\
"\n"\
"				if (N >= targetSum)\n"\
"				{\n"\
"					return;\n"\
"				}\n"\
"			}\n"\
"			++i;\n"\
"		} while ((i <= 3) && (targetSumHigh != 0));\n"\
"	}\n"\
"\n"\
"	// Collect divisibility data for first 15 odd primes\n"\
"	// No divergence here, compiler will use predicates for \"A ? B : C\" operator\n"\
"	uint pattern =\n"\
"		((N * primeInverses[ 1].x <= primeInverses[ 1].y) ? 16384 : 0) |\n"\
"		((N * primeInverses[ 2].x <= primeInverses[ 2].y) ? 8192  : 0) |\n"\
"		((N * primeInverses[ 3].x <= primeInverses[ 3].y) ? 4096  : 0) |\n"\
"		((N * primeInverses[ 4].x <= primeInverses[ 4].y) ? 2048  : 0) |\n"\
"		((N * primeInverses[ 5].x <= primeInverses[ 5].y) ? 1024  : 0) |\n"\
"		((N * primeInverses[ 6].x <= primeInverses[ 6].y) ? 512   : 0) |\n"\
"		((N * primeInverses[ 7].x <= primeInverses[ 7].y) ? 256   : 0) |\n"\
"		((N * primeInverses[ 8].x <= primeInverses[ 8].y) ? 128   : 0) |\n"\
"		((N * primeInverses[ 9].x <= primeInverses[ 9].y) ? 64    : 0) |\n"\
"		((N * primeInverses[10].x <= primeInverses[10].y) ? 32    : 0) |\n"\
"		((N * primeInverses[11].x <= primeInverses[11].y) ? 16    : 0) |\n"\
"		((N * primeInverses[12].x <= primeInverses[12].y) ? 8     : 0) |\n"\
"		((N * primeInverses[13].x <= primeInverses[13].y) ? 4     : 0) |\n"\
"		((N * primeInverses[14].x <= primeInverses[14].y) ? 2     : 0) |\n"\
"		((N * primeInverses[15].x <= primeInverses[15].y) ? 1     : 0);\n"\
"\n"\
"	i = 16;\n"\
"\n"\
"	// Do trial divisions by primes found\n"\
"	if (pattern)\n"\
"	{\n"\
"		ulong sumN = 1;\n"\
"		do\n"\
"		{\n"\
"			// Bit 0: divisibility for i == 15\n"\
"			// Bit 1: divisibility for i == 14\n"\
"			// Bit 2: divisibility for i == 13\n"\
"			// and so on ...\n"\
"			const uint num_leading_zeroes = clz(pattern);\n"\
"			pattern -= (1 << (31 - num_leading_zeroes));\n"\
"			const uint index = num_leading_zeroes - 16;\n"\
"\n"\
"			const ulong prevSumN = sumN;\n"\
"			const uint p1 = smallPrimes[index];\n"\
"			ulong q = N * primeInverses[index].x;\n"\
"			do\n"\
"			{\n"\
"				sumN = sumN * p1 + prevSumN;\n"\
"				N = q;\n"\
"				q *= primeInverses[index].x;\n"\
"			} while (q <= primeInverses[index].y);\n"\
"		} while (pattern);\n"\
"\n"\
"		if (N < 59 * 59)\n"\
"		{\n"\
"			if (N > 1)\n"\
"			{\n"\
"				sumN *= N + 1;\n"\
"			}\n"\
"			if (sumN == targetSum)\n"\
"			{\n"\
"				NumberFound(amicable_numbers_count, amicable_numbers_data, M, 0);\n"\
"			}\n"\
"			return;\n"\
"		}\n"\
"\n"\
"		const ulong oldTargetSum = targetSum;\n"\
"		targetSum /= sumN;\n"\
"		if (targetSum * sumN != oldTargetSum)\n"\
"		{\n"\
"			return;\n"\
"		}\n"\
"\n"\
"		if (N >= targetSum)\n"\
"		{\n"\
"			return;\n"\
"		}\n"\
"	}\n"\
"\n"\
"	// Check maximal possible sum for N. If it's too small, it can't be an amicable number\n"\
"	int k = 8;\n"\
"	while (N <= PQ[k].x)\n"\
"	{\n"\
"		--k;\n"\
"	}\n"\
"	++k;\n"\
"	if (mul_hi(N, PQ[k].y) < targetSum - N)\n"\
"	{\n"\
"		return;\n"\
"	}\n"\
"\n"\
"	do\n"\
"	{\n"\
"		// Collect divisibility data for the next 16 primes\n"\
"		// No divergence here, compiler will use predicates for \"A ? B : C\" operator\n"\
"		pattern =\n"\
"			((N * primeInverses[i +  0].x <= primeInverses[i +  0].y) ? 32768 : 0) |\n"\
"			((N * primeInverses[i +  1].x <= primeInverses[i +  1].y) ? 16384 : 0) |\n"\
"			((N * primeInverses[i +  2].x <= primeInverses[i +  2].y) ? 8192  : 0) |\n"\
"			((N * primeInverses[i +  3].x <= primeInverses[i +  3].y) ? 4096  : 0) |\n"\
"			((N * primeInverses[i +  4].x <= primeInverses[i +  4].y) ? 2048  : 0) |\n"\
"			((N * primeInverses[i +  5].x <= primeInverses[i +  5].y) ? 1024  : 0) |\n"\
"			((N * primeInverses[i +  6].x <= primeInverses[i +  6].y) ? 512   : 0) |\n"\
"			((N * primeInverses[i +  7].x <= primeInverses[i +  7].y) ? 256   : 0) |\n"\
"			((N * primeInverses[i +  8].x <= primeInverses[i +  8].y) ? 128   : 0) |\n"\
"			((N * primeInverses[i +  9].x <= primeInverses[i +  9].y) ? 64    : 0) |\n"\
"			((N * primeInverses[i + 10].x <= primeInverses[i + 10].y) ? 32    : 0) |\n"\
"			((N * primeInverses[i + 11].x <= primeInverses[i + 11].y) ? 16    : 0) |\n"\
"			((N * primeInverses[i + 12].x <= primeInverses[i + 12].y) ? 8     : 0) |\n"\
"			((N * primeInverses[i + 13].x <= primeInverses[i + 13].y) ? 4     : 0) |\n"\
"			((N * primeInverses[i + 14].x <= primeInverses[i + 14].y) ? 2     : 0) |\n"\
"			((N * primeInverses[i + 15].x <= primeInverses[i + 15].y) ? 1     : 0);\n"\
"\n"\
"		i += 16;\n"\
"\n"\
"		// Do trial divisions by prime factors found\n"\
"		if (pattern)\n"\
"		{\n"\
"			ulong sumN = 1;\n"\
"			do\n"\
"			{\n"\
"				// Bit 0: divisibility for i - 1\n"\
"				// Bit 1: divisibility for i - 2\n"\
"				// Bit 2: divisibility for i - 3\n"\
"				// and so on ...\n"\
"				const uint num_leading_zeroes = clz(pattern);\n"\
"				pattern -= (1 << (31 - num_leading_zeroes));\n"\
"				const uint index = i + num_leading_zeroes - 32;\n"\
"\n"\
"				const ulong prevSumN = sumN;\n"\
"				const uint p1 = smallPrimes[index];\n"\
"				ulong q = N * primeInverses[index].x;\n"\
"				do\n"\
"				{\n"\
"					sumN = sumN * p1 + prevSumN;\n"\
"					N = q;\n"\
"					q *= primeInverses[index].x;\n"\
"				} while (q <= primeInverses[index].y);\n"\
"			} while (pattern);\n"\
"\n"\
"			const ulong p2 = smallPrimes[i];\n"\
"			if (p2 * p2 > N)\n"\
"			{\n"\
"				if (N > 1)\n"\
"				{\n"\
"					sumN *= N + 1;\n"\
"				}\n"\
"				if (sumN == targetSum)\n"\
"				{\n"\
"					NumberFound(amicable_numbers_count, amicable_numbers_data, M, 0);\n"\
"				}\n"\
"				return;\n"\
"			}\n"\
"\n"\
"			const ulong oldTargetSum = targetSum;\n"\
"			targetSum /= sumN;\n"\
"			if (targetSum * sumN != oldTargetSum)\n"\
"			{\n"\
"				return;\n"\
"			}\n"\
"\n"\
"			if (N >= targetSum)\n"\
"			{\n"\
"				return;\n"\
"			}\n"\
"		}\n"\
"\n"\
"		// Check maximal possible sum for N. If it's too small, it can't be an amicable number\n"\
"		// Same i for all threads, no divergence\n"\
"		uint index = k * PQ_STRIDE_SIZE + (i >> 4) + 16;\n"\
"		while (N <= PQ[index].x)\n"\
"		{\n"\
"			--k;\n"\
"			index -= PQ_STRIDE_SIZE;\n"\
"		}\n"\
"		if (mul_hi(N, PQ[index].y) < targetSum - N)\n"\
"		{\n"\
"			return;\n"\
"		}\n"\
"	} while (i != 256);\n"\
"\n"\
"	// Data needed to resume from arbitrary loop iteration in phase 2:\n"\
"	// 1) M (8 bytes)\n"\
"	// 2) N (8 bytes)\n"\
"	// 3) targetSum (8 bytes)\n"\
"	// 4) i (<= 192724, can fit in 18 bits)\n"\
"	// 5) k (< 16, can fit in 4 bits)\n"\
"	// It can all fit in 32 bytes\n"\
"\n"\
"	const uint phase2_number_index = atom_inc(phase2_numbers_count);\n"\
"	if (phase2_number_index < PHASE2_MAX_COUNT)\n"\
"	{\n"\
"		const ulong i_to_save = i;\n"\
"		const ulong k_to_save = k;\n"\
"		phase2_numbers_data[phase2_number_index].x = M;\n"\
"		phase2_numbers_data[phase2_number_index].y = N;\n"\
"		phase2_numbers_data[phase2_number_index].z = targetSum;\n"\
"		phase2_numbers_data[phase2_number_index].w = i_to_save | (k_to_save << 32);\n"\
"	}\n"\
"	else if (phase2_number_index == PHASE2_MAX_COUNT)\n"\
"	{\n"\
"		*phase1_offset_to_resume_after_overflow = global_offset;\n"\
"	}\n"\
"}\n"\
"\n"\
"__kernel\n"\
"__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))\n"\
"void CheckPairPhase2(\n"\
"	__global const uint* smallPrimes,\n"\
"	__global const ulong4* numbersToCheck,\n"\
"	__global const ulong2* primeInverses,\n"\
"	__global const ulong2* PQ,\n"\
"	__global uint* filtered_numbers_count,\n"\
"	__global ulong4* filtered_numbers_data,\n"\
"	__global uint* amicable_numbers_count,\n"\
"	__global ulong2* amicable_numbers_data\n"\
")\n"\
"{\n"\
"	const ulong4 cur_number = numbersToCheck[get_global_id(0)];\n"\
"	if (cur_number.x == 0)\n"\
"	{\n"\
"		return;\n"\
"	}\n"\
"\n"\
"	const ulong M = cur_number.x;\n"\
"	ulong N = cur_number.y;\n"\
"	ulong targetSum = cur_number.z;\n"\
"	int i = cur_number.w & 0xFFFFFFFF;\n"\
"	int k = cur_number.w >> 32;\n"\
"\n"\
"	const uint root4_N_rounded_up = convert_uint_rtp(sqrt(sqrt(convert_float_rtp(N))));\n"\
"\n"\
"	// All threads start with the same value of \"i\" here, which is set to 256 at the end of phase 1,\n"\
"	// so there will be no divergence because of different values of \"i\"\n"\
"	do\n"\
"	{\n"\
"		// Collect divisibility data for the next 16 primes\n"\
"		// No divergence here, compiler will use predicates for \"A ? B : C\" operator\n"\
"		uint pattern =\n"\
"			((N * primeInverses[i +  0].x <= primeInverses[i +  0].y) ? 32768 : 0) |\n"\
"			((N * primeInverses[i +  1].x <= primeInverses[i +  1].y) ? 16384 : 0) |\n"\
"			((N * primeInverses[i +  2].x <= primeInverses[i +  2].y) ? 8192  : 0) |\n"\
"			((N * primeInverses[i +  3].x <= primeInverses[i +  3].y) ? 4096  : 0) |\n"\
"			((N * primeInverses[i +  4].x <= primeInverses[i +  4].y) ? 2048  : 0) |\n"\
"			((N * primeInverses[i +  5].x <= primeInverses[i +  5].y) ? 1024  : 0) |\n"\
"			((N * primeInverses[i +  6].x <= primeInverses[i +  6].y) ? 512   : 0) |\n"\
"			((N * primeInverses[i +  7].x <= primeInverses[i +  7].y) ? 256   : 0) |\n"\
"			((N * primeInverses[i +  8].x <= primeInverses[i +  8].y) ? 128   : 0) |\n"\
"			((N * primeInverses[i +  9].x <= primeInverses[i +  9].y) ? 64    : 0) |\n"\
"			((N * primeInverses[i + 10].x <= primeInverses[i + 10].y) ? 32    : 0) |\n"\
"			((N * primeInverses[i + 11].x <= primeInverses[i + 11].y) ? 16    : 0) |\n"\
"			((N * primeInverses[i + 12].x <= primeInverses[i + 12].y) ? 8     : 0) |\n"\
"			((N * primeInverses[i + 13].x <= primeInverses[i + 13].y) ? 4     : 0) |\n"\
"			((N * primeInverses[i + 14].x <= primeInverses[i + 14].y) ? 2     : 0) |\n"\
"			((N * primeInverses[i + 15].x <= primeInverses[i + 15].y) ? 1     : 0);\n"\
"\n"\
"		i += 16;\n"\
"\n"\
"		// Do trial divisions by prime factors found\n"\
"		if (pattern)\n"\
"		{\n"\
"			ulong sumN = 1;\n"\
"			do\n"\
"			{\n"\
"				// Bit 0: divisibility for i - 1\n"\
"				// Bit 1: divisibility for i - 2\n"\
"				// Bit 2: divisibility for i - 3\n"\
"				// and so on ...\n"\
"				const uint num_leading_zeroes = clz(pattern);\n"\
"				pattern -= (1 << (31 - num_leading_zeroes));\n"\
"				const uint index = i + num_leading_zeroes - 32;\n"\
"\n"\
"				const ulong prevSumN = sumN;\n"\
"				const uint p1 = smallPrimes[index];\n"\
"				ulong q = N * primeInverses[index].x;\n"\
"				do\n"\
"				{\n"\
"					sumN = sumN * p1 + prevSumN;\n"\
"					N = q;\n"\
"					q *= primeInverses[index].x;\n"\
"				} while (q <= primeInverses[index].y);\n"\
"			} while (pattern);\n"\
"\n"\
"			const ulong p2 = smallPrimes[i];\n"\
"			if (p2 * p2 > N)\n"\
"			{\n"\
"				if (N > 1)\n"\
"				{\n"\
"					sumN *= N + 1;\n"\
"				}\n"\
"				if (sumN == targetSum)\n"\
"				{\n"\
"					NumberFound(amicable_numbers_count, amicable_numbers_data, M, 0);\n"\
"				}\n"\
"				return;\n"\
"			}\n"\
"\n"\
"			const ulong oldTargetSum = targetSum;\n"\
"			targetSum /= sumN;\n"\
"			if (targetSum * sumN != oldTargetSum)\n"\
"			{\n"\
"				return;\n"\
"			}\n"\
"\n"\
"			if (N >= targetSum)\n"\
"			{\n"\
"				return;\n"\
"			}\n"\
"		}\n"\
"\n"\
"		uint index = k * PQ_STRIDE_SIZE + (i >> 4) + 16;\n"\
"		while (N <= PQ[index].x)\n"\
"		{\n"\
"			--k;\n"\
"			index -= PQ_STRIDE_SIZE;\n"\
"		}\n"\
"		if (mul_hi(N, PQ[index].y) < targetSum - N)\n"\
"		{\n"\
"			return;\n"\
"		}\n"\
"	} while (smallPrimes[i] <= root4_N_rounded_up);\n"\
"\n"\
"	// No more than 3 factors remain here\n"\
"	\n"\
"	// The case when N is prime (it has 1 factor)\n"\
"	if (N + 1 == targetSum)\n"\
"	{\n"\
"		NumberFound(amicable_numbers_count, amicable_numbers_data, M, N);\n"\
"		return;\n"\
"	}\n"\
"\n"\
"	const ulong root3_N_rounded_down = convert_ulong_rtn(cbrt(convert_float_rtn(N)));\n"\
"\n"\
"	// The case when N can only have 2 factors\n"\
"	if (N >= targetSum - root3_N_rounded_down * (root3_N_rounded_down + 1))\n"\
"	{\n"\
"		// The case when N is squared prime\n"\
"		const ulong sqrt_N = IntegerSquareRoot(N);\n"\
"		if ((sqrt_N * sqrt_N == N) && (N + sqrt_N + 1 == targetSum))\n"\
"		{\n"\
"			NumberFound(amicable_numbers_count, amicable_numbers_data, M, sqrt_N);\n"\
"			return;\n"\
"		}\n"\
"\n"\
"		// The case when N = p * q - a product of two different primes\n"\
"		// p + q == targetSum - N - 1\n"\
"		// p * q == N\n"\
"		const ulong B = targetSum - N - 1;\n"\
"		ulong2 B2;\n"\
"		B2.x = B * B;\n"\
"		B2.y = mul_hi(B, B);\n"\
"		const ulong carry = (B2.x < (N << 2)) ? 1 : 0;\n"\
"		B2.x -= (N << 2);\n"\
"		if (B2.y >= (N >> 62) + carry)\n"\
"		{\n"\
"			B2.y -= (N >> 62) + carry;\n"\
"			const ulong sqrt_D = IntegerSquareRoot128(B2);\n"\
"			const ulong p = (B - sqrt_D) >> 1;\n"\
"			const ulong q = (B + sqrt_D) >> 1;\n"\
"			if ((p < q) && (p + q == B) && (p * q == N))\n"\
"			{\n"\
"				NumberFound(amicable_numbers_count, amicable_numbers_data, M, q);\n"\
"				return;\n"\
"			}\n"\
"		}\n"\
"		return;\n"\
"	}\n"\
"\n"\
"	// N can have 3 factors\n"\
"	// Check maximal possible sum of divisors for numbers with 3 factors\n"\
"	uint p1 = smallPrimes[i];\n"\
"	ulong q = N / p1;\n"\
"\n"\
"	ulong p2 = p1;\n"\
"	p2 *= p1 + 1;\n"\
"\n"\
"	// If it's too small, N can't be amicable\n"\
"	if (N < targetSum - ((q + p2) << 1))\n"\
"	{\n"\
"		return;\n"\
"	}\n"\
"\n"\
"	// Data needed to resume from arbitrary loop iteration in phase 3:\n"\
"	// 1) M (8 bytes)\n"\
"	// 2) N (8 bytes)\n"\
"	// 3) targetSum (8 bytes)\n"\
"	// 4) i (<= 192724, can fit in 18 bits)\n"\
"	// It can all fit in 32 bytes\n"\
"\n"\
"	const uint filtered_numbers_data_id = atom_inc(filtered_numbers_count);\n"\
"	filtered_numbers_data[filtered_numbers_data_id].x = M;\n"\
"	filtered_numbers_data[filtered_numbers_data_id].y = N;\n"\
"	filtered_numbers_data[filtered_numbers_data_id].z = targetSum;\n"\
"	filtered_numbers_data[filtered_numbers_data_id].w = i;\n"\
"}\n"\
"\n"\
"__kernel\n"\
"__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))\n"\
"void CheckPairPhase3(\n"\
"	__global const ulong4* numbersToCheck,\n"\
"	const uint numbersToCheckOffset,\n"\
"	__global const ulong2* primeInverses,\n"\
"	__global const ulong2* primeReciprocals,\n"\
"	__global uint* filtered_numbers_count,\n"\
"	__global ulong4* filtered_numbers_data,\n"\
"	__global uint* amicable_numbers_count,\n"\
"	__global ulong2* amicable_numbers_data,\n"\
"	const int max_i\n"\
")\n"\
"{\n"\
"	const ulong4 cur_number = numbersToCheck[get_global_id(0) + numbersToCheckOffset];\n"\
"	if (cur_number.x == 0)\n"\
"	{\n"\
"		return;\n"\
"	}\n"\
"\n"\
"	const ulong M = cur_number.x;\n"\
"	ulong N = cur_number.y;\n"\
"	ulong targetSum = cur_number.z;\n"\
"	int i = cur_number.w;\n"\
"\n"\
"	// No more than 3 factors remain here\n"\
"	const ulong root3_N_rounded_up = convert_uint_rtp(cbrt(convert_float_rtp(N)));\n"\
"	const ulong MaximumSumOfDivisors3_target = targetSum - N - 1;\n"\
"\n"\
"	// All threads start with the same value of \"i\" here except for \"Phase 3, iteration 1\"\n"\
"	// so there will be no divergence because of different values of \"i\" on iteration 2, 3, etc.\n"\
"	// \"Phase 3, iteration 1\" is deliberately made shorter because of this\n"\
"	for (;;)\n"\
"	{\n"\
"		uint pattern =\n"\
"			((N * primeInverses[i +  0].x <= primeInverses[i +  0].y) ? (1U << 31) : 0) |\n"\
"			((N * primeInverses[i +  1].x <= primeInverses[i +  1].y) ? (1U << 30) : 0) |\n"\
"			((N * primeInverses[i +  2].x <= primeInverses[i +  2].y) ? (1U << 29) : 0) |\n"\
"			((N * primeInverses[i +  3].x <= primeInverses[i +  3].y) ? (1U << 28) : 0) |\n"\
"			((N * primeInverses[i +  4].x <= primeInverses[i +  4].y) ? (1U << 27) : 0) |\n"\
"			((N * primeInverses[i +  5].x <= primeInverses[i +  5].y) ? (1U << 26) : 0) |\n"\
"			((N * primeInverses[i +  6].x <= primeInverses[i +  6].y) ? (1U << 25) : 0) |\n"\
"			((N * primeInverses[i +  7].x <= primeInverses[i +  7].y) ? (1U << 24) : 0) |\n"\
"			((N * primeInverses[i +  8].x <= primeInverses[i +  8].y) ? (1U << 23) : 0) |\n"\
"			((N * primeInverses[i +  9].x <= primeInverses[i +  9].y) ? (1U << 22) : 0) |\n"\
"			((N * primeInverses[i + 10].x <= primeInverses[i + 10].y) ? (1U << 21) : 0) |\n"\
"			((N * primeInverses[i + 11].x <= primeInverses[i + 11].y) ? (1U << 20) : 0) |\n"\
"			((N * primeInverses[i + 12].x <= primeInverses[i + 12].y) ? (1U << 19) : 0) |\n"\
"			((N * primeInverses[i + 13].x <= primeInverses[i + 13].y) ? (1U << 18) : 0) |\n"\
"			((N * primeInverses[i + 14].x <= primeInverses[i + 14].y) ? (1U << 17) : 0) |\n"\
"			((N * primeInverses[i + 15].x <= primeInverses[i + 15].y) ? (1U << 16) : 0) |\n"\
"			((N * primeInverses[i + 16].x <= primeInverses[i + 16].y) ? (1U << 15) : 0) |\n"\
"			((N * primeInverses[i + 17].x <= primeInverses[i + 17].y) ? (1U << 14) : 0) |\n"\
"			((N * primeInverses[i + 18].x <= primeInverses[i + 18].y) ? (1U << 13) : 0) |\n"\
"			((N * primeInverses[i + 19].x <= primeInverses[i + 19].y) ? (1U << 12) : 0) |\n"\
"			((N * primeInverses[i + 20].x <= primeInverses[i + 20].y) ? (1U << 11) : 0) |\n"\
"			((N * primeInverses[i + 21].x <= primeInverses[i + 21].y) ? (1U << 10) : 0) |\n"\
"			((N * primeInverses[i + 22].x <= primeInverses[i + 22].y) ? (1U <<  9) : 0) |\n"\
"			((N * primeInverses[i + 23].x <= primeInverses[i + 23].y) ? (1U <<  8) : 0) |\n"\
"			((N * primeInverses[i + 24].x <= primeInverses[i + 24].y) ? (1U <<  7) : 0) |\n"\
"			((N * primeInverses[i + 25].x <= primeInverses[i + 25].y) ? (1U <<  6) : 0) |\n"\
"			((N * primeInverses[i + 26].x <= primeInverses[i + 26].y) ? (1U <<  5) : 0) |\n"\
"			((N * primeInverses[i + 27].x <= primeInverses[i + 27].y) ? (1U <<  4) : 0) |\n"\
"			((N * primeInverses[i + 28].x <= primeInverses[i + 28].y) ? (1U <<  3) : 0) |\n"\
"			((N * primeInverses[i + 29].x <= primeInverses[i + 29].y) ? (1U <<  2) : 0) |\n"\
"			((N * primeInverses[i + 30].x <= primeInverses[i + 30].y) ? (1U <<  1) : 0) |\n"\
"			((N * primeInverses[i + 31].x <= primeInverses[i + 31].y) ? (1U <<  0) : 0);\n"\
"\n"\
"		i += 32;\n"\
"\n"\
"		// Do trial divisions by prime factors found\n"\
"		if (pattern)\n"\
"		{\n"\
"			ulong sumN = 1;\n"\
"			do\n"\
"			{\n"\
"				// Bit 0: divisibility for i - 1\n"\
"				// Bit 1: divisibility for i - 2\n"\
"				// Bit 2: divisibility for i - 3\n"\
"				// and so on ...\n"\
"				const uint num_leading_zeroes = clz(pattern);\n"\
"				pattern -= (1 << (31 - num_leading_zeroes));\n"\
"				const uint index = i + num_leading_zeroes - 32;\n"\
"\n"\
"				const ulong prevSumN = sumN;\n"\
"				const uint p1 = primeReciprocals[index].y >> 32;\n"\
"				ulong q = N * primeInverses[index].x;\n"\
"				do\n"\
"				{\n"\
"					sumN = sumN * p1 + prevSumN;\n"\
"					N = q;\n"\
"					q *= primeInverses[index].x;\n"\
"				} while (q <= primeInverses[index].y);\n"\
"			} while (pattern);\n"\
"\n"\
"			const ulong p2 = primeReciprocals[i].y >> 32;\n"\
"			if (p2 * p2 > N)\n"\
"			{\n"\
"				if (N > 1)\n"\
"				{\n"\
"					sumN *= N + 1;\n"\
"				}\n"\
"				if (sumN == targetSum)\n"\
"				{\n"\
"					NumberFound(amicable_numbers_count, amicable_numbers_data, M, 0);\n"\
"				}\n"\
"				return;\n"\
"			}\n"\
"\n"\
"			const ulong oldTargetSum = targetSum;\n"\
"			targetSum /= sumN;\n"\
"			if (targetSum * sumN != oldTargetSum)\n"\
"			{\n"\
"				return;\n"\
"			}\n"\
"\n"\
"			if (N >= targetSum)\n"\
"			{\n"\
"				return;\n"\
"			}\n"\
"\n"\
"			// No more than 2 factors remain, we need to exit the loop now, otherwise MaximumSumOfDivisors3 will fail\n"\
"			break;\n"\
"		}\n"\
"\n"\
"		const ulong2 curReciprocal = primeReciprocals[i];\n"\
"		const ulong p1 = curReciprocal.y >> 32;\n"\
"		if (p1 > root3_N_rounded_up)\n"\
"		{\n"\
"			break;\n"\
"		}\n"\
"		const ulong q = mul_hi(N + ((curReciprocal.y >> 8) & 1), curReciprocal.x) >> (curReciprocal.y & 63);\n"\
"\n"\
"		ulong p2 = p1;\n"\
"		p2 *= p1 + 1;\n"\
"		if (((q + p2) << 1) <= MaximumSumOfDivisors3_target)\n"\
"		{\n"\
"			return;\n"\
"		}\n"\
"\n"\
"		if (i >= max_i)\n"\
"		{\n"\
"			const uint filtered_numbers_data_id = atom_inc(filtered_numbers_count);\n"\
"			filtered_numbers_data[filtered_numbers_data_id].x = M;\n"\
"			filtered_numbers_data[filtered_numbers_data_id].y = N;\n"\
"			filtered_numbers_data[filtered_numbers_data_id].z = targetSum;\n"\
"			filtered_numbers_data[filtered_numbers_data_id].w = max_i;\n"\
"			return;\n"\
"		}\n"\
"	}\n"\
"\n"\
"	// No more than 2 factors remain here\n"\
"\n"\
"	// 1) The case when N is prime\n"\
"	const ulong minSumN = N + 1;\n"\
"	if (minSumN >= targetSum)\n"\
"	{\n"\
"		if (minSumN == targetSum)\n"\
"		{\n"\
"			// M is an amicable number if N is prime\n"\
"			// We don't check for primality here, let CPU do it\n"\
"			NumberFound(amicable_numbers_count, amicable_numbers_data, M, N);\n"\
"		}\n"\
"		return;\n"\
"	}\n"\
"\n"\
"	// 2) The case when N is squared prime\n"\
"	const ulong sqrt_N = IntegerSquareRoot(N);\n"\
"	if ((sqrt_N * sqrt_N == N) && (N + sqrt_N + 1 == targetSum))\n"\
"	{\n"\
"		NumberFound(amicable_numbers_count, amicable_numbers_data, M, sqrt_N);\n"\
"		return;\n"\
"	}\n"\
"\n"\
"	// 3) The case when N = p * q - a product of two different primes\n"\
"	// p + q == targetSum - N - 1\n"\
"	// p * q == N\n"\
"	const ulong B = targetSum - N - 1;\n"\
"	ulong2 B2;\n"\
"	B2.x = B * B;\n"\
"	B2.y = mul_hi(B, B);\n"\
"	const ulong carry = (B2.x < (N << 2)) ? 1 : 0;\n"\
"	B2.x -= (N << 2);\n"\
"	if (B2.y >= (N >> 62) + carry)\n"\
"	{\n"\
"		B2.y -= (N >> 62) + carry;\n"\
"		const ulong sqrt_D = IntegerSquareRoot128(B2);\n"\
"		const ulong p = (B - sqrt_D) >> 1;\n"\
"		const ulong q = (B + sqrt_D) >> 1;\n"\
"		if ((p < q) && (p + q == B) && (p * q == N))\n"\
"		{\n"\
"			NumberFound(amicable_numbers_count, amicable_numbers_data, M, q);\n"\
"			return;\n"\
"		}\n"\
"	}\n"\
"}\n"\
"\n"\
"__kernel\n"\
"__attribute__((reqd_work_group_size(1, 1, 1)))\n"\
"void SaveCounter(__global uint* phase2_numbers_count)\n"\
"{\n"\
"	if (phase2_numbers_count[0] < PHASE2_MAX_COUNT)\n"\
"	{\n"\
"		phase2_numbers_count[1] = phase2_numbers_count[0];\n"\
"	}\n"\
"}\n"\
"\n"\
"// primes is an array of pointers to chunks of data\n"\
"// each chunk is 128 MB in size (2^24 uint2 elements) or bigger\n"\
"ulong GetNthPrime(uint n,\n"\
"	__global uint2* primes0\n"\
"#if NUM_PRIME_CHUNKS > 1\n"\
"	, __global uint2* primes1\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 2\n"\
"	, __global uint2* primes2\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 3\n"\
"	, __global uint2* primes3\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 4\n"\
"	, __global uint2* primes4\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 5\n"\
"	, __global uint2* primes5\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 6\n"\
"	, __global uint2* primes6\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 7\n"\
"	, __global uint2* primes7\n"\
"#endif\n"\
")\n"\
"{\n"\
"#if NUM_PRIME_CHUNKS == 1\n"\
"	uint2 data = primes0[n >> 2];\n"\
"#else\n"\
"	const uint global_offset = n >> 2;\n"\
"	__global const uint2* chunk;\n"\
"	switch (global_offset >> CHUNK_SIZE_SHIFT)\n"\
"	{\n"\
"	default: chunk = primes0; break;\n"\
"#if NUM_PRIME_CHUNKS > 1\n"\
"	case 1: chunk = primes1; break;\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 2\n"\
"	case 2: chunk = primes2; break;\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 3\n"\
"	case 3: chunk = primes3; break;\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 4\n"\
"	case 4: chunk = primes4; break;\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 5\n"\
"	case 5: chunk = primes5; break;\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 6\n"\
"	case 6: chunk = primes6; break;\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 7\n"\
"	case 7: chunk = primes7; break;\n"\
"#endif\n"\
"	}\n"\
"\n"\
"	const uint chunk_offset = global_offset & ((1 << CHUNK_SIZE_SHIFT) - 1);\n"\
"	uint2 data = chunk[chunk_offset];\n"\
"#endif\n"\
"\n"\
"	data.x += (data.y >> (30 - (n & 3) * 10)) & 1023;\n"\
"	const ulong result = data.x;\n"\
"	return result + result + 1;\n"\
"}\n"\
"\n"\
"__kernel\n"\
"__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))\n"\
"void SearchMultipleRanges(\n"\
"	__global const uint* smallPrimes,\n"\
"	__global const ulong2* primeInverses,\n"\
"	__global const ulong2* PQ,\n"\
"	__constant const ulong2* PowerOf2SumInverses,\n"\
"	__constant const ulong4* PowersOfP_128SumInverses,\n"\
"	__constant const ulong2* PrimeInverses_128,\n"\
"	__global const ulong4* RangesTable,\n"\
"	__global const uint2* RangeLookupTable,\n"\
"	const uint lookup_shift,\n"\
"	const uint global_offset,\n"\
"	const uint global_size,\n"\
"	__global uint* phase1_offset_to_resume_after_overflow,\n"\
"	__global uint* phase2_numbers_count,\n"\
"	__global ulong4* phase2_numbers_data,\n"\
"	__global uint* amicable_numbers_count,\n"\
"	__global ulong2* amicable_numbers_data,\n"\
"\n"\
"	__global uint2* primes0\n"\
"#if NUM_PRIME_CHUNKS > 1\n"\
"	, __global uint2* primes1\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 2\n"\
"	, __global uint2* primes2\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 3\n"\
"	, __global uint2* primes3\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 4\n"\
"	, __global uint2* primes4\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 5\n"\
"	, __global uint2* primes5\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 6\n"\
"	, __global uint2* primes6\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 7\n"\
"	, __global uint2* primes7\n"\
"#endif\n"\
")\n"\
"{\n"\
"	if (*phase2_numbers_count >= PHASE2_MAX_COUNT)\n"\
"	{\n"\
"		return;\n"\
"	}\n"\
"\n"\
"	const uint globalIndex = get_global_id(0) + global_offset;\n"\
"	if (globalIndex >= global_size)\n"\
"	{\n"\
"		return;\n"\
"	}\n"\
"\n"\
"	ulong M0, sumM0, p;\n"\
"\n"\
"	uint2 rangeIndexAndOffset = RangeLookupTable[globalIndex >> lookup_shift];\n"\
"	uint curIndex = globalIndex & (0xFFFFFFFFU << lookup_shift);\n"\
"	for (;;)\n"\
"	{\n"\
"		const uint numbersRemainingInCurRange = RangesTable[rangeIndexAndOffset.x].w - rangeIndexAndOffset.y;\n"\
"\n"\
"		// curIndex + numbersRemainingInCurRange now points exactly at the beginning of the next range\n"\
"		// If the beginning of the next range is > globalIndex, then globalIndex is within the current range\n"\
"		if (curIndex + numbersRemainingInCurRange > globalIndex)\n"\
"		{\n"\
"			M0 = RangesTable[rangeIndexAndOffset.x].x;\n"\
"			sumM0 = RangesTable[rangeIndexAndOffset.x].y;\n"\
"\n"\
"			//      primes[ start_prime_index                    + range_offset          + globalIndex - curIndex]\n"\
"			p = GetNthPrime(RangesTable[rangeIndexAndOffset.x].z + rangeIndexAndOffset.y + globalIndex - curIndex,\n"\
"				primes0\n"\
"#if NUM_PRIME_CHUNKS > 1\n"\
"				, primes1\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 2\n"\
"				, primes2\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 3\n"\
"				, primes3\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 4\n"\
"				, primes4\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 5\n"\
"				, primes5\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 6\n"\
"				, primes6\n"\
"#endif\n"\
"#if NUM_PRIME_CHUNKS > 7\n"\
"				, primes7\n"\
"#endif\n"\
"			);\n"\
"			break;\n"\
"		}\n"\
"\n"\
"		// Move curIndex to the beginning of the next range\n"\
"		curIndex += numbersRemainingInCurRange;\n"\
"		++rangeIndexAndOffset.x;\n"\
"		rangeIndexAndOffset.y = 0;\n"\
"	}\n"\
"\n"\
"	CheckPairPhase1(\n"\
"		smallPrimes,\n"\
"		primeInverses,\n"\
"		PQ,\n"\
"		PowerOf2SumInverses,\n"\
"		PowersOfP_128SumInverses,\n"\
"		PrimeInverses_128,\n"\
"#if LPP == 1\n"\
"		M0 * p, sumM0 * (p + 1), mul_hi(sumM0, p + 1),\n"\
"#elif LPP == 2\n"\
"		M0 * p * p, sumM0 * (p * p + p + 1), mul_hi(sumM0, p * p + p + 1),\n"\
"#elif LPP == 3\n"\
"		M0 * p * p * p, sumM0 * (p * p * p + p * p + p + 1), mul_hi(sumM0, p * p + p + 1),\n"\
"#endif\n"\
"		global_offset,\n"\
"		phase1_offset_to_resume_after_overflow,\n"\
"		phase2_numbers_count,\n"\
"		phase2_numbers_data,\n"\
"		amicable_numbers_count,\n"\
"		amicable_numbers_data\n"\
"	);\n"\
"}\n"\
"\n"\
"__kernel\n"\
"__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))\n"\
"void CheckPairs(\n"\
"	__global const uint* smallPrimes,\n"\
"	__global const ulong2* primeInverses,\n"\
"	__global const ulong2* PQ,\n"\
"	__constant const ulong2* PowerOf2SumInverses,\n"\
"	__constant const ulong4* PowersOfP_128SumInverses,\n"\
"	__constant const ulong2* PrimeInverses_128,\n"\
"	__global const ulong4* pairsToCheck,\n"\
"	__global uint* filtered_numbers_count,\n"\
"	__global ulong4* filtered_numbers_data,\n"\
"	__global uint* amicable_numbers_count,\n"\
"	__global ulong2* amicable_numbers_data\n"\
")\n"\
"{\n"\
"	const ulong4 cur_pair = pairsToCheck[get_global_id(0)];\n"\
"	if (cur_pair.x == 0)\n"\
"	{\n"\
"		return;\n"\
"	}\n"\
"	CheckPairPhase1(smallPrimes, primeInverses, PQ, PowerOf2SumInverses, PowersOfP_128SumInverses, PrimeInverses_128, cur_pair.x, cur_pair.y, cur_pair.z, 0, 0, filtered_numbers_count, filtered_numbers_data, amicable_numbers_count, amicable_numbers_data);\n"\
"}\n"\
"";

static const unsigned int kernel_cl_crc32 = 0xc412ade0;
