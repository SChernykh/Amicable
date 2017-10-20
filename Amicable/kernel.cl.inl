#pragma once

static const char* kernel_cl = "#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable\n"\
"\n"\
"#ifdef cl_clang_storage_class_specifiers\n"\
"#pragma OPENCL EXTENSION cl_clang_storage_class_specifiers : enable\n"\
"#endif\n"\
"\n"\
"#define PHASE1_DEPTH 256\n"\
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
"static void NumberFound(__global uint* amicable_numbers_count, __global ulong4* amicable_numbers_data, ulong A, ulong A_high, ulong B)\n"\
"{\n"\
"	const uint index = atom_inc(amicable_numbers_count);\n"\
"	amicable_numbers_data[index].x = A;\n"\
"	amicable_numbers_data[index].y = A_high;\n"\
"	amicable_numbers_data[index].z = B;\n"\
"}\n"\
"\n"\
"// 32-bit floats don't have enough precision\n"\
"// 64-bit floats are not supported by all OpenCL devices\n"\
"// So we have to use integer arithmetic to calculate square roots\n"\
"\n"\
"// Returns number x such that x^2 <= n < (x+1)^2\n"\
"static uint IntegerSquareRoot(const ulong n)\n"\
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
"static ulong IntegerSquareRoot128(const ulong2 n)\n"\
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
"// Returns number x such that x^3 <= n < (x+1)^3\n"\
"static ulong IntegerCubeRoot128(const ulong2 n)\n"\
"{\n"\
"	ulong result = 1;\n"\
"	const uint highest_bit_index = n.y ? (127 - clz(n.y)) : (63 - clz(n.x));\n"\
"	result <<= (highest_bit_index * ((65536 / 3) + 1)) >> 16;\n"\
"	for (ulong cur_bit = result >> 1; cur_bit > 0; cur_bit >>= 1)\n"\
"	{\n"\
"		const ulong k = result | cur_bit;\n"\
"		ulong2 k2;\n"\
"		k2.x = (k * k) * k;\n"\
"		k2.y = mul_hi(k * k, k);\n"\
"		if ((k2.y < n.y) || ((k2.y == n.y) && (k2.x <= n.x)))\n"\
"		{\n"\
"			result = k;\n"\
"		}\n"\
"	}\n"\
"	return result;\n"\
"}\n"\
"\n"\
"static void CheckPairPhase1(\n"\
"	__global const uint* smallPrimes,\n"\
"	__global const ulong2* primeInverses,\n"\
"	__global const ulong2* PQ,\n"\
"	__constant ulong2* PowerOf2SumInverses,\n"\
"	__global const uint* PowersOfP_128SumInverses_offsets,\n"\
"	__global const ulong4* PowersOfP_128SumInverses,\n"\
"	__global const ulong2* PrimeInverses_128,\n"\
"	__global const ulong* SumEstimates_128,\n"\
"	const ulong M,\n"\
"	const ulong M_high,\n"\
"	ulong targetSum,\n"\
"	ulong targetSumHigh,\n"\
"	const ulong global_offset,\n"\
"	__global ulong* phase1_offset_to_resume_after_overflow,\n"\
"	__global uint* phase2_numbers_count,\n"\
"	__global ulong4* phase2_numbers_data,\n"\
"	__global uint* amicable_numbers_count,\n"\
"	__global ulong4* amicable_numbers_data\n"\
")\n"\
"{\n"\
"	ulong N = targetSum - M;\n"\
"	ulong N_high = targetSumHigh - M_high - ((targetSum < M) ? 1 : 0);\n"\
"\n"\
"	// Divide out power of 2 using fast bitwise operations\n"\
"	// All numbers in the same range have the same parity, so threads won't diverge here\n"\
"	if ((N & 1) == 0)\n"\
"	{\n"\
"		if (N == 0)\n"\
"		{\n"\
"			// N = 2^64*k, k >= 0. It's not amicable for k = 0, 1, 2, ..., 5 and too large for k > 5.\n"\
"			return;\n"\
"		}\n"\
"\n"\
"		const uint powerOf2 = 63 - clz(N ^ (N - 1));\n"\
"\n"\
"		// M and targetSum (128-bit) come in ranges where they both grow monotonically,\n"\
"		// so there will be exactly one divergence in each range and one divergence on a border between ranges\n"\
"		// Threads don't diverge 99.9% of the time here\n"\
"		if (targetSumHigh)\n"\
"		{\n"\
"			const ulong2 inv = PowerOf2SumInverses[powerOf2 + 64];\n"\
"			const ulong q_high = mul_hi(targetSum, inv.x) + targetSum * inv.y + targetSumHigh * inv.x;\n"\
"\n"\
"			// The proper check would be \"if (128-bit value of q > (2^128 - 1) / divisor)\"\n"\
"			// \"if (q_high > targetSumHigh)\" doesn't give false positives because the result of division (by divisor > 1) can't be larger than the original value\n"\
"			// It can give very rare false negatives when the divisor is very large (larger than 2^64 / targetSumHigh), but false negatives don't jeopardize the correctness of the search here\n"\
"			// Checking only highest 64 bits and not having to store or calculate \"(2^128 - 1) / divisor\" is much faster on average, so it's worth it\n"\
"			if (q_high > targetSumHigh)\n"\
"			{\n"\
"				return;\n"\
"			}\n"\
"\n"\
"			targetSum *= inv.x;\n"\
"			targetSumHigh = q_high;\n"\
"\n"\
"			N = (N >> powerOf2) | (N_high << (64 - powerOf2));\n"\
"			N_high >>= powerOf2;\n"\
"\n"\
"			if ((N_high > targetSumHigh) || ((N_high == targetSumHigh) && (N >= targetSum)))\n"\
"			{\n"\
"				return;\n"\
"			}\n"\
"		}\n"\
"		else\n"\
"		{\n"\
"			targetSum *= PowerOf2SumInverses[powerOf2].x;\n"\
"			if (targetSum > PowerOf2SumInverses[powerOf2].y)\n"\
"			{\n"\
"				return;\n"\
"			}\n"\
"			N >>= powerOf2;\n"\
"\n"\
"			if (N >= targetSum)\n"\
"			{\n"\
"				return;\n"\
"			}\n"\
"		}\n"\
"	}\n"\
"\n"\
"	int i = 0;\n"\
"\n"\
"	// M and targetSum (128-bit) come in ranges where they both grow monotonically,\n"\
"	// so there will be exactly one divergence in each range and one divergence on a border between ranges\n"\
"	// Threads don't diverge 99.9% of the time here\n"\
"	if (targetSumHigh)\n"\
"	{\n"\
"		for (int cur_index = 0; cur_index < PHASE1_DEPTH; cur_index += 16)\n"\
"		{\n"\
"			const int e = cur_index + 16;\n"\
"			for (i = cur_index + 1; i <= e; ++i)\n"\
"			{\n"\
"				const ulong2 inv = PrimeInverses_128[i - 1];\n"\
"				ulong q_high = mul_hi(N, inv.x) + N * inv.y + N_high * inv.x;\n"\
"				if (q_high > N_high)\n"\
"				{\n"\
"					continue;\n"\
"				}\n"\
"\n"\
"				ulong q = N * inv.x;\n"\
"				int powerOfP = 0;\n"\
"				do\n"\
"				{\n"\
"					++powerOfP;\n"\
"\n"\
"					N = q;\n"\
"					N_high = q_high;\n"\
"\n"\
"					q = N * inv.x;\n"\
"					q_high = mul_hi(N, inv.x) + N * inv.y + N_high * inv.x;\n"\
"				} while (q_high <= N_high);\n"\
"\n"\
"				const ulong4 inv2 = PowersOfP_128SumInverses[PowersOfP_128SumInverses_offsets[i - 1] + powerOfP - 1];\n"\
"				if (inv2.z != 0)\n"\
"				{\n"\
"					if ((targetSum & inv2.w) != 0)\n"\
"					{\n"\
"						return;\n"\
"					}\n"\
"					targetSum = (targetSum >> inv2.z) | (targetSumHigh << (64 - inv2.z));\n"\
"					targetSumHigh >>= inv2.z;\n"\
"				}\n"\
"\n"\
"				q_high = mul_hi(targetSum, inv2.x) + targetSum * inv2.y + targetSumHigh * inv2.x;\n"\
"				if (q_high > targetSumHigh)\n"\
"				{\n"\
"					return;\n"\
"				}\n"\
"\n"\
"				targetSum *= inv2.x;\n"\
"				targetSumHigh = q_high;\n"\
"\n"\
"				const ulong p = smallPrimes[i];\n"\
"				if (p * p > N)\n"\
"				{\n"\
"					const ulong sumN = ((N > 1) ? N : 0) + 1;\n"\
"					if (sumN == targetSum)\n"\
"					{\n"\
"						NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, 0);\n"\
"					}\n"\
"					return;\n"\
"				}\n"\
"\n"\
"				if ((N_high > targetSumHigh) || ((N_high == targetSumHigh) && (N >= targetSum)))\n"\
"				{\n"\
"					return;\n"\
"				}\n"\
"\n"\
"				if (!targetSumHigh)\n"\
"				{\n"\
"					const int i1 = i + 1;\n"\
"					i = i1 - (i1 & 15);\n"\
"#ifdef DISABLE_GOTO\n"\
"					break;\n"\
"#else\n"\
"					goto small_numbers;\n"\
"#endif\n"\
"				}\n"\
"			}\n"\
"\n"\
"#ifdef DISABLE_GOTO\n"\
"			if (!targetSumHigh)\n"\
"			{\n"\
"				break;\n"\
"			}\n"\
"#endif\n"\
"\n"\
"			ulong2 max_sum;\n"\
"			const ulong max_sum_ratio = SumEstimates_128[cur_index >> 4];\n"\
"\n"\
"#if SUM_ESTIMATES_128_SHIFT == 64\n"\
"			max_sum.x = mul_hi(N, max_sum_ratio) + N_high * max_sum_ratio;\n"\
"			max_sum.y = 0;\n"\
"#else\n"\
"			max_sum.x = N * max_sum_ratio;\n"\
"			max_sum.y = mul_hi(N, max_sum_ratio) + N_high * max_sum_ratio;\n"\
"			max_sum.x = (max_sum.x >> SUM_ESTIMATES_128_SHIFT) | (max_sum.y << (64 - SUM_ESTIMATES_128_SHIFT));\n"\
"			max_sum.y >>= SUM_ESTIMATES_128_SHIFT;\n"\
"#endif\n"\
"\n"\
"			const ulong carry = (max_sum.x + N < max_sum.x) ? 1 : 0;\n"\
"			max_sum.x = max_sum.x + N;\n"\
"			max_sum.y = max_sum.y + N_high + carry;\n"\
"\n"\
"			if ((max_sum.y < targetSumHigh) || ((max_sum.y == targetSumHigh) && (max_sum.x < targetSum)))\n"\
"			{\n"\
"				return;\n"\
"			}\n"\
"		}\n"\
"\n"\
"#ifdef DISABLE_GOTO\n"\
"		if (targetSumHigh)\n"\
"#endif\n"\
"		{\n"\
"			const uint phase2_number_index = atom_inc(phase2_numbers_count);\n"\
"			if (phase2_number_index < PHASE2_MAX_COUNT)\n"\
"			{\n"\
"				const ulong i_to_save = i - (i & 15);\n"\
"				phase2_numbers_data[phase2_number_index].x = M;\n"\
"				phase2_numbers_data[phase2_number_index].y = N;\n"\
"				phase2_numbers_data[phase2_number_index].z = targetSum;\n"\
"				phase2_numbers_data[phase2_number_index].w = i_to_save | (M_high << 22) | (N_high << 36) | (targetSumHigh << 50);\n"\
"			}\n"\
"			else if (phase2_number_index == PHASE2_MAX_COUNT)\n"\
"			{\n"\
"				*phase1_offset_to_resume_after_overflow = global_offset;\n"\
"			}\n"\
"			return;\n"\
"		}\n"\
"	}\n"\
"\n"\
"#ifndef DISABLE_GOTO\n"\
"	small_numbers:\n"\
"#endif\n"\
"\n"\
"	if (i < 16)\n"\
"	{\n"\
"		// Collect divisibility data for first 15 odd primes\n"\
"		// No divergence here, compiler will use predicates for \"A ? B : C\" operator\n"\
"		uint pattern =\n"\
"			((N * primeInverses[ 1].x <= primeInverses[ 1].y) ? 16384 : 0) |\n"\
"			((N * primeInverses[ 2].x <= primeInverses[ 2].y) ? 8192  : 0) |\n"\
"			((N * primeInverses[ 3].x <= primeInverses[ 3].y) ? 4096  : 0) |\n"\
"			((N * primeInverses[ 4].x <= primeInverses[ 4].y) ? 2048  : 0) |\n"\
"			((N * primeInverses[ 5].x <= primeInverses[ 5].y) ? 1024  : 0) |\n"\
"			((N * primeInverses[ 6].x <= primeInverses[ 6].y) ? 512   : 0) |\n"\
"			((N * primeInverses[ 7].x <= primeInverses[ 7].y) ? 256   : 0) |\n"\
"			((N * primeInverses[ 8].x <= primeInverses[ 8].y) ? 128   : 0) |\n"\
"			((N * primeInverses[ 9].x <= primeInverses[ 9].y) ? 64    : 0) |\n"\
"			((N * primeInverses[10].x <= primeInverses[10].y) ? 32    : 0) |\n"\
"			((N * primeInverses[11].x <= primeInverses[11].y) ? 16    : 0) |\n"\
"			((N * primeInverses[12].x <= primeInverses[12].y) ? 8     : 0) |\n"\
"			((N * primeInverses[13].x <= primeInverses[13].y) ? 4     : 0) |\n"\
"			((N * primeInverses[14].x <= primeInverses[14].y) ? 2     : 0) |\n"\
"			((N * primeInverses[15].x <= primeInverses[15].y) ? 1     : 0);\n"\
"\n"\
"		i = 16;\n"\
"\n"\
"		// Do trial divisions by primes found\n"\
"		if (pattern)\n"\
"		{\n"\
"			ulong sumN = 1;\n"\
"			do\n"\
"			{\n"\
"				// Bit 0: divisibility for i == 15\n"\
"				// Bit 1: divisibility for i == 14\n"\
"				// Bit 2: divisibility for i == 13\n"\
"				// and so on ...\n"\
"				const uint num_leading_zeroes = clz(pattern);\n"\
"				pattern -= (1 << (31 - num_leading_zeroes));\n"\
"				const uint index = num_leading_zeroes - 16;\n"\
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
"			if (N < 59 * 59)\n"\
"			{\n"\
"				if (N > 1)\n"\
"				{\n"\
"					sumN *= N + 1;\n"\
"				}\n"\
"				if (sumN == targetSum)\n"\
"				{\n"\
"					NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, 0);\n"\
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
"	if (i < PHASE1_DEPTH)\n"\
"	{\n"\
"		do\n"\
"		{\n"\
"			// Collect divisibility data for the next 16 primes\n"\
"			// No divergence here, compiler will use predicates for \"A ? B : C\" operator\n"\
"			uint pattern =\n"\
"				((N * primeInverses[i +  0].x <= primeInverses[i +  0].y) ? 32768 : 0) |\n"\
"				((N * primeInverses[i +  1].x <= primeInverses[i +  1].y) ? 16384 : 0) |\n"\
"				((N * primeInverses[i +  2].x <= primeInverses[i +  2].y) ? 8192  : 0) |\n"\
"				((N * primeInverses[i +  3].x <= primeInverses[i +  3].y) ? 4096  : 0) |\n"\
"				((N * primeInverses[i +  4].x <= primeInverses[i +  4].y) ? 2048  : 0) |\n"\
"				((N * primeInverses[i +  5].x <= primeInverses[i +  5].y) ? 1024  : 0) |\n"\
"				((N * primeInverses[i +  6].x <= primeInverses[i +  6].y) ? 512   : 0) |\n"\
"				((N * primeInverses[i +  7].x <= primeInverses[i +  7].y) ? 256   : 0) |\n"\
"				((N * primeInverses[i +  8].x <= primeInverses[i +  8].y) ? 128   : 0) |\n"\
"				((N * primeInverses[i +  9].x <= primeInverses[i +  9].y) ? 64    : 0) |\n"\
"				((N * primeInverses[i + 10].x <= primeInverses[i + 10].y) ? 32    : 0) |\n"\
"				((N * primeInverses[i + 11].x <= primeInverses[i + 11].y) ? 16    : 0) |\n"\
"				((N * primeInverses[i + 12].x <= primeInverses[i + 12].y) ? 8     : 0) |\n"\
"				((N * primeInverses[i + 13].x <= primeInverses[i + 13].y) ? 4     : 0) |\n"\
"				((N * primeInverses[i + 14].x <= primeInverses[i + 14].y) ? 2     : 0) |\n"\
"				((N * primeInverses[i + 15].x <= primeInverses[i + 15].y) ? 1     : 0);\n"\
"\n"\
"			i += 16;\n"\
"\n"\
"			// Do trial divisions by prime factors found\n"\
"			if (pattern)\n"\
"			{\n"\
"				ulong sumN = 1;\n"\
"				do\n"\
"				{\n"\
"					// Bit 0: divisibility for i - 1\n"\
"					// Bit 1: divisibility for i - 2\n"\
"					// Bit 2: divisibility for i - 3\n"\
"					// and so on ...\n"\
"					const uint num_leading_zeroes = clz(pattern);\n"\
"					pattern -= (1 << (31 - num_leading_zeroes));\n"\
"					const uint index = i + num_leading_zeroes - 32;\n"\
"\n"\
"					const ulong prevSumN = sumN;\n"\
"					const uint p1 = smallPrimes[index];\n"\
"					ulong q = N * primeInverses[index].x;\n"\
"					do\n"\
"					{\n"\
"						sumN = sumN * p1 + prevSumN;\n"\
"						N = q;\n"\
"						q *= primeInverses[index].x;\n"\
"					} while (q <= primeInverses[index].y);\n"\
"				} while (pattern);\n"\
"\n"\
"				const ulong p = smallPrimes[i];\n"\
"				if (p * p > N)\n"\
"				{\n"\
"					if (N > 1)\n"\
"					{\n"\
"						sumN *= N + 1;\n"\
"					}\n"\
"					if (sumN == targetSum)\n"\
"					{\n"\
"						NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, 0);\n"\
"					}\n"\
"					return;\n"\
"				}\n"\
"\n"\
"				const ulong oldTargetSum = targetSum;\n"\
"				targetSum /= sumN;\n"\
"				if (targetSum * sumN != oldTargetSum)\n"\
"				{\n"\
"					return;\n"\
"				}\n"\
"\n"\
"				if (N >= targetSum)\n"\
"				{\n"\
"					return;\n"\
"				}\n"\
"			}\n"\
"\n"\
"			// Check maximal possible sum for N. If it's too small, it can't be an amicable number\n"\
"			// Same i for all threads, no divergence\n"\
"			uint index = k * PQ_STRIDE_SIZE + (i >> 4) + 16;\n"\
"			while (N <= PQ[index].x)\n"\
"			{\n"\
"				--k;\n"\
"				index -= PQ_STRIDE_SIZE;\n"\
"			}\n"\
"			if (mul_hi(N, PQ[index].y) < targetSum - N)\n"\
"			{\n"\
"				return;\n"\
"			}\n"\
"		} while (i != PHASE1_DEPTH);\n"\
"	}\n"\
"	else\n"\
"	{\n"\
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
"	}\n"\
"\n"\
"	// Data needed to resume from arbitrary loop iteration in phase 2:\n"\
"	// 1) M (8 bytes + 14 bits)\n"\
"	// 2) N (8 bytes + 14 bits)\n"\
"	// 3) targetSum (8 bytes + 14 bits)\n"\
"	// 4) i (<= 192724, can fit in 18 bits)\n"\
"	// 5) k (< 16, can fit in 4 bits)\n"\
"	// It can all fit in 32 bytes\n"\
"	// 8 bytes + 14 bits is enough to store numbers up to ~3.022 * 10^23, so 10^23 is the maximum search limit for this version\n"\
"\n"\
"	const uint phase2_number_index = atom_inc(phase2_numbers_count);\n"\
"	if (phase2_number_index < PHASE2_MAX_COUNT)\n"\
"	{\n"\
"		const ulong i_to_save = i;\n"\
"		const ulong k_to_save = k;\n"\
"		phase2_numbers_data[phase2_number_index].x = M;\n"\
"		phase2_numbers_data[phase2_number_index].y = N;\n"\
"		phase2_numbers_data[phase2_number_index].z = targetSum;\n"\
"		phase2_numbers_data[phase2_number_index].w = i_to_save | (k_to_save << 18) | (M_high << 22) | (N_high << 36) | (targetSumHigh << 50);\n"\
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
"	const uint numbersToCheckOffset,\n"\
"	__global const ulong2* primeInverses,\n"\
"	__global const uint* PowersOfP_128SumInverses_offsets,\n"\
"	__global const ulong4* PowersOfP_128SumInverses,\n"\
"	__global const ulong2* PrimeInverses_128,\n"\
"	__global const ulong* SumEstimates_128,\n"\
"	__global const ulong2* PQ,\n"\
"	__global uint* filtered_numbers_count,\n"\
"	__global ulong4* filtered_numbers_data,\n"\
"	__global uint* amicable_numbers_count,\n"\
"	__global ulong4* amicable_numbers_data\n"\
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
"	int i = cur_number.w & ((1 << 18) - 1);\n"\
"	int k = (cur_number.w >> 18) & ((1 << 4) - 1);\n"\
"\n"\
"	const ulong M_high = (cur_number.w >> 22) & ((1 << 14) - 1);\n"\
"	ulong N_high = (cur_number.w >> 36) & ((1 << 14) - 1);\n"\
"	ulong targetSumHigh = cur_number.w >> 50;\n"\
"\n"\
"	if (targetSumHigh)\n"\
"	{\n"\
"		ulong2 N_128;\n"\
"		N_128.x = N;\n"\
"		N_128.y = N_high;\n"\
"		uint root4_N_rounded_up = convert_uint_rtp(sqrt(convert_float_rtp(IntegerSquareRoot128(N_128))));\n"\
"\n"\
"		for (int cur_index = i; (smallPrimes[i] <= root4_N_rounded_up) && targetSumHigh; cur_index += 16)\n"\
"		{\n"\
"			const int e = cur_index + 16;\n"\
"			for (i = cur_index + 1; i <= e; ++i)\n"\
"			{\n"\
"				const ulong2 inv = PrimeInverses_128[i - 1];\n"\
"				ulong q_high = mul_hi(N, inv.x) + N * inv.y + N_high * inv.x;\n"\
"				if (q_high > N_high)\n"\
"				{\n"\
"					continue;\n"\
"				}\n"\
"\n"\
"				ulong q = N * inv.x;\n"\
"				int powerOfP = 0;\n"\
"				do\n"\
"				{\n"\
"					++powerOfP;\n"\
"\n"\
"					N = q;\n"\
"					N_high = q_high;\n"\
"\n"\
"					q = N * inv.x;\n"\
"					q_high = mul_hi(N, inv.x) + N * inv.y + N_high * inv.x;\n"\
"				} while (q_high <= N_high);\n"\
"\n"\
"				const ulong4 inv2 = PowersOfP_128SumInverses[PowersOfP_128SumInverses_offsets[i - 1] + powerOfP - 1];\n"\
"				if (inv2.z != 0)\n"\
"				{\n"\
"					if ((targetSum & inv2.w) != 0)\n"\
"					{\n"\
"						return;\n"\
"					}\n"\
"					targetSum = (targetSum >> inv2.z) | (targetSumHigh << (64 - inv2.z));\n"\
"					targetSumHigh >>= inv2.z;\n"\
"				}\n"\
"\n"\
"				q_high = mul_hi(targetSum, inv2.x) + targetSum * inv2.y + targetSumHigh * inv2.x;\n"\
"				if (q_high > targetSumHigh)\n"\
"				{\n"\
"					return;\n"\
"				}\n"\
"\n"\
"				targetSum *= inv2.x;\n"\
"				targetSumHigh = q_high;\n"\
"\n"\
"				const ulong p = smallPrimes[i];\n"\
"				if (p * p > N)\n"\
"				{\n"\
"					const ulong sumN = ((N > 1) ? N : 0) + 1;\n"\
"					if (sumN == targetSum)\n"\
"					{\n"\
"						NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, 0);\n"\
"					}\n"\
"					return;\n"\
"				}\n"\
"\n"\
"				if ((N_high > targetSumHigh) || ((N_high == targetSumHigh) && (N >= targetSum)))\n"\
"				{\n"\
"					return;\n"\
"				}\n"\
"\n"\
"				N_128.x = N;\n"\
"				N_128.y = N_high;\n"\
"				root4_N_rounded_up = convert_uint_rtp(sqrt(convert_float_rtp(IntegerSquareRoot128(N_128))));\n"\
"			}\n"\
"\n"\
"			ulong2 max_sum;\n"\
"			const ulong max_sum_ratio = SumEstimates_128[cur_index >> 4];\n"\
"\n"\
"#if SUM_ESTIMATES_128_SHIFT == 64\n"\
"			max_sum.x = mul_hi(N, max_sum_ratio) + N_high * max_sum_ratio;\n"\
"			max_sum.y = 0;\n"\
"#else\n"\
"			max_sum.x = N * max_sum_ratio;\n"\
"			max_sum.y = mul_hi(N, max_sum_ratio) + N_high * max_sum_ratio;\n"\
"			max_sum.x = (max_sum.x >> SUM_ESTIMATES_128_SHIFT) | (max_sum.y << (64 - SUM_ESTIMATES_128_SHIFT));\n"\
"			max_sum.y >>= SUM_ESTIMATES_128_SHIFT;\n"\
"#endif\n"\
"\n"\
"			const ulong carry = (max_sum.x + N < max_sum.x) ? 1 : 0;\n"\
"			max_sum.x = max_sum.x + N;\n"\
"			max_sum.y = max_sum.y + N_high + carry;\n"\
"\n"\
"			if ((max_sum.y < targetSumHigh) || ((max_sum.y == targetSumHigh) && (max_sum.x < targetSum)))\n"\
"			{\n"\
"				return;\n"\
"			}\n"\
"		}\n"\
"\n"\
"		if (targetSumHigh)\n"\
"		{\n"\
"			const uint filtered_numbers_data_id = atom_inc(filtered_numbers_count);\n"\
"			filtered_numbers_data[filtered_numbers_data_id].x = M;\n"\
"			filtered_numbers_data[filtered_numbers_data_id].y = N;\n"\
"			filtered_numbers_data[filtered_numbers_data_id].z = targetSum;\n"\
"			const ulong i_to_save = i;\n"\
"			filtered_numbers_data[filtered_numbers_data_id].w = i_to_save | (M_high << 22) | (N_high << 36) | (targetSumHigh << 50);\n"\
"			return;\n"\
"		}\n"\
"\n"\
"		i -= (i & 15);\n"\
"\n"\
"		k = 8;\n"\
"		while (N <= PQ[k].x)\n"\
"		{\n"\
"			--k;\n"\
"		}\n"\
"		++k;\n"\
"		uint index = k * PQ_STRIDE_SIZE + (i >> 4) + 16;\n"\
"		while (N <= PQ[index].x)\n"\
"		{\n"\
"			--k;\n"\
"			index -= PQ_STRIDE_SIZE;\n"\
"		}\n"\
"	}\n"\
"\n"\
"	uint root4_N_rounded_up = convert_uint_rtp(sqrt(sqrt(convert_float_rtp(N))));\n"\
"\n"\
"	// All threads start with the same value of \"i\" here, which is set to PHASE1_DEPTH at the end of phase 1,\n"\
"	// so there will be no divergence because of different values of \"i\"\n"\
"	if (smallPrimes[i] <= root4_N_rounded_up)\n"\
"	{\n"\
"		for (;;)\n"\
"		{\n"\
"			// Collect divisibility data for the next 16 primes\n"\
"			// No divergence here, compiler will use predicates for \"A ? B : C\" operator\n"\
"			uint pattern =\n"\
"				((N * primeInverses[i +  0].x <= primeInverses[i +  0].y) ? 32768 : 0) |\n"\
"				((N * primeInverses[i +  1].x <= primeInverses[i +  1].y) ? 16384 : 0) |\n"\
"				((N * primeInverses[i +  2].x <= primeInverses[i +  2].y) ? 8192  : 0) |\n"\
"				((N * primeInverses[i +  3].x <= primeInverses[i +  3].y) ? 4096  : 0) |\n"\
"				((N * primeInverses[i +  4].x <= primeInverses[i +  4].y) ? 2048  : 0) |\n"\
"				((N * primeInverses[i +  5].x <= primeInverses[i +  5].y) ? 1024  : 0) |\n"\
"				((N * primeInverses[i +  6].x <= primeInverses[i +  6].y) ? 512   : 0) |\n"\
"				((N * primeInverses[i +  7].x <= primeInverses[i +  7].y) ? 256   : 0) |\n"\
"				((N * primeInverses[i +  8].x <= primeInverses[i +  8].y) ? 128   : 0) |\n"\
"				((N * primeInverses[i +  9].x <= primeInverses[i +  9].y) ? 64    : 0) |\n"\
"				((N * primeInverses[i + 10].x <= primeInverses[i + 10].y) ? 32    : 0) |\n"\
"				((N * primeInverses[i + 11].x <= primeInverses[i + 11].y) ? 16    : 0) |\n"\
"				((N * primeInverses[i + 12].x <= primeInverses[i + 12].y) ? 8     : 0) |\n"\
"				((N * primeInverses[i + 13].x <= primeInverses[i + 13].y) ? 4     : 0) |\n"\
"				((N * primeInverses[i + 14].x <= primeInverses[i + 14].y) ? 2     : 0) |\n"\
"				((N * primeInverses[i + 15].x <= primeInverses[i + 15].y) ? 1     : 0);\n"\
"\n"\
"			i += 16;\n"\
"\n"\
"			// Do trial divisions by prime factors found\n"\
"			if (pattern)\n"\
"			{\n"\
"				ulong sumN = 1;\n"\
"				do\n"\
"				{\n"\
"					// Bit 0: divisibility for i - 1\n"\
"					// Bit 1: divisibility for i - 2\n"\
"					// Bit 2: divisibility for i - 3\n"\
"					// and so on ...\n"\
"					const uint num_leading_zeroes = clz(pattern);\n"\
"					pattern -= (1 << (31 - num_leading_zeroes));\n"\
"					const uint index = i + num_leading_zeroes - 32;\n"\
"\n"\
"					const ulong prevSumN = sumN;\n"\
"					const uint p1 = smallPrimes[index];\n"\
"					ulong q = N * primeInverses[index].x;\n"\
"					do\n"\
"					{\n"\
"						sumN = sumN * p1 + prevSumN;\n"\
"						N = q;\n"\
"						q *= primeInverses[index].x;\n"\
"					} while (q <= primeInverses[index].y);\n"\
"				} while (pattern);\n"\
"\n"\
"				const ulong p = smallPrimes[i];\n"\
"				if (p * p > N)\n"\
"				{\n"\
"					if (N > 1)\n"\
"					{\n"\
"						sumN *= N + 1;\n"\
"					}\n"\
"					if (sumN == targetSum)\n"\
"					{\n"\
"						NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, 0);\n"\
"					}\n"\
"					return;\n"\
"				}\n"\
"\n"\
"				const ulong oldTargetSum = targetSum;\n"\
"				targetSum /= sumN;\n"\
"				if (targetSum * sumN != oldTargetSum)\n"\
"				{\n"\
"					return;\n"\
"				}\n"\
"\n"\
"				if (N >= targetSum)\n"\
"				{\n"\
"					return;\n"\
"				}\n"\
"\n"\
"				root4_N_rounded_up = convert_uint_rtp(sqrt(sqrt(convert_float_rtp(N))));\n"\
"			}\n"\
"\n"\
"			if (smallPrimes[i] > root4_N_rounded_up)\n"\
"			{\n"\
"				break;\n"\
"			}\n"\
"\n"\
"			uint index = k * PQ_STRIDE_SIZE + (i >> 4) + 16;\n"\
"			while (N <= PQ[index].x)\n"\
"			{\n"\
"				--k;\n"\
"				index -= PQ_STRIDE_SIZE;\n"\
"			}\n"\
"			if (mul_hi(N, PQ[index].y) < targetSum - N)\n"\
"			{\n"\
"				return;\n"\
"			}\n"\
"		}\n"\
"	}\n"\
"\n"\
"	// No more than 3 factors remain here\n"\
"	\n"\
"	// The case when N is prime (it has 1 factor)\n"\
"	if (N + 1 == targetSum)\n"\
"	{\n"\
"		NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, N);\n"\
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
"			NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, sqrt_N);\n"\
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
"				NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, q);\n"\
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
"	// 1) M (8 bytes + 14 bits)\n"\
"	// 2) N (8 bytes + 14 bits)\n"\
"	// 3) targetSum (8 bytes + 14 bits)\n"\
"	// 4) i (<= 2798576, can fit in 22 bits)\n"\
"	// It can all fit in 32 bytes\n"\
"\n"\
"	const uint filtered_numbers_data_id = atom_inc(filtered_numbers_count);\n"\
"	filtered_numbers_data[filtered_numbers_data_id].x = M;\n"\
"	filtered_numbers_data[filtered_numbers_data_id].y = N;\n"\
"	filtered_numbers_data[filtered_numbers_data_id].z = targetSum;\n"\
"	const ulong i_to_save = i;\n"\
"	filtered_numbers_data[filtered_numbers_data_id].w = i_to_save | (M_high << 22) | (N_high << 36) | (targetSumHigh << 50);\n"\
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
"	__global ulong4* amicable_numbers_data,\n"\
"	const int max_i,\n"\
"	__global const uint* smallPrimes,\n"\
"	__global const uint* PowersOfP_128SumInverses_offsets,\n"\
"	__global const ulong4* PowersOfP_128SumInverses,\n"\
"	__global const ulong2* PrimeInverses_128\n"\
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
"\n"\
"	int i = cur_number.w & ((1 << 22) - 1);\n"\
"\n"\
"	const ulong M_high = (cur_number.w >> 22) & ((1 << 14) - 1);\n"\
"	ulong N_high = (cur_number.w >> 36) & ((1 << 14) - 1);\n"\
"	ulong targetSumHigh = cur_number.w >> 50;\n"\
"\n"\
"	if (targetSumHigh)\n"\
"	{\n"\
"		ulong2 N_128;\n"\
"		N_128.x = N;\n"\
"		N_128.y = N_high;\n"\
"		uint root3_N_rounded_up = IntegerCubeRoot128(N_128);\n"\
"\n"\
"		for (;(i & 31) || ((i < max_i) && (smallPrimes[i] <= root3_N_rounded_up) && targetSumHigh); ++i)\n"\
"		{\n"\
"			const ulong2 inv = PrimeInverses_128[i - 1];\n"\
"			ulong q_high = mul_hi(N, inv.x) + N * inv.y + N_high * inv.x;\n"\
"			if (q_high > N_high)\n"\
"			{\n"\
"				continue;\n"\
"			}\n"\
"\n"\
"			ulong q = N * inv.x;\n"\
"			int powerOfP = 0;\n"\
"			do\n"\
"			{\n"\
"				++powerOfP;\n"\
"\n"\
"				N = q;\n"\
"				N_high = q_high;\n"\
"\n"\
"				q = N * inv.x;\n"\
"				q_high = mul_hi(N, inv.x) + N * inv.y + N_high * inv.x;\n"\
"			} while (q_high <= N_high);\n"\
"\n"\
"			const ulong4 inv2 = PowersOfP_128SumInverses[PowersOfP_128SumInverses_offsets[i - 1] + powerOfP - 1];\n"\
"			if (inv2.z != 0)\n"\
"			{\n"\
"				if ((targetSum & inv2.w) != 0)\n"\
"				{\n"\
"					return;\n"\
"				}\n"\
"				targetSum = (targetSum >> inv2.z) | (targetSumHigh << (64 - inv2.z));\n"\
"				targetSumHigh >>= inv2.z;\n"\
"			}\n"\
"\n"\
"			q_high = mul_hi(targetSum, inv2.x) + targetSum * inv2.y + targetSumHigh * inv2.x;\n"\
"			if (q_high > targetSumHigh)\n"\
"			{\n"\
"				return;\n"\
"			}\n"\
"\n"\
"			targetSum *= inv2.x;\n"\
"			targetSumHigh = q_high;\n"\
"\n"\
"			const ulong p = smallPrimes[i];\n"\
"			if (p * p > N)\n"\
"			{\n"\
"				const ulong sumN = ((N > 1) ? N : 0) + 1;\n"\
"				if (sumN == targetSum)\n"\
"				{\n"\
"					NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, 0);\n"\
"				}\n"\
"				return;\n"\
"			}\n"\
"\n"\
"			if ((N_high > targetSumHigh) || ((N_high == targetSumHigh) && (N >= targetSum)))\n"\
"			{\n"\
"				return;\n"\
"			}\n"\
"\n"\
"			N_128.x = N;\n"\
"			N_128.y = N_high;\n"\
"			root3_N_rounded_up = IntegerCubeRoot128(N_128);\n"\
"		}\n"\
"\n"\
"		if (targetSumHigh)\n"\
"		{\n"\
"			if (smallPrimes[i] <= root3_N_rounded_up)\n"\
"			{\n"\
"				const uint filtered_numbers_data_id = atom_inc(filtered_numbers_count);\n"\
"				filtered_numbers_data[filtered_numbers_data_id].x = M;\n"\
"				filtered_numbers_data[filtered_numbers_data_id].y = N;\n"\
"				filtered_numbers_data[filtered_numbers_data_id].z = targetSum;\n"\
"				const ulong i_to_save = max_i;\n"\
"				filtered_numbers_data[filtered_numbers_data_id].w = i_to_save | (M_high << 22) | (N_high << 36) | (targetSumHigh << 50);\n"\
"			}\n"\
"			else\n"\
"			{\n"\
"				const ulong N2 = N + 1;\n"\
"				const ulong N2_high = N_high + ((N2 == 0) ? 1 : 0);\n"\
"				if (targetSum & 1)\n"\
"				{\n"\
"					const ulong carry = (targetSum < N2) ? 1 : 0;\n"\
"					const ulong p = targetSum - N2;\n"\
"					const ulong p_high = targetSumHigh - N2_high - carry;\n"\
"					if ((p_high == 0) && (p * p == N) && (mul_hi(p, p) == N_high))\n"\
"					{\n"\
"						NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, p);\n"\
"					}\n"\
"				}\n"\
"				else\n"\
"				{\n"\
"					if ((targetSum == N2) && (targetSumHigh == N2_high))\n"\
"					{\n"\
"						NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, 0);\n"\
"					}\n"\
"					else\n"\
"					{\n"\
"						const ulong carry = (targetSum < N2) ? 1 : 0;\n"\
"						ulong2 B;\n"\
"						B.x = targetSum - N2;\n"\
"						B.y = targetSumHigh - N2_high - carry;\n"\
"						B.x = (B.x >> 1) | (B.y << 63);\n"\
"						B.y >>= 1;\n"\
"						if (B.y == 0)\n"\
"						{\n"\
"							ulong2 D;\n"\
"							D.x = B.x * B.x;\n"\
"							D.y = mul_hi(B.x, B.x);\n"\
"							if ((D.y > N_high) || ((D.y == N_high) && (D.x > N)))\n"\
"							{\n"\
"								const ulong carry2 = (D.x < N) ? 1 : 0;\n"\
"								D.x -= N;\n"\
"								D.y -= N_high + carry2;\n"\
"								ulong sqrtD = IntegerSquareRoot128(D);\n"\
"								if ((sqrtD * sqrtD == D.x) && (mul_hi(sqrtD, sqrtD) == D.y))\n"\
"								{\n"\
"									const ulong p = B.x - sqrtD;\n"\
"									const ulong q = B.x + sqrtD;\n"\
"									if ((p * q == N) && (mul_hi(p, q) == N_high))\n"\
"									{\n"\
"										NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, p);\n"\
"									}\n"\
"								}\n"\
"							}\n"\
"						}\n"\
"					}\n"\
"				}\n"\
"			}\n"\
"			return;\n"\
"		}\n"\
"	}\n"\
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
"					NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, 0);\n"\
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
"			const ulong i_to_save = max_i;\n"\
"			filtered_numbers_data[filtered_numbers_data_id].w = i_to_save | (M_high << 22) | (N_high << 36) | (targetSumHigh << 50);\n"\
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
"			NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, N);\n"\
"		}\n"\
"		return;\n"\
"	}\n"\
"\n"\
"	// 2) The case when N is squared prime\n"\
"	const ulong sqrt_N = IntegerSquareRoot(N);\n"\
"	if ((sqrt_N * sqrt_N == N) && (N + sqrt_N + 1 == targetSum))\n"\
"	{\n"\
"		NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, sqrt_N);\n"\
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
"			NumberFound(amicable_numbers_count, amicable_numbers_data, M, M_high, q);\n"\
"			return;\n"\
"		}\n"\
"	}\n"\
"}\n"\
"\n"\
"__kernel\n"\
"__attribute__((reqd_work_group_size(1, 1, 1)))\n"\
"void SaveCounter(__global uint* phase2_numbers_count)\n"\
"{\n"\
"	if (phase2_numbers_count[0] <= PHASE2_MAX_COUNT)\n"\
"	{\n"\
"		phase2_numbers_count[1] = phase2_numbers_count[0];\n"\
"	}\n"\
"}\n"\
"\n"\
"// primes is an array of pointers to chunks of data\n"\
"// each chunk is 128 MB in size (2^24 uint2 elements) or bigger\n"\
"static ulong GetNthPrime(uint n,\n"\
"	__global uint2* primes0\n"\
"#if NUM_DATA_CHUNKS > 1\n"\
"	, __global uint2* primes1\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 2\n"\
"	, __global uint2* primes2\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 3\n"\
"	, __global uint2* primes3\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 4\n"\
"	, __global uint2* primes4\n"\
"#endif\n"\
")\n"\
"{\n"\
"#if NUM_DATA_CHUNKS == 1\n"\
"	const uint2 data = primes0[n >> 2];\n"\
"#else\n"\
"	const uint global_offset = n >> 2;\n"\
"	__global const uint2* chunk;\n"\
"	switch (global_offset >> (CHUNK_SIZE_SHIFT - 3))\n"\
"	{\n"\
"	default: chunk = primes0; break;\n"\
"#if NUM_DATA_CHUNKS > 1\n"\
"	case 1: chunk = primes1; break;\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 2\n"\
"	case 2: chunk = primes2; break;\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 3\n"\
"	case 3: chunk = primes3; break;\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 4\n"\
"	case 4: chunk = primes4; break;\n"\
"#endif\n"\
"	}\n"\
"\n"\
"	const uint chunk_offset = global_offset & ((1 << (CHUNK_SIZE_SHIFT - 3)) - 1);\n"\
"	const uint2 data = chunk[chunk_offset];\n"\
"#endif\n"\
"\n"\
"	ulong base = data.y & 31;\n"\
"	base = (base << 32) + data.x;\n"\
"\n"\
"	const ulong offsets = data.y >> 5;\n"\
"\n"\
"	return base + ((offsets >> ((~n & 3) * 9)) & 511);\n"\
"}\n"\
"\n"\
"__kernel\n"\
"__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))\n"\
"void SearchMultipleRanges(\n"\
"	__global const uint* smallPrimes,\n"\
"	__global const ulong2* primeInverses,\n"\
"	__global const ulong2* PQ,\n"\
"	__constant ulong2* PowerOf2SumInverses,\n"\
"	__global const uint* PowersOfP_128SumInverses_offsets,\n"\
"	__global const ulong4* PowersOfP_128SumInverses,\n"\
"	__global const ulong2* PrimeInverses_128,\n"\
"	__global ulong* SumEstimates_128,\n"\
"	__global const ulong* RangesTable,\n"\
"	__global const uint2* RangeLookupTable,\n"\
"	const uint lookup_shift,\n"\
"	const ulong global_offset,\n"\
"	const uint global_size,\n"\
"	__global ulong* phase1_offset_to_resume_after_overflow,\n"\
"	__global uint* phase2_numbers_count,\n"\
"	__global ulong4* phase2_numbers_data,\n"\
"	__global uint* amicable_numbers_count,\n"\
"	__global ulong4* amicable_numbers_data,\n"\
"\n"\
"	__global uint2* primes0\n"\
"#if NUM_DATA_CHUNKS > 1\n"\
"	, __global uint2* primes1\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 2\n"\
"	, __global uint2* primes2\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 3\n"\
"	, __global uint2* primes3\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 4\n"\
"	, __global uint2* primes4\n"\
"#endif\n"\
")\n"\
"{\n"\
"	if (*phase2_numbers_count > PHASE2_MAX_COUNT)\n"\
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
"	ulong2 M0, sumM0;\n"\
"	ulong p;\n"\
"\n"\
"	uint2 rangeIndexAndOffset = RangeLookupTable[globalIndex >> lookup_shift];\n"\
"	uint curIndex = globalIndex & (0xFFFFFFFFU << lookup_shift);\n"\
"	for (;;)\n"\
"	{\n"\
"		__global const ulong* curRange = RangesTable + rangeIndexAndOffset.x * 5;\n"\
"		const uint numbersRemainingInCurRange = (curRange[4] >> 32) - rangeIndexAndOffset.y;\n"\
"\n"\
"		// curIndex + numbersRemainingInCurRange now points exactly at the beginning of the next range\n"\
"		// If the beginning of the next range is > globalIndex, then globalIndex is within the current range\n"\
"		if (curIndex + numbersRemainingInCurRange > globalIndex)\n"\
"		{\n"\
"			M0.x = curRange[0];\n"\
"			M0.y = curRange[1];\n"\
"			sumM0.x = curRange[2];\n"\
"			sumM0.y = curRange[3];\n"\
"\n"\
"			//      primes[ start_prime_index            + range_offset          + globalIndex - curIndex]\n"\
"			p = GetNthPrime((curRange[4] & 0xFFFFFFFFUL) + rangeIndexAndOffset.y + globalIndex - curIndex,\n"\
"				primes0\n"\
"#if NUM_DATA_CHUNKS > 1\n"\
"				, primes1\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 2\n"\
"				, primes2\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 3\n"\
"				, primes3\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 4\n"\
"				, primes4\n"\
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
"#if LPP == 1\n"\
"	ulong value = p;\n"\
"	ulong value_sum = p + 1;\n"\
"#elif LPP == 2\n"\
"	const ulong value = p * p;\n"\
"	const ulong value_sum = p * (p + 1) + 1;\n"\
"#elif LPP == 3\n"\
"	const ulong value = p * p * p;\n"\
"	const ulong value_sum = p * (p * (p + 1) + 1) + 1;\n"\
"#endif\n"\
"\n"\
"	CheckPairPhase1(\n"\
"		smallPrimes,\n"\
"		primeInverses,\n"\
"		PQ,\n"\
"		PowerOf2SumInverses,\n"\
"		PowersOfP_128SumInverses_offsets,\n"\
"		PowersOfP_128SumInverses,\n"\
"		PrimeInverses_128,\n"\
"		SumEstimates_128,\n"\
"		M0.x * value, mul_hi(M0.x, value) + M0.y * value, sumM0.x * value_sum, mul_hi(sumM0.x, value_sum) + sumM0.y * value_sum,\n"\
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
"void SearchLargePrimes(\n"\
"	__global const uint* smallPrimes,\n"\
"	__global const ulong2* primeInverses,\n"\
"	__global const ulong2* PQ,\n"\
"	__constant ulong2* PowerOf2SumInverses,\n"\
"	__global const uint* PowersOfP_128SumInverses_offsets,\n"\
"	__global const ulong4* PowersOfP_128SumInverses,\n"\
"	__global const ulong2* PrimeInverses_128,\n"\
"	__global const ulong* SumEstimates_128,\n"\
"	__global const ulong* largePrimes,\n"\
"	const uint largePrimesCount,\n"\
"	const ulong largePrimesCountReciprocal,\n"\
"	const uint largePrimesCountIncrementAndShift,\n"\
"	const ulong global_offset,\n"\
"	const ulong global_size,\n"\
"	__global ulong* phase1_offset_to_resume_after_overflow,\n"\
"	__global uint* phase2_numbers_count,\n"\
"	__global ulong4* phase2_numbers_data,\n"\
"	__global uint* amicable_numbers_count,\n"\
"	__global ulong4* amicable_numbers_data,\n"\
"\n"\
"	__global ulong2* amicableCandidates0\n"\
"#if NUM_DATA_CHUNKS > 1\n"\
"	, __global ulong2* amicableCandidates1\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 2\n"\
"	, __global ulong2* amicableCandidates2\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 3\n"\
"	, __global ulong2* amicableCandidates3\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 4\n"\
"	, __global ulong2* amicableCandidates4\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 5\n"\
"	, __global ulong2* amicableCandidates5\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 6\n"\
"	, __global ulong2* amicableCandidates6\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 7\n"\
"	, __global ulong2* amicableCandidates7\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 8\n"\
"	, __global ulong2* amicableCandidates8\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 9\n"\
"	, __global ulong2* amicableCandidates9\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 10\n"\
"	, __global ulong2* amicableCandidates10\n"\
"#endif\n"\
")\n"\
"{\n"\
"	if (*phase2_numbers_count > PHASE2_MAX_COUNT)\n"\
"	{\n"\
"		return;\n"\
"	}\n"\
"\n"\
"	const ulong globalIndex = global_offset + get_global_id(0);\n"\
"	if (globalIndex >= global_size)\n"\
"	{\n"\
"		return;\n"\
"	}\n"\
"\n"\
"	ulong amicableCandidateIndex;\n"\
"	uint largePrimeIndex;\n"\
"\n"\
"	if (largePrimesCountReciprocal != 0)\n"\
"	{\n"\
"		amicableCandidateIndex = mul_hi(globalIndex + (largePrimesCountIncrementAndShift & 1), largePrimesCountReciprocal);\n"\
"		amicableCandidateIndex >>= (largePrimesCountIncrementAndShift >> 1);\n"\
"		largePrimeIndex = globalIndex - amicableCandidateIndex * largePrimesCount;\n"\
"	}\n"\
"	else\n"\
"	{\n"\
"		amicableCandidateIndex = globalIndex >> largePrimesCountIncrementAndShift;\n"\
"		largePrimeIndex = globalIndex - (amicableCandidateIndex << largePrimesCountIncrementAndShift);\n"\
"	}\n"\
"\n"\
"#if NUM_DATA_CHUNKS == 1\n"\
"	const ulong2 value = amicableCandidates0[amicableCandidateIndex];\n"\
"#else\n"\
"	__global const ulong2* chunk;\n"\
"	switch (amicableCandidateIndex >> (CHUNK_SIZE_SHIFT - 4))\n"\
"	{\n"\
"	default: chunk = amicableCandidates0; break;\n"\
"#if NUM_DATA_CHUNKS > 1\n"\
"	case 1: chunk = amicableCandidates1; break;\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 2\n"\
"	case 2: chunk = amicableCandidates2; break;\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 3\n"\
"	case 3: chunk = amicableCandidates3; break;\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 4\n"\
"	case 4: chunk = amicableCandidates4; break;\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 5\n"\
"	case 5: chunk = amicableCandidates5; break;\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 6\n"\
"	case 6: chunk = amicableCandidates6; break;\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 7\n"\
"	case 7: chunk = amicableCandidates7; break;\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 8\n"\
"	case 8: chunk = amicableCandidates8; break;\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 9\n"\
"	case 9: chunk = amicableCandidates9; break;\n"\
"#endif\n"\
"#if NUM_DATA_CHUNKS > 10\n"\
"	case 10: chunk = amicableCandidates10; break;\n"\
"#endif\n"\
"	}\n"\
"	const uint chunk_offset = amicableCandidateIndex & ((1 << (CHUNK_SIZE_SHIFT - 4)) - 1);\n"\
"	const ulong2 value = chunk[chunk_offset];\n"\
"#endif\n"\
"\n"\
"	const ulong p = largePrimes[largePrimeIndex];\n"\
"\n"\
"	CheckPairPhase1(\n"\
"		smallPrimes,\n"\
"		primeInverses,\n"\
"		PQ,\n"\
"		PowerOf2SumInverses,\n"\
"		PowersOfP_128SumInverses_offsets,\n"\
"		PowersOfP_128SumInverses,\n"\
"		PrimeInverses_128,\n"\
"		SumEstimates_128,\n"\
"		value.x * p, mul_hi(value.x, p), value.y * (p + 1), mul_hi(value.y, p + 1),\n"\
"		global_offset,\n"\
"		phase1_offset_to_resume_after_overflow,\n"\
"		phase2_numbers_count,\n"\
"		phase2_numbers_data,\n"\
"		amicable_numbers_count,\n"\
"		amicable_numbers_data\n"\
"	);\n"\
"}\n"\
"\n"\
"\n"\
"__kernel\n"\
"__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))\n"\
"void CheckPairs(\n"\
"	__global const uint* smallPrimes,\n"\
"	__global const ulong2* primeInverses,\n"\
"	__global const ulong2* PQ,\n"\
"	__constant ulong2* PowerOf2SumInverses,\n"\
"	__global const uint* PowersOfP_128SumInverses_offsets,\n"\
"	__global const ulong4* PowersOfP_128SumInverses,\n"\
"	__global const ulong2* PrimeInverses_128,\n"\
"	__global const ulong* SumEstimates_128,\n"\
"	__global const ulong4* pairsToCheck,\n"\
"	__global uint* filtered_numbers_count,\n"\
"	__global ulong4* filtered_numbers_data,\n"\
"	__global uint* amicable_numbers_count,\n"\
"	__global ulong4* amicable_numbers_data\n"\
")\n"\
"{\n"\
"	const ulong4 cur_pair = pairsToCheck[get_global_id(0)];\n"\
"	if (cur_pair.x == 0)\n"\
"	{\n"\
"		return;\n"\
"	}\n"\
"	CheckPairPhase1(smallPrimes, primeInverses, PQ, PowerOf2SumInverses, PowersOfP_128SumInverses_offsets, PowersOfP_128SumInverses, PrimeInverses_128, SumEstimates_128, cur_pair.x, cur_pair.y, cur_pair.z, cur_pair.w, 0, 0, filtered_numbers_count, filtered_numbers_data, amicable_numbers_count, amicable_numbers_data);\n"\
"}\n"\
"";

static const unsigned int kernel_cl_crc32 = 0xb2dfcdac;
