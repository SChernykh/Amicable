#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Check amicability for a range of numbers of the form "M0*p" where p is prime
// p runs through all values in primes[p_offset]...primes[p_offset + get_global_size(0) - 1] (inclusive)
// primes stores prime numbers in compact form: 4 primes per each 8 bytes
// primeInverses is an array of multiplicative inverses of first 192725 prime numbers, this is used for checking divisibility by doing just one multiplication
// primeReciprocals is an array of reciprocals of first 192725 prime numbers, this is used for calculating N / p (not just checking divisibility)
// PQ corresponds to tables P and Q in lemma 2.1 from "Computation of All the Amicable Pairs Below 10^10 By H.J.J.te Riele": http://www.ams.org/journals/mcom/1986-47-175/S0025-5718-1986-0842142-3/S0025-5718-1986-0842142-3.pdf
// The only difference is that we calculate exact (hence better) upper bounds for S(m)/m instead of inexact estimates
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// A: smaller member of an amicable pair
// B: if it's not zero, it must be prime in order for A to be amicable number
void NumberFound(__global uint* amicable_numbers_count, __global ulong2* amicable_numbers_data, ulong A, ulong B)
{
	const uint index = atom_inc(amicable_numbers_count);
	amicable_numbers_data[index].x = A; 
	amicable_numbers_data[index].y = B;
}

// 32-bit floats don't have enough precision
// 64-bit floats are not supported by all OpenCL devices
// So we have to use integer arithmetic to calculate square roots

// Returns number x such that x^2 <= n < (x+1)^2
uint IntegerSquareRoot(const ulong n)
{
	uint result = 1;
	result <<= ((63 - clz(n)) >> 1);
	for (uint cur_bit = result >> 1; cur_bit > 0; cur_bit >>= 1)
	{
		const uint result1 = result | cur_bit;
		const ulong k = result1;
		if (k * k <= n)
		{
			result = result1;
		}
	}
	return result;
}

// Returns number x such that x^2 <= n < (x+1)^2
ulong IntegerSquareRoot128(const ulong2 n)
{
	ulong result = 1;
	const ulong highest_bit_index = n.y ? (127 - clz(n.y)) : (63 - clz(n.x));
	result <<= (highest_bit_index >> 1);
	for (ulong cur_bit = result >> 1; cur_bit > 0; cur_bit >>= 1)
	{
		const ulong k = result | cur_bit;
		ulong2 k2;
		k2.x = k * k;
		k2.y = mul_hi(k, k);
		if ((k2.y < n.y) || ((k2.y == n.y) && (k2.x <= n.x)))
		{
			result = k;
		}
	}
	return result;
}

void CheckPairPhase1(
	__global const uint* smallPrimes,
	__global const ulong2* primeInverses,
	__global const ulong2* PQ,
	__constant const ulong2* PowerOf2SumInverses,
	__constant const ulong4* PowersOfP_128SumInverses,
	__constant const ulong2* PrimeInverses_128,
	const ulong M,
	ulong targetSum,
	ulong targetSumHigh,
	const uint global_offset,
	__global uint* phase1_offset_to_resume_after_overflow,
	__global uint* phase2_numbers_count,
	__global ulong4* phase2_numbers_data,
	__global uint* amicable_numbers_count,
	__global ulong2* amicable_numbers_data
)
{
	ulong N = targetSum - M;
	ulong N_high = targetSumHigh - ((targetSum < M) ? 1 : 0);

	// Divide out power of 2 using fast bitwise operations
	// All numbers in the same range have the same parity, so threads won't diverge here
	if ((N & 1) == 0)
	{
		if (N == 0)
		{
			// N = 2^64*k, k >= 0. It's not amicable for k = 0, 1, 2 and too large for k > 2.
			return;
		}

		const uint powerOf2 = 63 - clz(N ^ (N - 1));

		// M and targetSum (128-bit) come in ranges where they both grow monotonically,
		// so there will be exactly one divergence in each range and one divergence on a border between ranges
		// Threads don't diverge 99.9% of the time here
		if (targetSumHigh > 0)
		{
			const ulong2 inv = PowerOf2SumInverses[powerOf2 + 64];
			if (mul_hi(targetSum, inv.x) + targetSum * inv.y + targetSumHigh * inv.x != 0)
			{
				return;
			}

			targetSum = targetSum * inv.x;
			targetSumHigh = 0;

			N = (N >> powerOf2) | (N_high << (64 - powerOf2));
		}
		else
		{
			targetSum *= PowerOf2SumInverses[powerOf2].x;
			if (targetSum > PowerOf2SumInverses[powerOf2].y)
			{
				return;
			}
			N >>= powerOf2;
		}
	}

	int i = 1;

	// M and targetSum (128-bit) come in ranges where they both grow monotonically,
	// so there will be exactly one divergence in each range and one divergence on a border between ranges
	// Threads don't diverge 99.9% of the time here
	if (targetSumHigh > 0)
	{
		do
		{
			const ulong2 inv = PrimeInverses_128[i - 1];
			if (mul_hi(N, inv.x) + N * inv.y + N_high * inv.x == 0)
			{
				ulong q = N * inv.x;
				int powerOfP = 0;
				do
				{
					++powerOfP;
					N = q;
					q *= primeInverses[i].x;
				} while (q <= primeInverses[i].y);

				const ulong4 inv2 = PowersOfP_128SumInverses[(i - 1) * 64 + powerOfP];
				if (inv2.z != 0)
				{
					if (targetSum & ((1 << inv2.z) - 1))
					{
						return;
					}
					targetSum = (targetSum >> inv2.z) | (targetSumHigh << (64 - inv2.z));
					targetSumHigh >>= inv2.z;
				}

				if (mul_hi(targetSum, inv2.x) + targetSum * inv2.y + targetSumHigh * inv2.x != 0)
				{
					return;
				}

				targetSum *= inv2.x;
				targetSumHigh = 0;

				const ulong p2 = (0x31190900 >> (i << 3)) & 0xFF;
				// i = 1: (0x31190900 >>  8) & 0xFF =  9 = 3^2
				// i = 2: (0x31190900 >> 16) & 0xFF = 25 = 5^2
				// i = 3: (0x31190900 >> 24) & 0xFF = 49 = 7^2
				if (p2 > N)
				{
					const ulong sumN = ((N > 1) ? N : 0) + 1;
					if (sumN == targetSum)
					{
						NumberFound(amicable_numbers_count, amicable_numbers_data, M, 0);
					}
					return;
				}

				if (N >= targetSum)
				{
					return;
				}
			}
			++i;
		} while ((i <= 3) && (targetSumHigh != 0));
	}

	// Collect divisibility data for first 15 odd primes
	// No divergence here, compiler will use predicates for "A ? B : C" operator
	uint pattern =
		((N * primeInverses[ 1].x <= primeInverses[ 1].y) ? 16384 : 0) |
		((N * primeInverses[ 2].x <= primeInverses[ 2].y) ? 8192  : 0) |
		((N * primeInverses[ 3].x <= primeInverses[ 3].y) ? 4096  : 0) |
		((N * primeInverses[ 4].x <= primeInverses[ 4].y) ? 2048  : 0) |
		((N * primeInverses[ 5].x <= primeInverses[ 5].y) ? 1024  : 0) |
		((N * primeInverses[ 6].x <= primeInverses[ 6].y) ? 512   : 0) |
		((N * primeInverses[ 7].x <= primeInverses[ 7].y) ? 256   : 0) |
		((N * primeInverses[ 8].x <= primeInverses[ 8].y) ? 128   : 0) |
		((N * primeInverses[ 9].x <= primeInverses[ 9].y) ? 64    : 0) |
		((N * primeInverses[10].x <= primeInverses[10].y) ? 32    : 0) |
		((N * primeInverses[11].x <= primeInverses[11].y) ? 16    : 0) |
		((N * primeInverses[12].x <= primeInverses[12].y) ? 8     : 0) |
		((N * primeInverses[13].x <= primeInverses[13].y) ? 4     : 0) |
		((N * primeInverses[14].x <= primeInverses[14].y) ? 2     : 0) |
		((N * primeInverses[15].x <= primeInverses[15].y) ? 1     : 0);

	i = 16;

	// Do trial divisions by primes found
	if (pattern)
	{
		ulong sumN = 1;
		do
		{
			// Bit 0: divisibility for i == 15
			// Bit 1: divisibility for i == 14
			// Bit 2: divisibility for i == 13
			// and so on ...
			const uint num_leading_zeroes = clz(pattern);
			pattern -= (1 << (31 - num_leading_zeroes));
			const uint index = num_leading_zeroes - 16;

			const ulong prevSumN = sumN;
			const uint p1 = smallPrimes[index];
			ulong q = N * primeInverses[index].x;
			do
			{
				sumN = sumN * p1 + prevSumN;
				N = q;
				q *= primeInverses[index].x;
			} while (q <= primeInverses[index].y);
		} while (pattern);

		if (N < 59 * 59)
		{
			if (N > 1)
			{
				sumN *= N + 1;
			}
			if (sumN == targetSum)
			{
				NumberFound(amicable_numbers_count, amicable_numbers_data, M, 0);
			}
			return;
		}

		const ulong oldTargetSum = targetSum;
		targetSum /= sumN;
		if (targetSum * sumN != oldTargetSum)
		{
			return;
		}

		if (N >= targetSum)
		{
			return;
		}
	}

	// Check maximal possible sum for N. If it's too small, it can't be an amicable number
	int k = 8;
	while (N <= PQ[k].x)
	{
		--k;
	}
	++k;
	if (mul_hi(N, PQ[k].y) < targetSum - N)
	{
		return;
	}

	do
	{
		// Collect divisibility data for the next 16 primes
		// No divergence here, compiler will use predicates for "A ? B : C" operator
		pattern =
			((N * primeInverses[i +  0].x <= primeInverses[i +  0].y) ? 32768 : 0) |
			((N * primeInverses[i +  1].x <= primeInverses[i +  1].y) ? 16384 : 0) |
			((N * primeInverses[i +  2].x <= primeInverses[i +  2].y) ? 8192  : 0) |
			((N * primeInverses[i +  3].x <= primeInverses[i +  3].y) ? 4096  : 0) |
			((N * primeInverses[i +  4].x <= primeInverses[i +  4].y) ? 2048  : 0) |
			((N * primeInverses[i +  5].x <= primeInverses[i +  5].y) ? 1024  : 0) |
			((N * primeInverses[i +  6].x <= primeInverses[i +  6].y) ? 512   : 0) |
			((N * primeInverses[i +  7].x <= primeInverses[i +  7].y) ? 256   : 0) |
			((N * primeInverses[i +  8].x <= primeInverses[i +  8].y) ? 128   : 0) |
			((N * primeInverses[i +  9].x <= primeInverses[i +  9].y) ? 64    : 0) |
			((N * primeInverses[i + 10].x <= primeInverses[i + 10].y) ? 32    : 0) |
			((N * primeInverses[i + 11].x <= primeInverses[i + 11].y) ? 16    : 0) |
			((N * primeInverses[i + 12].x <= primeInverses[i + 12].y) ? 8     : 0) |
			((N * primeInverses[i + 13].x <= primeInverses[i + 13].y) ? 4     : 0) |
			((N * primeInverses[i + 14].x <= primeInverses[i + 14].y) ? 2     : 0) |
			((N * primeInverses[i + 15].x <= primeInverses[i + 15].y) ? 1     : 0);

		i += 16;

		// Do trial divisions by prime factors found
		if (pattern)
		{
			ulong sumN = 1;
			do
			{
				// Bit 0: divisibility for i - 1
				// Bit 1: divisibility for i - 2
				// Bit 2: divisibility for i - 3
				// and so on ...
				const uint num_leading_zeroes = clz(pattern);
				pattern -= (1 << (31 - num_leading_zeroes));
				const uint index = i + num_leading_zeroes - 32;

				const ulong prevSumN = sumN;
				const uint p1 = smallPrimes[index];
				ulong q = N * primeInverses[index].x;
				do
				{
					sumN = sumN * p1 + prevSumN;
					N = q;
					q *= primeInverses[index].x;
				} while (q <= primeInverses[index].y);
			} while (pattern);

			const ulong p2 = smallPrimes[i];
			if (p2 * p2 > N)
			{
				if (N > 1)
				{
					sumN *= N + 1;
				}
				if (sumN == targetSum)
				{
					NumberFound(amicable_numbers_count, amicable_numbers_data, M, 0);
				}
				return;
			}

			const ulong oldTargetSum = targetSum;
			targetSum /= sumN;
			if (targetSum * sumN != oldTargetSum)
			{
				return;
			}

			if (N >= targetSum)
			{
				return;
			}
		}

		// Check maximal possible sum for N. If it's too small, it can't be an amicable number
		// Same i for all threads, no divergence
		uint index = k * PQ_STRIDE_SIZE + (i >> 4) + 16;
		while (N <= PQ[index].x)
		{
			--k;
			index -= PQ_STRIDE_SIZE;
		}
		if (mul_hi(N, PQ[index].y) < targetSum - N)
		{
			return;
		}
	} while (i != 256);

	// Data needed to resume from arbitrary loop iteration in phase 2:
	// 1) M (8 bytes)
	// 2) N (8 bytes)
	// 3) targetSum (8 bytes)
	// 4) i (<= 192724, can fit in 18 bits)
	// 5) k (< 16, can fit in 4 bits)
	// It can all fit in 32 bytes

	const uint phase2_number_index = atom_inc(phase2_numbers_count);
	if (phase2_number_index < PHASE2_MAX_COUNT)
	{
		const ulong i_to_save = i;
		const ulong k_to_save = k;
		phase2_numbers_data[phase2_number_index].x = M;
		phase2_numbers_data[phase2_number_index].y = N;
		phase2_numbers_data[phase2_number_index].z = targetSum;
		phase2_numbers_data[phase2_number_index].w = i_to_save | (k_to_save << 32);
	}
	else if (phase2_number_index == PHASE2_MAX_COUNT)
	{
		*phase1_offset_to_resume_after_overflow = global_offset;
	}
}

__kernel
__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
void CheckPairPhase2(
	__global const uint* smallPrimes,
	__global const ulong4* numbersToCheck,
	__global const ulong2* primeInverses,
	__global const ulong2* PQ,
	__global uint* filtered_numbers_count,
	__global ulong4* filtered_numbers_data,
	__global uint* amicable_numbers_count,
	__global ulong2* amicable_numbers_data
)
{
	const ulong4 cur_number = numbersToCheck[get_global_id(0)];
	if (cur_number.x == 0)
	{
		return;
	}

	const ulong M = cur_number.x;
	ulong N = cur_number.y;
	ulong targetSum = cur_number.z;
	int i = cur_number.w & 0xFFFFFFFF;
	int k = cur_number.w >> 32;

	const uint root4_N_rounded_up = convert_uint_rtp(sqrt(sqrt(convert_float_rtp(N))));

	// All threads start with the same value of "i" here, which is set to 256 at the end of phase 1,
	// so there will be no divergence because of different values of "i"
	do
	{
		// Collect divisibility data for the next 16 primes
		// No divergence here, compiler will use predicates for "A ? B : C" operator
		uint pattern =
			((N * primeInverses[i +  0].x <= primeInverses[i +  0].y) ? 32768 : 0) |
			((N * primeInverses[i +  1].x <= primeInverses[i +  1].y) ? 16384 : 0) |
			((N * primeInverses[i +  2].x <= primeInverses[i +  2].y) ? 8192  : 0) |
			((N * primeInverses[i +  3].x <= primeInverses[i +  3].y) ? 4096  : 0) |
			((N * primeInverses[i +  4].x <= primeInverses[i +  4].y) ? 2048  : 0) |
			((N * primeInverses[i +  5].x <= primeInverses[i +  5].y) ? 1024  : 0) |
			((N * primeInverses[i +  6].x <= primeInverses[i +  6].y) ? 512   : 0) |
			((N * primeInverses[i +  7].x <= primeInverses[i +  7].y) ? 256   : 0) |
			((N * primeInverses[i +  8].x <= primeInverses[i +  8].y) ? 128   : 0) |
			((N * primeInverses[i +  9].x <= primeInverses[i +  9].y) ? 64    : 0) |
			((N * primeInverses[i + 10].x <= primeInverses[i + 10].y) ? 32    : 0) |
			((N * primeInverses[i + 11].x <= primeInverses[i + 11].y) ? 16    : 0) |
			((N * primeInverses[i + 12].x <= primeInverses[i + 12].y) ? 8     : 0) |
			((N * primeInverses[i + 13].x <= primeInverses[i + 13].y) ? 4     : 0) |
			((N * primeInverses[i + 14].x <= primeInverses[i + 14].y) ? 2     : 0) |
			((N * primeInverses[i + 15].x <= primeInverses[i + 15].y) ? 1     : 0);

		i += 16;

		// Do trial divisions by prime factors found
		if (pattern)
		{
			ulong sumN = 1;
			do
			{
				// Bit 0: divisibility for i - 1
				// Bit 1: divisibility for i - 2
				// Bit 2: divisibility for i - 3
				// and so on ...
				const uint num_leading_zeroes = clz(pattern);
				pattern -= (1 << (31 - num_leading_zeroes));
				const uint index = i + num_leading_zeroes - 32;

				const ulong prevSumN = sumN;
				const uint p1 = smallPrimes[index];
				ulong q = N * primeInverses[index].x;
				do
				{
					sumN = sumN * p1 + prevSumN;
					N = q;
					q *= primeInverses[index].x;
				} while (q <= primeInverses[index].y);
			} while (pattern);

			const ulong p2 = smallPrimes[i];
			if (p2 * p2 > N)
			{
				if (N > 1)
				{
					sumN *= N + 1;
				}
				if (sumN == targetSum)
				{
					NumberFound(amicable_numbers_count, amicable_numbers_data, M, 0);
				}
				return;
			}

			const ulong oldTargetSum = targetSum;
			targetSum /= sumN;
			if (targetSum * sumN != oldTargetSum)
			{
				return;
			}

			if (N >= targetSum)
			{
				return;
			}
		}

		uint index = k * PQ_STRIDE_SIZE + (i >> 4) + 16;
		while (N <= PQ[index].x)
		{
			--k;
			index -= PQ_STRIDE_SIZE;
		}
		if (mul_hi(N, PQ[index].y) < targetSum - N)
		{
			return;
		}
	} while (smallPrimes[i] <= root4_N_rounded_up);

	// No more than 3 factors remain here
	
	// The case when N is prime (it has 1 factor)
	if (N + 1 == targetSum)
	{
		NumberFound(amicable_numbers_count, amicable_numbers_data, M, N);
		return;
	}

	const ulong root3_N_rounded_down = convert_ulong_rtn(cbrt(convert_float_rtn(N)));

	// The case when N can only have 2 factors
	if (N >= targetSum - root3_N_rounded_down * (root3_N_rounded_down + 1))
	{
		// The case when N is squared prime
		const ulong sqrt_N = IntegerSquareRoot(N);
		if ((sqrt_N * sqrt_N == N) && (N + sqrt_N + 1 == targetSum))
		{
			NumberFound(amicable_numbers_count, amicable_numbers_data, M, sqrt_N);
			return;
		}

		// The case when N = p * q - a product of two different primes
		// p + q == targetSum - N - 1
		// p * q == N
		const ulong B = targetSum - N - 1;
		ulong2 B2;
		B2.x = B * B;
		B2.y = mul_hi(B, B);
		const ulong carry = (B2.x < (N << 2)) ? 1 : 0;
		B2.x -= (N << 2);
		if (B2.y >= (N >> 62) + carry)
		{
			B2.y -= (N >> 62) + carry;
			const ulong sqrt_D = IntegerSquareRoot128(B2);
			const ulong p = (B - sqrt_D) >> 1;
			const ulong q = (B + sqrt_D) >> 1;
			if ((p < q) && (p + q == B) && (p * q == N))
			{
				NumberFound(amicable_numbers_count, amicable_numbers_data, M, q);
				return;
			}
		}
		return;
	}

	// N can have 3 factors
	// Check maximal possible sum of divisors for numbers with 3 factors
	uint p1 = smallPrimes[i];
	ulong q = N / p1;

	ulong p2 = p1;
	p2 *= p1 + 1;

	// If it's too small, N can't be amicable
	if (N < targetSum - ((q + p2) << 1))
	{
		return;
	}

	// Data needed to resume from arbitrary loop iteration in phase 3:
	// 1) M (8 bytes)
	// 2) N (8 bytes)
	// 3) targetSum (8 bytes)
	// 4) i (<= 192724, can fit in 18 bits)
	// It can all fit in 32 bytes

	const uint filtered_numbers_data_id = atom_inc(filtered_numbers_count);
	filtered_numbers_data[filtered_numbers_data_id].x = M;
	filtered_numbers_data[filtered_numbers_data_id].y = N;
	filtered_numbers_data[filtered_numbers_data_id].z = targetSum;
	filtered_numbers_data[filtered_numbers_data_id].w = i;
}

__kernel
__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
void CheckPairPhase3(
	__global const ulong4* numbersToCheck,
	const uint numbersToCheckOffset,
	__global const ulong2* primeInverses,
	__global const ulong2* primeReciprocals,
	__global uint* filtered_numbers_count,
	__global ulong4* filtered_numbers_data,
	__global uint* amicable_numbers_count,
	__global ulong2* amicable_numbers_data,
	const int max_i
)
{
	const ulong4 cur_number = numbersToCheck[get_global_id(0) + numbersToCheckOffset];
	if (cur_number.x == 0)
	{
		return;
	}

	const ulong M = cur_number.x;
	ulong N = cur_number.y;
	ulong targetSum = cur_number.z;
	int i = cur_number.w;

	// No more than 3 factors remain here
	const ulong root3_N_rounded_up = convert_uint_rtp(cbrt(convert_float_rtp(N)));
	const ulong MaximumSumOfDivisors3_target = targetSum - N - 1;

	// All threads start with the same value of "i" here except for "Phase 3, iteration 1"
	// so there will be no divergence because of different values of "i" on iteration 2, 3, etc.
	// "Phase 3, iteration 1" is deliberately made shorter because of this
	for (;;)
	{
		uint pattern =
			((N * primeInverses[i +  0].x <= primeInverses[i +  0].y) ? (1U << 31) : 0) |
			((N * primeInverses[i +  1].x <= primeInverses[i +  1].y) ? (1U << 30) : 0) |
			((N * primeInverses[i +  2].x <= primeInverses[i +  2].y) ? (1U << 29) : 0) |
			((N * primeInverses[i +  3].x <= primeInverses[i +  3].y) ? (1U << 28) : 0) |
			((N * primeInverses[i +  4].x <= primeInverses[i +  4].y) ? (1U << 27) : 0) |
			((N * primeInverses[i +  5].x <= primeInverses[i +  5].y) ? (1U << 26) : 0) |
			((N * primeInverses[i +  6].x <= primeInverses[i +  6].y) ? (1U << 25) : 0) |
			((N * primeInverses[i +  7].x <= primeInverses[i +  7].y) ? (1U << 24) : 0) |
			((N * primeInverses[i +  8].x <= primeInverses[i +  8].y) ? (1U << 23) : 0) |
			((N * primeInverses[i +  9].x <= primeInverses[i +  9].y) ? (1U << 22) : 0) |
			((N * primeInverses[i + 10].x <= primeInverses[i + 10].y) ? (1U << 21) : 0) |
			((N * primeInverses[i + 11].x <= primeInverses[i + 11].y) ? (1U << 20) : 0) |
			((N * primeInverses[i + 12].x <= primeInverses[i + 12].y) ? (1U << 19) : 0) |
			((N * primeInverses[i + 13].x <= primeInverses[i + 13].y) ? (1U << 18) : 0) |
			((N * primeInverses[i + 14].x <= primeInverses[i + 14].y) ? (1U << 17) : 0) |
			((N * primeInverses[i + 15].x <= primeInverses[i + 15].y) ? (1U << 16) : 0) |
			((N * primeInverses[i + 16].x <= primeInverses[i + 16].y) ? (1U << 15) : 0) |
			((N * primeInverses[i + 17].x <= primeInverses[i + 17].y) ? (1U << 14) : 0) |
			((N * primeInverses[i + 18].x <= primeInverses[i + 18].y) ? (1U << 13) : 0) |
			((N * primeInverses[i + 19].x <= primeInverses[i + 19].y) ? (1U << 12) : 0) |
			((N * primeInverses[i + 20].x <= primeInverses[i + 20].y) ? (1U << 11) : 0) |
			((N * primeInverses[i + 21].x <= primeInverses[i + 21].y) ? (1U << 10) : 0) |
			((N * primeInverses[i + 22].x <= primeInverses[i + 22].y) ? (1U <<  9) : 0) |
			((N * primeInverses[i + 23].x <= primeInverses[i + 23].y) ? (1U <<  8) : 0) |
			((N * primeInverses[i + 24].x <= primeInverses[i + 24].y) ? (1U <<  7) : 0) |
			((N * primeInverses[i + 25].x <= primeInverses[i + 25].y) ? (1U <<  6) : 0) |
			((N * primeInverses[i + 26].x <= primeInverses[i + 26].y) ? (1U <<  5) : 0) |
			((N * primeInverses[i + 27].x <= primeInverses[i + 27].y) ? (1U <<  4) : 0) |
			((N * primeInverses[i + 28].x <= primeInverses[i + 28].y) ? (1U <<  3) : 0) |
			((N * primeInverses[i + 29].x <= primeInverses[i + 29].y) ? (1U <<  2) : 0) |
			((N * primeInverses[i + 30].x <= primeInverses[i + 30].y) ? (1U <<  1) : 0) |
			((N * primeInverses[i + 31].x <= primeInverses[i + 31].y) ? (1U <<  0) : 0);

		i += 32;

		// Do trial divisions by prime factors found
		if (pattern)
		{
			ulong sumN = 1;
			do
			{
				// Bit 0: divisibility for i - 1
				// Bit 1: divisibility for i - 2
				// Bit 2: divisibility for i - 3
				// and so on ...
				const uint num_leading_zeroes = clz(pattern);
				pattern -= (1 << (31 - num_leading_zeroes));
				const uint index = i + num_leading_zeroes - 32;

				const ulong prevSumN = sumN;
				const uint p1 = primeReciprocals[index].y >> 32;
				ulong q = N * primeInverses[index].x;
				do
				{
					sumN = sumN * p1 + prevSumN;
					N = q;
					q *= primeInverses[index].x;
				} while (q <= primeInverses[index].y);
			} while (pattern);

			const ulong p2 = primeReciprocals[i].y >> 32;
			if (p2 * p2 > N)
			{
				if (N > 1)
				{
					sumN *= N + 1;
				}
				if (sumN == targetSum)
				{
					NumberFound(amicable_numbers_count, amicable_numbers_data, M, 0);
				}
				return;
			}

			const ulong oldTargetSum = targetSum;
			targetSum /= sumN;
			if (targetSum * sumN != oldTargetSum)
			{
				return;
			}

			if (N >= targetSum)
			{
				return;
			}

			// No more than 2 factors remain, we need to exit the loop now, otherwise MaximumSumOfDivisors3 will fail
			break;
		}

		const ulong2 curReciprocal = primeReciprocals[i];
		const ulong p1 = curReciprocal.y >> 32;
		if (p1 > root3_N_rounded_up)
		{
			break;
		}
		const ulong q = mul_hi(N + ((curReciprocal.y >> 8) & 1), curReciprocal.x) >> (curReciprocal.y & 63);

		ulong p2 = p1;
		p2 *= p1 + 1;
		if (((q + p2) << 1) <= MaximumSumOfDivisors3_target)
		{
			return;
		}

		if (i >= max_i)
		{
			const uint filtered_numbers_data_id = atom_inc(filtered_numbers_count);
			filtered_numbers_data[filtered_numbers_data_id].x = M;
			filtered_numbers_data[filtered_numbers_data_id].y = N;
			filtered_numbers_data[filtered_numbers_data_id].z = targetSum;
			filtered_numbers_data[filtered_numbers_data_id].w = max_i;
			return;
		}
	}

	// No more than 2 factors remain here

	// 1) The case when N is prime
	const ulong minSumN = N + 1;
	if (minSumN >= targetSum)
	{
		if (minSumN == targetSum)
		{
			// M is an amicable number if N is prime
			// We don't check for primality here, let CPU do it
			NumberFound(amicable_numbers_count, amicable_numbers_data, M, N);
		}
		return;
	}

	// 2) The case when N is squared prime
	const ulong sqrt_N = IntegerSquareRoot(N);
	if ((sqrt_N * sqrt_N == N) && (N + sqrt_N + 1 == targetSum))
	{
		NumberFound(amicable_numbers_count, amicable_numbers_data, M, sqrt_N);
		return;
	}

	// 3) The case when N = p * q - a product of two different primes
	// p + q == targetSum - N - 1
	// p * q == N
	const ulong B = targetSum - N - 1;
	ulong2 B2;
	B2.x = B * B;
	B2.y = mul_hi(B, B);
	const ulong carry = (B2.x < (N << 2)) ? 1 : 0;
	B2.x -= (N << 2);
	if (B2.y >= (N >> 62) + carry)
	{
		B2.y -= (N >> 62) + carry;
		const ulong sqrt_D = IntegerSquareRoot128(B2);
		const ulong p = (B - sqrt_D) >> 1;
		const ulong q = (B + sqrt_D) >> 1;
		if ((p < q) && (p + q == B) && (p * q == N))
		{
			NumberFound(amicable_numbers_count, amicable_numbers_data, M, q);
			return;
		}
	}
}

__kernel
__attribute__((reqd_work_group_size(1, 1, 1)))
void SaveCounter(__global uint* phase2_numbers_count)
{
	if (phase2_numbers_count[0] < PHASE2_MAX_COUNT)
	{
		phase2_numbers_count[1] = phase2_numbers_count[0];
	}
}

// primes is an array of pointers to chunks of data
// each chunk is 128 MB in size (2^24 uint2 elements) or bigger
ulong GetNthPrime(uint n,
	__global uint2* primes0
#if NUM_PRIME_CHUNKS > 1
	, __global uint2* primes1
#endif
#if NUM_PRIME_CHUNKS > 2
	, __global uint2* primes2
#endif
#if NUM_PRIME_CHUNKS > 3
	, __global uint2* primes3
#endif
#if NUM_PRIME_CHUNKS > 4
	, __global uint2* primes4
#endif
#if NUM_PRIME_CHUNKS > 5
	, __global uint2* primes5
#endif
#if NUM_PRIME_CHUNKS > 6
	, __global uint2* primes6
#endif
#if NUM_PRIME_CHUNKS > 7
	, __global uint2* primes7
#endif
)
{
#if NUM_PRIME_CHUNKS == 1
	uint2 data = primes0[n >> 2];
#else
	const uint global_offset = n >> 2;
	__global const uint2* chunk;
	switch (global_offset >> CHUNK_SIZE_SHIFT)
	{
	default: chunk = primes0; break;
#if NUM_PRIME_CHUNKS > 1
	case 1: chunk = primes1; break;
#endif
#if NUM_PRIME_CHUNKS > 2
	case 2: chunk = primes2; break;
#endif
#if NUM_PRIME_CHUNKS > 3
	case 3: chunk = primes3; break;
#endif
#if NUM_PRIME_CHUNKS > 4
	case 4: chunk = primes4; break;
#endif
#if NUM_PRIME_CHUNKS > 5
	case 5: chunk = primes5; break;
#endif
#if NUM_PRIME_CHUNKS > 6
	case 6: chunk = primes6; break;
#endif
#if NUM_PRIME_CHUNKS > 7
	case 7: chunk = primes7; break;
#endif
	}

	const uint chunk_offset = global_offset & ((1 << CHUNK_SIZE_SHIFT) - 1);
	uint2 data = chunk[chunk_offset];
#endif

	data.x += (data.y >> (30 - (n & 3) * 10)) & 1023;
	const ulong result = data.x;
	return result + result + 1;
}

__kernel
__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
void SearchMultipleRanges(
	__global const uint* smallPrimes,
	__global const ulong2* primeInverses,
	__global const ulong2* PQ,
	__constant const ulong2* PowerOf2SumInverses,
	__constant const ulong4* PowersOfP_128SumInverses,
	__constant const ulong2* PrimeInverses_128,
	__global const ulong4* RangesTable,
	__global const uint2* RangeLookupTable,
	const uint lookup_shift,
	const uint global_offset,
	const uint global_size,
	__global uint* phase1_offset_to_resume_after_overflow,
	__global uint* phase2_numbers_count,
	__global ulong4* phase2_numbers_data,
	__global uint* amicable_numbers_count,
	__global ulong2* amicable_numbers_data,

	__global uint2* primes0
#if NUM_PRIME_CHUNKS > 1
	, __global uint2* primes1
#endif
#if NUM_PRIME_CHUNKS > 2
	, __global uint2* primes2
#endif
#if NUM_PRIME_CHUNKS > 3
	, __global uint2* primes3
#endif
#if NUM_PRIME_CHUNKS > 4
	, __global uint2* primes4
#endif
#if NUM_PRIME_CHUNKS > 5
	, __global uint2* primes5
#endif
#if NUM_PRIME_CHUNKS > 6
	, __global uint2* primes6
#endif
#if NUM_PRIME_CHUNKS > 7
	, __global uint2* primes7
#endif
)
{
	if (*phase2_numbers_count >= PHASE2_MAX_COUNT)
	{
		return;
	}

	const uint globalIndex = get_global_id(0) + global_offset;
	if (globalIndex >= global_size)
	{
		return;
	}

	ulong M0, sumM0, p;

	uint2 rangeIndexAndOffset = RangeLookupTable[globalIndex >> lookup_shift];
	uint curIndex = globalIndex & (0xFFFFFFFFU << lookup_shift);
	for (;;)
	{
		const uint numbersRemainingInCurRange = RangesTable[rangeIndexAndOffset.x].w - rangeIndexAndOffset.y;

		// curIndex + numbersRemainingInCurRange now points exactly at the beginning of the next range
		// If the beginning of the next range is > globalIndex, then globalIndex is within the current range
		if (curIndex + numbersRemainingInCurRange > globalIndex)
		{
			M0 = RangesTable[rangeIndexAndOffset.x].x;
			sumM0 = RangesTable[rangeIndexAndOffset.x].y;

			//      primes[ start_prime_index                    + range_offset          + globalIndex - curIndex]
			p = GetNthPrime(RangesTable[rangeIndexAndOffset.x].z + rangeIndexAndOffset.y + globalIndex - curIndex,
				primes0
#if NUM_PRIME_CHUNKS > 1
				, primes1
#endif
#if NUM_PRIME_CHUNKS > 2
				, primes2
#endif
#if NUM_PRIME_CHUNKS > 3
				, primes3
#endif
#if NUM_PRIME_CHUNKS > 4
				, primes4
#endif
#if NUM_PRIME_CHUNKS > 5
				, primes5
#endif
#if NUM_PRIME_CHUNKS > 6
				, primes6
#endif
#if NUM_PRIME_CHUNKS > 7
				, primes7
#endif
			);
			break;
		}

		// Move curIndex to the beginning of the next range
		curIndex += numbersRemainingInCurRange;
		++rangeIndexAndOffset.x;
		rangeIndexAndOffset.y = 0;
	}

	CheckPairPhase1(
		smallPrimes,
		primeInverses,
		PQ,
		PowerOf2SumInverses,
		PowersOfP_128SumInverses,
		PrimeInverses_128,
#if LPP == 1
		M0 * p, sumM0 * (p + 1), mul_hi(sumM0, p + 1),
#elif LPP == 2
		M0 * p * p, sumM0 * (p * p + p + 1), mul_hi(sumM0, p * p + p + 1),
#elif LPP == 3
		M0 * p * p * p, sumM0 * (p * p * p + p * p + p + 1), mul_hi(sumM0, p * p + p + 1),
#endif
		global_offset,
		phase1_offset_to_resume_after_overflow,
		phase2_numbers_count,
		phase2_numbers_data,
		amicable_numbers_count,
		amicable_numbers_data
	);
}

__kernel
__attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
void CheckPairs(
	__global const uint* smallPrimes,
	__global const ulong2* primeInverses,
	__global const ulong2* PQ,
	__constant const ulong2* PowerOf2SumInverses,
	__constant const ulong4* PowersOfP_128SumInverses,
	__constant const ulong2* PrimeInverses_128,
	__global const ulong4* pairsToCheck,
	__global uint* filtered_numbers_count,
	__global ulong4* filtered_numbers_data,
	__global uint* amicable_numbers_count,
	__global ulong2* amicable_numbers_data
)
{
	const ulong4 cur_pair = pairsToCheck[get_global_id(0)];
	if (cur_pair.x == 0)
	{
		return;
	}
	CheckPairPhase1(smallPrimes, primeInverses, PQ, PowerOf2SumInverses, PowersOfP_128SumInverses, PrimeInverses_128, cur_pair.x, cur_pair.y, cur_pair.z, 0, 0, filtered_numbers_count, filtered_numbers_data, amicable_numbers_count, amicable_numbers_data);
}
