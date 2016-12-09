#pragma once

#include "sprp_bases.inl"

#pragma warning(push, 1)
#pragma warning(disable : 4146 4838)

static FORCEINLINE number mont_prod64(const number a, const number b, const number n, const number npi)
{
	number t_hi, mn_hi, u;

	// t_hi * 2^64 + t_lo = a*b
	number t_lo = _umul128(a, b, &t_hi);

	const number m = t_lo * npi;

	// mn_hi * 2^64 + mn_lo = m*n
	const number mn_lo = _umul128(m, n, &mn_hi);

	add128(t_lo, t_hi, mn_lo, mn_hi, &t_lo, &u);

	number n1 = 0;
	if (u < t_hi) n1 = n;
	if (u >= n) n1 = n;
	return u - n1;
}

static FORCEINLINE number mont_prod_add64(const number a, const number b, const number c, const number n, const number npi)
{
	number t_hi, mn_hi, u;

	// t_hi * 2^64 + t_lo = a*b
	number t_lo = _umul128(a, b, &t_hi);

	// t_hi * 2^64 + t_lo = a*b+c
	add128(t_lo, t_hi, c, 0, &t_lo, &t_hi);

	const number m = t_lo * npi;

	// mn_hi * 2^64 + mn_lo = m*n
	const number mn_lo = _umul128(m, n, &mn_hi);

	add128(t_lo, t_hi, mn_lo, mn_hi, &t_lo, &u);

	number n1 = 0;
	if (u < t_hi) n1 = n;
	if (u >= n) n1 = n;
	return u - n1;
}

// WARNING: a must be odd
// returns -a^-1 mod 2^64
// method from B. Arazi 'On Primality Testing Using Purely Divisionless Operations'
// The Computer Journal (1994) 37 (3): 219-222, Procedure 5.
// modified to process 4 or 8 bits at a time
static FORCEINLINE number modular_inverse64(const number a)
{
	static const unsigned char mask[128] = { 255,85,51,73,199,93,59,17,15,229,195,89,215,237,203,33,31,117,83,105,231,125,91,49,47,5,227,121,247,13,235,65,63,149,115,137,7,157,123,81,79,37,3,153,23,45,11,97,95,181,147,169,39,189,155,113,111,69,35,185,55,77,43,129,127,213,179,201,71,221,187,145,143,101,67,217,87,109,75,161,159,245,211,233,103,253,219,177,175,133,99,249,119,141,107,193,191,21,243,9,135,29,251,209,207,165,131,25,151,173,139,225,223,53,19,41,167,61,27,241,239,197,163,57,183,205,171,1 };

	const unsigned char amask = mask[static_cast<unsigned char>(a) >> 1];

	union
	{
		unsigned char resultBytes[8];
		number result;
	};

	number S = 1;
	for (number i = 0; i < 8; ++i)
	{
		resultBytes[i] = amask * S;
		S = (S + a * resultBytes[i]) >> 8;
	}
	return result;
}

static FORCEINLINE void modular_inverse64_2(const number a1, const number a2, number* results)
{
	static const unsigned char mask[128] = { 255,85,51,73,199,93,59,17,15,229,195,89,215,237,203,33,31,117,83,105,231,125,91,49,47,5,227,121,247,13,235,65,63,149,115,137,7,157,123,81,79,37,3,153,23,45,11,97,95,181,147,169,39,189,155,113,111,69,35,185,55,77,43,129,127,213,179,201,71,221,187,145,143,101,67,217,87,109,75,161,159,245,211,233,103,253,219,177,175,133,99,249,119,141,107,193,191,21,243,9,135,29,251,209,207,165,131,25,151,173,139,225,223,53,19,41,167,61,27,241,239,197,163,57,183,205,171,1 };

	const unsigned char amask1 = mask[static_cast<unsigned char>(a1) >> 1];
	const unsigned char amask2 = mask[static_cast<unsigned char>(a2) >> 1];

	number S1 = 1;
	number S2 = 1;

	for (number i = 0; i < 8; ++i)
	{
		reinterpret_cast<unsigned char*>(results + 0)[i] = amask1 * S1;
		reinterpret_cast<unsigned char*>(results + 1)[i] = amask2 * S2;
		S1 = (S1 + a1 * reinterpret_cast<unsigned char*>(results + 0)[i]) >> 8;
		S2 = (S2 + a2 * reinterpret_cast<unsigned char*>(results + 1)[i]) >> 8;
	}
}

// returns 2^64 mod n
static FORCEINLINE number compute_modn64(const number n)
{
	if (n <= (1ULL << 63)) {
		number res = ((1ULL << 63) % n) << 1;

		return res < n ? res : res-n;
	} else
		return -n;
}

static FORCEINLINE number compute_a_times_2_64_mod_n(const number a, const number n, const number r)
{
	return mulmod64(a, r, n);
}

static FORCEINLINE bool efficient_mr64(const number N)
{
	if (bitOffset[N % PrimeTableParameters::Modulo] >= PrimeTableParameters::NumOffsets)
	{
		return false;
	}

	const number npi = modular_inverse64(N);
	const number r = compute_modn64(N);

	number u=N-1;
	const number nr = N-r;	

	int t;
	_BitScanForward64((DWORD*) &t, u);
	u >>= t;

	// These 7 bases are enough to make Miller-Rabin test deterministic for N < 2^64
	//static const number bases[][2] = {{2, 0}, {325, 2047}, {9375, 4033}, {28178, 284301751}, {450775, 3874471147}, {9780504, 3874471147}, {1795265022, 107528788110061}};

	number bases[3][2] = { {2, 0}, {0, 0}, {0, 0} };
	{
		number h = N;
		h = ((h >> 32) ^ h) * 0x45d9f3b3335b369;
		h = ((h >> 32) ^ h) * 0x3335b36945d9f3b;
		h = ((h >> 32) ^ h);
		const unsigned int b = sprp_bases[h & 16383];
		bases[1][0] = b & 4095;
		bases[2][0] = b >> 12;
	}

	for (const number (&a)[2] : bases)
	{
		if (N < a[1])
		{
			return true;
		}

		number A = compute_a_times_2_64_mod_n(a[0], N, r); // a * 2^64 mod n

		if (A)
		{
			number d = r, u_copy = u;

			// compute a^u mod n
			do
			{
				if (u_copy & 1) d = mont_prod64(d, A, N, npi);
				A = mont_prod64(A, A, N, npi);
			} while (u_copy >>= 1);

			if (d != r && d != nr)
			{
				int i;
				for (i = 1; i < t; i++)
				{
					d = mont_prod64(d, d, N, npi);
					if (d == r) return false;
					if (d == nr) break; // PRIME in subtest
				}

				if (i == t)
				{
					return false;
				}
			}
		}
	}

	return true;
}

#pragma warning(pop)
