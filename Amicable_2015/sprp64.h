#pragma once

#pragma warning(push, 1)
#pragma warning(disable : 4146 4838)

#if _MSC_VER >= 1900
static FORCEINLINE number mont_prod64(const number a, const number b, const number n, const number npi)
{
	number t_hi, mn_hi, u;

	// t_hi * 2^64 + t_lo = a*b
	number t_lo = _umul128(a, b, &t_hi);

	const number m = t_lo * npi;

	// mn_hi * 2^64 + mn_lo = m*n
	const number mn_lo = _umul128(m, n, &mn_hi);

	const unsigned char carry = _addcarry_u64(0, t_lo, mn_lo, &t_lo);
	_addcarry_u64(carry, t_hi, mn_hi, &u);

	if (u < t_hi) return u - n;

	return u >= n ? u - n : u;
}

static FORCEINLINE number mont_prod_add64(const number a, const number b, const number c, const number n, const number npi)
{
	number t_hi, mn_hi, u;

	// t_hi * 2^64 + t_lo = a*b
	number t_lo = _umul128(a, b, &t_hi);

	// t_hi * 2^64 + t_lo = a*b+c
	const unsigned char carry1 = _addcarry_u64(0, t_lo, c, &t_lo);
	_addcarry_u64(carry1, t_hi, 0, &t_hi);

	const number m = t_lo * npi;

	// mn_hi * 2^64 + mn_lo = m*n
	const number mn_lo = _umul128(m, n, &mn_hi);

	const unsigned char carry2 = _addcarry_u64(0, t_lo, mn_lo, &t_lo);
	_addcarry_u64(carry2, t_hi, mn_hi, &u);

	if (u < t_hi) return u - n;

	return u >= n ? u - n : u;
}

#else

static FORCEINLINE number mont_prod64(const number a, const number b, const number n, const number npi)
{
	number t_hi, t_lo, m, mn_hi, mn_lo, u;
	int carry;

	// t_hi * 2^64 + t_lo = a*b
	t_lo = _umul128(a, b, &t_hi);

	m = t_lo * npi;

	// mn_hi * 2^64 + mn_lo = m*n
	mn_lo = _umul128(m, n, &mn_hi);

	carry = t_lo + mn_lo < t_lo ? 1 : 0;

	u = t_hi + mn_hi + carry;
	if (u < t_hi) return u - n;

	return u >= n ? u - n : u;
}

static FORCEINLINE number mont_prod_add64(const number a, const number b, const number c, const number n, const number npi)
{
	number t_hi, t_lo, m, mn_hi, mn_lo, u;

	// t_hi * 2^64 + t_lo = a*b
	t_lo = _umul128(a, b, &t_hi);

	// t_hi * 2^64 + t_lo = a*b+c
	int carry = (t_lo + c < t_lo) ? 1 : 0;
	t_lo += c;
	t_hi += carry;

	m = t_lo * npi;

	// mn_hi * 2^64 + mn_lo = m*n
	mn_lo = _umul128(m, n, &mn_hi);

	carry = (t_lo + mn_lo < t_lo) ? 1 : 0;

	u = t_hi + mn_hi + carry;
	if (u < t_hi) return u - n;

	return u >= n ? u - n : u;
}
#endif

// WARNING: a must be odd
// returns -a^-1 mod 2^64
// method from B. Arazi 'On Primality Testing Using Purely Divisionless Operations'
// The Computer Journal (1994) 37 (3): 219-222, Procedure 5.
// modified to process 4 or 8 bits at a time
static FORCEINLINE number modular_inverse64(const number a)
{
	number S = 1;

	static const char mask[128] = {255,85,51,73,199,93,59,17,15,229,195,89,215,237,203,33,31,117,83,105,231,125,91,49,47,5,227,121,247,13,235,65,63,149,115,137,7,157,123,81,79,37,3,153,23,45,11,97,95,181,147,169,39,189,155,113,111,69,35,185,55,77,43,129,127,213,179,201,71,221,187,145,143,101,67,217,87,109,75,161,159,245,211,233,103,253,219,177,175,133,99,249,119,141,107,193,191,21,243,9,135,29,251,209,207,165,131,25,151,173,139,225,223,53,19,41,167,61,27,241,239,197,163,57,183,205,171,1};

	const char amask = mask[(a >> 1) & 127];

	int index = (amask * (S & 255)) & 255;
	number J = index;
	S = (S + a * index) >> 8;

	index = (amask * (S & 255)) & 255;
	J |= (number)index << 8;
	S = (S + a * index) >> 8;

	index = (amask * (S & 255)) & 255;
	J |= (number)index << 16;
	S = (S + a * index) >> 8;

	index = (amask * (S & 255)) & 255;
	J |= (number)index << 24;
	unsigned int S2 = (S + a * index) >> 8;

	index = (amask * (S2 & 255)) & 255;
	J |= (number)index << 32;
	S2 = (S2 + a * index) >> 8;

	index = (amask * (S2 & 255)) & 255;
	J |= (number)index << 40;
	S2 = (S2 + a * index) >> 8;

	index = (amask * (S2 & 255)) & 255;
	J |= (number)index << 48;
	S2 = (S2 + a * index) >> 8;

	index = (amask * (S2 & 255)) & 255;
	J |= (number)index << 56;

	return J;
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
	if ((N & 1) == 0)
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
	static const number bases[][2] = {{2, 0}, {325, 2047}, {9375, 4033}, {28178, 284301751}, {450775, 3874471147}, {9780504, 3874471147}, {1795265022, 107528788110061}};
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
