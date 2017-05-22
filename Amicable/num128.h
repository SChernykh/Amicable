#pragma once

#ifndef __GNUG__

struct num128
{
	// Default constructor is empty for performance reasons
	FORCEINLINE num128() {}

	FORCEINLINE num128(const num128& a) : lo(a.lo), hi(a.hi) {}

	FORCEINLINE num128(int a) : lo(static_cast<unsigned int>(a)), hi(0) {}
	FORCEINLINE num128(num64 a) : lo(a), hi(0) {}

	FORCEINLINE num128& operator=(num64 a)
	{
		lo = a;
		hi = 0;
		return *this;
	}

	FORCEINLINE num128& operator=(const num128& a)
	{
		lo = a.lo;
		hi = a.hi;
		return *this;
	}

	FORCEINLINE num128 operator+(num64 a) const
	{
		num128 result;
		add128(lo, hi, a, 0, &result.lo, &result.hi);
		return result;
	}

	FORCEINLINE num128& operator+=(num64 a)
	{
		add128(lo, hi, a, 0, &lo, &hi);
		return *this;
	}

	FORCEINLINE num128 operator-() const
	{
		num128 result;
		add128(~lo, ~hi, 1, 0, &result.lo, &result.hi);
		return result;
	}

	FORCEINLINE num128 operator-(num64 a) const
	{
		num128 result;
		sub128(lo, hi, a, 0, &result.lo, &result.hi);
		return result;
	}

	FORCEINLINE num128& operator-=(num64 a)
	{
		sub128(lo, hi, a, 0, &lo, &hi);
		return *this;
	}

	FORCEINLINE num128 operator*(num64 a) const
	{
		num128 result;
		result.lo = _umul128(lo, a, &result.hi);
		result.hi += hi * a;
		return result;
	}

	FORCEINLINE num128& operator*=(num64 a)
	{
		const num64 old_hi = hi;
		lo = _umul128(lo, a, &hi);
		hi += old_hi * a;
		return *this;
	}

	FORCEINLINE num128 operator/(num64 a) const
	{
		const num128 t = a;
		return *this / t;
	}

	FORCEINLINE num128& operator/=(num64 a)
	{
		const num128 t = a;
		return operator/=(t);
	}

	FORCEINLINE bool operator==(num64 a) const { return (lo == a) && (hi == 0); }
	FORCEINLINE bool operator!=(num64 a) const { return (lo != a) || (hi != 0); }
	FORCEINLINE byte operator<=(num64 a) const { return leq128(lo, hi, a, 0); }
	FORCEINLINE byte operator>=(num64 a) const { return leq128(a, 0, lo, hi); }
	FORCEINLINE byte operator<(num64 a) const { return less128(lo, hi, a, 0); }
	FORCEINLINE byte operator>(num64 a) const { return less128(a, 0, lo, hi); }

	FORCEINLINE num128 operator+(const num128& a) const
	{
		num128 result;
		add128(lo, hi, a.lo, a.hi, &result.lo, &result.hi);
		return result;
	}

	FORCEINLINE num128& operator+=(const num128& a)
	{
		add128(lo, hi, a.lo, a.hi, &lo, &hi);
		return *this;
	}

	FORCEINLINE num128 operator-(const num128& a) const
	{
		num128 result;
		sub128(lo, hi, a.lo, a.hi, &result.lo, &result.hi);
		return result;
	}

	FORCEINLINE num128& operator-=(const num128& a)
	{
		sub128(lo, hi, a.lo, a.hi, &lo, &hi);
		return *this;
	}

	FORCEINLINE num128 operator<<(unsigned int count) const
	{
		num128 result = *this;
		shl128(result.lo, result.hi, static_cast<unsigned char>(count));
		return result;
	}

	FORCEINLINE num128& operator<<=(unsigned int count)
	{
		shl128(lo, hi, static_cast<unsigned char>(count));
		return *this;
	}

	FORCEINLINE num128 operator>>(unsigned int count) const
	{
		num128 result = *this;
		shr128(result.lo, result.hi, static_cast<unsigned char>(count));
		return result;
	}

	FORCEINLINE num128& operator>>=(unsigned int count)
	{
		shr128(lo, hi, static_cast<unsigned char>(count));
		return *this;
	}

	FORCEINLINE num128 operator*(const num128& a) const
	{
		num128 result;
		result.lo = _umul128(lo, a.lo, &result.hi);
		result.hi += lo * a.hi + hi * a.lo;
		return result;
	}

	FORCEINLINE num128& operator*=(const num128& a)
	{
		num128 t;
		t.lo = _umul128(lo, a.lo, &t.hi);
		t.hi += lo * a.hi + hi * a.lo;
		*this = t;
		return *this;
	}

	NOINLINE num128 operator/(const num128& a) const;

	FORCEINLINE num128 operator%(const num128& a) const
	{
		return *this - a * operator/(a);
	}

	FORCEINLINE num128& operator/=(const num128& a)
	{
		*this = operator/(a);
		return *this;
	}

	FORCEINLINE num128& operator%=(const num128& a)
	{
		return operator-=(a * operator/(a));
	}

	FORCEINLINE bool operator==(const num128& a) const { return (lo == a.lo) && (hi == a.hi); }
	FORCEINLINE bool operator!=(const num128& a) const { return (lo != a.lo) || (hi != a.hi); }
	FORCEINLINE byte operator<=(const num128& a) const { return leq128(lo, hi, a.lo, a.hi); }
	FORCEINLINE byte operator>=(const num128& a) const { return leq128(a.lo, a.hi, lo, hi); }
	FORCEINLINE byte operator<(const num128& a) const { return less128(lo, hi, a.lo, a.hi); }
	FORCEINLINE byte operator>(const num128& a) const { return less128(a.lo, a.hi, lo, hi); }

	num64 lo;
	num64 hi;
};

#endif

num128 atoi128(const char* s);
char* itoa128(num128 a, char* buf, size_t bufSize);
std::ostream& operator<<(std::ostream& s, num128 a);
num64 IntegerSquareRoot(const num128& n);

extern const num128 NUM128_MAX;
