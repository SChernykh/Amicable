#pragma once

template<number Pred, number A, number B> struct Predicate;

template<number A, number B>
struct Predicate<1, A, B>
{
	enum Value : number
	{
		value = A,
	};
};

template<number A, number B>
struct Predicate<0, A, B>
{
	enum Value : number
	{
		value = B,
	};
};

template<number N, number A, number B>
struct CompileTimeSQRTImpl
{
	enum Value : number
	{
		C = (A / 2) + (B / 2) + ((A & B) & 1) + 1,
		Pred = (C * C > N) ? 1 : 0,
		value = CompileTimeSQRTImpl<N, Predicate<Pred, A, C>::value, Predicate<Pred, C - 1, B>::value>::value,
	};
};

template<number N, number C>
struct CompileTimeSQRTImpl<N, C, C>
{
	enum Value : number
	{
		value = C,
	};
};

template<number N>
struct CompileTimeSQRT
{
	enum Value : number
	{
		value = CompileTimeSQRTImpl<N, 1, 4294967295>::value,
	};
};

template<number A, number B>
struct Min
{
	enum X : number
	{
		value = (A < B) ? A : B,
	};
};

template<number A, number B>
struct Max
{
	enum X : number
	{
		value = (A > B) ? A : B,
	};
};

#pragma warning(push)

// Disable this warning here because number overflows are the essence of this algorithm: it computes everything modulo 2^64.
#pragma warning(disable : 4307) // warning C4307: '*' : integral constant overflow

template<number x1, number x2, number v1, number v2>
struct ExtendedEuclidCompileTime
{
	enum Variables : number
	{
		q = v1 / v2,
		x3 = x1 - q * x2,
		v3 = v1 - q * v2,
		value = ExtendedEuclidCompileTime<x2, x3, v2, v3>::value,
	};
};

template<number x1, number x2, number v1>
struct ExtendedEuclidCompileTime<x1, x2, v1, 1>
{
	enum Variables : number
	{
		value = x2,
	};
};

template<number N>
struct MultiplicativeInverse
{
	enum Variables : number
	{
		value = ExtendedEuclidCompileTime<number(-1), 1, ~N + 1, N>::value,
	};
};

// For even numbers
template<number N, number K, bool IsEven>
struct MultiplicativeInverseEvenImpl
{
	enum Variables : number
	{
		value = MultiplicativeInverseEvenImpl<N / 2, K + 1, ((N / 2) % 2) == 0>::value,
		shift = MultiplicativeInverseEvenImpl<N / 2, K + 1, ((N / 2) % 2) == 0>::shift,
	};
};

template<number N, number K>
struct MultiplicativeInverseEvenImpl<N, K, false>
{
	enum Variables : number
	{
		value = MultiplicativeInverse<N>::value,
		shift = K,
	};
};

template<number N>
struct MultiplicativeInverseEven
{
	enum Variables : number
	{
		value = MultiplicativeInverseEvenImpl<N, 0, (N % 2) == 0>::value,
		shift = MultiplicativeInverseEvenImpl<N, 0, (N % 2) == 0>::shift,
	};
};

#pragma warning(pop)

template<number N> struct CompileTimePrimes;

template<> struct CompileTimePrimes<0>  { enum { value =  2 }; };
template<> struct CompileTimePrimes<1>  { enum { value =  3 }; };
template<> struct CompileTimePrimes<2>  { enum { value =  5 }; };
template<> struct CompileTimePrimes<3>  { enum { value =  7 }; };
template<> struct CompileTimePrimes<4>  { enum { value = 11 }; };
template<> struct CompileTimePrimes<5>  { enum { value = 13 }; };
template<> struct CompileTimePrimes<6>  { enum { value = 17 }; };
template<> struct CompileTimePrimes<7>  { enum { value = 19 }; };
template<> struct CompileTimePrimes<8>  { enum { value = 23 }; };
template<> struct CompileTimePrimes<9>  { enum { value = 29 }; };
template<> struct CompileTimePrimes<10> { enum { value = 31 }; };
template<> struct CompileTimePrimes<11> { enum { value = 37 }; };
template<> struct CompileTimePrimes<12> { enum { value = 41 }; };
template<> struct CompileTimePrimes<13> { enum { value = 43 }; };
template<> struct CompileTimePrimes<14> { enum { value = 47 }; };
template<> struct CompileTimePrimes<15> { enum { value = 53 }; };
template<> struct CompileTimePrimes<16> { enum { value = 59 }; };
template<> struct CompileTimePrimes<17> { enum { value = 61 }; };
template<> struct CompileTimePrimes<18> { enum { value = 67 }; };
template<> struct CompileTimePrimes<19> { enum { value = 71 }; };
template<> struct CompileTimePrimes<20> { enum { value = 73 }; };
template<> struct CompileTimePrimes<21> { enum { value = 79 }; };
template<> struct CompileTimePrimes<22> { enum { value = 83 }; };
template<> struct CompileTimePrimes<23> { enum { value = 89 }; };
template<> struct CompileTimePrimes<24> { enum { value = 97 }; };
template<> struct CompileTimePrimes<25> { enum { value = 101 }; };
template<> struct CompileTimePrimes<26> { enum { value = 103 }; };
template<> struct CompileTimePrimes<27> { enum { value = 107 }; };
template<> struct CompileTimePrimes<28> { enum { value = 109 }; };
template<> struct CompileTimePrimes<29> { enum { value = 113 }; };
template<> struct CompileTimePrimes<30> { enum { value = 127 }; };
template<> struct CompileTimePrimes<31> { enum { value = 131 }; };
template<> struct CompileTimePrimes<32> { enum { value = 137 }; };
template<> struct CompileTimePrimes<33> { enum { value = 139 }; };
template<> struct CompileTimePrimes<34> { enum { value = 149 }; };
template<> struct CompileTimePrimes<35> { enum { value = 151 }; };
template<> struct CompileTimePrimes<36> { enum { value = 157 }; };
template<> struct CompileTimePrimes<37> { enum { value = 163 }; };
template<> struct CompileTimePrimes<38> { enum { value = 167 }; };
template<> struct CompileTimePrimes<39> { enum { value = 173 }; };
template<> struct CompileTimePrimes<40> { enum { value = 179 }; };
template<> struct CompileTimePrimes<41> { enum { value = 181 }; };
template<> struct CompileTimePrimes<42> { enum { value = 191 }; };
template<> struct CompileTimePrimes<43> { enum { value = 193 }; };
template<> struct CompileTimePrimes<44> { enum { value = 197 }; };
template<> struct CompileTimePrimes<45> { enum { value = 199 }; };
template<> struct CompileTimePrimes<46> { enum { value = 211 }; };
template<> struct CompileTimePrimes<47> { enum { value = 223 }; };
template<> struct CompileTimePrimes<48> { enum { value = 227 }; };
template<> struct CompileTimePrimes<49> { enum { value = 229 }; };
template<> struct CompileTimePrimes<50> { enum { value = 233 }; };
template<> struct CompileTimePrimes<51> { enum { value = 239 }; };
template<> struct CompileTimePrimes<52> { enum { value = 241 }; };
template<> struct CompileTimePrimes<53> { enum { value = 251 }; };
template<> struct CompileTimePrimes<54> { enum { value = 257 }; };
template<> struct CompileTimePrimes<55> { enum { value = 263 }; };
template<> struct CompileTimePrimes<56> { enum { value = 269 }; };
template<> struct CompileTimePrimes<57> { enum { value = 271 }; };
template<> struct CompileTimePrimes<58> { enum { value = 277 }; };
template<> struct CompileTimePrimes<59> { enum { value = 281 }; };
template<> struct CompileTimePrimes<60> { enum { value = 283 }; };
template<> struct CompileTimePrimes<61> { enum { value = 293 }; };
template<> struct CompileTimePrimes<62> { enum { value = 307 }; };
template<> struct CompileTimePrimes<63> { enum { value = 311 }; };
template<> struct CompileTimePrimes<64> { enum { value = 313 }; };
template<> struct CompileTimePrimes<65> { enum { value = 317 }; };
template<> struct CompileTimePrimes<66> { enum { value = 331 }; };
template<> struct CompileTimePrimes<67> { enum { value = 337 }; };
template<> struct CompileTimePrimes<68> { enum { value = 347 }; };
template<> struct CompileTimePrimes<69> { enum { value = 349 }; };
template<> struct CompileTimePrimes<70> { enum { value = 353 }; };
template<> struct CompileTimePrimes<71> { enum { value = 359 }; };
template<> struct CompileTimePrimes<72> { enum { value = 367 }; };
template<> struct CompileTimePrimes<73> { enum { value = 373 }; };
template<> struct CompileTimePrimes<74> { enum { value = 379 }; };
template<> struct CompileTimePrimes<75> { enum { value = 383 }; };
template<> struct CompileTimePrimes<76> { enum { value = 389 }; };
template<> struct CompileTimePrimes<77> { enum { value = 397 }; };
template<> struct CompileTimePrimes<78> { enum { value = 401 }; };
template<> struct CompileTimePrimes<79> { enum { value = 409 }; };
template<> struct CompileTimePrimes<80> { enum { value = 419 }; };
template<> struct CompileTimePrimes<81> { enum { value = 421 }; };
template<> struct CompileTimePrimes<82> { enum { value = 431 }; };
template<> struct CompileTimePrimes<83> { enum { value = 433 }; };
template<> struct CompileTimePrimes<84> { enum { value = 439 }; };
template<> struct CompileTimePrimes<85> { enum { value = 443 }; };
template<> struct CompileTimePrimes<86> { enum { value = 449 }; };
template<> struct CompileTimePrimes<87> { enum { value = 457 }; };
template<> struct CompileTimePrimes<88> { enum { value = 461 }; };
template<> struct CompileTimePrimes<89> { enum { value = 463 }; };
template<> struct CompileTimePrimes<90> { enum { value = 467 }; };
template<> struct CompileTimePrimes<91> { enum { value = 479 }; };
template<> struct CompileTimePrimes<92> { enum { value = 487 }; };
template<> struct CompileTimePrimes<93> { enum { value = 491 }; };
template<> struct CompileTimePrimes<94> { enum { value = 499 }; };
template<> struct CompileTimePrimes<95> { enum { value = 503 }; };
template<> struct CompileTimePrimes<96> { enum { value = 509 }; };
template<> struct CompileTimePrimes<97> { enum { value = 521 }; };

enum { CompileTimePrimesCount = 96 };
