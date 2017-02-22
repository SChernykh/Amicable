#pragma once

// Ugly setKernelArguments implementation for compilers without variadic templates

template<typename T0>
FORCEINLINE cl_int setKernelArguments(cl_kernel kernel, const T0& arg0)
{
	cl_int ciErrNum = clSetKernelArg(kernel, 0, sizeof(T0), &arg0);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	return CL_SUCCESS;
}

template<typename T0, typename T1>
FORCEINLINE cl_int setKernelArguments(cl_kernel kernel, const T0& arg0, const T1& arg1)
{
	cl_int ciErrNum = clSetKernelArg(kernel, 0, sizeof(T0), &arg0);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 1, sizeof(T1), &arg1);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	return CL_SUCCESS;
}

template<typename T0, typename T1, typename T2>
FORCEINLINE cl_int setKernelArguments(cl_kernel kernel, const T0& arg0, const T1& arg1, const T2& arg2)
{
	cl_int ciErrNum = clSetKernelArg(kernel, 0, sizeof(T0), &arg0);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 1, sizeof(T1), &arg1);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 2, sizeof(T2), &arg2);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	return CL_SUCCESS;
}

template<typename T0, typename T1, typename T2, typename T3>
FORCEINLINE cl_int setKernelArguments(cl_kernel kernel, const T0& arg0, const T1& arg1, const T2& arg2, const T3& arg3)
{
	cl_int ciErrNum = clSetKernelArg(kernel, 0, sizeof(T0), &arg0);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 1, sizeof(T1), &arg1);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 2, sizeof(T2), &arg2);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 3, sizeof(T3), &arg3);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	return CL_SUCCESS;
}

template<typename T0, typename T1, typename T2, typename T3, typename T4>
FORCEINLINE cl_int setKernelArguments(cl_kernel kernel, const T0& arg0, const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4)
{
	cl_int ciErrNum = clSetKernelArg(kernel, 0, sizeof(T0), &arg0);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 1, sizeof(T1), &arg1);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 2, sizeof(T2), &arg2);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 3, sizeof(T3), &arg3);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 4, sizeof(T4), &arg4);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	return CL_SUCCESS;
}

template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
FORCEINLINE cl_int setKernelArguments(cl_kernel kernel, const T0& arg0, const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5)
{
	cl_int ciErrNum = clSetKernelArg(kernel, 0, sizeof(T0), &arg0);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 1, sizeof(T1), &arg1);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 2, sizeof(T2), &arg2);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 3, sizeof(T3), &arg3);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 4, sizeof(T4), &arg4);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 5, sizeof(T5), &arg5);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	return CL_SUCCESS;
}

template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
FORCEINLINE cl_int setKernelArguments(cl_kernel kernel, const T0& arg0, const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6)
{
	cl_int ciErrNum = clSetKernelArg(kernel, 0, sizeof(T0), &arg0);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 1, sizeof(T1), &arg1);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 2, sizeof(T2), &arg2);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 3, sizeof(T3), &arg3);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 4, sizeof(T4), &arg4);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 5, sizeof(T5), &arg5);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 6, sizeof(T6), &arg6);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	return CL_SUCCESS;
}

template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
FORCEINLINE cl_int setKernelArguments(cl_kernel kernel, const T0& arg0, const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7)
{
	cl_int ciErrNum = clSetKernelArg(kernel, 0, sizeof(T0), &arg0);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 1, sizeof(T1), &arg1);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 2, sizeof(T2), &arg2);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 3, sizeof(T3), &arg3);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 4, sizeof(T4), &arg4);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 5, sizeof(T5), &arg5);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 6, sizeof(T6), &arg6);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 7, sizeof(T7), &arg7);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	return CL_SUCCESS;
}

template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
FORCEINLINE cl_int setKernelArguments(cl_kernel kernel, const T0& arg0, const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8)
{
	cl_int ciErrNum = clSetKernelArg(kernel, 0, sizeof(T0), &arg0);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 1, sizeof(T1), &arg1);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 2, sizeof(T2), &arg2);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 3, sizeof(T3), &arg3);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 4, sizeof(T4), &arg4);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 5, sizeof(T5), &arg5);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 6, sizeof(T6), &arg6);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 7, sizeof(T7), &arg7);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 8, sizeof(T8), &arg8);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	return CL_SUCCESS;
}

template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
FORCEINLINE cl_int setKernelArguments(cl_kernel kernel, const T0& arg0, const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9)
{
	cl_int ciErrNum = clSetKernelArg(kernel, 0, sizeof(T0), &arg0);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 1, sizeof(T1), &arg1);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 2, sizeof(T2), &arg2);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 3, sizeof(T3), &arg3);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 4, sizeof(T4), &arg4);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 5, sizeof(T5), &arg5);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 6, sizeof(T6), &arg6);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 7, sizeof(T7), &arg7);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 8, sizeof(T8), &arg8);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 9, sizeof(T9), &arg9);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	return CL_SUCCESS;
}

template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
FORCEINLINE cl_int setKernelArguments(cl_kernel kernel, const T0& arg0, const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10)
{
	cl_int ciErrNum = clSetKernelArg(kernel, 0, sizeof(T0), &arg0);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 1, sizeof(T1), &arg1);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 2, sizeof(T2), &arg2);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 3, sizeof(T3), &arg3);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 4, sizeof(T4), &arg4);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 5, sizeof(T5), &arg5);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 6, sizeof(T6), &arg6);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 7, sizeof(T7), &arg7);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 8, sizeof(T8), &arg8);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 9, sizeof(T9), &arg9);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 10, sizeof(T10), &arg10);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	return CL_SUCCESS;
}

template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
FORCEINLINE cl_int setKernelArguments(cl_kernel kernel, const T0& arg0, const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, const T11& arg11)
{
	cl_int ciErrNum = clSetKernelArg(kernel, 0, sizeof(T0), &arg0);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 1, sizeof(T1), &arg1);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 2, sizeof(T2), &arg2);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 3, sizeof(T3), &arg3);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 4, sizeof(T4), &arg4);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 5, sizeof(T5), &arg5);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 6, sizeof(T6), &arg6);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 7, sizeof(T7), &arg7);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 8, sizeof(T8), &arg8);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 9, sizeof(T9), &arg9);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 10, sizeof(T10), &arg10);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 11, sizeof(T11), &arg11);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	return CL_SUCCESS;
}

template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12>
FORCEINLINE cl_int setKernelArguments(cl_kernel kernel, const T0& arg0, const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, const T11& arg11, const T12& arg12)
{
	cl_int ciErrNum = clSetKernelArg(kernel, 0, sizeof(T0), &arg0);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 1, sizeof(T1), &arg1);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 2, sizeof(T2), &arg2);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 3, sizeof(T3), &arg3);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 4, sizeof(T4), &arg4);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 5, sizeof(T5), &arg5);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 6, sizeof(T6), &arg6);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 7, sizeof(T7), &arg7);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 8, sizeof(T8), &arg8);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 9, sizeof(T9), &arg9);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 10, sizeof(T10), &arg10);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 11, sizeof(T11), &arg11);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 12, sizeof(T12), &arg12);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	return CL_SUCCESS;
}

template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13>
FORCEINLINE cl_int setKernelArguments(cl_kernel kernel, const T0& arg0, const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, const T11& arg11, const T12& arg12, const T13& arg13)
{
	cl_int ciErrNum = clSetKernelArg(kernel, 0, sizeof(T0), &arg0);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 1, sizeof(T1), &arg1);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 2, sizeof(T2), &arg2);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 3, sizeof(T3), &arg3);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 4, sizeof(T4), &arg4);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 5, sizeof(T5), &arg5);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 6, sizeof(T6), &arg6);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 7, sizeof(T7), &arg7);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 8, sizeof(T8), &arg8);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 9, sizeof(T9), &arg9);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 10, sizeof(T10), &arg10);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 11, sizeof(T11), &arg11);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 12, sizeof(T12), &arg12);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 13, sizeof(T13), &arg13);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	return CL_SUCCESS;
}

template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14>
FORCEINLINE cl_int setKernelArguments(cl_kernel kernel, const T0& arg0, const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14)
{
	cl_int ciErrNum = clSetKernelArg(kernel, 0, sizeof(T0), &arg0);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 1, sizeof(T1), &arg1);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 2, sizeof(T2), &arg2);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 3, sizeof(T3), &arg3);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 4, sizeof(T4), &arg4);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 5, sizeof(T5), &arg5);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 6, sizeof(T6), &arg6);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 7, sizeof(T7), &arg7);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 8, sizeof(T8), &arg8);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 9, sizeof(T9), &arg9);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 10, sizeof(T10), &arg10);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 11, sizeof(T11), &arg11);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 12, sizeof(T12), &arg12);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 13, sizeof(T13), &arg13);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 14, sizeof(T14), &arg14);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	return CL_SUCCESS;
}

template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15>
FORCEINLINE cl_int setKernelArguments(cl_kernel kernel, const T0& arg0, const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15)
{
	cl_int ciErrNum = clSetKernelArg(kernel, 0, sizeof(T0), &arg0);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 1, sizeof(T1), &arg1);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 2, sizeof(T2), &arg2);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 3, sizeof(T3), &arg3);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 4, sizeof(T4), &arg4);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 5, sizeof(T5), &arg5);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 6, sizeof(T6), &arg6);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 7, sizeof(T7), &arg7);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 8, sizeof(T8), &arg8);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 9, sizeof(T9), &arg9);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 10, sizeof(T10), &arg10);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 11, sizeof(T11), &arg11);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 12, sizeof(T12), &arg12);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 13, sizeof(T13), &arg13);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 14, sizeof(T14), &arg14);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 15, sizeof(T15), &arg15);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	return CL_SUCCESS;
}

template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16>
FORCEINLINE cl_int setKernelArguments(cl_kernel kernel, const T0& arg0, const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16)
{
	cl_int ciErrNum = clSetKernelArg(kernel, 0, sizeof(T0), &arg0);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 1, sizeof(T1), &arg1);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 2, sizeof(T2), &arg2);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 3, sizeof(T3), &arg3);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 4, sizeof(T4), &arg4);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 5, sizeof(T5), &arg5);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 6, sizeof(T6), &arg6);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 7, sizeof(T7), &arg7);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 8, sizeof(T8), &arg8);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 9, sizeof(T9), &arg9);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 10, sizeof(T10), &arg10);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 11, sizeof(T11), &arg11);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 12, sizeof(T12), &arg12);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 13, sizeof(T13), &arg13);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 14, sizeof(T14), &arg14);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 15, sizeof(T15), &arg15);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	ciErrNum = clSetKernelArg(kernel, 16, sizeof(T16), &arg16);
	if (ciErrNum != CL_SUCCESS) return ciErrNum;

	return CL_SUCCESS;
}
