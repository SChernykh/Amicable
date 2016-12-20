#include "stdafx.h"
#include "Definitions.h"

#if DYNAMIC_SEARCH_LIMIT

namespace SearchLimit
{
	number value;
	number LinearLimit;
	number MainPrimeTableBound;
	number PrimeInversesBound;
	number SafeLimit;
}

#endif

void atoi128(const char* s, number &numlo, number &numhi)
{
	numlo = 0;
	numhi = 0;
	for (;;)
	{
		const unsigned char d = static_cast<unsigned char>(*(s++) - '0');
		if (d > 9)
		{
			break;
		}

		number k[2];
		k[0] = _umul128(numlo, 10, &k[1]);
		k[1] += numhi * 10;
		add128(k[0], k[1], d, 0, k, k + 1);
		numlo = k[0];
		numhi = k[1];
	}
}
