#pragma once

struct InverseData128
{
	unsigned char shift;
	number inverse[2];
	number max_value[2];
};

CACHE_ALIGNED extern const number PowersOf2_128DivisibilityData[64][2][2];
CACHE_ALIGNED extern const InverseData128 PowersOfP_128DivisibilityData[3][128];
CACHE_ALIGNED extern const number PrimeInverses_128[3][2][2];
