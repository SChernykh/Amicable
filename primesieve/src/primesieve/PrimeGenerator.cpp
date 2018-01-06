///
/// @file   PrimeGenerator.cpp
/// @brief  After a segment has been sieved PrimeGenerator is
///         used to reconstruct primes and prime k-tuplets from
///         1 bits of the sieve array.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <stdafx.h>

#include <primesieve/config.hpp>
#include <primesieve/littleendian_cast.hpp>
#include <primesieve/pmath.hpp>
#include <primesieve/PrimeGenerator.hpp>
#include <primesieve/PrimeSieve.hpp>
#include <primesieve/SieveOfEratosthenes.hpp>
#include <primesieve/StorePrimes.hpp>

#include <stdint.h>
#include <algorithm>
#include <iostream>

using namespace std;

namespace primesieve {

const uint64_t PrimeGenerator::bitmasks_[6][5] =
{
  { END },
  { 0x06, 0x18, 0xc0, END },       // Twin primes:       b00000110, b00011000, b11000000
  { 0x07, 0x0e, 0x1c, 0x38, END }, // Prime triplets:    b00000111, b00001110, ...
  { 0x1e, END },                   // Prime quadruplets: b00011110
  { 0x1f, 0x3e, END },             // Prime quintuplets
  { 0x3f, END }                    // Prime sextuplets
};

PrimeGenerator::PrimeGenerator(PrimeSieve& ps, const PreSieve& preSieve) :
  SieveOfEratosthenes(max<uint64_t>(7, ps.getStart()),
                      ps.getStop(),
                      static_cast<unsigned int>(ps.getSieveSize()),
                      preSieve),
  ps_(ps),
  counts_(ps_.getCounts())
{
  if (ps_.isFlag(ps_.COUNT_TWINS, ps_.COUNT_SEXTUPLETS))
    init_kCounts();
}

/// Calculate the number of twins, triplets, ...
/// for each possible byte value
///
void PrimeGenerator::init_kCounts()
{
  for (uint_t i = 1; i < counts_.size(); i++)
  {
    if (ps_.isCount(static_cast<int>(i)))
    {
      kCounts_[i].resize(256);

      for (uint64_t j = 0; j < 256; j++)
      {
        byte_t count = 0;
        for (const uint64_t* b = bitmasks_[i]; *b <= j; b++)
        {
          if ((j & *b) == *b)
            count++;
        }
        kCounts_[i][j] = count;
      }
    }
  }
}

/// Executed after each sieved segment.
/// @see sieveSegment() in SieveOfEratosthenes.cpp
///
void PrimeGenerator::generatePrimes(const byte_t* sieve, uint64_t sieveSize)
{
  if (ps_.isStore())
    storePrimes(ps_.getStore(), sieve, sieveSize);
  if (ps_.isStatus())
    ps_.updateStatus(sieveSize * NUMBERS_PER_BYTE);
}

void PrimeGenerator::storePrimes(Store& store, const byte_t* sieve, uint64_t sieveSize) const
{
  uint64_t low = getSegmentLow();

  for (uint64_t i = 0; i < sieveSize; i += 8)
  {
    uint64_t bits = littleendian_cast<uint64_t>(&sieve[i]); 
    while (bits)
      store(nextPrime(&bits, low));

    low += NUMBERS_PER_BYTE * 8;
  }
}

} // namespace
