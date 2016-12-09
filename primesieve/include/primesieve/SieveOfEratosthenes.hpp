///
/// @file  SieveOfEratosthenes.hpp
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef SIEVEOFERATOSTHENES_HPP
#define SIEVEOFERATOSTHENES_HPP

#include "config.hpp"
#include <stdint.h>
#include <stdafx.h>

namespace primesieve {

class PreSieve;
class EratSmall;
class EratMedium;
class EratBig;

/// The abstract SieveOfEratosthenes class sieves primes using the
/// segmented sieve of Eratosthenes. It uses a bit array for sieving,
/// the bit array uses 8 flags for 30 numbers. SieveOfEratosthenes
/// uses three different sieve of Eratosthenes algorithms optimized
/// for small, medium and big sieving primes to cross-off multiples.
///
class SieveOfEratosthenes {
public:
  uint64_t getStart() const;
  uint64_t getStop() const;
  uint_t getSqrtStop() const;
  uint_t getSieveSize() const;
  void addSievingPrime(uint_t);
  void sieve();
protected:
  SieveOfEratosthenes(uint64_t, uint64_t, uint_t, const PreSieve&);
  virtual ~SieveOfEratosthenes();
  virtual void segmentFinished(const byte_t*, uint_t) = 0;

  static FORCEINLINE uint64_t getNextPrime(uint64_t* bits, uint64_t base)
  {
	  // calculate bitValues_[ bitScanForward(*bits) ]
	  // using a custom De Bruijn bitscan
	  //uint64_t debruijn64 = UINT64_C(0x3F08A4C6ACB9DBD);
	  uint64_t mask = *bits - 1;
	  //uint64_t bitValue = bruijnBitValues_[((*bits ^ mask) * debruijn64) >> 58];
	  //uint64_t prime = base + bitValue;
	  unsigned long index;
	  _BitScanForward64(&index, *bits);
	  uint64_t prime = base + bitValuesRaw_[index];
	  *bits &= mask;
	  return prime;
  }

  FORCEINLINE uint64_t getSegmentLow() const { return segmentLow_; }
private:
  static const uint_t bitValues_[8];
  static const uint_t bruijnBitValues_[64];
  static CACHE_ALIGNED const uint8_t bitValuesRaw_[64];
  /// Lower bound of the current segment
  uint64_t segmentLow_;
  /// Upper bound of the current segment
  uint64_t segmentHigh_;
  /// Sieve primes >= start_
  const uint64_t start_;
  /// Sieve primes <= stop_
  const uint64_t stop_;
  const PreSieve& preSieve_;
  /// sqrt(stop_)
  uint_t sqrtStop_;
  uint_t limitPreSieve_;
  uint_t limitEratSmall_;
  uint_t limitEratMedium_;
  /// Size of sieve_ in bytes (power of 2)
  uint_t sieveSize_;
  /// Sieve of Eratosthenes array
  byte_t* sieve_;
  EratSmall* eratSmall_;
  EratMedium* eratMedium_;
  EratBig* eratBig_;
  static uint64_t getByteRemainder(uint64_t);
  void init();
  void cleanUp();
  void preSieve();
  void crossOffMultiples();
  void sieveSegment();
  DISALLOW_COPY_AND_ASSIGN(SieveOfEratosthenes);
};

} // namespace

#endif
