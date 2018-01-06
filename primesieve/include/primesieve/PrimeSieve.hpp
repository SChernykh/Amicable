///
/// @file   PrimeSieve.hpp
/// @brief  The PrimeSieve class is a high level class that
///         manages prime sieving using the PreSieve, SievingPrimes
///         and PrimeGenerator classes.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PRIMESIEVE_CLASS_HPP
#define PRIMESIEVE_CLASS_HPP

#include <stdint.h>
#include <vector>
#include <array>
#include <string>

#include <primesieve/SievingPrimes.hpp>
#include <primesieve/PrimeGenerator.hpp>

namespace primesieve {

class Store;

struct SmallPrime
{
	uint64_t first;
	uint64_t last;
	int index;
	std::string str;
};

extern const std::array<SmallPrime, 8> smallPrimes;

class PrimeSieve
{
public:
  enum
  {
    COUNT_PRIMES      = 1 << 0,
    COUNT_TWINS       = 1 << 1,
    COUNT_TRIPLETS    = 1 << 2,
    COUNT_QUADRUPLETS = 1 << 3,
    COUNT_QUINTUPLETS = 1 << 4,
    COUNT_SEXTUPLETS  = 1 << 5,
    PRINT_PRIMES      = 1 << 6,
    PRINT_TWINS       = 1 << 7,
    PRINT_TRIPLETS    = 1 << 8,
    PRINT_QUADRUPLETS = 1 << 9,
    PRINT_QUINTUPLETS = 1 << 10,
    PRINT_SEXTUPLETS  = 1 << 11,
    PRINT_STATUS      = 1 << 12,
    CALCULATE_STATUS  = 1 << 13
  };
  PrimeSieve();
  PrimeSieve(PrimeSieve*);
  virtual ~PrimeSieve();
  // Getters
  uint64_t getStart() const;
  uint64_t getStop() const;
  int getSieveSize() const;
  double getStatus() const;
  double getSeconds() const;
  Store& getStore();
  // Setters
  void setStart(uint64_t);
  void setStop(uint64_t);
  void setSieveSize(int);
  void setFlags(int);
  void addFlags(int);
  // Bool is*
  bool isCount() const;
  bool isCount(int) const;
  bool isPrint() const;
  bool isPrint(int) const;
  bool isFlag(int) const;
  bool isFlag(int, int) const;
  bool isStatus() const;
  bool isStore() const;
  // Sieve
  virtual void sieve();

  template<typename T>
  void sieveTemplated(uint64_t start, uint64_t stop, T&& callback)
  {
	  setStart(start);
	  setStop(stop);

	  reset();
	  if (start_ > stop_)
		  return;

	  // small primes and k-tuplets (first prime <= 5)
	  if (start_ <= 5)
	  {
		  for (auto& p : smallPrimes)
		  {
			  if (p.index != 0)
			  {
				  break;
			  }
			  if (p.first >= start_ && p.last <= stop_)
			  {
				  callback(p.first);
			  }
		  }
	  }

	  if (stop_ >= 7)
	  {
		  PreSieve preSieve(start_, stop_);
		  PrimeGeneratorTemplated<T> primeGen(*this, preSieve, static_cast<T&&>(callback));

		  if (primeGen.getSqrtStop() > preSieve.getMaxPrime())
		  {
			  // generate the sieving primes <= sqrt(stop)
			  SievingPrimes sp(primeGen, preSieve);
			  sp.generate();
		  }

		  // sieve the primes within [start, stop]
		  primeGen.sieve();
	  }
  }

  void sieve(uint64_t, uint64_t);
  void sieve(uint64_t, uint64_t, int);
  void storePrimes(uint64_t, uint64_t, Store*);
  // nth prime
  uint64_t nthPrime(uint64_t);
  uint64_t nthPrime(int64_t, uint64_t);
  // Print
  void printPrimes(uint64_t, uint64_t);
  void printTwins(uint64_t, uint64_t);
  void printTriplets(uint64_t, uint64_t);
  void printQuadruplets(uint64_t, uint64_t);
  void printQuintuplets(uint64_t, uint64_t);
  void printSextuplets(uint64_t, uint64_t);
  // Count
  uint64_t countPrimes(uint64_t, uint64_t);
  uint64_t countTwins(uint64_t, uint64_t);
  uint64_t countTriplets(uint64_t, uint64_t);
  uint64_t countQuadruplets(uint64_t, uint64_t);
  uint64_t countQuintuplets(uint64_t, uint64_t);
  uint64_t countSextuplets(uint64_t, uint64_t);
  // Count getters
  std::vector<uint64_t>& getCounts();
  uint64_t getPrimeCount() const;
  uint64_t getTwinCount() const;
  uint64_t getTripletCount() const;
  uint64_t getQuadrupletCount() const;
  uint64_t getQuintupletCount() const;
  uint64_t getSextupletCount() const;
  uint64_t getCount(size_t) const;
  virtual bool updateStatus(uint64_t, bool tryLock = true);
protected:
  /// Sieve primes >= start_
  uint64_t start_;
  /// Sieve primes <= stop_
  uint64_t stop_;
  /// Prime number and prime k-tuplet counts
  std::vector<uint64_t> counts_;
  /// Time elapsed of sieve()
  double seconds_;
  uint64_t getDistance() const;
  void reset();
private:
  /// Sum of all processed segments
  uint64_t processed_;
  /// Sum of processed segments to update
  uint64_t toUpdate_;
  /// Status of sieve() in percent
  double percent_;
  /// Sieve size in kilobytes
  int sieveSize_;
  /// Setter methods set flags e.g. COUNT_PRIMES
  int flags_;
  /// parent ParallelPrimeSieve object
  PrimeSieve* parent_;
  Store* store_;
  static void printStatus(double, double);
  bool isParallelPrimeSieve() const;
  void processSmallPrimes();
};

} // namespace

#endif
