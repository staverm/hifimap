// Copyright (c) 2021 Mauro Staver
// Some of the code is taken from https://github.com/lbcb-sci/ram

#ifndef HIFIMAP_MINIMIZER_ENGINE_HPP_
#define HIFIMAP_MINIMIZER_ENGINE_HPP_

#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "thread_pool/thread_pool.hpp"

namespace hifimap {

class MinimizerEngine {
public:
  double minimize_query = 0;
  double index_lookup = 0;
  double find_hits = 0;
  double minimize_region = 0;
  double find_matches = 0;
  double chain_matches = 0;
  double best_regions_ = 0;

  MinimizerEngine(
      std::shared_ptr<thread_pool::ThreadPool> thread_pool = nullptr,
      std::uint32_t k = 15, // element of [1, 31]
      std::uint32_t w = 5, std::uint32_t bandwidth = 500,
      std::uint32_t chain = 4, std::uint32_t matches = 100,
      std::uint32_t gap = 10000, std::uint32_t avg_len = 25000,
      std::uint32_t max_regions = 2, std::uint32_t min_hits = 10);

  MinimizerEngine(const MinimizerEngine &) = delete;
  MinimizerEngine &operator=(const MinimizerEngine &) = delete;

  MinimizerEngine(MinimizerEngine &&) = default;
  MinimizerEngine &operator=(MinimizerEngine &&) = default;

  ~MinimizerEngine() = default;

  void CreateIndex(
      std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator first,
      std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator last,
      double tfidf_offset, bool minhash = false, bool fast = true);

  std::vector<biosoup::Overlap>
  HifiMap(const std::unique_ptr<biosoup::NucleicAcid> &sequence, bool minhash,
          bool approx = false);

private:
  struct Kmer {
  public:
    Kmer() = default;
    Kmer(std::uint64_t value, std::uint64_t origin)
        : value(value), origin(origin) {}

    std::uint32_t id() const {
      return static_cast<std::uint32_t>(origin >> 32);
    }

    std::uint32_t position() const {
      return static_cast<std::uint32_t>(origin) >> 1;
    }

    bool strand() const { return origin & 1; }

    static std::uint64_t SortByValue(const Kmer &kmer) { return kmer.value; }

    static std::uint64_t SortByPosition(const Kmer &kmer) {
      return kmer.position();
    }

    std::uint64_t value;
    std::uint64_t origin; // 32 bit id, 31 bit position, 1 bit strand
  };

  struct Match {
  public:
    Match() = default;
    Match(std::uint64_t group, std::uint64_t positions)
        : group(group), positions(positions) {}

    std::uint32_t rhs_id() const {
      return static_cast<std::uint32_t>(group >> 33);
    }

    bool strand() const { return (group >> 32) & 1; }

    std::uint32_t diagonal() const { return static_cast<std::uint32_t>(group); }

    std::uint32_t lhs_position() const {
      return static_cast<std::uint32_t>(positions >> 32);
    }

    std::uint32_t rhs_position() const {
      return static_cast<std::uint32_t>(positions);
    }

    static std::uint64_t SortByGroup(const Match &match) { return match.group; }
    static std::uint64_t SortByPositions(const Match &match) {
      return match.positions;
    }

    std::uint64_t group;
    std::uint64_t positions;
  };

  class Index {
  public:
    Index() = default;

    std::uint32_t Find(std::uint64_t key, std::uint64_t *dst) const;

    struct Hash {
      std::size_t operator()(std::uint64_t key) const {
        return std::hash<std::uint64_t>()(key >> 1);
      }
    };
    struct KeyEqual {
      bool operator()(std::uint64_t lhs, std::uint64_t rhs) const {
        return (lhs >> 1) == (rhs >> 1);
      }
    };

    std::vector<std::pair<std::uint64_t, std::uint32_t>>
        regions; // id_region, cnt

    std::unordered_map<std::uint64_t, std::uint64_t, Hash, KeyEqual> locator;
  };

  std::vector<Kmer>
  Minimize(const std::unique_ptr<biosoup::NucleicAcid> &sequence,
           bool minhash = false) const;

  std::vector<Kmer>
  Minimize(const std::unique_ptr<biosoup::NucleicAcid> &sequence,
           std::uint32_t pos, std::uint32_t len, bool minhash = false) const;

  std::vector<biosoup::Overlap> Chain(std::uint64_t lhs_id,
                                      std::vector<Match> &&matches) const;

  template <typename RandomAccessIterator, typename Compare>
  static void RadixSort(RandomAccessIterator first, RandomAccessIterator last,
                        std::uint8_t max_bits,
                        Compare comp); //  unary comparison function

  template <typename Compare>
  static std::vector<std::uint64_t>
  LongestSubsequence(std::vector<Match>::const_iterator first,
                     std::vector<Match>::const_iterator last,
                     Compare comp); // binary comparison function

  std::pair<std::uint64_t, std::uint64_t>
  RegionsOf(std::uint32_t pos, std::uint32_t num_regions) const;

  std::vector<std::uint64_t>
  BestRegions(std::vector<std::vector<std::uint32_t>> &hits) const;

  std::uint32_t k_;
  std::uint32_t w_;
  std::uint32_t bandwidth_;
  std::uint32_t chain_;
  std::uint32_t matches_;
  std::uint64_t gap_;
  std::shared_ptr<thread_pool::ThreadPool> thread_pool_;

  // hifimap members
  std::uint32_t avg_len_; // avg read len
  std::uint32_t region_len_;
  std::uint32_t max_regions_;
  std::uint32_t min_hits_;
  std::unordered_map<std::uint32_t,
                     const std::unique_ptr<biosoup::NucleicAcid> &>
      targets_;
  std::vector<std::uint32_t> num_regions_; // index = id
  std::uint32_t total_regions_ = 0;
  std::vector<Index> index_;
  bool fast_ = false;
  std::vector<std::vector<Kmer>> fast_minimizers;
  std::vector<std::vector<std::uint64_t>>
      fast_locator; // index = id, region   ,  value = pos | cnt
};

} // namespace hifimap

#endif // HIFIMAP_MINIMIZER_ENGINE_HPP_
