// Copyright (c) 2021 Mauro Staver
// Some of the code is taken from https://github.com/lbcb-sci/ram

#include "hifimap/minimizer_engine.hpp"

#include "biosoup/timer.hpp"
#include <cmath>
#include <deque>
#include <unordered_set>

namespace hifimap {

MinimizerEngine::MinimizerEngine(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, std::uint32_t k,
    std::uint32_t w, std::uint32_t bandwidth, std::uint32_t chain,
    std::uint32_t matches, std::uint32_t gap, std::uint32_t avg_len,
    std::uint32_t max_regions, std::uint32_t min_hits)
    : k_(std::min(std::max(k, 1U), 31U)), w_(w), bandwidth_(bandwidth),
      chain_(chain), matches_(matches), gap_(gap),
      thread_pool_(thread_pool ? thread_pool
                               : std::make_shared<thread_pool::ThreadPool>(1)),
      avg_len_(avg_len), region_len_(avg_len * 2), max_regions_(max_regions),
      min_hits_(min_hits), index_(1U << std::min(14U, 2 * k_)) {}

std::pair<std::uint64_t, std::uint64_t>
MinimizerEngine::RegionsOf(std::uint32_t pos, std::uint32_t id) const {
  std::uint32_t region_a = pos / avg_len_;
  std::uint32_t region_b = 0;

  // if 0 -> (0,0)
  // if num_regions -> (num_regions, 0)
  // else (region_a - 1, region_a)

  std::uint64_t id_ = static_cast<std::uint64_t>(id) << 32;

  if (region_a == num_regions_[id])
    return {id_ | region_a, id_ | region_b};

  region_b = region_a;

  if (region_a != 0) {
    region_a--;
  }

  // if region_b is zero it should be ignored
  return {id_ | region_a, id_ | region_b};
}

void MinimizerEngine::CreateIndex(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator first,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator last,
    double tfidf_offset, bool minhash, bool fast) {
  for (auto &it : index_) {
    it.locator.clear();
  }

  if (first >= last) {
    return;
  }

  fast_ = fast;

  auto update_flocator = [&](std::uint32_t id, std::uint32_t region,
                             std::uint32_t pos) -> void {
    if (fast_locator[id][region] == 0) {
      fast_locator[id][region] = (static_cast<std::uint64_t>(pos) << 32) | 1;
    } else {
      fast_locator[id][region]++;
    }
  };

  std::vector<std::vector<Kmer>> minimizers(index_.size());
  { // computes minimizers
    std::uint64_t mask = index_.size() - 1;

    while (first != last) {
      std::size_t batch_size = 0;
      std::vector<std::future<std::vector<Kmer>>> futures;
      for (; first != last && batch_size < 50000000; ++first) {
        targets_.emplace((*first)->id, *first); // save id -> target
        if (num_regions_.size() != (*first)->id)
          continue;
        num_regions_.emplace_back((*first)->inflated_len / avg_len_);
        total_regions_ += (*first)->inflated_len / avg_len_;

        if (fast_) {
          fast_locator.emplace_back(num_regions_[(*first)->id], 0);
          fast_minimizers.emplace_back(std::vector<Kmer>());
        }

        batch_size += (*first)->inflated_len;
        futures.emplace_back(thread_pool_->Submit(
            [&](decltype(first) it) -> std::vector<Kmer> {
              auto minimzs = Minimize(*it, minhash);
              if (minhash) {
                RadixSort(minimzs.begin(), minimzs.end(), k_ * 2,
                          Kmer::SortByPosition);
              }
              return minimzs;
            },
            first));
      }

      for (auto &it : futures) {
        for (const auto &jt : it.get()) { // jt = kmer
          auto &m = minimizers[jt.value & mask];
          if (m.capacity() == m.size()) {
            m.reserve(m.capacity() * 1.5);
          }
          m.emplace_back(jt);

          if (fast_) {
            auto jt_id = jt.id();
            auto &fm = fast_minimizers[jt_id];
            if (fm.capacity() == fm.size()) {
              fm.reserve(fm.capacity() * 1.5);
            }
            fm.emplace_back(jt);

            auto regions = RegionsOf(jt.position(), jt_id);

            if (static_cast<std::uint32_t>(regions.first) ==
                num_regions_[jt_id]) {
              regions.first--;
            }
            update_flocator(jt_id, static_cast<std::uint32_t>(regions.first),
                            fm.size() - 1);

            if (static_cast<std::uint32_t>(regions.second) != 0) {
              update_flocator(jt_id, static_cast<std::uint32_t>(regions.second),
                              fm.size() - 1);
            }
          }
        }
      }
    }
  }

  {
    std::vector<std::future<std::pair<std::size_t, std::size_t>>> futures;
    for (std::uint32_t i = 0; i < minimizers.size(); ++i) {
      futures.emplace_back(thread_pool_->Submit(
          [&](std::uint32_t i) -> std::pair<std::size_t, std::size_t> {
            if (minimizers[i].empty()) {
              return std::make_pair(0, 0);
            }

            RadixSort(minimizers[i].begin(), minimizers[i].end(), k_ * 2,
                      Kmer::SortByValue);

            minimizers[i].emplace_back(-1, -1); // stop dummy

            std::size_t num_origins = 0;
            std::size_t num_keys = 0;

            for (std::uint64_t j = 1, c = 1; j < minimizers[i].size();
                 ++j, ++c) { // NOLINT
              if (minimizers[i][j - 1].value != minimizers[i][j].value) {
                if (c > 1) {
                  num_origins += c;
                }
                ++num_keys;
                c = 0;
              }
            }

            return std::make_pair(num_origins, num_keys);
          },
          i));
    }

    auto idf = [](std::uint32_t num_regions, size_t region_size) -> double {
      return log(num_regions / (double)region_size);
    };

    for (std::uint32_t i = 0; i < minimizers.size(); ++i) {
      auto num_entries = futures[i].get();
      if (minimizers[i].empty()) {
        continue;
      }

      index_[i].regions.reserve(num_entries.first * 2);
      index_[i].locator.reserve(num_entries.second);

      for (std::uint64_t j = 1, c = 1; j < minimizers[i].size(); ++j, ++c) {
        if (minimizers[i][j - 1].value != minimizers[i][j].value) {
          if (c == 1) { // singleton optimization
            auto regions = RegionsOf(minimizers[i][j - 1].position(),
                                     minimizers[i][j - 1].id());
            if (regions.second == 0)
              index_[i].locator.emplace(minimizers[i][j - 1].value << 1 | 1,
                                        regions.first);
            else
              index_[i].locator.emplace(minimizers[i][j - 1].value << 1 | 1,
                                        regions.second);
          } else {

            auto prev = RegionsOf(minimizers[i][j - c].position(),
                                  minimizers[i][j - c].id());

            int tf1 = 1;
            int tf2 = 1;
            int cnt = 0;

            std::vector<std::pair<std::uint64_t, std::uint32_t>> buffer;
            // emplace (region, tf) pairs of last c minimizers
            for (std::uint64_t k = j - c + 1; k < j; ++k) {
              auto curr =
                  RegionsOf(minimizers[i][k].position(), minimizers[i][k].id());

              // 2nd half of last region
              if (static_cast<std::uint32_t>(prev.first) ==
                  num_regions_[minimizers[i][k].id()])
                prev.first = prev.first - 1;
              if (static_cast<std::uint32_t>(curr.first) ==
                  num_regions_[minimizers[i][k].id()])
                curr.first = curr.first - 1;

              if (prev.first == curr.first) {
                tf1++;
                if (static_cast<std::uint32_t>(prev.second) != 0)
                  tf2++;
              } else if (prev.second == curr.first) {
                buffer.emplace_back(prev.first, tf1);
                cnt++;

                tf1 = tf2;
                tf1++;
                tf2 = 1;
              } else {
                buffer.emplace_back(prev.first, tf1);
                cnt++;
                if (static_cast<std::uint32_t>(prev.second) != 0) {
                  buffer.emplace_back(prev.second, tf2);
                  cnt++;
                }
                tf1 = 1;
                tf2 = 1;
              }
              if (k == j - 1) {
                buffer.emplace_back(curr.first, tf1);
                cnt++;
                if (static_cast<std::uint32_t>(curr.second) != 0) {
                  buffer.emplace_back(curr.second, tf2);
                  cnt++;
                }
              }
              prev = curr;
            }

            std::uint32_t emplaced = 0;
            for (auto &[id_region, tf] : buffer) {
              auto num_regions =
                  num_regions_[static_cast<std::uint32_t>(id_region >> 32)];
              if (tf * idf(num_regions, cnt) >=
                  idf(num_regions, 1) + tfidf_offset) {
                index_[i].regions.emplace_back(id_region, tf);
                emplaced++;
              }
            }

            index_[i].locator.emplace(
                minimizers[i][j - 1].value << 1,
                (index_[i].regions.size() - emplaced) << 32 | emplaced);
          }
          c = 0;
        }
      }
      std::vector<Kmer>().swap(minimizers[i]);
    }
  }
}

std::uint32_t MinimizerEngine::Index::Find(std::uint64_t key,
                                           std::uint64_t *dst) const {
  auto it = locator.find(key << 1);
  if (it == locator.end()) {
    return 0;
  }
  if (it->first & 1) {
    *dst = it->second; // id_region
    return 1;          // singleton
  }
  *dst = (it->second >> 32);                     // pos
  return static_cast<std::uint32_t>(it->second); // cnt
}

std::vector<std::uint64_t> MinimizerEngine::BestRegions(
    std::vector<std::vector<std::uint32_t>> &hits) const {
  std::unordered_set<std::uint64_t> best;

  for (std::uint64_t id = 0; id < hits.size(); id++) {
    std::vector<std::uint32_t> idx(hits[id].size());
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(),
                     [&](std::uint32_t i1, std::uint32_t i2) {
                       return hits[id][i1] > hits[id][i2];
                     });

    best.emplace((id << 32) | idx[0]);
    std::uint32_t cnt = 1;
    std::uint32_t j = 1;
    while (cnt < max_regions_) {
      if (j >= idx.size())
        break;
      if (hits[id][idx[j]] < min_hits_)
        break;
      if (!best.count((id << 32) | (idx[j] - 1)) &&
          !best.count((id << 32) | (idx[j] + 1))) {
        best.emplace((id << 32) | idx[j]);
        cnt++;
      }
      j++;
    }
  }

  return std::vector<std::uint64_t>(best.begin(), best.end());
}

std::vector<biosoup::Overlap>
MinimizerEngine::HifiMap(const std::unique_ptr<biosoup::NucleicAcid> &sequence,
                         bool minhash, bool approx) {
  biosoup::Timer timer{};

  timer.Start();
  auto sketch = Minimize(sequence, minhash);

  if (sketch.empty()) {
    this->minimize_query += timer.Stop();
    return std::vector<biosoup::Overlap>{};
  }

  std::uint64_t mask = index_.size() - 1;
  std::vector<Match> matches;
  std::vector<biosoup::Overlap> overlaps;
  std::vector<std::vector<std::uint32_t>> hits;
  std::vector<std::uint64_t> max_hits(
      num_regions_.size(), 0); // index = id, value = region | max_hits

  hits.reserve(num_regions_.size());
  for (std::uint32_t i = 0; i < num_regions_.size(); i++) {
    hits.emplace_back(num_regions_[i], 0);
  }

  if (!minhash) {
    RadixSort(sketch.begin(), sketch.end(), k_ * 2, Kmer::SortByValue);
  }
  this->minimize_query += timer.Stop();

  for (std::uint64_t j = 1, c = 1; j < sketch.size(); ++j, ++c) {
    if (sketch[j - 1].value != sketch[j].value) {

      timer.Start();
      auto &it = sketch[j - 1];
      std::uint32_t i = it.value & mask;

      std::uint64_t
          region; // id_region if singleton, otherwise pos in regions list
      auto n = index_[i].Find(it.value, &region);
      if (n == 0) {
        c = 0;
        continue;
      }
      this->index_lookup += timer.Stop();

      timer.Start();

      // if singleton
      // if 0 -> only 0
      // else if num_regions -> only num_regions - 1
      // else if x: x and x-1

      auto update_hits = [&](std::uint64_t id_region,
                             std::uint32_t tf) -> void {
        std::uint32_t id = id_region >> 32;
        std::uint32_t reg = static_cast<std::uint32_t>(id_region);

        hits[id][reg] += std::min(tf, static_cast<std::uint32_t>(c));
        if (hits[id][reg] > static_cast<std::uint32_t>(max_hits[id])) {
          max_hits[id] =
              (static_cast<std::uint64_t>(reg) << 32) | hits[id][reg];
        }
      };

      // calculate region hits
      for (std::uint32_t k = 0; k < n; ++k, ++region) {
        if (n == 1) {
          if (region == 0) {
            update_hits(region, 1);
          } else if (region == num_regions_[region >> 32]) {
            update_hits(region - 1, 1);
          } else {
            update_hits(region, 1);
            update_hits(region - 1, 1);
          }
        } else {
          update_hits(index_[i].regions[region].first,
                      index_[i].regions[region].second);
        }
      }
      this->find_hits += timer.Stop();
      c = 0;
    }
  }

  if (hits.size() < 1)
    return std::vector<biosoup::Overlap>{};

  if (approx) {
    for (std::uint32_t id = 0; id < max_hits.size(); id++) {
      std::uint32_t reg = max_hits[id] >> 32;
      std::uint32_t num_hits = static_cast<std::uint32_t>(max_hits[id]);
      std::uint32_t reg_midpos = reg * avg_len_ + avg_len_;

      if (reg == 0) {
        double right = hits[id][reg + 1] / (double)num_hits;
        overlaps.emplace_back(sequence->id, 0, sequence->inflated_len - 1, id,
                              reg_midpos + right * sequence->inflated_len -
                                  sequence->inflated_len,
                              reg_midpos + right * sequence->inflated_len, 0);
      } else {
        double left = hits[id][reg - 1] / (double)num_hits;
        overlaps.emplace_back(sequence->id, 0, sequence->inflated_len - 1, id,
                              reg_midpos - left * sequence->inflated_len,
                              reg_midpos - left * sequence->inflated_len +
                                  sequence->inflated_len,
                              0);
      }
    }
    return overlaps;
  }

  timer.Start();
  std::vector<std::uint64_t> best_regions;
  if (max_regions_ == 1) {
    for (std::uint64_t i = 0; i < max_hits.size(); i++) {
      if (static_cast<std::uint32_t>(max_hits[i]) < min_hits_)
        continue;
      best_regions.emplace_back((i << 32) | (max_hits[i] >> 32));
    }
  } else {
    best_regions = BestRegions(hits);
  }
  this->best_regions_ += timer.Stop();

  if (best_regions.size() < 1)
    return std::vector<biosoup::Overlap>{};

  // map query to each chosen region
  for (const auto &id_region : best_regions) {
    auto target = targets_.find(id_region >> 32);
    if (target == targets_.end()) {
      continue;
    }

    timer.Start();
    std::vector<Kmer> region_sketch;
    if (fast_) {
      std::uint64_t &pos_cnt =
          fast_locator[id_region >> 32][static_cast<std::uint32_t>(id_region)];
      region_sketch.reserve(static_cast<std::uint32_t>(pos_cnt));

      for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(pos_cnt); i++) {
        region_sketch.push_back(
            fast_minimizers[id_region >> 32][(pos_cnt >> 32) + i]);
      }
    } else {
      region_sketch = Minimize(target->second,
                               static_cast<std::uint32_t>(id_region) * avg_len_,
                               region_len_, minhash);
    }

    this->minimize_region += timer.Stop();

    timer.Start();

    if (fast_ || !minhash) {
      RadixSort(region_sketch.begin(), region_sketch.end(), k_ * 2,
                Kmer::SortByValue);
    }

    std::uint64_t rhs_id = id_region >> 32;

    // find matches
    for (std::uint32_t q = 0, r = 0; q < sketch.size(); ++q) {
      while (r < region_sketch.size()) {
        if (sketch[q].value < region_sketch[r].value) {
          break;
        } else if (sketch[q].value == region_sketch[r].value) {
          for (std::uint32_t k = r; k < region_sketch.size(); ++k) {
            if (sketch[q].value != region_sketch[k].value) {
              break;
            }

            std::uint64_t strand =
                (sketch[q].strand() & 1) == (region_sketch[k].strand() & 1);
            std::uint64_t lhs_pos = sketch[q].position();
            std::uint64_t rhs_pos = region_sketch[k].position();
            std::uint64_t diagonal =
                !strand ? rhs_pos + lhs_pos : rhs_pos - lhs_pos + (3ULL << 30);

            matches.emplace_back((((rhs_id << 1) | strand) << 32) | diagonal,
                                 (lhs_pos << 32) | rhs_pos);
          }
          break;
        } else {
          ++r;
        }
      }
    }

    this->find_matches += timer.Stop();
  }

  timer.Start();
  overlaps = Chain(sequence->id, std::move(matches));
  this->chain_matches += timer.Stop();
  return overlaps;
}

std::vector<biosoup::Overlap>
MinimizerEngine::Chain(std::uint64_t lhs_id,
                       std::vector<Match> &&matches) const {

  RadixSort(matches.begin(), matches.end(), 64, Match::SortByGroup);
  matches.emplace_back(-1, -1); // stop dummy

  std::vector<std::pair<std::uint64_t, std::uint64_t>> intervals;
  for (std::uint64_t i = 1, j = 0; i < matches.size(); ++i) { // NOLINT
    if (matches[i].group - matches[j].group > bandwidth_) {
      if (i - j >= 4) {
        if (!intervals.empty() && intervals.back().second > j) { // extend
          intervals.back().second = i;
        } else { // new
          intervals.push_back({j, i});
        }
      }
      ++j;
      while (j < i && matches[i].group - matches[j].group > bandwidth_) {
        ++j;
      }
    }
  }

  std::vector<biosoup::Overlap> dst;

  for (const auto &it : intervals) {
    std::uint64_t j = it.first;
    std::uint64_t i = it.second;

    if (i - j < chain_) {
      continue;
    }

    RadixSort(matches.begin() + j, matches.begin() + i, 64,
              Match::SortByPositions);

    std::uint64_t strand = matches[j].strand();

    std::vector<std::uint64_t> indices;
    if (strand) {                   // same strand
      indices = LongestSubsequence( // increasing
          matches.begin() + j, matches.begin() + i, std::less<std::uint64_t>());
    } else {                        // different strand
      indices = LongestSubsequence( // decreasing
          matches.begin() + j, matches.begin() + i,
          std::greater<std::uint64_t>());
    }

    if (indices.size() < chain_) {
      continue;
    }

    indices.emplace_back(matches.size() - 1 - j); // stop dummy from above
    for (std::uint64_t k = 1, l = 0; k < indices.size(); ++k) {
      if (matches[j + indices[k]].lhs_position() -
              matches[j + indices[k - 1]].lhs_position() >
          gap_) {
        if (k - l < chain_) {
          l = k;
          continue;
        }

        std::uint32_t lhs_matches = 0;
        std::uint32_t lhs_begin = 0;
        std::uint32_t lhs_end = 0;
        std::uint32_t rhs_matches = 0;
        std::uint32_t rhs_begin = 0;
        std::uint32_t rhs_end = 0;

        for (std::uint64_t m = l; m < k; ++m) {
          std::uint32_t lhs_pos = matches[j + indices[m]].lhs_position();
          if (lhs_pos > lhs_end) {
            lhs_matches += lhs_end - lhs_begin;
            lhs_begin = lhs_pos;
          }
          lhs_end = lhs_pos + k_;

          std::uint32_t rhs_pos = matches[j + indices[m]].rhs_position();
          rhs_pos = strand ? rhs_pos : (1U << 31) - (rhs_pos + k_ - 1);
          if (rhs_pos > rhs_end) {
            rhs_matches += rhs_end - rhs_begin;
            rhs_begin = rhs_pos;
          }
          rhs_end = rhs_pos + k_;
        }
        lhs_matches += lhs_end - lhs_begin;
        rhs_matches += rhs_end - rhs_begin;
        if (std::min(lhs_matches, rhs_matches) < matches_) {
          l = k;
          continue;
        }

        dst.emplace_back(lhs_id, matches[j + indices[l]].lhs_position(),
                         k_ + matches[j + indices[k - 1]].lhs_position(),
                         matches[j].rhs_id(),
                         strand ? matches[j + indices[l]].rhs_position()
                                : matches[j + indices[k - 1]].rhs_position(),
                         k_ + (strand
                                   ? matches[j + indices[k - 1]].rhs_position()
                                   : matches[j + indices[l]].rhs_position()),
                         std::min(lhs_matches, rhs_matches), strand);

        l = k;
      }
    }
  }
  return dst;
}

std::vector<MinimizerEngine::Kmer>
MinimizerEngine::Minimize(const std::unique_ptr<biosoup::NucleicAcid> &sequence,
                          bool minhash) const {
  return Minimize(sequence, 0, sequence->inflated_len, minhash);
}

std::vector<MinimizerEngine::Kmer>
MinimizerEngine::Minimize(const std::unique_ptr<biosoup::NucleicAcid> &sequence,
                          std::uint32_t pos, std::uint32_t len,
                          bool minhash) const {

  if (len < k_) {
    return std::vector<Kmer>{};
  }

  std::uint64_t mask = (1ULL << (k_ * 2)) - 1;

  auto hash = [&](std::uint64_t key) -> std::uint64_t {
    key = ((~key) + (key << 21)) & mask;
    key = key ^ (key >> 24);
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ (key >> 14);
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ (key >> 28);
    key = (key + (key << 31)) & mask;
    return key;
  };

  std::deque<Kmer> window;
  auto window_add = [&](std::uint64_t value, std::uint64_t location) -> void {
    while (!window.empty() && window.back().value > value) {
      window.pop_back();
    }
    window.emplace_back(value, location);
  };
  auto window_update = [&](std::uint32_t position) -> void {
    while (!window.empty() && (window.front().position()) < position) {
      window.pop_front();
    }
  };

  std::uint64_t shift = (k_ - 1) * 2;
  std::uint64_t minimizer = 0;
  std::uint64_t reverse_minimizer = 0;
  std::uint64_t id = static_cast<std::uint64_t>(sequence->id) << 32;
  std::uint64_t is_stored = 1ULL << 63;

  std::vector<Kmer> dst;

  for (std::uint32_t i = 0; i < len; ++i) {
    std::uint64_t c = sequence->Code(pos + i);
    minimizer = ((minimizer << 2) | c) & mask;
    reverse_minimizer = (reverse_minimizer >> 2) | ((c ^ 3) << shift);
    if (i >= k_ - 1U) {
      if (minimizer < reverse_minimizer) {
        window_add(hash(minimizer), (pos + i - (k_ - 1U)) << 1 | 0);
      } else if (minimizer > reverse_minimizer) {
        window_add(hash(reverse_minimizer), (pos + i - (k_ - 1U)) << 1 | 1);
      }
    }
    if (i >= (k_ - 1U) + (w_ - 1U)) {
      for (auto it = window.begin(); it != window.end(); ++it) {
        if (it->value != window.front().value) {
          break;
        }
        if (it->origin & is_stored) {
          continue;
        }
        dst.emplace_back(it->value, id | it->origin);
        it->origin |= is_stored;
      }
      window_update(pos + i - (k_ - 1U) - (w_ - 1U) + 1);
    }
  }

  if (minhash) {
    RadixSort(dst.begin(), dst.end(), k_ * 2, Kmer::SortByValue);
    dst.resize(len / k_);
  }

  return dst;
}

template <typename RandomAccessIterator, typename Compare>
void MinimizerEngine::RadixSort(RandomAccessIterator first,
                                RandomAccessIterator last,
                                std::uint8_t max_bits,
                                Compare comp) { //  unary comparison function
  if (first >= last) {
    return;
  }

  std::vector<typename std::iterator_traits<RandomAccessIterator>::value_type>
      tmp(last - first); // NOLINT
  auto begin = tmp.begin();
  auto end = tmp.end();

  std::uint64_t buckets[0x100]{}; // 256 b
  std::uint8_t shift = 0;
  for (; shift < max_bits; shift += 8) {
    std::uint64_t counts[0x100]{};
    for (auto it = first; it != last; ++it) {
      ++counts[comp(*it) >> shift & 0xFF];
    }
    for (std::uint64_t i = 0, j = 0; i < 0x100; j += counts[i++]) {
      buckets[i] = j;
    }
    for (auto it = first; it != last; ++it) {
      *(begin + buckets[comp(*it) >> shift & 0xFF]++) = *it;
    }
    std::swap(begin, first);
    std::swap(end, last);
  }

  if (shift / 8 & 1) { // copy the sorted array for odd cases
    for (; first != last; ++first, ++begin) {
      *begin = *first;
    }
  }
}

template <typename Compare>
std::vector<std::uint64_t> MinimizerEngine::LongestSubsequence(
    std::vector<Match>::const_iterator first,
    std::vector<Match>::const_iterator last,
    Compare comp) { // binary comparison function
  if (first >= last) {
    return std::vector<std::uint64_t>{};
  }

  std::vector<std::uint64_t> minimal(last - first + 1, 0);
  std::vector<std::uint64_t> predecessor(last - first, 0);

  std::uint64_t longest = 0;
  for (auto it = first; it != last; ++it) {
    std::uint64_t lo = 1, hi = longest;
    while (lo <= hi) {
      std::uint64_t mid = lo + (hi - lo) / 2;
      if ((first + minimal[mid])->lhs_position() < it->lhs_position() &&
          comp((first + minimal[mid])->rhs_position(), it->rhs_position())) {
        lo = mid + 1;
      } else {
        hi = mid - 1;
      }
    }

    predecessor[it - first] = minimal[lo - 1];
    minimal[lo] = it - first;
    longest = std::max(longest, lo);
  }

  std::vector<std::uint64_t> dst;
  for (std::uint64_t i = 0, j = minimal[longest]; i < longest; ++i) {
    dst.emplace_back(j);
    j = predecessor[j];
  }
  std::reverse(dst.begin(), dst.end());

  return dst;
}

} // namespace hifimap
