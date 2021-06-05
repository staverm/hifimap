// Copyright (c) 2021 Mauro Staver
// Some of the code is taken from https://github.com/lbcb-sci/ram

#include <getopt.h>

#include <bitset>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/progress_bar.hpp"
#include "biosoup/timer.hpp"

#include "hifimap/minimizer_engine.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

#if defined(_WIN32)
#include <psapi.h>
#include <windows.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) ||                 \
    (defined(__APPLE__) && defined(__MACH__))
#include <sys/resource.h>
#include <unistd.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) ||                              \
    (defined(__sun__) || defined(__sun) ||                                     \
     defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) ||              \
    defined(__gnu_linux__)
#include <stdio.h>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif

/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t getPeakRSS() {
#if defined(_WIN32)
  /* Windows -------------------------------------------------- */
  PROCESS_MEMORY_COUNTERS info;
  GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
  return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) ||                              \
    (defined(__sun__) || defined(__sun) ||                                     \
     defined(sun) && (defined(__SVR4) || defined(__svr4__)))
  /* AIX and Solaris ------------------------------------------ */
  struct psinfo psinfo;
  int fd = -1;
  if ((fd = open("/proc/self/psinfo", O_RDONLY)) == -1)
    return (size_t)0L; /* Can't open? */
  if (read(fd, &psinfo, sizeof(psinfo)) != sizeof(psinfo)) {
    close(fd);
    return (size_t)0L; /* Can't read? */
  }
  close(fd);
  return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) ||                 \
    (defined(__APPLE__) && defined(__MACH__))
  /* BSD, Linux, and OSX -------------------------------------- */
  struct rusage rusage;
  getrusage(RUSAGE_SELF, &rusage);
#if defined(__APPLE__) && defined(__MACH__)
  return (size_t)rusage.ru_maxrss;
#else
  return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
  /* Unknown OS ----------------------------------------------- */
  return (size_t)0L; /* Unsupported. */
#endif
}

namespace {

static struct option options[] = {
    {"kmer-length", required_argument, nullptr, 'k'},
    {"window-length", required_argument, nullptr, 'w'},
    {"tfidf-offset", required_argument, nullptr, 'f'},
    {"avg-len", required_argument, nullptr, 'l'},
    {"max-regions", required_argument, nullptr, 'r'},
    {"min-hits", required_argument, nullptr, 'H'},
    {"approx", no_argument, nullptr, 'a'},
    {"overhead", no_argument, nullptr, 'o'},
    {"bandwidth", required_argument, nullptr, 'b'},
    {"chain", required_argument, nullptr, 'c'},
    {"matches", required_argument, nullptr, 'm'},
    {"gap", required_argument, nullptr, 'g'},
    {"minhash", no_argument, nullptr, 'M'},
    {"threads", required_argument, nullptr, 't'},
    {"version", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0}};

std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>>
CreateParser(const std::string &path) {
  auto is_suffix = [](const std::string &s, const std::string &suff) {
    return s.size() < suff.size()
               ? false
               : s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
  };

  if (is_suffix(path, ".fasta") || is_suffix(path, ".fasta.gz") ||
      is_suffix(path, ".fna") || is_suffix(path, ".fna.gz") ||
      is_suffix(path, ".fa") || is_suffix(path, ".fa.gz")) {
    try {
      return bioparser::Parser<biosoup::NucleicAcid>::Create<
          bioparser::FastaParser>(path); // NOLINT
    } catch (const std::invalid_argument &exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }
  if (is_suffix(path, ".fastq") || is_suffix(path, ".fastq.gz") ||
      is_suffix(path, ".fq") || is_suffix(path, ".fq.gz")) {
    try {
      return bioparser::Parser<biosoup::NucleicAcid>::Create<
          bioparser::FastqParser>(path); // NOLINT
    } catch (const std::invalid_argument &exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }

  std::cerr << "[hifimap::CreateParser] error: file " << path
            << " has unsupported format extension (valid extensions: .fasta, "
            << ".fasta.gz, .fna, .fna.gz, .fa, .fa.gz, .fastq, .fastq.gz, "
            << ".fq, .fq.gz)" << std::endl;
  return nullptr;
}

void Help() {
  std::cout
      << "usage: hifimap [options ...] <target> [<sequences>]\n"
         "\n"
         "  # default output is stdout\n"
         "  <target>/<sequences> \n"
         "    input file in FASTA/FASTQ format (can be compressed with gzip)\n"
         "\n"
         "  options:\n"
         "    -k, --kmer-length <int>\n"
         "      default: 15\n"
         "      length of minimizers\n"
         "    -w, --window-length <int>\n"
         "      default: 5\n"
         "      length of sliding window from which minimizers are sampled\n"
         "    -f, --tfidf-offset <float>\n"
         "      default: -0.0001\n"
         "      add offset to tfidf limit\n"
         "    -l, --avg-len <int>\n"
         "      default: 20000\n"
         "      average read length or half of region length\n"
         "    -r, --max-regions <int>\n"
         "      default: 1\n"
         "      maximal number of regions on which to perform chaining\n"
         "    -H, --min-hits <int>\n"
         "      default: 10\n"
         "      minimal number of region hits\n"
         "    -a, --approx\n"
         "      approximate mapping positions based on region hits\n"
         "    -o, --overhead\n"
         "      increase time overhead but save space\n"
         "    --bandwidth <int>\n"
         "      default: 500\n"
         "      size of bandwidth in which minimizer hits can be chained\n"
         "    --chain <int>\n"
         "      default: 4\n"
         "      minimal number of chained minimizer hits in overlap\n"
         "    --matches <int>\n"
         "      default: 100\n"
         "      minimal number of matching bases in overlap\n"
         "    --gap <int>\n"
         "      default: 10000\n"
         "      maximal gap between minimizer hits in a chain\n"
         "    --minhash\n"
         "      use only a portion of all minimizers\n"
         "    -t, --threads <int>\n"
         "      default: 4\n"
         "      number of threads\n"
         "    --version\n"
         "      prints the version number\n"
         "    -h, --help\n"
         "      prints the usage\n";
}

} // namespace

void printStats(hifimap::MinimizerEngine &minimizer_engine) {
  double sum = minimizer_engine.minimize_query + minimizer_engine.index_lookup +
               minimizer_engine.find_hits + minimizer_engine.minimize_region +
               minimizer_engine.find_matches + minimizer_engine.chain_matches +
               minimizer_engine.best_regions_;
  std::cerr << " Minimize query: " << std::fixed << std::setprecision(2)
            << minimizer_engine.minimize_query / sum * 100
            << "% Index lookup: " << minimizer_engine.index_lookup / sum * 100
            << "% Find hits: " << minimizer_engine.find_hits / sum * 100
            << "% Best regions: " << minimizer_engine.best_regions_ / sum * 100
            << "% Minimize region: "
            << minimizer_engine.minimize_region / sum * 100
            << "% Find matches: " << minimizer_engine.find_matches / sum * 100
            << "% Chain matches: " << minimizer_engine.chain_matches / sum * 100
            << '%' << std::endl;
}

int main(int argc, char **argv) {
  std::uint32_t k = 15;
  std::uint32_t w = 5;
  std::uint32_t bandwidth = 500;
  std::uint32_t chain = 4;
  std::uint32_t matches = 100;
  std::uint32_t gap = 10000;
  bool minhash = false;
  std::uint32_t num_threads = 4;
  std::uint32_t avg_len = 20000;
  std::uint32_t max_regions = 1;
  std::uint32_t min_hits = 10;
  bool approx = false;
  std::uint32_t fast = true;
  double tfidf_offset = -0.0001;

  std::vector<std::string> input_paths;

  // clang-format off
  const char *optstr = "k:w:l:r:H:f:b:c:g:t:Maovh";
  char arg;
  while ((arg = getopt_long(argc, argv, optstr, options, nullptr)) != -1) {
    switch (arg) {
    case 'k': k = std::atoi(optarg); break;
    case 'w': w = std::atoi(optarg); break;
    case 'l': avg_len = std::atoi(optarg); break;
    case 'r': max_regions = std::atoi(optarg); break;
    case 'b': bandwidth = std::atoi(optarg); break;
    case 'c': chain = std::atoi(optarg); break;
    case 'm': matches = std::atoi(optarg); break;
    case 'g': gap = std::atoi(optarg); break;
    case 'H': min_hits = std::atoi(optarg); break;
    case 'f': tfidf_offset = std::atof(optarg); break;
    case 'M': minhash = true; break;
    case 'a': approx = true; break;
    case 'o': fast = false; break;
    case 't': num_threads = std::atoi(optarg); break;
    case 'v': std::cout << VERSION << std::endl; return 0;
    case 'h': Help(); return 0;
    default: return 1;
    }
  }
  // clang-format on

  if (argc == 1) {
    Help();
    return 0;
  }

  if (approx)
    fast = false;

  for (auto i = optind; i < argc; ++i) {
    input_paths.emplace_back(argv[i]);
  }

  if (input_paths.empty()) {
    std::cerr << "[hifimap::] error: missing target file" << std::endl;
    return 1;
  }

  auto tparser = CreateParser(input_paths[0]);
  if (tparser == nullptr) {
    return 1;
  }

  bool is_ava = false;
  std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> sparser = nullptr;
  if (input_paths.size() > 1) {
    sparser = CreateParser(input_paths[1]);
    if (sparser == nullptr) {
      return 1;
    }
    is_ava = input_paths[0] == input_paths[1];
  } else {
    sparser = CreateParser(input_paths[0]);
    is_ava = true;
  }

  if (is_ava) {
    std::cerr << "[hifimap::] error: ava not supported" << std::endl;
    return 1;
  }

  auto thread_pool = std::make_shared<thread_pool::ThreadPool>(num_threads);
  hifimap::MinimizerEngine minimizer_engine{
      thread_pool, k,   w,       bandwidth,   chain,
      matches,     gap, avg_len, max_regions, min_hits};

  biosoup::Timer timer{};

  while (true) {
    timer.Start();

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> targets;
    try {
      targets = tparser->Parse(1ULL << 32);
    } catch (std::invalid_argument &exception) {
      std::cerr << exception.what() << std::endl;
      return 1;
    }

    if (targets.empty()) {
      break;
    }

    std::cerr << "[hifimap::] parsed " << targets.size() << " targets "
              << std::fixed << timer.Stop() << "s" << std::endl;

    timer.Start();

    minimizer_engine.CreateIndex(targets.begin(), targets.end(), tfidf_offset,
                                 minhash, fast);

    std::cerr << "[hifimap::] minimized targets " << std::fixed << timer.Stop()
              << "s" << std::endl;

    std::uint64_t num_targets = biosoup::NucleicAcid::num_objects;
    biosoup::NucleicAcid::num_objects = 0;

    while (true) {
      timer.Start();

      std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
      try {
        sequences = sparser->Parse(1U << 30);
      } catch (std::invalid_argument &exception) {
        std::cerr << exception.what() << std::endl;
        return 1;
      }

      if (sequences.empty()) {
        break;
      }

      std::vector<std::future<std::vector<biosoup::Overlap>>> futures;
      for (const auto &it : sequences) {
        if (is_ava && it->id >= num_targets) {
          continue;
        }
        futures.emplace_back(thread_pool->Submit(
            [&](const std::unique_ptr<biosoup::NucleicAcid> &sequence)
                -> std::vector<biosoup::Overlap> {
              return minimizer_engine.HifiMap(sequence, minhash, approx);
            },
            std::ref(it)));
      }

      biosoup::ProgressBar bar{static_cast<std::uint32_t>(futures.size()), 16};

      std::uint64_t rhs_offset = targets.front()->id;
      std::uint64_t lhs_offset = sequences.front()->id;
      for (auto &it : futures) {
        for (const auto &jt : it.get()) {
          std::cout << sequences[jt.lhs_id - lhs_offset]->name << "\t"
                    << sequences[jt.lhs_id - lhs_offset]->inflated_len << "\t"
                    << jt.lhs_begin << "\t" << jt.lhs_end << "\t"
                    << (jt.strand ? "+" : "-") << "\t"
                    << targets[jt.rhs_id - rhs_offset]->name << "\t"
                    << targets[jt.rhs_id - rhs_offset]->inflated_len << "\t"
                    << jt.rhs_begin << "\t" << jt.rhs_end << "\t" << jt.score
                    << "\t"
                    << std::max(jt.lhs_end - jt.lhs_begin,
                                jt.rhs_end - jt.rhs_begin)
                    << "\t" << 255 << std::endl;
        }

        if (++bar) {
          std::cerr << "[hifimap::] mapped " << bar.event_counter()
                    << " sequences "
                    << "[" << bar << "] " << std::fixed << timer.Lap() << "s"
                    << "\r";
        }
      }
      std::cerr << std::endl;
      timer.Stop();

      printStats(minimizer_engine);
      if (is_ava && biosoup::NucleicAcid::num_objects >= num_targets) {
        break;
      }
    }

    sparser->Reset();
    biosoup::NucleicAcid::num_objects = num_targets;
  }

  std::cerr << "[hifimap::] " << timer.elapsed_time() << "s" << std::endl;
  size_t peakSize = getPeakRSS();
  std::cerr << "[hifimap::] peak RSS: " << peakSize / 1000000000.0 << " GB"
            << std::endl;

  printStats(minimizer_engine);

  return 0;
}
