#include <bits/stdc++.h>
#include <regex>
#include <sys/time.h>

void generate_reads(const std::string &name, const std::string &reference,
                    int TARGET_SIZE, int NO_OF_READS, float PROBABILITY,
                    int GENERATE_COMPLEMENTED, int SEED) {
  std::default_random_engine gene;
  gene.seed(SEED);

  std::uniform_int_distribution<int> position_distribution(
      0, reference.size() - TARGET_SIZE - 1);
  std::uniform_real_distribution<double> error_distribution(0.0, 1.0);
  std::uniform_int_distribution<int> bases_distribution(0, 3);
  std::vector<char> bases = {'A', 'C', 'G', 'T'};

  std::vector<std::string> queries;
  std::vector<int> positions;

  for (int i = 0; i < NO_OF_READS; i++) {
    int position = position_distribution(gene);
    char text[TARGET_SIZE + 1];

    for (int j = 0; j < TARGET_SIZE; j++) {
      if (error_distribution(gene) <= PROBABILITY) {
        text[j] = bases[bases_distribution(gene)];
      } else {
        text[j] = reference[position + j];
      }
    }
    if (GENERATE_COMPLEMENTED && error_distribution(gene) < 0.5) {
      for (int j = 0; j < TARGET_SIZE; j++) {
        switch (text[j]) {
        case 'A':
          text[j] = 'T';
          break;
        case 'C':
          text[j] = 'G';
          break;
        case 'G':
          text[j] = 'C';
          break;
        case 'T':
          text[j] = 'A';
          break;
        }
      }
      for (int j = 0; j < (TARGET_SIZE + 1) / 2; j++) {
        std::swap(text[j], text[TARGET_SIZE - 1 - j]);
      }
    }

    text[TARGET_SIZE] = '\0';
    std::cout << name;
    printf("%d_TARGET_SIZE_%d_POSITION_%09d\n%s\n\n", i, TARGET_SIZE, position,
           text);
  }
}

int main(int argc, char **argv) {
  if (argc != 6 + 1) {
    std::cout << "wrong args, usage: TARGET_SIZE NO_OF_READS "
                 "REVCOMP(0/1) PROBABILITY "
                 "reference_path SEED"
              << std::endl;
    exit(1);
  }

  int TARGET_SIZE = atoi(argv[1]);
  int NO_OF_READS = atoi(argv[2]);
  int GENERATE_COMPLEMENTED = atoi(argv[3]);
  float PROBABILITY = atof(argv[4]);
  char *reference_path = argv[5];
  int SEED = atof(argv[6]);

  std::ifstream ref{reference_path};

  std::string NAME;
  std::string REFERENCE;
  std::getline(ref, NAME);
  std::string str;
  std::regex newlines_re("\n+");

  while (std::getline(ref, str)) {
    if (str[0] == '>')
      break;
    REFERENCE += std::regex_replace(str, newlines_re, "");
  }

  // std::getline(ref, REFERENCE);

  generate_reads(NAME, REFERENCE, TARGET_SIZE, NO_OF_READS, PROBABILITY,
                 GENERATE_COMPLEMENTED, SEED);
}
