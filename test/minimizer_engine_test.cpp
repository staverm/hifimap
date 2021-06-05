// Copyright (c) 2021 Mauro Staver

#include "hifimap/minimizer_engine.hpp"

#include "bioparser/fasta_parser.hpp"
#include "gtest/gtest.h"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace hifimap {
namespace test {

class HifimapMinimizerEngineTest : public ::testing::Test {
public:
  void SetUp() override {
    biosoup::NucleicAcid::num_objects = 0;
    auto p =
        bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastaParser>(
            TEST_DATA); // NOLINT
    s = p->Parse(-1);
    EXPECT_EQ(2, s.size());
  }

  std::vector<std::unique_ptr<biosoup::NucleicAcid>> s;
};

// TODO: write tests here
// TEST_F(HifimapMinimizerEngineTest, HifiMap) {}

} // namespace test
} // namespace hifimap
