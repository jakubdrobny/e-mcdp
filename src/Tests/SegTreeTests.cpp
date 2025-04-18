#include "../SegTree/SegTree.hpp"
#include <gtest/gtest.h>

TEST(SegTreeBasicOps, Sum) {
  MarkovChain mc;
  SegTree<int> t([](int x, int y, const MarkovChain &mc) { return x + y; }, 0,
                 std::vector<int>{1, 5, 10, -1, 0, 2, -8, -3, 6}, mc);
  ASSERT_EQ(5, t.query(1, 2));
  ASSERT_EQ(-9, t.query(4, 8));
}
