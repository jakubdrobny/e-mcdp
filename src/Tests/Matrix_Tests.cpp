#include "../Matrix/Matrix.h"
#include <gtest/gtest.h>

TEST(MatrixTest, Constructor) {
  Matrix<int> m1(3, 4);
  EXPECT_EQ(m1.getRowCount(), 3);
  EXPECT_EQ(m1.getColCount(), 4);

  std::vector<std::vector<int>> data = {{1, 2, 3}};
  Matrix<int> m2(data);
  EXPECT_EQ(m2.getRowCount(), 1);
  EXPECT_EQ(m2.getColCount(), 3);

  std::vector<std::vector<int>> data2 = {{1, 2, 3}, {4, 5, 6}};
  Matrix<int> m3(data2);
  EXPECT_EQ(m3.getRowCount(), 2);
  EXPECT_EQ(m3.getColCount(), 3);
}

TEST(MatrixTest, OperatorParenthesesOneIndex) {
  Matrix<int> m(2, 2);
  m(0) = std::vector<int>{1, 2};
  m(1) = std::vector<int>{1, 2};

  EXPECT_EQ(m(0, 0), 1);
  EXPECT_EQ(m(0, 1), 2);
  EXPECT_EQ(m(1, 0), 1);
  EXPECT_EQ(m(1, 1), 2);

  const Matrix<int> &cm = m;
  EXPECT_EQ(cm(0, 0), 1);
  EXPECT_EQ(cm(0, 1), 2);
  EXPECT_EQ(cm(1, 0), 1);
  EXPECT_EQ(cm(1, 1), 2);
}

TEST(MatrixTest, OperatorParenthesesTwoIndices) {
  Matrix<int> m(2, 2);
  m(0, 0) = 1;
  m(0, 1) = 2;
  m(1, 0) = 3;
  m(1, 1) = 4;

  EXPECT_EQ(m(0, 0), 1);
  EXPECT_EQ(m(0, 1), 2);
  EXPECT_EQ(m(1, 0), 3);
  EXPECT_EQ(m(1, 1), 4);

  const Matrix<int> &cm = m;
  EXPECT_EQ(cm(0, 0), 1);
  EXPECT_EQ(cm(0, 1), 2);
  EXPECT_EQ(cm(1, 0), 3);
  EXPECT_EQ(cm(1, 1), 4);
}

TEST(MatrixTest, OperatorMultiplication) {
  std::vector<std::vector<int>> data1 = {{1, 2}, {3, 4}};
  std::vector<std::vector<int>> data2 = {{5, 6}, {7, 8}};
  Matrix<int> m1(data1);
  Matrix<int> m2(data2);

  Matrix<int> result = m1 * m2;
  EXPECT_EQ(result(0, 0), 19);
  EXPECT_EQ(result(0, 1), 22);
  EXPECT_EQ(result(1, 0), 43);
  EXPECT_EQ(result(1, 1), 50);
}

TEST(MatrixTest, OperatorExponentiation) {
  std::vector<std::vector<int>> data = {{1, 1}, {1, 0}};
  Matrix<int> m(data);

  Matrix<int> result = m ^ 3;
  EXPECT_EQ(result(0, 0), 3);
  EXPECT_EQ(result(0, 1), 2);
  EXPECT_EQ(result(1, 0), 2);
  EXPECT_EQ(result(1, 1), 1);
}

TEST(MatrixTest, GetRowAndColCount) {
  Matrix<int> m(5, 6);
  EXPECT_EQ(m.getRowCount(), 5);
  EXPECT_EQ(m.getColCount(), 6);
}

TEST(MatrixTest, OperatorAssignmentVector) {
  Matrix<int> m(2, 2);
  std::vector<int> row = {5, 6};
  m = row;
  EXPECT_EQ(m(0, 0), 5);
  EXPECT_EQ(m(0, 1), 6);
}

TEST(MatrixTest, OperatorAssignmentMatrix) {
  Matrix<int> m1(2, 2);
  std::vector<std::vector<int>> data = {{5, 6}};
  Matrix<int> m2(data);
  m1 = m2;
  EXPECT_EQ(m1(0, 0), 5);
  EXPECT_EQ(m1(0, 1), 6);
}
