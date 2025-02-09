#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

template <typename T> class Matrix {
public:
  Matrix(int rows, int cols);

  T &operator[](int i);
  const T &operator[](int i) const;

  std::vector<T> operator[](const std::vector<int> &index) const;

  T &operator()(int i, int j);
  const T &operator()(int i, int j) const;

  Matrix<T> operator*(const Matrix<T> &other) const;

  Matrix<T> operator^(int n) const;

private:
  int rows;
  int cols;
  std::vector<T> data;
};

#endif
