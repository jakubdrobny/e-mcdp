#ifndef MATRIX_H
#define MATRIX_H

#include <stdexcept>
#include <vector>

template <typename T> class Matrix {
public:
  Matrix(size_t rows, size_t cols);
  Matrix(const std::vector<T> &data);
  Matrix(const std::vector<std::vector<T>> &data);

  Matrix<T> &operator()(size_t row);
  const Matrix<T> &operator()(size_t row) const;

  T &operator()(size_t row, size_t col);
  const T &operator()(size_t row, size_t col) const;

  Matrix<T> operator*(const Matrix<T> &other) const;
  Matrix<T> operator^(long long power) const;

  size_t getRowCount() const;
  size_t getColCount() const;

  Matrix<T> &operator=(const std::vector<T> &rowData);
  Matrix<T> &operator=(const Matrix<T> &rowMatrix);

  void print() const;

private:
  size_t rows, cols;
  std::vector<std::vector<T>> data;
};

template class Matrix<int>;
template class Matrix<long double>;

#endif // MATRIX_H
