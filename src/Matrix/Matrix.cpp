#include "Matrix.h"
#include "../Logger/Logger.h"

template <typename T>
Matrix<T>::Matrix(int rows, int cols) : rows(rows), cols(cols) {
  data.resize(rows * cols);
}

template <typename T> T &Matrix<T>::operator[](int i) { return data[i]; }

template <typename T> const T &Matrix<T>::operator[](int i) const {
  return data[i];
}

template <typename T>
std::vector<T> Matrix<T>::operator[](const std::vector<int> &index) const {
  if (index.size() != 1) {
    logger.error("Invalid index size");
    exit(1);
  }
  int i = index[0];
  if (i < 0 || i >= rows) {
    logger.error("Row index out of range");
    exit(1);
  }
  std::vector<T> row(cols);
  for (int j = 0; j < cols; ++j) {
    row[j] = data[i * cols + j];
  }
  return row;
}

template <typename T> T &Matrix<T>::operator()(int i, int j) {
  return data[i * cols + j];
}

template <typename T> const T &Matrix<T>::operator()(int i, int j) const {
  return data[i * cols + j];
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &other) const {
  if (cols != other.rows) {
    logger.error("Incompatible matrix dimensions for multiplication");
    exit(1);
  }
  Matrix<T> result(rows, other.cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < other.cols; ++j) {
      for (int k = 0; k < cols; ++k) {
        result(i, j) += (*this)(i, k) * other(k, j);
      }
    }
  }
  return result;
}

template <typename T> Matrix<T> Matrix<T>::operator^(int n) const {
  if (rows != cols) {
    logger.error("Exponentiation is only defined for square matrices");
    exit(1);
  }
  Matrix<T> result(rows, cols);
  for (int i = 0; i < rows; ++i) {
    result(i, i) = 1;
  }
  Matrix<T> base = *this;
  while (n > 0) {
    if (n % 2 == 1) {
      result = result * base;
    }
    base = base * base;
    n /= 2;
  }
  return result;
}
