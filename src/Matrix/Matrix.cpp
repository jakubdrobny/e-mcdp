#include "Matrix.h"
#include "../Logger/Logger.h"
#include <iostream>

template <typename T>
Matrix<T>::Matrix(size_t rows, size_t cols) : rows(rows), cols(cols) {
  data.resize(rows, std::vector<T>(cols, T()));
}

template <typename T> Matrix<T>::Matrix(const std::vector<T> &data) {
  if (data.empty()) {
    logger.error("Matrix data cannot be empty.");
    exit(1);
  }
  this->rows = 1;
  this->cols = data.size();
  this->data = std::vector<std::vector<T>>(1, data);
}

template <typename T>
Matrix<T>::Matrix(const std::vector<std::vector<T>> &data) {
  if (data.empty() || data[0].empty()) {
    logger.error("Matrix data cannot be empty.");
    exit(1);
  }
  this->rows = data.size();
  this->cols = data[0].size();
  this->data = data;
}

template <typename T> Matrix<T> &Matrix<T>::operator()(size_t row) {
  if (row >= rows) {
    logger.error("Matrix indices out of range.");
    exit(1);
  }
  return Matrix(data[row]);
}

template <typename T> const Matrix<T> &Matrix<T>::operator()(size_t row) const {
  if (row >= rows) {
    logger.error("Matrix indices out of range.");
    exit(1);
  }
  return Matrix(data[row]);
}

template <typename T> T &Matrix<T>::operator()(size_t row, size_t col) {
  if (row >= rows || col >= cols) {
    logger.error("Matrix indices out of range.");
    exit(1);
  }
  return data[row][col];
}

template <typename T>
const T &Matrix<T>::operator()(size_t row, size_t col) const {
  if (row >= rows || col >= cols) {
    logger.error("Matrix indices out of range.");
    exit(1);
  }
  return data[row][col];
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &other) const {
  if (cols != other.rows) {
    logger.error("Matrix dimensions do not match for multiplication.");
    exit(1);
  }

  Matrix<T> result(rows, other.cols);
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < other.cols; ++j) {
      for (size_t k = 0; k < cols; ++k) {
        result(i, j) += data[i][k] * other(k, j);
      }
    }
  }
  return result;
}

template <typename T> Matrix<T> Matrix<T>::operator^(long long power) const {
  if (rows != cols) {
    logger.error("Matrix must be square for exponentiation.");
    exit(1);
  }

  Matrix<T> result(rows, cols);
  for (size_t i = 0; i < rows; ++i) {
    result(i, i) = 1;
  }

  Matrix<T> base = *this;

  while (power > 0) {
    if (power % 2 == 1)
      result = result * base;
    power /= 2;
    base = base * base;
  }

  return result;
}

template <typename T> size_t Matrix<T>::getRowCount() const { return rows; }
template <typename T> size_t Matrix<T>::getColCount() const { return cols; }

template <typename T>
Matrix<T> &Matrix<T>::operator=(const std::vector<T> &rowData) {
  if (rowData.size() != cols) {
    logger.error("Row data size does not match matrix columns.");
    exit(1);
  }
  if (rows == 0) {
    rows = 1;
    data.push_back(rowData);
  } else {
    data[0] = rowData;
  }
  return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &rowMatrix) {
  if (rowMatrix.getRowCount() != 1 || rowMatrix.getColCount() != cols) {
    logger.error("Row matrix dimensions do not match.");
    exit(1);
  }
  if (rows == 0) {
    rows = 1;
    data.push_back(rowMatrix(0));
  } else {
    data[0] = rowMatrix(0);
  }
  return *this;
}

template <typename T> void Matrix<T>::print() const {
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      std::cout << data[i][j] << " ";
    }
    std::cout << std::endl;
  }
}
