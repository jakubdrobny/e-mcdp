#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

template <typename T> class Matrix {
public:
  Matrix(std::vector<int> dimensions);

  T &operator[](const std::vector<int> &index);
  const T &operator[](const std::vector<int> &index) const;

  void print() const;

private:
  std::vector<int> dimensions;
  std::vector<T> data;

  int get_index(const std::vector<int> &index) const;
};

#endif
