#include "Matrix.h"

template <typename T>
Matrix<T>::Matrix(std::vector<int> dimensions) : dimensions(dimensions) {
  int size = 1;
  for (int dim : dimensions) {
    size *= dim;
  }
  data.resize(size);
}

template <typename T> T &Matrix<T>::operator[](const std::vector<int> &index) {
  return data[get_index(index)];
}

template <typename T>
const T &Matrix<T>::operator[](const std::vector<int> &index) const {
  return data[get_index(index)];
}

template <typename T>
int Matrix<T>::get_index(const std::vector<int> &index) const {
  int index_1d = 0;
  int multiplier = 1;
  for (int i = dimensions.size() - 1; i >= 0; --i) {
    index_1d += index[i] * multiplier;
    multiplier *= dimensions[i];
  }
  return index_1d;
}
