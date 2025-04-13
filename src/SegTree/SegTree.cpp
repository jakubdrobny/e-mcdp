#include "SegTree.hpp"

template <class T>
SegTree<T>::SegTree(int n, SegTreeOperation operation, T neutral_element)
    : N(n), op(operation), neutral_element(neutral_element) {
  this->t.resize(2 * N);
}

template <class T>
SegTree<T>::SegTree(SegTreeOperation operation, T neutral_element, const std::vector<T> &values)
    : N(values.size()), op(operation), neutral_element(neutral_element) {
  this->t.resize(2 * this->N);

  for (int i = 0; i < this->N; i++) {
    this->t[i + this->N] = values[i];
  }

  for (int i = this->N - 1; i > 0; --i) {
    this->t[i] = this->op(this->t[i << 1], this->t[i << 1 | 1]);
  }
}

template <class T> void SegTree<T>::set(int pos, T el) {
  for (this->t[pos += this->N] = el; pos > 1; pos >>= 1) {
    this->t[pos >> 1] = this->op(this->t[pos], this->t[pos ^ 1]);
  }
}

template <class T> T SegTree<T>::query(int l, int r) {
  T la = this->neutral_element, ra = this->neutral_element;
  for (l += this->N, r += this->N; l < r; l /= 2, r /= 2) {
    if (l & 1)
      la = this->op(la, this->t[l++]);
    if (r & 1)
      ra = this->op(ra, this->t[--r]);
  }
  return this->op(la, ra);
}

template class SegTree<int>;
