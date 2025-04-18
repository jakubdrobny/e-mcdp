#include "SegTree.hpp"
#include "../Interval/Section.hpp"

template <class T> SegTree<T>::SegTree() {}

template <class T>
SegTree<T>::SegTree(int n, SegTreeOperation operation, T neutral_element, const MarkovChain &mc)
    : op(operation), neutral_element(neutral_element), markov_chain(mc) {
  this->init(n);
}

template <class T> void SegTree<T>::init(int n) {
  this->N = 1;
  while (this->N < n)
    this->N <<= 1;
  this->t.assign(this->N << 1, this->neutral_element);
}

template <class T> void SegTree<T>::_build(const std::vector<T> &values, int x, int lx, int rx) {
  if (rx - lx == 1) {
    if (lx < (int)values.size())
      this->t[x] = values[lx];
    return;
  }

  int m = (lx + rx) / 2;
  _build(values, 2 * x + 1, lx, m);
  _build(values, 2 * x + 2, m, rx);
  t[x] = this->op(this->t[2 * x + 1], this->t[2 * x + 2], this->markov_chain);
}

template <class T>
SegTree<T>::SegTree(SegTreeOperation operation, T neutral_element, const std::vector<T> &values, const MarkovChain &mc)
    : op(operation), neutral_element(neutral_element), markov_chain(mc) {
  this->init(values.size());

  this->_build(values, 0, 0, this->N);
}

template <class T> T SegTree<T>::_query(int l, int r, int x, int lx, int rx) {
  if (rx <= l || lx >= r)
    return this->neutral_element;
  if (l <= lx && rx <= r)
    return this->t[x];

  int m = (lx + rx) / 2;
  T e1 = this->_query(l, r, 2 * x + 1, lx, m);
  T e2 = this->_query(l, r, 2 * x + 2, m, rx);
  return this->op(e1, e2, this->markov_chain);
}

template <class T> T SegTree<T>::query(int l, int r) { return this->_query(l, r, 0, 0, this->N); }

template class SegTree<int>;
template class SegTree<Section>;
