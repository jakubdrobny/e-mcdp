#ifndef SEGTREE_H
#define SEGTREE_H

#include <functional>
#include <vector>

template <class T> class SegTree {
public:
  using SegTreeOperation = std::function<T(T, T)>;

  SegTree(int n, SegTreeOperation operation, T neutral_element);
  SegTree(SegTreeOperation operation, T neutral_element, const std::vector<T> &values);

  void set(int idx, T el);
  T query(int l, int r); // returns combined elements from interval [l, r) with the operation provided
private:
  int N;
  SegTreeOperation op;
  std::vector<T> t;
  T neutral_element;
};

#endif // SEGTREE_H
