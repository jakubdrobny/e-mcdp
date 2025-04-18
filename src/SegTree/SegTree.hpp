#ifndef SEGTREE_H
#define SEGTREE_H

#include "../MarkovChain/MarkovChain.hpp"

#include <functional>
#include <vector>

template <class T> class SegTree {
public:
  using SegTreeOperation = std::function<T(T, T, const MarkovChain &)>;

  SegTree();
  SegTree(int n, SegTreeOperation operation, T neutral_element, const MarkovChain &mc);
  SegTree(SegTreeOperation operation, T neutral_element, const std::vector<T> &values, const MarkovChain &mc);

  void set(int idx, T el);
  T query(int l, int r); // returns combined elements from interval [l, r) with the operation provided
private:
  int N;
  SegTreeOperation op;
  std::vector<T> t;
  T neutral_element;
  MarkovChain markov_chain;

  void init(int n);
  void _build(const std::vector<T> &values, int x, int lx, int rx);
  T _query(int l, int r, int x, int lx, int rx);
};

#endif // SEGTREE_H
