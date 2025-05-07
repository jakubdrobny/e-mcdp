#include "SegTree.hpp"
#include "../Interval/Section.hpp"
#include <iostream>
#include <queue>

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

template <class T> void SegTree<T>::dump() {
  if (this->t.empty()) {
    std::cout << "DUMPING SEGTREE (EMPTY)\n=========\n";
    std::cout << "DUMP FINISHED\n";
    return;
  }

  std::queue<std::tuple<int, int, int>> q_bfs;
  q_bfs.push({0, 0, this->N});

  std::cout << "DUMPING SEGTREE\n=========\n";

  int nodes_in_current_level = 1;
  int nodes_in_next_level = 0;

  while (!q_bfs.empty()) {
    auto [x, lx, rx] = q_bfs.front();
    q_bfs.pop();

    if (x >= (int)this->t.size()) {
      std::cerr << "SegTree::dump() Error: Attempted to access invalid index " << x << " (t.size() = " << this->t.size()
                << "). Skipping." << std::endl;
      std::cout << "[INVALID_IDX:" << x << "] ";
      nodes_in_current_level--;
      if (nodes_in_current_level == 0) {
        std::cout << "\n";
        nodes_in_current_level = nodes_in_next_level;
        nodes_in_next_level = 0;
      }
      continue;
    }

    std::cout << this->t[x] << " (idx:" << x << " [" << lx << "," << rx << ")) ";

    nodes_in_current_level--;

    if (rx - lx > 1) {
      int m = lx + (rx - lx) / 2;
      int left_child_idx = 2 * x + 1;
      int right_child_idx = 2 * x + 2;

      if (right_child_idx < (int)this->t.size()) {
        q_bfs.push({left_child_idx, lx, m});
        nodes_in_next_level++;
        q_bfs.push({right_child_idx, m, rx});
        nodes_in_next_level++;
      } else {
        if (left_child_idx < (int)this->t.size()) {
          q_bfs.push({left_child_idx, lx, m});
          nodes_in_next_level++;
          std::cerr << "SegTree::dump() Warning: Node " << x << " (range [" << lx << "," << rx << ")) has left child "
                    << left_child_idx << " but right child " << right_child_idx
                    << " is out of bounds (t.size()=" << this->t.size() << ").\n";
        } else {
          std::cerr << "SegTree::dump() Warning: Internal node " << x << " (range [" << lx << "," << rx
                    << ")) has children indices " << left_child_idx << ", " << right_child_idx
                    << " which are out of bounds (t.size()=" << this->t.size() << ").\n";
        }
      }
    }

    if (nodes_in_current_level == 0) {
      std::cout << "\n";
      nodes_in_current_level = nodes_in_next_level;
      nodes_in_next_level = 0;
    }
  }

  std::cout << "DUMP FINISHED\n";
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
