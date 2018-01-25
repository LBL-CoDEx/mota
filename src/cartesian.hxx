#ifndef _0ea7eae1_1ec6_4917_a76e_8c69a177f5df
#define _0ea7eae1_1ec6_4917_a76e_8c69a177f5df

#include <cassert>
#include <array>
#include <ostream>
#include <sstream>

namespace mota {

template <class T, long unsigned N>
class Cartesian {
public:
  Cartesian() = default;
  Cartesian(const Cartesian &) = default;
  Cartesian(Cartesian &&) = default;
  Cartesian &operator=(const Cartesian &other) = default;
  Cartesian &operator=(Cartesian &&other) = default;

  explicit Cartesian(T x) {
    for (unsigned i = 0; i < N; ++i) val_[i] = x;
  }
  Cartesian(T shape_id, Cartesian shape) {
    for (unsigned i = 0; i < N; ++i) {
      val_[i] = shape_id % shape[i];
      shape_id /= shape[i]; // integer division drops remainder
    }
  }

  // implicit conversions

  Cartesian(std::array<T,N> ar) : val_(std::move(ar)) {}
  operator std::array<T,N> () const { return val_; }

  // access

  T& operator[] (unsigned i) { return val_[i]; }
  const T& operator[] (unsigned i) const { return val_[i]; }

  // arithmetic

  Cartesian &operator+=(Cartesian b) {
    for (unsigned i = 0; i < N; ++i) { val_[i] += b[i]; }
    return *this;
  }
  Cartesian &operator-=(Cartesian b) {
    for (unsigned i = 0; i < N; ++i) { val_[i] -= b[i]; }
    return *this;
  }
  Cartesian &operator*=(Cartesian b) {
    for (unsigned i = 0; i < N; ++i) { val_[i] *= b[i]; }
    return *this;
  }
  Cartesian &operator/=(Cartesian b) {
    for (unsigned i = 0; i < N; ++i) { val_[i] /= b[i]; }
    return *this;
  }
  Cartesian &operator%=(Cartesian b) {
    for (unsigned i = 0; i < N; ++i) {
      val_[i] %= b[i];
      if (val_[i] < 0) { val_[i] += b[i]; } // ensure non-negative
    }
    return *this;
  }

  Cartesian operator-() const {
    Cartesian c;
    for (unsigned i = 0; i < N; ++i) { c[i] = -val_[i]; }
    return c;
  }
  Cartesian operator+(Cartesian b) const { Cartesian c{*this}; return (c+=b); }
  Cartesian operator-(Cartesian b) const { Cartesian c{*this}; return (c-=b); }
  Cartesian operator*(Cartesian b) const { Cartesian c{*this}; return (c*=b); }
  Cartesian operator/(Cartesian b) const { Cartesian c{*this}; return (c/=b); }
  Cartesian operator%(Cartesian b) const { Cartesian c{*this}; return (c%=b); }

  // IO

  std::string str() const {
    std::ostringstream oss;
    oss << "(";
    for (unsigned i = 0; i < N-1; ++i) oss << val_[i] << ",";
    oss << val_[N-1] << ")";
    return oss.str();
  }

  friend std::ostream &operator<<(std::ostream &os, Cartesian c) {
    os << c.str();
    return os;
  }

  // comparison

  bool operator==(Cartesian b) const {
    for (unsigned i = 0; i < N; ++i) { if (val_[i] != b[i]) return false; }
    return true;
  }
  bool operator!=(Cartesian b) const {
    return !(*this == b);
  }

  // these have "all" semantics (not "any" semantics)
  //   with these semantics, some pairs of objects are incomparable
  //   also b/c of this, cannot implement e.g. <= as (!>).
  bool operator<(Cartesian b) const {
    for (unsigned i = 0; i < N; ++i) { if (val_[i] >= b[i]) return false; }
    return true;
  }
  bool operator<=(Cartesian b) const {
    for (unsigned i = 0; i < N; ++i) { if (val_[i] > b[i]) return false; }
    return true;
  }
  bool operator>(Cartesian b) const {
    for (unsigned i = 0; i < N; ++i) { if (val_[i] <= b[i]) return false; }
    return true;
  }
  bool operator>=(Cartesian b) const {
    for (unsigned i = 0; i < N; ++i) { if (val_[i] < b[i]) return false; }
    return true;
  }

  // element-wise max and min
  friend Cartesian max(Cartesian a, Cartesian b) {
    Cartesian result {a};
    for (unsigned i = 0; i < N; ++i) { if (a[i] < b[i]) result[i] = b[i]; }
    return result;
  }
  friend Cartesian min(Cartesian a, Cartesian b) {
    Cartesian result {a};
    for (unsigned i = 0; i < N; ++i) { if (a[i] > b[i]) result[i] = b[i]; }
    return result;
  }

  // get the ID of the point within a shape

  T shape_id(Cartesian shape) const {
    assert(val_[N-1] < shape[N-1]);
    T id = val_[N-1];
    for (int i = N-2; i >= 0; --i) {
      assert(val_[i] < shape[i]);
      id = id * shape[i] + val_[i];
    }
    return id;
  }

private:
  std::array<T, N> val_;
};

template <class T, long unsigned N>
std::string str(Cartesian<T,N> x) {
  std::stringstream ss;
  bool first = true;
  ss << "(";
  for (unsigned i = 0; i < N; ++i) {
    if (!first) ss << ",";
    ss << x[i];
    first = false;
  }
  ss << ")";
  return ss.str();
}
template <class T, long unsigned N>
T prod(Cartesian<T,N> x) {
  T result = 1;
  for (unsigned i = 0; i < N; ++i) result *= x[i];
  return result;
}

}

namespace std {
  template <class T, long unsigned N>
  struct hash<mota::Cartesian<T,N>> {
    inline size_t operator()(const mota::Cartesian<T,N> &x) const {
      return hash<std::array<T,N>>()((std::array<T,N>) x);
    }
  };
}

namespace mota {

template <class T, long unsigned N>
class CartesianProd {
public:
  class iterator {
  public:
    iterator(Cartesian<T,N> s, Cartesian<T,N> c, bool v) : shape_(s), cur_(c), valid(v) {}
    iterator &operator++() {
      int i = 0;
      while (i < (int) N && cur_[i] >= shape_[i]-1) ++i;
      if (i < (int) N) {
        ++cur_[i--];
        while (i >= 0) cur_[i--] = 0;
      } else {
        valid = false;
      }
      return *this;
    }
    bool operator==(const iterator &other) {
      if (!valid && !other.valid) return true; // invalid iterators are equal
      if (valid != other.valid) return false;
      int i = N-1;
      while (i >= 0 && cur_[i] == other.cur_[i]) --i;
      return i < 0;
    }
    bool operator!=(const iterator &other) { return !(*this == other); }
    Cartesian<T,N> operator*() { return cur_; }

  private:
    Cartesian<T, N> shape_;
    Cartesian<T, N> cur_;
    bool valid = true;
  };

  CartesianProd (Cartesian<T, N> shape) : shape_(shape) { }
  iterator begin() { return iterator{shape_, Cartesian<T,N>(0), true}; }
  iterator end() { return iterator{shape_, Cartesian<T,N>(0), false}; }

private:
  Cartesian<T, N> shape_;
};

using Car3 = Cartesian<int,3>;
using Car4 = Cartesian<int,4>;

}

#endif
