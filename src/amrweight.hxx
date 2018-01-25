#ifndef AMRWEIGHT_HXX
#define AMRWEIGHT_HXX

#include <array>
#include <vector>
#include <functional>
#include <iomanip>
#include <cmath>
#include <ostream>
#include <limits>

namespace mota {

// N is the number of AMR levels
// weight/capacity is a N+1 tuple, where the last element is the aggregate
// aggregate represents memory constraint
// per-level weights represent compute load balance constraint
template <unsigned N>
struct AMRWeight {
  // array of weights/capacities listed per-level
  // extra last element is the aggregate weight/capacity (not necessarily the sum)
  std::array<double, N+1> weights;

  AMRWeight() = default;
  AMRWeight(const AMRWeight &) = default;
  AMRWeight(AMRWeight &&) = default;
  AMRWeight &operator=(const AMRWeight &) = default;
  AMRWeight &operator=(AMRWeight &&) = default;

  // set one level and the aggregate to specified value
  AMRWeight(double x, int level) : AMRWeight(0.0) { weights[level] = weights[N] = x; }

  // implicit conversion from double
  AMRWeight(double x) { for (auto &w : weights) { w = x; } }

  // casting to double just returns aggregate entry
  explicit operator double() const { return weights[N]; }

  AMRWeight &operator+=(AMRWeight other) {
    for (unsigned i = 0; i < N+1; ++i) { weights[i] += other.weights[i]; }
    return *this;
  }
  AMRWeight operator-() const { // unary operator
    AMRWeight result;
    for (unsigned i = 0; i < N+1; ++i) { result.weights[i] = -weights[i]; }
    return result;
  }
  AMRWeight &operator-=(AMRWeight other) {
    (*this) += (-other);
    return *this;
  }
  AMRWeight &operator*=(AMRWeight other) {
    for (unsigned i = 0; i < N+1; ++i) { weights[i] *= other.weights[i]; }
    return *this;
  }
  AMRWeight &operator/=(AMRWeight other) {
    for (unsigned i = 0; i < N+1; ++i) { weights[i] /= other.weights[i]; }
    return *this;
  }

  AMRWeight operator+(AMRWeight other) const {
    return AMRWeight{*this} += other;
  }
  AMRWeight operator-(AMRWeight other) const { // binary operator
    return AMRWeight{*this} -= other;
  }
  AMRWeight operator*(AMRWeight other) const {
    return AMRWeight{*this} *= other;
  }
  AMRWeight operator/(AMRWeight other) const {
    return AMRWeight{*this} /= other;
  }

  // comparison operators
  bool operator==(AMRWeight other) const {
    bool result {true};
    for (unsigned i = 0; i < N+1; ++i) {
      if (weights[i] != other.weights[i]) result = false;
    }
    return result;
  }
  bool operator<(AMRWeight other) const {
    bool result {true};
    for (unsigned i = 0; i < N+1; ++i) {
      if (weights[i] >= other.weights[i]) result = false;
    }
    return result;
  }
  bool operator>(AMRWeight other) const {
    bool result {true};
    for (unsigned i = 0; i < N+1; ++i) {
      if (weights[i] <= other.weights[i]) result = false;
    }
    return result;
  }
  bool operator<=(AMRWeight other) const {
    bool result {true};
    for (unsigned i = 0; i < N+1; ++i) {
      if (weights[i] > other.weights[i]) result = false;
    }
    return result;
  }
  bool operator>=(AMRWeight other) const {
    bool result {true};
    for (unsigned i = 0; i < N+1; ++i) {
      if (weights[i] < other.weights[i]) result = false;
    }
    return result;
  }

  AMRWeight sqrt() const {
    AMRWeight result;
    for (unsigned i = 0; i < N+1; ++i) { result.weights[i] = std::sqrt(weights[i]); }
    return result;
  }
  friend AMRWeight sqrt(AMRWeight x) {
    return x.sqrt();
  }

  friend AMRWeight min(AMRWeight a, AMRWeight b) {
    AMRWeight result;
    for (unsigned i = 0; i < N+1; ++i) {
      result.weights[i] = std::min(a.weights[i], b.weights[i]);
    }
    return result;
  }

  friend AMRWeight max(AMRWeight a, AMRWeight b) {
    AMRWeight result;
    for (unsigned i = 0; i < N+1; ++i) {
      result.weights[i] = std::max(a.weights[i], b.weights[i]);
    }
    return result;
  }

  friend std::ostream &operator<<(std::ostream &os, AMRWeight w) {
    bool first = true;
    os << "(";
    for (unsigned i = 0; i < N+1; ++i) {
      if (!first) os << ",";
      os << std::setw(12) << w.weights[i];
      first = false;
    }
    os << ")";
    return os;
  }

  friend AMRWeight operator*(double d, AMRWeight w) {
    return AMRWeight{w} * d;
  }
};

template <class W>
struct MetisCons;

template <>
template <unsigned N>
struct MetisCons<AMRWeight<N>> {
  // NOTE: use this to invoke METIS multi-constraint partitioner
  // set of N+1 constraints, one for each AMR level and one for aggregate
  static std::vector<std::function<double(AMRWeight<N>)>> cons(bool flag_multicons) {
    if (flag_multicons) {
      std::vector<std::function<double(AMRWeight<N>)>> result;
      result.reserve(N+1);
      for (unsigned i = 0; i < N+1; ++i) {
        result.push_back( [=] (AMRWeight<N> x) { return x.weights[i]; } );
      }
      return result;
    } else {
      return std::vector<std::function<double(AMRWeight<N>)>> {
        [] (AMRWeight<N> x) { return static_cast<double>(x); }
      };
    }
  }
};

}

namespace std {
  // std specializations/overloads for AMRWeight
//template <unsigned N>
//AMRWeight<N> numeric_limits<AMRWeight<N>>::min() {
//  return AMRWeight<N> {std::numeric_limits<double>::min()};
//}
  template <>
  inline mota::AMRWeight<3> numeric_limits<mota::AMRWeight<3>>::min() {
    return mota::AMRWeight<3> {std::numeric_limits<double>::min()};
  }
  template <>
  inline mota::AMRWeight<3> numeric_limits<mota::AMRWeight<3>>::max() {
    return mota::AMRWeight<3> {std::numeric_limits<double>::max()};
  }
}

#endif // AMRWEIGHT_HXX
