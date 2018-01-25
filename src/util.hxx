#ifndef _81ba2418_3919_44e6_8ee3_5b23d22881c2
#define _81ba2418_3919_44e6_8ee3_5b23d22881c2

#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <sys/stat.h>

#ifdef __PROGRAMR__
#include "diagnostic.hxx"
#include "env.hxx"
#include "typeclass.hxx"
namespace mota {
  using programr::Say;
  using programr::env;
}
#else
namespace mota {
  // simplified Say grabbed from jdbachan's diagnostic.hxx
  class Say {
    std::stringstream ss;
  public:
    Say() { }
    ~Say() {
      ss << '\n';
      std::cerr << ss.str();
      std::cerr.flush();
    }
    template<class T>
    Say& operator<<(const T &x) {
      ss << x;
      return *this;
    }
  };

  template<class T>
  T c_str_to(const char *str) {
    T val; {
      std::istringstream ss{str};
      ss >> val;
    }
    return val;
  }
  template<>
  inline std::string c_str_to<std::string>(const char *str) {
    return static_cast<std::string>(str);
  }

  template<class T>
  T env(std::string envkey, T defval) {
    const char *str = std::getenv(envkey.c_str());
    return str ? c_str_to<T>(str) : defval;
  }
}

namespace std {
  // grabbed from jdbachan's typeclass.hxx
  template<class A, class B>
  struct hash<pair<A,B>> {
    inline size_t operator()(const pair<A,B> &x) const {
      size_t h = hash<A>()(x.first);
      h ^= h >> 13;
      h *= 41;
      h += hash<B>()(x.second);
      return h;
    }
  };
  template<class T, std::size_t n>
  struct hash<array<T,n>> {
    inline size_t operator()(const array<T,n> &x) const {
      size_t h = 0;
      for(std::size_t i=0; i < n; i++) {
        h ^= h >> 13;
        h *= 41;
        h += hash<T>()(x[i]);
      }
      return h;
    }
  };
}

#endif

namespace mota {

using Int3 = std::array<int,3>;
using Int4 = std::array<int,4>;

const int kNone = -1;

inline int wrap(int x, int m) {
  int r = x % m;
  return r < 0 ? r + m : r;
}
inline int div_ceil(int x, int y) {
  return (x + y-1) / y;
}

// creates map from each element in xs to its index
// U is a iterable container of hashable Ts
template <class T, class U>
std::unordered_map<T,int> make_reverse_map(const U &xs) {
  std::unordered_map<T,int> result;
  int ix = 0;
  for (const auto &x : xs) {
    result[x] = ix++;
  }
  return result;
}

template <class T, size_t N>
std::ostream &operator<<(std::ostream &os, const std::array<T,N> &xs) {
  bool first = true;
  os << "[";
  for (const auto &x : xs) {
    if (!first) os << ", ";
    os << x;
    first = false;
  }
  os << "]";
  return os;
}

template <class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &xs) {
  bool first = true;
  os << "[";
  for (const auto &x : xs) {
    if (!first) os << ", ";
    os << x;
    first = false;
  }
  os << "]";
  return os;
}

template <class K, class V>
std::ostream &operator<<(std::ostream &os, const std::unordered_map<K,V> &xs) {
  bool first = true;
  os << "{";
  for (const auto &p : xs) {
    if (!first) os << ", ";
    os << p.first << " -> " << p.second;
    first = false;
  }
  os << "}";
  return os;
}

template <class T>
std::string to_str(const T &x) {
  std::ostringstream oss;
  oss << x;
  return oss.str();
}

inline bool dir_exists(std::string pathname) {
  struct stat sb;
  return stat(pathname.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode);
}

template <class T>
std::string seq_to_str(const T &xs) {
  std::ostringstream ss;
  bool first = true;
  ss << "(";
  for (auto it = xs.begin(); it != xs.end(); ++it) {
    if (!first) { ss << ","; }
    first = false;
    ss << *it;
  }
  ss << ")";
  return ss.str();
}

}

#endif
