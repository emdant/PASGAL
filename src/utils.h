
#ifndef UTILS_H
#define UTILS_H

#include <climits>
#include <cstring>
#include <iostream>
#include <random>

template <typename ET>
inline bool atomic_compare_and_swap(ET *a, ET oldval, ET newval) {
  static_assert(sizeof(ET) <= 8, "Bad CAS length");
  if constexpr (sizeof(ET) == 1) {
    uint8_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<uint8_t *>(a), r_oval,
                                        r_nval);
  } else if constexpr (sizeof(ET) == 4) {
    uint32_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<uint32_t *>(a), r_oval,
                                        r_nval);
  } else if constexpr (sizeof(ET) == 8) {
    uint64_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<uint64_t *>(a), r_oval,
                                        r_nval);
  } else {
    std::cerr << "Bad CAS length" << std::endl;
  }
}

template <class ET>
inline bool compare_and_swap(std::atomic<ET> *a, ET oldval, ET newval) {
  return a->load(std::memory_order_relaxed) == oldval &&
         atomic_compare_exchange_weak(a, &oldval, newval);
}

template <class ET> inline bool compare_and_swap(ET *a, ET oldval, ET newval) {
  return (*a) == oldval && atomic_compare_and_swap(a, oldval, newval);
}

template <typename E, typename EV> inline E fetch_and_add(E *a, EV b) {
  volatile E newV, oldV;
  do {
    oldV = *a;
    newV = oldV + b;
  } while (!atomic_compare_and_swap(a, oldV, newV));
  return oldV;
}

template <typename E, typename EV> inline void write_add(E *a, EV b) {
  // volatile E newV, oldV;
  E newV, oldV;
  do {
    oldV = *a;
    newV = oldV + b;
  } while (!atomic_compare_and_swap(a, oldV, newV));
}

template <typename ET, typename F = std::less<ET>>
inline bool write_min(ET *a, ET b, F less = {}) {
  ET c;
  bool r = 0;
  do
    c = *a;
  while (less(b, c) && !(r = atomic_compare_and_swap(a, c, b)));
  return r;
}

template <typename ET, typename F = std::less<ET>>
inline bool write_max(ET *a, ET b, F less = {}) {
  ET c;
  bool r = 0;
  do
    c = *a;
  while (less(c, b) && !(r = atomic_compare_and_swap(a, c, b)));
  return r;
}

template <typename NodeID_, typename rng_t_,
          typename uNodeID_ = typename std::make_unsigned<NodeID_>::type>
class UniDist {
public:
  UniDist(NodeID_ max_value, rng_t_ &rng) : rng_(rng) {
    no_mod_ = rng_.max() == static_cast<uNodeID_>(max_value);
    mod_ = max_value + 1;
    uNodeID_ remainder_sub_1 = rng_.max() % mod_;
    if (remainder_sub_1 == mod_ - 1)
      cutoff_ = 0;
    else
      cutoff_ = rng_.max() - remainder_sub_1;
  }

  NodeID_ operator()() {
    uNodeID_ rand_num = rng_();
    if (no_mod_)
      return rand_num;
    if (cutoff_ != 0) {
      while (rand_num >= cutoff_)
        rand_num = rng_();
    }
    return rand_num % mod_;
  }

private:
  rng_t_ &rng_;
  bool no_mod_;
  uNodeID_ mod_;
  uNodeID_ cutoff_;
};

template <typename ValueT_> class VectorReader {
  std::string filename_;

public:
  explicit VectorReader(std::string filename) : filename_(filename) {
    if (filename == "") {
      std::cout << "No sources filename given (Use -h for help)" << std::endl;
      std::exit(-8);
    }
  }

  std::vector<ValueT_> Read() {
    std::ifstream file(filename_, std::ios::binary);
    if (!file.is_open()) {
      std::cout << "Couldn't open file " << filename_ << std::endl;
      std::exit(-2);
    }

    std::vector<ValueT_> sources;
    while (!file.eof()) {
      ValueT_ source;
      file >> source;
      sources.push_back(source);
    }
    file.close();

    return sources;
  }

  std::vector<ValueT_> ReadSerialized() {
    std::ifstream file(filename_, std::ios::binary);
    if (!file.is_open()) {
      std::cout << "Couldn't open file " << filename_ << std::endl;
      std::exit(-2);
    }

    std::vector<ValueT_> values;
    int64_t num_values; // must be 64-bit value
    file.read(reinterpret_cast<char *>(&num_values), sizeof(num_values));

    values.resize(num_values);
    file.read(reinterpret_cast<char *>(values.data()),
              num_values * sizeof(ValueT_));
    file.close();

    return values;
  }
};

template <typename GraphT_, typename NodeID_> class SourcePicker {
public:
  explicit SourcePicker(const GraphT_ &g, std::string filename = "",
                        NodeID_ source = UINT_MAX)
      : g_(g), given_source_(source), rng_(27491095), udist_(g.n - 1, rng_) {
    if (filename != "") {
      VectorReader<NodeID_> reader(filename);
      file_sources_ = reader.Read();
    }
  }

  NodeID_ PickNext() {
    // Fixed source
    if (given_source_ != UINT_MAX)
      return given_source_;

    // File sources
    if (!file_sources_.empty()) {
      static size_t current = 0;
      return file_sources_[current++];
    }

    // Random sources
    NodeID_ s, deg;
    do {
      s = udist_();
      deg = g_.offsets[s + 1] - g_.offsets[s];
    } while (deg == 0);

    return s;
  }

private:
  const GraphT_ &g_;
  NodeID_ given_source_;
  std::vector<NodeID_> file_sources_;
  std::mt19937_64 rng_;
  UniDist<NodeID_, std::mt19937_64> udist_;
};

#endif // UTILS_H
