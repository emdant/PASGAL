#include "sssp.h"
#include "utils.h"

#include <iostream>
#include <random>

#include "dijkstra.h"
#include "graph.h"

typedef uint32_t NodeId;
typedef uint64_t EdgeId;
#ifdef FLOAT
typedef float EdgeTy;
#else
typedef uint32_t EdgeTy;
#endif
constexpr int NUM_SRC = 1;
constexpr int NUM_ROUND = 22;
constexpr int LOG2_WEIGHT = 18;
constexpr int WEIGHT_RANGE = 1 << LOG2_WEIGHT;

template <class Algo, class Graph, class NodeId = typename Graph::NodeId>
void run(Algo &algo, [[maybe_unused]] const Graph &G, NodeId s, int rounds,
         bool verify, bool dump) {
  double total_time = 0;
  sequence<EdgeTy> dist;
  for (int i = 0; i < rounds; i++) {
    internal::timer t;
    dist = algo.sssp(s);
    t.stop();
    printf("Round %d: %f\n", i, t.total_time());
    total_time += t.total_time();
  }
  double average_time = total_time / NUM_ROUND;
  printf("Average time: %f\n", average_time);

  auto not_max_cmp = [&](EdgeTy a, EdgeTy b) {
    if (b == Algo::DIST_MAX)
      return false;
    if (a == Algo::DIST_MAX)
      return true;
    return a < b;
  };

  auto not_max = [&](EdgeTy e) { return e != Algo::DIST_MAX; };

  auto longest_distance = *parlay::max_element(dist, not_max_cmp);
  auto reached = parlay::count_if(dist, not_max);
  std::cout << "Longest distance: " << longest_distance << std::endl;
  std::cout << "Nodes reached: " << reached << std::endl;

  if (verify) {
    printf("Running verifier...\n");
    Dijkstra verifier(G);
    auto exp_dist = verifier.dijkstra(s);
    assert(dist == exp_dist);
    printf("Passed!\n");
  }
  if (dump) {
    ofstream ofs("sssp.out");
    for (size_t i = 0; i < dist.size(); i++) {
      ofs << dist[i] << '\n';
    }
    ofs.close();
  }
  printf("\n");
}

template <class Algo, class Graph>
void run(Algo &algo, const Graph &G, int sources, int rounds, bool verify,
         bool dump) {
  std::mt19937_64 rng(27491095);
  UniDist<NodeId, std::mt19937_64> udist(G.n - 1, rng);

  for (int v = 0; v < sources; v++) {
    NodeId s, deg;
    do {
      s = udist();
      deg = G.offsets[s + 1] - G.offsets[s];
    } while (deg == 0);

    printf("source %d: %-10d\n", v, s);
    run(algo, G, s, rounds, verify, dump);
  }
}

int main(int argc, char *argv[]) {
  if (argc == 1) {
    fprintf(stderr,
            "Usage: %s [-i input_file] [-a algorithm] [-p parameter] [-s] [-v] "
            "[-d] [-S sources] [-n trials]\n"
            "Options:\n"
            "\t-i,\tinput file path\n"
            "\t-a,\talgorithm: [rho-stepping] [delta-stepping] [bellman-ford]\n"
            "\t-p,\tparameter(e.g. delta, rho)\n"
            "\t-s,\tsymmetrized input graph\n"
            "\t-v,\tverify result\n"
            "\t-d,\tdump distances to file\n"
            "\t-S,\tnumber of sources\n"
            "\t-n,\tnumber of trials per source\n",
            argv[0]);
    return 0;
  }
  char c;
  char const *input_path = nullptr;
  int algorithm = rho_stepping;
  string parameter;
  uint32_t source = UINT_MAX;
  bool symmetrized = false;
  bool verify = false;
  bool dump = false;
  int rounds = NUM_ROUND;
  int sources = NUM_SRC;

  while ((c = getopt(argc, argv, "i:a:p:r:svdS:n:")) != -1) {
    switch (c) {
    case 'i':
      input_path = optarg;
      break;
    case 'a':
      if (!strcmp(optarg, "rho-stepping")) {
        algorithm = rho_stepping;
      } else if (!strcmp(optarg, "delta-stepping")) {
        algorithm = delta_stepping;
      } else if (!strcmp(optarg, "bellman-ford")) {
        algorithm = bellman_ford;
      } else {
        std::cerr << "Error: Unknown algorithm " << optarg << std::endl;
        abort();
      }
      break;
    case 'p':
      parameter = string(optarg);
      break;
    case 'r':
      source = atol(optarg);
      break;
    case 's':
      symmetrized = true;
      break;
    case 'v':
      verify = true;
      break;
    case 'd':
      dump = true;
      break;
    case 'n':
      rounds = atoi(optarg);
      break;
    case 'S':
      sources = atoi(optarg);
      break;
    default:
      std::cerr << "Error: Unknown option " << optopt << std::endl;
      abort();
    }
  }

  printf("Reading graph...\n");
  Graph<NodeId, EdgeId, EdgeTy> G;
  G.symmetrized = symmetrized;
  if (!strcmp(input_path, "random")) {
    G.generate_random_graph();
  } else {
    G.read_graph(input_path);
  }
  if (!G.weighted) {
    printf("Generating edge weights...\n");
    G.generate_random_weight(1, WEIGHT_RANGE);
  }

  fprintf(stdout,
          "Running on %s: |V|=%zu, |E|=%zu, num_src=%d, num_round=%d\n\n",
          input_path, G.n, G.m, sources, rounds);

  if (algorithm == rho_stepping) {
    size_t rho = 1 << 20;
    if (!parameter.empty()) {
      Rho_Stepping solver(G);
      rho = stoull(parameter);
    }
    Rho_Stepping solver(G, rho);
    if (source == UINT_MAX) {
      run(solver, G, sources, rounds, verify, dump);
    } else {
      run(solver, G, source, rounds, verify, dump);
    }
  } else if (algorithm == delta_stepping) {
    EdgeTy delta = 1 << 15;
    if (!parameter.empty()) {
      if constexpr (is_integral_v<EdgeTy>) {
        delta = stoull(parameter);
      } else {
        delta = stod(parameter);
      }
    }
    Delta_Stepping solver(G, delta);
    if (source == UINT_MAX) {
      run(solver, G, sources, rounds, verify, dump);
    } else {
      run(solver, G, source, rounds, verify, dump);
    }
  } else if (algorithm == bellman_ford) {
    Bellman_Ford solver(G);
    if (source == UINT_MAX) {
      run(solver, G, sources, rounds, verify, dump);
    } else {
      run(solver, G, source, rounds, verify, dump);
    }
  }
  return 0;
}
