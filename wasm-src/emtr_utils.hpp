#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <stdexcept>


namespace emtr {

    // ---------- Dirichlet sampler ----------
    using Md = std::array<std::array<double,4>,4>;

    inline std::array<double, 4> sample_dirichlet(const std::array<double, 4>& alpha,
                                                std::mt19937_64& gen) {
        std::array<double, 4> x{};
        double sum = 0.0;
        for (std::size_t i = 0; i < 4; ++i) {
            std::gamma_distribution<double> gamma(alpha[i], 1.0);
            x[i] = gamma(gen);
            sum += x[i];
        }
        for (auto& v : x) v /= sum;
        return x;
    }

    inline constexpr std::array<double,4> D_pi_default  {100,100,100,100};
    inline constexpr std::array<double,4> D_M_row_default {100,2,2,2};


    inline std::mt19937_64& rng() {
        static thread_local std::mt19937_64 g{0xC0FFEEULL};
        return g;
    }

    inline std::array<double, 4> sample_pi(const std::array<double,4>& alpha = D_pi_default) {
    return sample_dirichlet(alpha, rng());
    }

    inline std::array<double, 4> sample_M_row(const std::array<double,4>& alpha = D_M_row_default) {
    return sample_dirichlet(alpha, rng());
    }

    
    inline bool starts_with(const std::string& str, const std::string& prefix) {
        return str.size() >= prefix.size() && str.compare(0, prefix.size(), prefix) == 0;
    }

    inline std::pair<int,int> ord_pair(int a, int b) {
    return (a < b) ? std::make_pair(a,b) : std::make_pair(b,a);
    }

    inline std::vector<std::string> split_ws(const std::string& s) {
    std::istringstream iss(s);
    std::vector<std::string> out;
    for (std::string tok; iss >> tok;) out.push_back(tok);
    return out;
    }

   
    inline Md MT(const Md& P) {
    Md Pt{};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
        Pt[j][i] = P[i][j];
    return Pt;
    }

    inline Md MM(const Md& A, const Md& B) {
    Md R{};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
        for (int k = 0; k < 4; ++k)
            R[i][j] += A[i][k] * B[k][j];
    return R;
    }

    struct prim_graph {
    int n;
    std::vector<double> W;

    prim_graph(int n_, const std::pair<int,int>* edges, const double* weights, int m)
    : n(n_), W(static_cast<std::size_t>(n_) * static_cast<std::size_t>(n_),
                std::numeric_limits<double>::infinity()) {
        for (int i = 0; i < n; ++i) W[static_cast<std::size_t>(i)*n + i] = 0.0;
        for (int k = 0; k < m; ++k) {
        int u = edges[k].first;
        int v = edges[k].second;
        double w = weights[k];
        if (!(0 <= u && u < n && 0 <= v && v < n))
            throw mt_error("check input for prim's algorithm");
        W[static_cast<std::size_t>(u)*n + v] = w;
        W[static_cast<std::size_t>(v)*n + u] = w;
        }
    }

    int num_vertices() const { return n; }
    double weight(int u, int v) const { return W[static_cast<std::size_t>(u)*n + v]; }
    };

    inline void prim(const prim_graph& g, int* parent_out) {
    const int n = g.num_vertices();
    const double INF = std::numeric_limits<double>::infinity();
    std::vector<double> key(n, INF);
    std::vector<char>   inMST(n, 0);

    key[0] = 0.0;
    parent_out[0] = 0;

    for (int it = 0; it < n; ++it) {
        int u = -1;
        double best = INF;
        for (int v = 0; v < n; ++v) {
        if (!inMST[v] && key[v] < best) { best = key[v]; u = v; }
        }
        if (u == -1) break;
        inMST[u] = 1;

        for (int v = 0; v < n; ++v) {
        if (inMST[v] || v == u) continue;
        double w = g.weight(u, v);
        if (w < key[v]) { key[v] = w; parent_out[v] = u; }
        }
    }
    }
}
