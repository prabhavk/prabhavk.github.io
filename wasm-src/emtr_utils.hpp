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
#include <unordered_map>
#include <stdexcept>
#include <cstdio>
#include <tuple>


namespace emtr {

    
    
    
    
    // Escape a string for JSON (handles quotes, backslashes, control chars)
    inline std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 8);
    for (unsigned char c : s) {
        switch (c) {
        case '\"': out += "\\\""; break;
        case '\\': out += "\\\\"; break;
        case '\b': out += "\\b";  break;
        case '\f': out += "\\f";  break;
        case '\n': out += "\\n";  break;
        case '\r': out += "\\r";  break;
        case '\t': out += "\\t";  break;
        default:
            if (c < 0x20) {
            char buf[7];
            std::snprintf(buf, sizeof(buf), "\\u%04x", static_cast<unsigned>(c));
            out += buf;
            } else {
            out.push_back(static_cast<char>(c));
            }
        }
    }
    return out;
    }

    // Alias that matches your EMTR_results element type
    using Row = std::tuple<
    std::string, // method ("Parsimony" | "Dirichlet" | "SSH")
    std::string, // root (one of 17 names)
    int,         // repetition (1..N)
    int,         // iter (max iterations taken by EM)
    double,      // ll_initial
    double,      // ecd_ll_first
    double,      // ecd_ll_final
    double       // ll_final
    >;

    inline void debug_counts(const std::vector<Row>& v, const char* where) {
    std::unordered_map<std::string,size_t> by;
    for (const auto& t : v) by[std::get<0>(t)]++;
    std::printf("[COUNT]\t%s size=%zu :: ", where, v.size());
    bool first=true;
    for (auto& kv : by) { if(!first) std::printf(", "); std::printf("%s:%zu", kv.first.c_str(), kv.second); first=false; }
    if (first) std::printf("(empty)");
    std::printf("\n"); std::fflush(stdout);
    }

    // Emit one row as JSON line prefixed with a tag (default "[ROW]\t")
    inline void emit_row_json(
        const std::string& method,
        const std::string& root,
        int repetition,
        int iter,
        double ll_initial,
        double ecd_ll_first,
        double ecd_ll_final,
        double ll_final,
        const char* tag = "[ROW]\t") {            
            const std::string im = json_escape(method);
            const std::string rt = json_escape(root);

            // %.17g keeps double precision compact; fflush for prompt delivery to worker
            std::printf(
                "%s"
                "{\"method\":\"%s\",\"root\":\"%s\",\"repetition\":%d,"
                "\"iter\":%d,\"ll_initial\":%.17g,"
                "\"ecd_ll_first\":%.17g,\"ecd_ll_final\":%.17g,\"ll_final\":%.17g}\n",
                tag,
                im.c_str(), rt.c_str(), repetition, iter,
                ll_initial, ecd_ll_first, ecd_ll_final, ll_final);
            std::fflush(stdout);}
    

    // Convenience: push one row into your vector (keeps call sites tidy)
    inline void push_result(
        std::vector<Row>& results,
        const std::string& method,
        const std::string& root,
        int repetition,
        int iter,
        double ll_initial,
        double ecd_ll_first,
        double ecd_ll_final,
        double ll_final)
    {
    results.emplace_back(method, root, repetition, iter,
                        ll_initial, ecd_ll_first, ecd_ll_final, ll_final);
    }

    // Flush & clear a vector of rows, emitting each as JSON
    inline void flush_rows_json(std::vector<Row>& results, const char* tag = "[ROW]\t") {
    for (auto& t : results) {
        emit_row_json(
        std::get<0>(t), std::get<1>(t), std::get<2>(t), std::get<3>(t),
        std::get<4>(t), std::get<5>(t), std::get<6>(t), std::get<7>(t),
        tag
        );
    }
    results.clear();
    }

    // ---------- Dirichlet sampler ----------
    using Md = std::array<std::array<double,4>,4>;

    inline std::mt19937_64& rng() {
        static thread_local std::mt19937_64 g{0xC0FFEEULL};
        return g;
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
