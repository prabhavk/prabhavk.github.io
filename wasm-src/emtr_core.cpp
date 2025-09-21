#include "emtr_core.hpp"
#include "emtr_utils.hpp"

using namespace std; 
using namespace emtr;


int ll_precision = 14;

///...///...///...///...///...///...///... Family Joining ...///...///...///...///...///...///...///
#include <algorithm>
#include <cctype>
#include <sstream>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include "third_party/eigen3/Eigen/Dense"
#include "third_party/eigen3/unsupported/Eigen/MatrixFunctions"
#include <tuple>


static inline void emit_progress(const char* kind,
                                 const char* method = nullptr,
                                 int rep = -1,
                                 int nodeIndex = -1,
                                 int nodeTotal = -1,
                                 int layer = -1,
                                 const char* rootName = nullptr,
                                 int totalReps = -1) {
  nlohmann::json j{
    {"type", "progress"},
    {"kind", kind}
  };
  if (method)     j["method"]     = method;
  if (rep >= 0)   j["rep"]        = rep;
  if (nodeIndex>=0) j["nodeIndex"]= nodeIndex;
  if (nodeTotal>=0) j["nodeTotal"]= nodeTotal;
  if (layer >= 0) j["layer"]      = layer;
  if (rootName)   j["rootName"]   = rootName;
  if (totalReps>=0) j["totalReps"]= totalReps;

  std::cout << "[EM_Method] " << j.dump() << "\n";
  std::fflush(stdout);
}

static inline void emit_layer_done(const char* method,
                                   int rep1,
                                   int nodeIndex1,
                                   int nodeTotal,
                                   int layerIdx,
                                   const std::string& rootName,
                                   double ll = NAN) {
  nlohmann::json j{
    {"type", "progress"},
    {"kind", "layer:done"},
    {"method", method},
    {"rep", rep1},
    {"nodeIndex", nodeIndex1},
    {"nodeTotal", nodeTotal},
    {"layer", layerIdx},
    {"rootName", rootName}
  };
  if (!std::isnan(ll)) j["ll"] = ll;
  std::cout << "[EM_Method] " << j.dump() << "\n";
  std::fflush(stdout);
}



using Matrix20d  = Eigen::Matrix<double,20,20>;
using Matrix20cd = Eigen::Matrix<std::complex<double>,20,20>;
using Vector20cd = Eigen::Matrix<std::complex<double>,20,1>;

struct Spectrum {
    Matrix20cd V;
    Matrix20cd V_inv;
    Vector20cd lambda;
    bool ok = false;
};


namespace fj {

// generic hasher for std::pair<T1,T2>
struct PairHash {
  template <class T1, class T2>
  std::size_t operator()(const std::pair<T1,T2>& p) const noexcept {
    std::size_t h1 = std::hash<T1>{}(p.first);
    std::size_t h2 = std::hash<T2>{}(p.second);
    return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1<<6) + (h1>>2));
  }
};

struct Graph {
    struct Edge { int u, v; double length; Edge():u(-1),v(-1),length(0.0){} Edge(int uu,int vv,double L):u(uu),v(vv),length(L){} };
    std::vector<std::string> name;
    std::vector<std::vector<int>> adj;
    std::vector<Edge> edges;
    std::unordered_map<std::pair<int,int>, int, PairHash> eid;

    int vcount() const { return (int)name.size(); }
    int ecount() const { return (int)edges.size(); }
    void reserve_vertices(std::size_t n) { name.resize(n); adj.resize(n); }

    int add_vertex(const std::string& nm) {
      int id = (int)name.size();
      name.push_back(nm);
      adj.emplace_back();
      return id;
    }
    int add_edge(int u, int v, double len=0.0) {
      if (u==v) throw std::runtime_error("Graph: self-edge not allowed");
      int id = (int)edges.size();
      edges.emplace_back(u,v,len);
      adj[u].push_back(v);
      adj[v].push_back(u);
      eid[std::make_pair(std::min(u,v), std::max(u,v))] = id;
      return id;
    }
    int get_eid(int u, int v) const {
      auto it = eid.find(std::make_pair(std::min(u,v), std::max(u,v)));
      if (it == eid.end()) throw std::runtime_error("Graph: edge not found");
      return it->second;
    }
    double get_length(int u, int v) const { return edges[get_eid(u,v)].length; }
    void   set_length(int u, int v, double L) { edges[get_eid(u,v)].length = L; }
    const std::vector<int>& neighbors(int v) const { return adj[v]; }
    int degree(int v) const { return (int)adj[v].size(); }

    void rebuild() {
      eid.clear();
      for (auto& L : adj) L.clear();
      adj.resize(name.size());
      for (int id=0; id<(int)edges.size(); ++id) {
        const auto& e = edges[id];
        adj[e.u].push_back(e.v);
        adj[e.v].push_back(e.u);
        eid[std::make_pair(std::min(e.u,e.v), std::max(e.u,e.v))] = id;
      }
    }
    void remove_edge(int u, int v) {
      int id = get_eid(u,v);
      std::vector<Edge> kept; kept.reserve(edges.size());
      for (int i=0;i<(int)edges.size();++i) if (i!=id) kept.push_back(edges[i]);
      edges.swap(kept);
      rebuild();
    }
    void delete_vertex(int v) {
      int last = vcount()-1;
      if (v != last) {
        name[v] = std::move(name[last]);
        adj[v]  = std::move(adj[last]);
        for (auto& e : edges) { if (e.u==last) e.u=v; if (e.v==last) e.v=v; }
      }
      name.pop_back();
      adj.pop_back();
      std::vector<Edge> kept; kept.reserve(edges.size());
      for (auto& e : edges) if (e.u!=v && e.v!=v) kept.push_back(e);
      edges.swap(kept);
      rebuild();
    }
  };

// distance map shared across the project
using DistanceMap =
  std::unordered_map<std::pair<std::string,std::string>, double, PairHash>;

// optional helper (inline to avoid multiple definition issues)
inline std::pair<std::string,std::string>
SortedKey(const std::string& a, const std::string& b) {
  return (a < b) ? std::make_pair(a,b) : std::make_pair(b,a);
}

}

static inline std::string json_escape_(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 16);
    static const char* HEX = "0123456789ABCDEF";

    for (unsigned char uc : s) {          // important: unsigned char
        switch (uc) {
            case '\"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\b': out += "\\b";  break;
            case '\f': out += "\\f";  break;
            case '\n': out += "\\n";  break;
            case '\r': out += "\\r";  break;
            case '\t': out += "\\t";  break;
            default:
                if (uc < 0x20 || uc == 0x7F) {
                    out += "\\u00";
                    out += HEX[(uc >> 4) & 0x0F];
                    out += HEX[(uc      ) & 0x0F];
                } else {
                    out += static_cast<char>(uc); // pass through UTF-8 bytes as-is
                }
        }
    }
    return out;
}

static inline std::string num3(double x) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(3) << x;
    return ss.str();
}

static inline std::string estr(const std::string& s) {
    std::ostringstream os;
    os << '"' << json_escape_(s) << '"';
    return os.str();
}

static inline std::string hstring(double x) { return num3(x); }

static inline bool is_vec4_finite(const std::array<double,4>& v) {
    for (double x : v) if (!std::isfinite(x)) return false;
    return true;
}
static inline void dump_vec4_fixed3(std::ostringstream& os, const std::array<double,4>& v) {
    os << '[' << num3(v[0]) << ',' << num3(v[1]) << ','
              << num3(v[2]) << ',' << num3(v[3]) << ']';
}

static inline bool is_md4_finite(const Md& m) {
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            if (!std::isfinite(m[i][j])) return false;
    return true;
}
static inline void dump_md4_fixed3(std::ostringstream& os, const Md& m) {
    os << '[';
    for (int i = 0; i < 4; ++i) {
        if (i) os << ',';
        os << '['
           << num3(m[i][0]) << ',' << num3(m[i][1]) << ','
           << num3(m[i][2]) << ',' << num3(m[i][3]) << ']';
    }
    os << ']';
}

static inline std::string pvec4(const std::array<double,4>& v) {
    std::ostringstream os;
    if (is_vec4_finite(v)) {
        dump_vec4_fixed3(os, v);
    } else {
        os << "[]";
    }
    return os.str();
}

// map<string, Md> -> {"A":[[...],[...],[...],[...]], ...}
static inline std::string hmap_mat4(const std::map<std::string, Md>& mp) {
    std::ostringstream os;
    os << '{';
    bool first = true;
    for (const auto& kv : mp) {
        if (!is_md4_finite(kv.second)) continue;
        if (!first) os << ',';
        first = false;
        os << '"' << json_escape_(kv.first) << "\":";
        dump_md4_fixed3(os, kv.second);
    }
    os << '}';
    return os.str();
}

static inline std::string dump_iter_ll_obj(const std::map<int,double>& m) {
    std::ostringstream os;
    os << '{';
    bool first = true;
    for (const auto& kv : m) {
        if (!std::isfinite(kv.second)) continue;
        if (!first) os << ',';
        first = false;
        os << '"' << kv.first << "\":" << num3(kv.second);
    }
    os << '}';
    return os.str();
}


using namespace fj;

class FamilyJoining {
public:
	void SetNewick(bool striphidden_node = true) {
		if (T.vcount() == 0) {
			throw std::runtime_error("SetNewick: tree is empty. Call SetFJTree() first.");
		}
		Graph copy = T;
		Graph leaf = ConvertGenerallyLabeledTreeToLeafLabeledTree(copy);
		Graph bin  = ConvertMultifurcatingTreeToBifurcatingTree(leaf);
		std::string nwk = GetNewickLabelOfLeafLabeledTree(bin);
		// cout << "newick label with all node names" << nwk;
		// vector <tuple <string, string, double>> edge_vector = this->GetEdgeVector();
		if (striphidden_node) nwk = strip_hidden_node_labels_(nwk);
		newick_ = std::move(nwk);
	}	
  // ------------------------ Public data types ------------------------
  const Graph& GetTree() const { return T; }

  // ------------------------ Construction ------------------------
  FamilyJoining(const std::string& alignment_path, double epsilon) : epsilon_(epsilon) {
    ReadAlignment(alignment_path, alignment_);

    // dm_ = ComputeNHDistances(alignment_);
	// cout << "here 1 \n" << endl;
	int k = 5;
	int repeats = 2;
	dm_ = ComputeLogDetDistances(alignment_);
	SetFJTree(epsilon);	
  }

  // Build FJ tree from distances computed in constructor
  void SetFJTree(double epsilon) {
    // collect names from distances
    std::unordered_set<std::string> setnames;
    for (const auto& kv : dm_) {
      setnames.insert(kv.first.first);
      setnames.insert(kv.first.second);
    }
    std::vector<std::string> vertList(setnames.begin(), setnames.end());
    std::sort(vertList.begin(), vertList.end());

    // build topology, then fit lengths, then clamp negatives
    DistanceMap dm_copy = dm_;
    T = GetFJTreeTopology(std::move(dm_copy), vertList, epsilon);		
		AssignOLSLengthsToTree(T, dm_);
		SetNegativeLengthsToZero(T); // edges with negative length cause plotting problems
  	}

  // Optional getter if you need the string elsewhere
  const std::string& newick() const { return newick_; }

  // ------------------------ Tree Output (file) ------------------------
  // Kept for convenience; not used by JSON path.
  void WriteTree(const Graph& tree,
                 const std::string& outputFileName,
                 const std::string& fileFormat = "edgeList") const {
    if (fileFormat == "edgeList") {
      std::ofstream out(outputFileName.c_str());
      if (!out) throw std::runtime_error("Cannot open output: " + outputFileName);
      for (const auto& e : tree.edges) {
        out << tree.name[e.u] << '\t' << tree.name[e.v] << '\t' << e.length << '\n';
      }
      return;
    }
    if (fileFormat == "newick") {
      Graph copy = tree;
      Graph leaf = ConvertGenerallyLabeledTreeToLeafLabeledTree(copy);
      Graph bin  = ConvertMultifurcatingTreeToBifurcatingTree(leaf);
      std::string nwk = GetNewickLabelOfLeafLabeledTree(bin);
      // strip hidden_node labels to match Python re.sub
      nwk = strip_hidden_node_labels_(nwk);
      std::ofstream out(outputFileName.c_str());
      if (!out) throw std::runtime_error("Cannot open output: " + outputFileName);
      out << nwk;
      return;
    }
    throw std::runtime_error("Unknown format (use 'edgeList' or 'newick'): " + fileFormat);
  }

  Graph T;

  vector <tuple <string, string, double>> GetEdgeVector() {
	vector <tuple <string, string, double>> edge_vector;
	for (auto& e : T.edges) {
		tuple <string,string,double> edge;
		get<0>(edge)  = T.name[e.u];
		get<1>(edge)  = T.name[e.v];
		get<2>(edge)  = e.length;
		edge_vector.push_back(edge);
	}
	return (edge_vector);
  }

  static void AssignOLSLengthsToTree(Graph& T, const DistanceMap& DM) {
    const std::size_t M = DM.size();
    int N = (int)((1.0 + std::sqrt(1.0 + 8.0 * (double)M)) / 2.0);
    IntSet vertSet; for (int i=0;i<N;++i) vertSet.insert(i);
    auto splits = ComputeSplits(T, vertSet);
    auto A      = ComputeAtd(T, splits, vertSet, DM);
    auto bl     = ComputeEdges(T, splits, vertSet, A);
    for (const auto& kv : bl) {
      int u = kv.first.first, v = kv.first.second;
      T.set_length(u,v, kv.second);
    }
  }

  static void SetNegativeLengthsToZero(Graph& T, double minLen = 0.0) {
    for (auto& e : T.edges) if (e.length <= 0.0) e.length = minLen;
  }

private: 	
  	string newick_;
	static std::string strip_hidden_node_labels_(const std::string& s) {
	std::string out; out.reserve(s.size());
	std::size_t i = 0, n = s.size();		
	// get distance_matrices for test/train

	while (i < n) {
		if (i + 2 <= n && s.compare(i, 2, "h_") == 0) {
		std::size_t j = i + 12;
		if (j < n && s[j] == 'T') ++j;
		if (j < n && s[j] == '_') ++j;
		while (j < n && std::isdigit((int)s[j])) ++j;
		if (j < n && s[j] == '_') {
			++j;
			while (j < n && std::isdigit((int)s[j])) ++j;
		}		
		i = j;
		continue;
		}
		out.push_back(s[i++]);
	}
		return out;
	}
	
	static std::string json_escape_(const std::string& s) {
	std::string out;
	out.reserve(s.size() + 16);
	static const char* HEX = "0123456789ABCDEF";

	for (char c : s) {
		int uc = static_cast<int>(c);
		switch (c) {
		case '\"': out += "\\\""; break;
		case '\\': out += "\\\\"; break;
		case '\b': out += "\\b";  break;
		case '\f': out += "\\f";  break;
		case '\n': out += "\\n";  break;
		case '\r': out += "\\r";  break;
		case '\t': out += "\\t";  break;
		default:
			if (uc < 0x20) {
			out += "\\u00";
			out += HEX[(uc >> 4) & 0x0F];
			out += HEX[(uc      ) & 0x0F];
			} else {
			out += c;  // pass through UTF-8 bytes as-is
			}
		}
	}
	return out;
	}


  // ------------------------ Internal State ------------------------
  double epsilon_;
  std::unordered_map<std::string, std::string> alignment_;
  

//   using DistanceMap = std::unordered_map<std::pair<std::string,std::string>, double, PairHash>;
  DistanceMap dm_;

  // ------------------------ String helpers ------------------------
  static std::pair<std::string,std::string> SortedKey(const std::string& a, const std::string& b) {
    return (a<b) ? std::make_pair(a,b) : std::make_pair(b,a);
  }
  static void Trim(std::string& s) {
    auto sp = [](int ch){ return std::isspace(ch); };
    while (!s.empty() && sp((int)s.front())) s.erase(s.begin());
    while (!s.empty() && sp((int)s.back()))  s.pop_back();
  }
  static void ToUpper(std::string& s) {
    for (char& c : s) if (c>='a' && c<='z') c = char(c - 'a' + 'A');
  }
  static bool IsHiddenNode(const std::string& nm) {
    return nm.size()>=2 && (nm.compare(0,2,"h_")==0);
  }
 
  double GetEpsilonUpper(const DistanceMap& distanceMatrixOriginal) {
    // Local copy
    DistanceMap distanceMatrix = distanceMatrixOriginal;

    // Helper to canonicalize keys (sorted pair)
    auto keyOf = [](const std::string& a, const std::string& b) -> std::pair<std::string,std::string> {
        return (a < b) ? std::make_pair(a,b) : std::make_pair(b,a);
    };

    // ---- Build vertex list (unique, sorted) ----
    std::vector<std::string> fullvertList;
    fullvertList.reserve(distanceMatrix.size() * 2);
    for (const auto& kv : distanceMatrix) {
        fullvertList.push_back(kv.first.first);
        fullvertList.push_back(kv.first.second);
    }
    std::sort(fullvertList.begin(), fullvertList.end());
    fullvertList.erase(std::unique(fullvertList.begin(), fullvertList.end()), fullvertList.end());
    std::vector<std::string> vertList = fullvertList;

    int n = static_cast<int>(vertList.size());
    if (n < 3) return 0.0; 

    // ---- uList initialization ----
    std::unordered_map<std::string, double> uList;
    uList.reserve(vertList.size());
    for (const auto& v : vertList) uList[v] = 0.0;

    // Accumulate sums of pairwise distances
    for (int ii = 0; ii < n; ++ii) {
        const auto& i = vertList[ii];
        for (int jj = ii + 1; jj < n; ++jj) {
            const auto& j = vertList[jj];
            const auto it = distanceMatrix.find(keyOf(i,j));
            if (it == distanceMatrix.end()) {
                throw std::runtime_error("GetEpsilonUpper: missing distance for pair (" + i + "," + j + ")");
            }
            const double dist = it->second;
            uList[i] += dist;
            uList[j] += dist;
        }
    }
    // Divide by (n - 2)
    if (n > 2) {
        for (const auto& v : vertList) {
            uList[v] /= static_cast<double>(n - 2);
        }
    }

    double epsilon_upper = 0.0;

    // ---- Main reduction loop ----
    while (static_cast<int>(uList.size()) > 3) {
        // Find pair (i_selected, j_selected) minimizing neighborDist_current = d(i,j) - u[i] - u[j]
        double neighborDist = std::numeric_limits<double>::infinity();
        std::string i_selected, j_selected;

        for (int ii = 0; ii < n; ++ii) {
            const auto& i = vertList[ii];
            for (int jj = ii + 1; jj < n; ++jj) {
                const auto& j = vertList[jj];
                const auto it = distanceMatrix.find(keyOf(i,j));
                if (it == distanceMatrix.end()) continue;
                const double neighborDist_current = it->second - uList[i] - uList[j];
                if (neighborDist_current < neighborDist) {
                    neighborDist = neighborDist_current;
                    i_selected = i;
                    j_selected = j;
                }
            }
        }

        // Distances & deltas for the selected pair
        const double d_ij = distanceMatrix.at(keyOf(i_selected, j_selected));
        const double delta_ij = std::fabs(d_ij + uList[i_selected] - uList[j_selected]);
        const double delta_ji = std::fabs(d_ij - uList[i_selected] + uList[j_selected]);
        const double delta_pcvss = std::min(delta_ij, delta_ji) / 2.0;

        double delta_ps2 = delta_pcvss + 1.0;

        // Consider all other k
        for (const auto& k : vertList) {
            if (k == i_selected || k == j_selected) continue;
            const double dik = distanceMatrix.at(keyOf(k, i_selected));
            const double djk = distanceMatrix.at(keyOf(k, j_selected));
            const double delta = dik + djk - d_ij;
            const double delta_sMax = std::fabs(delta) / 4.0;
            if (delta_sMax < delta_ps2) delta_ps2 = delta_sMax;
        }

        epsilon_upper = std::max(epsilon_upper, std::min(delta_pcvss, delta_ps2));

        // ----- Case 1: remove 1 node (child) -----
        if (delta_pcvss < delta_ps2) {
            // Choose child per Python logic
            const std::string child = (delta_ij < delta_ji) ? j_selected : i_selected;

            // remove child from uList and update others
            // Python: n -= 1, then for each vert != child:
            // u[vert] = (u[vert]*(n-1) - d(vert,child)) / (n-2)
            // and erase all distances that involve child
            // Also remove from vertList.

            // Decrement n first (to match Python's order)
            n -= 1;

            // Update u for all remaining vertices (excluding child)
            for (const auto& vert : vertList) {
                if (vert == child) continue;
                const auto it = distanceMatrix.find(keyOf(vert, child));
                if (it == distanceMatrix.end()) {
                    // If missing, assume 0 or throw; Python assumes present
                    throw std::runtime_error("GetEpsilonUpper: missing distance involving child (" + vert + "," + child + ")");
                }
                const double d_vc = it->second;
                // (n here is already decremented, matching Python)
                uList[vert] = (uList[vert] * static_cast<double>(n) - d_vc) / static_cast<double>(n - 1);
            }

            // Erase child from uList
            uList.erase(child);

            // Erase all distances involving child
            for (const auto& vert : vertList) {
                if (vert == child) continue;
                distanceMatrix.erase(keyOf(vert, child));
            }

            // Remove child from vertList
            vertList.erase(std::remove(vertList.begin(), vertList.end(), child), vertList.end());
        }
        // ----- Case 2: remove 2 nodes (both i_selected and j_selected) -----
        else {
            // vertList.remove(i_selected)
            // vertList.remove(j_selected)
            // n -= 2
            // del u[i_selected], u[j_selected]
            // for each remaining vert:
            // u[vert] = (u[vert]*n - d(vert,i) - d(vert,j)) / (n-2)
            // and erase distances to i and j.
            // If n <= 2, return epsilon_upper.

            // Decrement n by 2 (to match Python's order)
            n -= 2;

            // Update u for all remaining vertices
            for (const auto& vert : vertList) {
                if (vert == i_selected || vert == j_selected) continue;
                const double d_vi = distanceMatrix.at(keyOf(vert, i_selected));
                const double d_vj = distanceMatrix.at(keyOf(vert, j_selected));
                // (n here is already decremented by 2)
                if (n > 2) {
                    uList[vert] = (uList[vert] * static_cast<double>(n) - d_vi - d_vj) / static_cast<double>(n - 2);
                }
            }

            // Erase i_selected, j_selected from uList
            uList.erase(i_selected);
            uList.erase(j_selected);

            // Erase all distances involving i_selected or j_selected
            for (const auto& vert : vertList) {
                if (vert != i_selected) distanceMatrix.erase(keyOf(vert, i_selected));
                if (vert != j_selected) distanceMatrix.erase(keyOf(vert, j_selected));
            }

            // Remove i_selected and j_selected from vertList
            vertList.erase(std::remove(vertList.begin(), vertList.end(), i_selected), vertList.end());
            vertList.erase(std::remove(vertList.begin(), vertList.end(), j_selected), vertList.end());

            if (n <= 2) {
                return epsilon_upper;
            }
        }
	}
	if (vertList.size() == 3) {
		const std::string& v0 = vertList[0];
		const std::string& v1 = vertList[1];
		const std::string& v2 = vertList[2];

		const double d01 = distanceMatrix.at(keyOf(v0, v1));
		const double d02 = distanceMatrix.at(keyOf(v0, v2));
		const double d12 = distanceMatrix.at(keyOf(v1, v2));

		const double m = std::min(
			{ std::fabs(d12 - d01 - d02),
			std::fabs(d02 - d01 - d12),
			std::fabs(d01 - d02 - d12) }
		) / 2.0;

		epsilon_upper = std::max(epsilon_upper, m);
		}

		return epsilon_upper + 2.0e-10;
	}

  // ------------------------ Alignment + Distances ------------------------
  static void ReadAlignment(const std::string& path,
                            std::unordered_map<std::string,std::string>& aln) {
    std::ifstream fin(path.c_str());
    if (!fin) throw std::runtime_error("Cannot open alignment: " + path);

    std::string first;
    while (std::getline(fin, first)) if (!first.empty()) break;
    if (!fin && first.empty()) throw std::runtime_error("Empty alignment: " + path);

    bool isFasta = !first.empty() && first[0] == '>';
    fin.clear();
    fin.seekg(0, std::ios::beg);

    if (isFasta) {
      std::string name, seq, line;
      while (std::getline(fin, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
          if (!name.empty()) { ToUpper(seq); aln[name] = seq; seq.clear(); }
          name = line.substr(1); Trim(name);
        } else { Trim(line); seq += line; }
      }
      if (!name.empty()) { ToUpper(seq); aln[name] = seq; }
    } else {
      std::string line;
      while (std::getline(fin, line)) {
        if (line.empty()) continue;
        auto tab = line.find('\t');
        if (tab == std::string::npos)
          throw std::runtime_error("Expected NAME<TAB>SEQ; got: " + line);
        std::string name = line.substr(0, tab);
        std::string seq  = line.substr(tab+1);
        Trim(name); Trim(seq); ToUpper(seq);
        if (name.empty() || seq.empty())
          throw std::runtime_error("Bad alignment line: " + line);
        aln[name] = seq;
      }
    }
    if (aln.empty()) throw std::runtime_error("No sequences parsed");
    const std::size_t L = aln.begin()->second.size();
    for (const auto& kv : aln) {
      if (kv.second.size() != L)
        throw std::runtime_error("Sequences must be aligned (equal length): " + kv.first);
    }
  }

  struct CVFold {
    	std::vector<std::size_t> train_pos;
    	std::vector<std::size_t> test_pos;
	};

  inline std::array<CVFold, 5> make5FoldCV(std::size_t Alen, std::uint64_t seed = 22, bool shuffle = true) {
    std::array<CVFold, 5> folds;

    // 1) Build index list [0, 1, 2, ..., Alen-1]
    std::vector<std::size_t> idx(Alen);
    std::iota(idx.begin(), idx.end(), 0);

    // 2) Optional shuffle (for reproducibility via seed)
    if (shuffle) {
        std::mt19937_64 rng(seed);
        std::shuffle(idx.begin(), idx.end(), rng);
    }

    // 3) Compute fold sizes as even as possible
    const std::size_t base = Alen / 5;
    const std::size_t rem  = Alen % 5; // first 'rem' folds get one extra

    std::size_t offset = 0;
    for (std::size_t f = 0; f < 5; ++f) {
        const std::size_t take = base + (f < rem ? 1 : 0);

        // test = this slice
        std::vector<std::size_t> test;
        test.reserve(take);
        for (std::size_t k = 0; k < take; ++k)
            test.push_back(idx[offset + k]);

        // train = everything except this slice
        std::vector<std::size_t> train;
        train.reserve(Alen - take);
        // before slice
        for (std::size_t k = 0; k < offset; ++k)
            train.push_back(idx[k]);
        // after slice
        for (std::size_t k = offset + take; k < idx.size(); ++k)
            train.push_back(idx[k]);

        folds[f] = CVFold{ std::move(train), std::move(test) };
        offset += take;
    }

    return folds;
}		

  void compute_distance_matrices_for_train_test_pairs(const std::unordered_map<std::string,std::string>& aln, int num_replicates=100, int k=5) {
	// only compute log-det distances
	std::vector<std::string> names; names.reserve(aln.size());
    for (const auto& kv : aln) names.push_back(kv.first);
    std::sort(names.begin(), names.end());
	const std::size_t Alen = aln.at(names[0]).size();
    if (names.size() < 2) throw std::runtime_error("Need >= 2 sequences.");
	// generate positions for train/test alignment pairs
	// generate distance_matrices for train/test alignment pairs
  	map <pair<int,int>,pair<DistanceMap,DistanceMap>> kfoldcv_dm;
	int seq_len;
  }
  
  static DistanceMap ComputeNHDistances(const std::unordered_map<std::string,std::string>& aln) {
    std::vector<std::string> names; names.reserve(aln.size());
    for (const auto& kv : aln) names.push_back(kv.first);
    std::sort(names.begin(), names.end());
    if (names.size() < 2) throw std::runtime_error("Need >= 2 sequences.");

    const std::size_t L = aln.at(names[0]).size();
    DistanceMap D;
    for (std::size_t i=0;i<names.size();++i) {
      const std::string& A = names[i]; const std::string& sA = aln.at(A);
      for (std::size_t j=i+1;j<names.size();++j) {
        const std::string& B = names[j]; const std::string& sB = aln.at(B);
        double diff = 0.0; for (std::size_t k=0;k<L;++k) if (sA[k]!=sB[k]) diff += 1.0;
        D[SortedKey(A,B)] = diff / (double)L;
      }
    }
    return D;
  }

// Compute LogDet distances for gappy sequences

//   Md DivergenceMatrix;
  	static int ntIndex(char c) {
		switch (c) {
			case 'A': case 'a': return 0;
			case 'C': case 'c': return 1;
			case 'G': case 'g': return 2;
			case 'T': case 't': case 'U': case 'u': return 3;
			default: return -1;
		}
	}

	// Build 4x4 divergence/count matrix D(i,j) for aligned sequences sA vs sB.
	static Eigen::Matrix4d GetDivergenceMatrix(const std::string& sA, const std::string& sB) {
		Eigen::Matrix4d D = Eigen::Matrix4d::Zero();
		const std::size_t L = std::min(sA.size(), sB.size());

		for (std::size_t k = 0; k < L; ++k) {
			const int ia = ntIndex(sA[k]);
			const int ib = ntIndex(sB[k]);
			if (ia >= 0 && ib >= 0) {
				D(ia, ib) += 1.0;
			}
		}
		return D;
	}  

	// LogDet (paralinear) distance for a 4x4 divergence matrix D.
	// D(i,j) = count of aligned sites with pair (i,j).
	// alpha: Laplace pseudocount added to every cell (0.5 is common).
	static double LogDet_LD(const Eigen::Matrix4d& D, double alpha = 0.5)
	{
		// Smoothed copy (Matrix4d)
		Eigen::Matrix4d M = D;
		if (alpha != 0.0) M.array() += alpha;  // elementwise add, still Matrix4d

		// Row/column sums (Vector4d)
		Eigen::Vector4d rowSums = M.rowwise().sum();
		Eigen::Vector4d colSums = M.colwise().sum();

		// Guards
		if ((rowSums.array() <= 0.0).any() || (colSums.array() <= 0.0).any()) {
			return std::numeric_limits<double>::infinity();
		}

		const double detM = M.determinant();
		if (!(detM > 0.0) || !std::isfinite(detM)) {
			return std::numeric_limits<double>::infinity();
		}

		// Sum of logs of row/column sums
		const double sumLogRows = rowSums.array().log().sum();
		const double sumLogCols = colSums.array().log().sum();

		// N = 4 for DNA
		constexpr double N = 4.0;
		return -(1.0 / N) * ( std::log(detM) - 0.5 * (sumLogCols + sumLogRows) );
	}

   static double GetLogDetDistance(const std::string& sA, const std::string& sB) {
	Eigen::Matrix4d D = GetDivergenceMatrix(sA,sB);
	double ldd = LogDet_LD(D);
	return (ldd);
   }

  static DistanceMap ComputeLogDetDistances(const std::unordered_map<std::string,std::string>& aln) {
    std::vector<std::string> names; names.reserve(aln.size());
    for (const auto& kv : aln) names.push_back(kv.first);
    std::sort(names.begin(), names.end());
    if (names.size() < 2) throw std::runtime_error("Need >= 2 sequences.");

    const std::size_t L = aln.at(names[0]).size(); // alignment length
    DistanceMap D;
    for (std::size_t i=0;i<names.size();++i) {
      const std::string& A = names[i]; const std::string& sA = aln.at(A);	  
      for (std::size_t j=i+1;j<names.size();++j) {
        const std::string& B = names[j]; const std::string& sB = aln.at(B);

        D[SortedKey(A,B)] = GetLogDetDistance(sA,sB);
      }
    }
    return D;
  }

  // ----------------------------- k fold CV -----------------------------

  inline Eigen::Matrix4d DivergenceMatrixAt(const std::string& sA, const std::string& sB, const std::vector<std::size_t>& cols) {
    const std::size_t Lmin = std::min(sA.size(), sB.size());
    Eigen::Matrix4d D = Eigen::Matrix4d::Zero();
    for (std::size_t p : cols) {
        if (p >= Lmin) continue;
        const int ia = ntIndex(sA[p]);
        const int ib = ntIndex(sB[p]);
        if (ia >= 0 && ib >= 0) D(ia, ib) += 1.0;
    	}
    	return D;
	}

	inline double GetLogDetDistanceAt(const std::string& sA, const std::string& sB, const std::vector<std::size_t>& cols) {
		if (cols.empty()) return std::numeric_limits<double>::quiet_NaN();
		Eigen::Matrix4d D = DivergenceMatrixAt(sA, sB, cols);
		return LogDet_LD(D);
	}

// Split a shuffled index vector into K contiguous folds (balanced sizes)
	inline std::vector<std::vector<std::size_t>> makeKFoldsFromIndex(const std::vector<std::size_t>& idx, std::size_t K) {
    std::vector<std::vector<std::size_t>> folds;
    folds.reserve(K);
    const std::size_t N = idx.size();
    const std::size_t base = K ? (N / K) : 0;
    const std::size_t rem  = K ? (N % K) : 0;

    std::size_t off = 0;
    for (std::size_t k = 0; k < K; ++k) {
        const std::size_t take = base + (k < rem ? 1 : 0);
        std::vector<std::size_t> part;
        part.reserve(take);
        for (std::size_t t = 0; t < take; ++t) part.push_back(idx[off + t]);
        folds.push_back(std::move(part));
        off += take;
    }
    	return folds;
	}

// --- : repeated K-fold CV distance matrices (train/test) ---
	inline std::map<std::pair<int,int>, std::pair<DistanceMap,DistanceMap>> 
ComputeKFoldCVDistances(const std::unordered_map<std::string, std::string>& aln, int K = 5, int R = 32, std::uint64_t seed = 22) {
    // Collect and sort names
    std::vector<std::string> names;
    names.reserve(aln.size());
    for (const auto& kv : aln) names.push_back(kv.first);
    std::sort(names.begin(), names.end());
    if (names.size() < 2) throw std::runtime_error("Need >= 2 sequences.");

    // Use min length across all sequences (robust if some lengths differ)
    std::size_t NAC = std::numeric_limits<std::size_t>::max();
    for (const auto& nm : names) NAC = std::min(NAC, aln.at(nm).size());
	// cout << "alignment length is " << NAC << endl;
    if (NAC == 0) throw std::runtime_error("Alignment length is zero.");

    // Base index 0..NAC-1
    std::vector<std::size_t> base(NAC);
    std::iota(base.begin(), base.end(), 0);

    std::map<std::pair<int,int>, std::pair<DistanceMap,DistanceMap>> out;

    // Repeats
    for (int r = 0; r < R; ++r) {
        // Shuffle indices deterministically per repeat
        std::vector<std::size_t> idx = base;
        std::mt19937_64 rng(seed + 0x9E3779B97F4A7C15ULL * static_cast<std::uint64_t>(r));
        std::shuffle(idx.begin(), idx.end(), rng);

        // Make K (contiguous) folds on this shuffled order
        const auto testFolds = makeKFoldsFromIndex(idx, static_cast<std::size_t>(K));

        // Precompute train positions for each fold as complement of test
        std::vector<std::vector<std::size_t>> trainFolds;
        trainFolds.reserve(testFolds.size());
        for (int k = 0; k < K; ++k) {
            std::vector<std::size_t> train;
            train.reserve(NAC - testFolds[k].size());
            for (int kk = 0; kk < K; ++kk) {
                if (kk == k) continue;
                train.insert(train.end(), testFolds[kk].begin(), testFolds[kk].end());
            }
            trainFolds.push_back(std::move(train));
        }

        // Debug print only for first repetition
        if (r == 0) {
            for (int k = 0; k < K; ++k) {
                // std::cout << "Fold " << k 
                        //   << " -> test size: " << testFolds[k].size()
                        //   << ", train size: " << trainFolds[k].size() << "\n";
            }
            // Print positions in test data for first fold
            if (!testFolds.empty()) {
                // std::cout << "Test positions (fold 0): ";
                // for (auto pos : testFolds[0]) std::cout << pos << " ";
                // std::cout << "\n";
				// std::cout << "Train positions (fold 0): ";
                // for (auto pos : trainFolds[0]) std::cout << pos << " ";
                // std::cout << "\n";
            }
            // Print positions in train data for fifth fold (index 4)
            if (K >= 5) {
				// std::cout << "Test positions (fold 4): ";
                // for (auto pos : testFolds[4]) std::cout << pos << " ";
                // std::cout << "\n";
                // std::cout << "Train positions (fold 4): ";
                // for (auto pos : trainFolds[4]) std::cout << pos << " ";
                // std::cout << "\n";
            }
        }

        // For each fold: compute train/test distance maps
        for (int k = 0; k < K; ++k) {
            DistanceMap dmTrain;
            DistanceMap dmTest;

            const auto& test_pos  = testFolds[k];
            const auto& train_pos = trainFolds[k];

            for (std::size_t i = 0; i < names.size(); ++i) {
                const std::string& A  = names[i];
                const std::string& sA = aln.at(A);
                for (std::size_t j = i + 1; j < names.size(); ++j) {
                    const std::string& B  = names[j];
                    const std::string& sB = aln.at(B);                    

                    const double dTrain = GetLogDetDistanceAt(sA, sB, train_pos);
                    dmTrain[SortedKey(A, B)] = dTrain;

                    const double dTest  = GetLogDetDistanceAt(sA, sB, test_pos);
                    dmTest[SortedKey(A, B)]  = dTest;
                }
            }

            out[{r, k}] = std::make_pair(std::move(dmTrain), std::move(dmTest));
        }
    }
    return out;
}

Graph GetFJTree(DistanceMap DM, double epsilon) {
    // collect names from distances
    std::unordered_set<std::string> setnames;
    for (const auto& kv : DM) {
      setnames.insert(kv.first.first);
      setnames.insert(kv.first.second);
    }
    std::vector<std::string> vertList(setnames.begin(), setnames.end());
    std::sort(vertList.begin(), vertList.end());

    // build topology, then fit lengths, then clamp negatives
    
    Graph fj = GetFJTreeTopology(std::move(DM), vertList, epsilon);
	AssignOLSLengthsToTree(fj, DM);
	SetNegativeLengthsToZero(fj); // edges with negative length cause plotting problems
	return (fj);
  }

 void GetErrorForFold(Graph fj_train, DistanceMap dm_test) {

 }
 
 void GetTestErrorForKfoldCVatEpsilon(int K, int R, const std::map<std::pair<int, int>, std::pair<DistanceMap, DistanceMap>>&kfoldDM, double epsilon) {	
	DistanceMap dm_train;
	DistanceMap dm_test;
	Graph fj_train;
	for (int r = 0; r < R; r++) {
		for (int k = 0; k < K; k++) {
			auto it = kfoldDM.find({r, k});
            dm_train = it->second.first;
            dm_test  = it->second.second;
			fj_train = GetFJTree(dm_train,epsilon);
		}
	}
}

  // ------------------------ FJ (topology + OLS) ------------------------
  using IntSet  = std::unordered_set<int>;
  using EdgeKey = std::pair<int,int>;
  struct EdgeKeyHash { std::size_t operator()(const EdgeKey& p) const noexcept { return PairHash{}(p); } };
  using SplitMap = std::unordered_map<EdgeKey, std::pair<IntSet,IntSet>, EdgeKeyHash>;
  using AtdMap   = std::unordered_map<EdgeKey, double, EdgeKeyHash>;

  static int GetInsertIndex(const std::vector<std::string>& sortedList, int lengthOfSortedList, const std::string& item) {
    int i = (int)std::floor(lengthOfSortedList/2.0), lower=0, upper=lengthOfSortedList-1;
    while (true) {
      if (sortedList[i] > item) {
        if (i==0) return i;
        if (sortedList[i-1] > item) { upper=i-1; i=(int)std::floor((lower+i)/2.0); }
        else return i;
      } else {
        if (i==lengthOfSortedList-1) return i+1;
        if (sortedList[i+1] <= item) { lower=i+1; i=(int)std::ceil((upper+i)/2.0); }
        else return i+1;
      }
    }
  }

  
  static SplitMap ComputeSplits(const Graph& T, const IntSet& vertSet) {
    SplitMap splits;
    std::vector<int> degrees(T.vcount());
    for (int v=0; v<T.vcount(); ++v) degrees[v] = T.degree(v);

    std::vector<EdgeKey> activeEdges;
    std::unordered_set<int> passive;
    std::unordered_map<int,int> countIncident;

    for (int v=0; v<T.vcount(); ++v)
      if (degrees[v]==1)
        activeEdges.emplace_back(v, T.neighbors(v)[0]);

    auto addCandidate = [&](int from, int to) {
      EdgeKey e(from,to), rev(to,from);
      if (std::find(activeEdges.begin(), activeEdges.end(), e) == activeEdges.end() &&
          std::find(activeEdges.begin(), activeEdges.end(), rev) == activeEdges.end()) {
        activeEdges.push_back(e);
      }
    };

    while (!activeEdges.empty()) {
      EdgeKey e = activeEdges.front(); activeEdges.erase(activeEdges.begin());
      int beta = e.first, alpha = e.second;
      passive.insert(beta);

      if (degrees[beta]==1) {
        IntSet L, R; L.insert(beta);
        for (int x : vertSet) if (x!=beta) R.insert(x);
        splits[e] = std::make_pair(L,R);

        countIncident[alpha] += 1;
        if (countIncident[alpha] == T.degree(alpha)-1) {
          int nxt=-1; for (int nb : T.neighbors(alpha)) if (passive.find(nb)==passive.end()) { nxt=nb; break; }
          if (nxt!=-1) addCandidate(alpha,nxt);
        }
      } else {
        std::vector<int> visited;
        for (int nb : T.neighbors(beta)) if (nb!=alpha) visited.push_back(nb);

        IntSet Cleft = splits[EdgeKey(visited[0], beta)].first;
        for (std::size_t i=1;i<visited.size();++i) {
          const IntSet& addL = splits[EdgeKey(visited[i], beta)].first;
          Cleft.insert(addL.begin(), addL.end());
        }
        IntSet Cright; for (int x : vertSet) if (Cleft.find(x)==Cleft.end()) Cright.insert(x);
        if (Cright.find(beta)!=Cright.end()) { Cright.erase(beta); Cleft.insert(beta); }
        splits[e] = std::make_pair(Cleft, Cright);

        if (passive.find(alpha)==passive.end()) {
          countIncident[alpha] += 1;
          if (countIncident[alpha] == T.degree(alpha)-1) {
            int nxt=-1; for (int nb : T.neighbors(alpha)) if (passive.find(nb)==passive.end()) { nxt=nb; break; }
            if (nxt!=-1) addCandidate(alpha,nxt);
          }
        }
      }
    }
    return splits;
  }

  static AtdMap ComputeAtd(const Graph& T, const SplitMap& splits,
                           const IntSet& vertSet, const DistanceMap& DM) {
    int maxIdx=-1; for (int x : vertSet) maxIdx = std::max(maxIdx, x);
    std::vector<std::string> vertexNames(maxIdx+1);
    for (int i=0;i<=maxIdx;++i) vertexNames[i] = T.name[i];

    std::vector<std::string> sortedNames = vertexNames;
    std::sort(sortedNames.begin(), sortedNames.end());

    std::vector<int> degrees(T.vcount());
    for (int v=0; v<T.vcount(); ++v) degrees[v] = T.degree(v);

    std::vector<EdgeKey> activeEdges;
    std::unordered_set<int> passive;
    std::unordered_map<int,int> countIncident;
    AtdMap A;

    for (int v=0; v<T.vcount(); ++v)
      if (degrees[v]==1)
        activeEdges.emplace_back(v, T.neighbors(v)[0]);

    auto addCandidate = [&](int from, int to) {
      EdgeKey e(from,to), rev(to,from);
      if (std::find(activeEdges.begin(), activeEdges.end(), e) == activeEdges.end() &&
          std::find(activeEdges.begin(), activeEdges.end(), rev) == activeEdges.end()) {
        activeEdges.push_back(e);
      }
    };

    while (!activeEdges.empty()) {
      EdgeKey e = activeEdges.front(); activeEdges.erase(activeEdges.begin());
      int beta = e.first, alpha = e.second;
      passive.insert(beta);

      if (degrees[beta]==1) {
        double s = 0.0;
        const std::string& nm = vertexNames[beta];
        int idx = (int)(std::lower_bound(sortedNames.begin(), sortedNames.end(), nm) - sortedNames.begin());
        for (int j=0;j<idx;++j) s += DM.at(SortedKey(sortedNames[j], nm));
        for (int j=idx+1;j<(int)sortedNames.size();++j) s += DM.at(SortedKey(nm, sortedNames[j]));
        A[e] = s;

        countIncident[alpha] += 1;
        if (countIncident[alpha] == T.degree(alpha)-1) {
          int nxt=-1; for (int nb : T.neighbors(alpha)) if (passive.find(nb)==passive.end()) { nxt=nb; break; }
          if (nxt!=-1) addCandidate(alpha,nxt);
        }
      } else {
        double s = 0.0;
        std::vector<int> visited;
        for (int nb : T.neighbors(beta)) if (nb!=alpha) visited.push_back(nb);

        for (int k : visited) s += A.at(EdgeKey(k,beta));

        for (std::size_t i=0;i<visited.size();++i) {
          int k = visited[i];
          const IntSet& CkL = splits.at(EdgeKey(k,beta)).first;
          for (std::size_t j=i+1;j<visited.size();++j) {
            int l = visited[j];
            const IntSet& ClL = splits.at(EdgeKey(l,beta)).first;
            for (int a : CkL) for (int b : ClL) s -= 2.0 * DM.at(SortedKey(vertexNames[a], vertexNames[b]));
          }
        }

        if (vertSet.find(beta)!=vertSet.end()) {
          for (int k : visited) {
            const IntSet& CkL = splits.at(EdgeKey(k,beta)).first;
            for (int b : CkL) s -= DM.at(SortedKey(vertexNames[beta], vertexNames[b]));
          }
          const auto& Ci = splits.at(EdgeKey(beta,alpha));
          for (int b : Ci.second) s += DM.at(SortedKey(vertexNames[beta], vertexNames[b]));
        }

        A[e] = s;

        countIncident[alpha] += 1;
        if (countIncident[alpha] == T.degree(alpha)-1) {
          int nxt=-1; for (int nb : T.neighbors(alpha)) if (passive.find(nb)==passive.end()) { nxt=nb; break; }
          if (nxt!=-1) addCandidate(alpha,nxt);
        }
      }
    }
    return A;
  }

  static std::unordered_map<EdgeKey,double,EdgeKeyHash>
  ComputeEdges(const Graph& T, const SplitMap& splits, const IntSet& vertSet, const AtdMap& A) {
    const double n = (double)vertSet.size();
    std::unordered_map<EdgeKey,double,EdgeKeyHash> bl;

    for (const auto& kv : splits) {
      EdgeKey edge = kv.first;
      const auto& C_beta = kv.second.first;
      const auto& C_alpha = kv.second.second;
      const double P0 = A.at(edge);
      double numerator = P0, denominator = 0.0;

      std::vector<double> sList, wList, PList, nList;
      bool caseNby2 = false; int kIdx = -1;

      int beta = edge.first, alpha = edge.second;

      if ((int)C_beta.size() == 1) {
        denominator = n - 1.0;

        std::vector<int> nbs;
        for (int nb : T.neighbors(alpha)) if (nb!=beta) nbs.push_back(nb);

        for (int nb : nbs) {
          const auto* Ci = [&]()->const std::unordered_set<int>*{
            auto it = splits.find(EdgeKey(alpha,nb));
            return (it!=splits.end()) ? &it->second.second : &splits.at(EdgeKey(nb,alpha)).first;
          }();
          auto itA = A.find(EdgeKey(alpha,nb));
          PList.push_back(itA!=A.end() ? itA->second : A.at(EdgeKey(nb,alpha)));

          double ni = (double)Ci->size();
          nList.push_back(ni);
          if (std::abs(ni - n/2.0) < 1e-12) { sList.push_back(0.0); kIdx=(int)sList.size()-1; caseNby2=true; }
          else sList.push_back(ni / (n - 2.0*ni));
        }

        if (caseNby2) {
          numerator -= PList[kIdx] / nList[kIdx];
          denominator -= 1.0;
        } else {
          double kappa = 1.0; for (double ni : nList) kappa += ni / (n - 2.0*ni);
          double gamma = 0.0; for (double ni : nList) gamma -= (1.0/kappa) / ((n/ni) - 2.0);
          double sumW = 0.0;
          for (std::size_t i=0;i<nList.size();++i) {
            double wi = (gamma + 1.0) / ((n/nList[i]) - 2.0);
            wList.push_back(wi);
            numerator -= wi * PList[i] / nList[i];
            sumW += wi;
          }
          denominator -= sumW;
        }

      } else {
        double n_alpha = (double)C_alpha.size();
        double n_beta  = (double)C_beta.size();
        denominator = n_alpha * n_beta;

        auto pushNeighbors = [&](int vertex, int other){
          for (int nb : T.neighbors(vertex)) if (!(nb==other)) {
            const auto* Ci = [&]()->const std::unordered_set<int>*{
              auto it = splits.find(EdgeKey(vertex,nb));
              return (it!=splits.end()) ? &it->second.second : &splits.at(EdgeKey(nb,vertex)).first;
            }();
            auto itA = A.find(EdgeKey(vertex,nb));
            PList.push_back(itA!=A.end() ? itA->second : A.at(EdgeKey(nb,vertex)));
            double ni = (double)Ci->size();
            nList.push_back(ni);
            if (std::abs(ni - n/2.0) < 1e-12) { sList.push_back(0.0); kIdx=(int)sList.size()-1; caseNby2=true; }
            else sList.push_back(ni / (n - 2.0*ni));
          }
        };
        pushNeighbors(alpha, beta);
        int mMinusk = (int)T.neighbors(alpha).size() - 1;
        pushNeighbors(beta,  alpha);

        if (caseNby2) {
          std::vector<double> w(sList.size(), 0.0);
          if (kIdx < (int)sList.size() - mMinusk) {
            w[kIdx] = n_beta;
            for (int i=(int)sList.size()-mMinusk; i<(int)sList.size(); ++i) {
              w[i] = (n_alpha - n_beta) / ((n/nList[i]) - 2.0);
              w[kIdx] += (nList[i] * (n_beta - n_alpha)) / (n - 2.0*nList[i]);
              numerator -= w[i] * PList[i] / nList[i];
            }
            numerator -= w[kIdx] * PList[kIdx] / nList[kIdx];
            double sumW = std::accumulate(w.begin(), w.end(), 0.0);
            denominator -= sumW * n_alpha + w[kIdx] * (n_beta - n_alpha);
          } else {
            w[kIdx] = n_alpha;
            for (int i=0; i<(int)sList.size()-mMinusk; ++i) {
              w[i] = (n_beta - n_alpha) / ((n/nList[i]) - 2.0);
              w[kIdx] += (nList[i] * (n_alpha - n_beta)) / (n - 2.0*nList[i]);
              numerator -= w[i] * PList[i] / nList[i];
            }
            numerator -= w[kIdx] * PList[kIdx] / nList[kIdx];
            double sumW = std::accumulate(w.begin(), w.end(), 0.0);
            denominator -= sumW * n_beta + w[kIdx] * (n_alpha - n_beta);
          }
        } else {
          double kappa = 1.0; for (double ni : nList) kappa += ni / (n - 2.0*ni);
          double gamma = 0.0;
          for (int i=0; i<(int)sList.size()-mMinusk; ++i) gamma -= (1.0/kappa) * n_beta  / ((n/nList[i]) - 2.0);
          for (int i=(int)sList.size()-mMinusk; i<(int)sList.size(); ++i) gamma -= (1.0/kappa) * n_alpha / ((n/nList[i]) - 2.0);
          double sumL=0.0, sumR=0.0;
          for (int i=0; i<(int)sList.size()-mMinusk; ++i) {
            double wi = (gamma + n_beta) / ((n/nList[i]) - 2.0);
            wList.push_back(wi); numerator -= wi * PList[i] / nList[i]; sumL += wi;
          }
          for (int i=(int)sList.size()-mMinusk; i<(int)sList.size(); ++i) {
            double wi = (gamma + n_alpha) / ((n/nList[i]) - 2.0);
            wList.push_back(wi); numerator -= wi * PList[i] / nList[i]; sumR += wi;
          }
          denominator -= sumL * n_beta + sumR * n_alpha;
        }
      }

      bl[edge] = numerator / denominator;
    }
    return bl;
  }

  

  static void GetTreeWithLSLengthsLargerThanThreshold(Graph& T, const DistanceMap& DM, double thresh) {
    AssignOLSLengthsToTree(T, DM);
    bool removed = true;
    while (removed) {
      removed = false;
      if (T.edges.empty()) return;
      double mn = std::numeric_limits<double>::infinity();
      for (const auto& e : T.edges) mn = std::min(mn, e.length);
      if (mn > thresh) return;

      std::vector<int> order(T.ecount());
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), [&](int a, int b){ return T.edges[a].length < T.edges[b].length; });

      for (int ei : order) {
        if (ei >= T.ecount()) break;
        if (T.edges[ei].length >= thresh) break;
        int i = T.edges[ei].u, j = T.edges[ei].v;
        bool iLat = IsHiddenNode(T.name[i]), jLat = IsHiddenNode(T.name[j]);
        if (iLat || jLat) {
          removed = true;
          int rm  = iLat ? i : j;
          int keep= iLat ? j : i;

          auto nbs = T.neighbors(rm);
          for (int nb : nbs) if (nb != keep) {
            double edge_length = T.get_length(rm, nb);
            bool exists=false; for (int x : T.neighbors(keep)) if (x==nb) { exists=true; break; }
            if (!exists) T.add_edge(keep, nb, edge_length);
          }
          T.delete_vertex(rm);
          break;
        }
      }
    }
  }

  static Graph GetFJTree(DistanceMap DM, std::vector<std::string> vertList, double epsilon) {
	Graph T = GetFJTreeTopology(std::move(DM), vertList, epsilon);
	AssignOLSLengthsToTree(T, DM);
	SetNegativeLengthsToZero(T);
	return (T);
  }

  static Graph GetFJTreeTopology(DistanceMap DM, std::vector<std::string> vertList, double epsilon) {
    epsilon += 1e-10;
    Graph T;
    T.reserve_vertices(vertList.size());
    for (std::size_t i=0;i<vertList.size();++i) T.name[i] = vertList[i];
    int n = (int)vertList.size();

    std::unordered_map<std::string,double> R;
    for (const auto& nm : vertList) R[nm] = 0.0;
    for (int i=0;i<n;++i) for (int j=i+1;j<n;++j) {
      double d = DM.at(SortedKey(vertList[i], vertList[j]));
      R[vertList[i]] += d; R[vertList[j]] += d;
    }
    for (const auto& nm : vertList) R[nm] /= (n-2.0);

    int nodeCounter = 1;

    auto ensure_idx = [&](const std::string& nm)->int {
      for (int i=0;i<T.vcount();++i) if (T.name[i]==nm) return i;
      return T.add_vertex(nm);
    };

    while ((int)R.size() > 3) {
      double best = std::numeric_limits<double>::infinity();
      std::string isel, jsel;
      for (int i=0;i<(int)vertList.size();++i)
        for (int j=i+1;j<(int)vertList.size();++j) {
          const std::string& a = vertList[i];
          const std::string& b = vertList[j];
          double cand = DM.at(SortedKey(a,b)) - R[a] - R[b];
          if (cand < best) { best=cand; isel=a; jsel=b; }
        }

      double d_ij = DM.at(SortedKey(isel,jsel));
      double delta_ij = std::abs(d_ij + R[isel] - R[jsel]) / 2.0;
      double delta_ji = std::abs(d_ij - R[isel] + R[jsel]) / 2.0;

      if (std::min(delta_ij, delta_ji) <= epsilon) {
        const std::string& parent = (delta_ij < delta_ji) ? isel : jsel;
        const std::string& child  = (delta_ij < delta_ji) ? jsel : isel;
        // if (verbose) std::cerr << "parent-child: " << parent << " <- " << child << "\n";

        int childIndex = (int)(std::find(vertList.begin(), vertList.end(), child) - vertList.begin());
        n -= 1;
        for (int t=0; t<childIndex; ++t) {
          const std::string& v = vertList[t];
          R[v] = (R[v]*(n-1) - DM.at(SortedKey(v,child))) / (n-2.0);
          DM.erase(SortedKey(v,child));
        }
        for (int t=childIndex+1; t<(int)vertList.size(); ++t) {
          const std::string& v = vertList[t];
          R[v] = (R[v]*(n-1) - DM.at(SortedKey(child,v))) / (n-2.0);
          DM.erase(SortedKey(child,v));
        }
        vertList.erase(vertList.begin()+childIndex);
        int ip = ensure_idx(isel), jp = ensure_idx(jsel);
        bool exists = false; for (int x : T.neighbors(ip)) if (x==jp) { exists=true; break; }
        if (!exists) T.add_edge(ip,jp,0.0);

      } else {
        bool ps2 = false; std::string ksel;
        for (const auto& k : vertList) if (k!=isel && k!=jsel) {
          double delta = std::abs(DM.at(SortedKey(k,isel)) + DM.at(SortedKey(k,jsel)) - d_ij) / 2.0;
          if (delta <= epsilon) { ps2 = true; ksel = k; break; }
        }
        if (ps2) {
        //   if (verbose) std::cerr << "parent-siblings: parent " << ksel << " with ("<<isel<<","<<jsel<<")\n";
          vertList.erase(std::remove(vertList.begin(), vertList.end(), isel), vertList.end());
          vertList.erase(std::remove(vertList.begin(), vertList.end(), jsel), vertList.end());
          n -= 2;

          int ki = ensure_idx(ksel), ii = ensure_idx(isel), jj = ensure_idx(jsel);
          bool h1=false,h2=false;
          for (int x : T.neighbors(ki)) if (x==ii) { h1=true; break; }
          if (!h1) T.add_edge(ki,ii,0.0);
          for (int x : T.neighbors(ki)) if (x==jj) { h2=true; break; }
          if (!h2) T.add_edge(ki,jj,0.0);

          R.erase(isel); R.erase(jsel);
          if (n > 2) {
            for (auto& v : vertList) if (v!=ksel) {
              double newD = 0.5 * (DM.at(SortedKey(v, isel)) + DM.at(SortedKey(v, jsel)) - d_ij);
              R[v] = (R[v]*n - DM.at(SortedKey(v,isel)) - DM.at(SortedKey(v,jsel)) + newD) / (n-2.0);
              DM.erase(SortedKey(v, isel));
              DM.erase(SortedKey(v, jsel));
              DM[SortedKey(v, ksel)] = newD;
              R[ksel] += newD;
            }
            R[ksel] /= (n-2.0);
          } else {
            int a = ensure_idx(vertList[0]), b = ensure_idx(vertList[1]);
            T.add_edge(a,b,0.0);
            return T;
          }

        } else {
        //   if (verbose) std::cerr << "siblings: ("<<isel<<","<<jsel<<")\n";
          std::string newNode = "h_" + std::to_string(nodeCounter++);
          int iNew = T.add_vertex(newNode);

          vertList.erase(std::remove(vertList.begin(), vertList.end(), isel), vertList.end());
          vertList.erase(std::remove(vertList.begin(), vertList.end(), jsel), vertList.end());
          int idx = GetInsertIndex(vertList, (int)vertList.size(), newNode);
          n -= 1;
          vertList.insert(vertList.begin()+idx, newNode);

          int ii = ensure_idx(isel), jj = ensure_idx(jsel);
          T.add_edge(iNew, ii, 0.0);
          T.add_edge(iNew, jj, 0.0);

          R.erase(isel); R.erase(jsel);
          R[newNode] = 0.0;

          for (int t=0; t<(int)vertList.size(); ++t) if (t!=idx) {
            const std::string& v = vertList[t];
            double newD = 0.5 * (DM.at(SortedKey(v, isel)) + DM.at(SortedKey(v, jsel)) - DM.at(SortedKey(isel, jsel)));
            R[v] = (R[v]*(n-1) - DM.at(SortedKey(v, isel)) - DM.at(SortedKey(v, jsel)) + newD) / (n-2.0);
            DM.erase(SortedKey(v, isel));
            DM.erase(SortedKey(v, jsel));
            DM[SortedKey(v, newNode)] = newD;
            R[newNode] += newD;
          }
          R[newNode] /= (n-2.0);
        }
        DM.erase(SortedKey(isel, jsel));
      }
    }

    if ((int)vertList.size() == 3) {
      const std::string& v0 = vertList[0];
      const std::string& v1 = vertList[1];
      const std::string& v2 = vertList[2];
      double d01 = DM.at(SortedKey(v0,v1));
      double d02 = DM.at(SortedKey(v0,v2));
      double d12 = DM.at(SortedKey(v1,v2));
      double D0 = std::abs(d01 + d02 - d12) / 2.0;
      double D1 = std::abs(d01 + d12 - d02) / 2.0;
      double D2 = std::abs(d02 + d12 - d01) / 2.0;
      int i0 = -1, i1=-1, i2=-1;
      for (int i=0;i<T.vcount();++i) { if (T.name[i]==v0) i0=i; if (T.name[i]==v1) i1=i; if (T.name[i]==v2) i2=i; }
      if (std::min({D0,D1,D2}) <= epsilon) {
        if (D0 < D1 && D0 < D2) { T.add_edge(i0,i1,0.0); T.add_edge(i0,i2,0.0); }
        else if (D1 < D0 && D1 < D2) { T.add_edge(i0,i1,0.0); T.add_edge(i1,i2,0.0); }
        else { T.add_edge(i0,i2,0.0); T.add_edge(i1,i2,0.0); }
      } else {
        std::string newNode = "h_" + std::to_string(nodeCounter++);
        int x = T.add_vertex(newNode);
        T.add_edge(i0,x,0.0); T.add_edge(i1,x,0.0); T.add_edge(i2,x,0.0);
      }
    } else if ((int)vertList.size() == 2) {
      int a=-1,b=-1; for (int i=0;i<T.vcount();++i) { if (T.name[i]==vertList[0]) a=i; if (T.name[i]==vertList[1]) b=i; }
      T.add_edge(a,b,0.0);
    }
    return T;
  }

  // ------------------------ Newick helpers ------------------------	

  Graph ConvertGenerallyLabeledTreeToLeafLabeledTree(Graph tree) const {
    std::vector<int> deg(tree.vcount());
    for (int v=0; v<tree.vcount(); ++v) deg[v] = tree.degree(v);

    std::vector<int> internalLabels; internalLabels.reserve(tree.vcount());
    int hidden_node_count = 1;
    for (int v=0; v<tree.vcount(); ++v) {
      if (IsHiddenNode(tree.name[v])) { ++hidden_node_count; continue; }
      if (deg[v] > 1) internalLabels.push_back(v);
    }

    for (int vidx=0; vidx<(int)internalLabels.size(); ++vidx) {
      std::string vertexName = tree.name[internalLabels[vidx]];
      int v = -1; for (int i=0;i<tree.vcount();++i) if (tree.name[i]==vertexName) { v=i; break; }

      std::string newNode = "h_" + std::to_string(hidden_node_count++);
      int t = tree.add_vertex(newNode);

      auto nbs = tree.neighbors(v);
      for (int nb : nbs) {
        double edge_length = tree.get_length(v, nb);
        tree.add_edge(nb, t, edge_length);
        tree.remove_edge(v, nb);
      }
      tree.add_edge(v, t, 0.0);
    }
    return tree;
  }

  Graph ConvertMultifurcatingTreeToBifurcatingTree(Graph tree) const {
    std::vector<int> initialDeg(tree.vcount());
    for (int v=0; v<tree.vcount(); ++v) initialDeg[v] = tree.degree(v);

    int hidden_node_count = 1;
    for (int v=0; v<tree.vcount(); ++v)
      if (IsHiddenNode(tree.name[v])) ++hidden_node_count;

    std::vector<std::string> needResolve;
    for (int v=0; v<tree.vcount(); ++v)
      if (initialDeg[v] > 3) needResolve.push_back(tree.name[v]);

    for (const auto& vname : needResolve) {
      while (true) {
        int v=-1; for (int i=0;i<tree.vcount();++i) if (tree.name[i]==vname) { v=i; break; }
        if (v==-1) break;
        if (tree.degree(v) <= 3) break;

        auto nbs = tree.neighbors(v);
        int v1 = nbs[0], v2 = nbs[1];
        std::string newNode = "h_" + std::to_string(hidden_node_count++);
        int t = tree.add_vertex(newNode);

        double L1 = tree.get_length(v, v1);
        double L2 = tree.get_length(v, v2);
        tree.add_edge(v1, t, L1);
        tree.add_edge(v2, t, L2);
        tree.add_edge(v,  t, 0.0);
        tree.remove_edge(v, v1);
        tree.remove_edge(v, v2);
      }
    }
    return tree;
  }

  std::string GetNewickLabelOfLeafLabeledTree(Graph tree) const {
    int leaf=-1; for (int i=0;i<tree.vcount();++i) if (tree.degree(i)==1) { leaf=i; break; }
    if (leaf==-1) throw std::runtime_error("Cannot find a leaf to root.");

    int parent = tree.neighbors(leaf)[0];
    double edge_length = tree.get_length(leaf, parent);
    int root = tree.add_vertex("h_root");
    tree.add_edge(leaf, root, edge_length/2.0);
    tree.add_edge(parent, root, edge_length/2.0);
    tree.remove_edge(leaf, parent);

    std::function<std::string(int,int)> dfs = [&](int u, int p)->std::string {
      std::vector<int> children;
      for (int nb : tree.neighbors(u)) if (nb!=p) children.push_back(nb);

      auto labelFor = [&](int node)->std::string {
        const std::string& nm = tree.name[node];
        if (IsHiddenNode(nm) || nm=="h_root") return std::string();
        return nm;
      };

      if (p == -1) {
        if (children.empty()) {
          return (labelFor(u).empty() ? std::string("") : labelFor(u)) + ";";
        }
        std::string inside;
        for (std::size_t i=0;i<children.size();++i) {
          if (i) inside += ",";
          inside += dfs(children[i], u);
        }
        return "(" + inside + ");";
      } else if (children.empty()) {
        std::string nm = labelFor(u);
        double len = tree.get_length(u, p);
        return (nm.empty() ? std::string("") : nm) + ":" + to_string_trim(len);
      } else {
        std::string inside;
        for (std::size_t i=0;i<children.size();++i) {
          if (i) inside += ",";
          inside += dfs(children[i], u);
        }
        std::string nm = labelFor(u);
        double len = tree.get_length(u, p);
        return "(" + inside + ")" + (nm.empty() ? std::string("") : nm) + ":" + to_string_trim(len);
      }
    };

    return dfs(root, -1);
  }

  // Convert a Graph to Newick and strip hidden_node labels; ensure trailing ';'.
  std::string BuildNewickStringFromGraph(const Graph& tree) const {
    Graph copy = tree;
    Graph leaf = ConvertGenerallyLabeledTreeToLeafLabeledTree(copy);
    Graph bin  = ConvertMultifurcatingTreeToBifurcatingTree(leaf);
    std::string nwk = GetNewickLabelOfLeafLabeledTree(bin);    
    nwk = this->strip_hidden_node_labels_(nwk);
    if (nwk.empty() || nwk.back()!=';') nwk.push_back(';');
    return nwk;
  }

  static std::string to_string_trim(double x) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed); oss.precision(10);
    oss << x;
    std::string s = oss.str();
    while (!s.empty() && s.back()=='0') s.pop_back();
    if (!s.empty() && s.back()=='.') s.push_back('0');
    return s;
  }
};


///...///...///...///...///...///...///... mst vertex ...///...///...///...///...///...///...///

class MST_vertex {
public:
	string name;
	int degree = 0;
	int numberOfLargeEdgesInSubtree = 0;
	int timesVisited = 0;
	int id;
	int idOfExternalVertex;
	int rank = 0;	
	vector <int> sequence;
	vector <int> globallyCompressedSequence;
	vector <int> idsOfVerticesInSubtree;
	vector <MST_vertex *> neighbors;
	vector <string> dupl_seq_names;
	void AddNeighbor(MST_vertex * v_ptr);		
	MST_vertex(int idToAdd, string nameToAdd, vector <int> sequenceToAdd) {
		id = idToAdd;
		sequence = sequenceToAdd;
		name = nameToAdd;
		idsOfVerticesInSubtree.push_back(idToAdd);
		idOfExternalVertex = -1;
	}
	~MST_vertex() {
		
	}
};

void MST_vertex::AddNeighbor(MST_vertex* v_ptr) {
	degree += 1;
	neighbors.push_back(v_ptr);
}

///...///...///...///...///...///...///...///... mst ...///...///...///...///...///...///...///...///

class MST {
private:	
	int largestVertexIndex;
	int edgeWeightThreshold = 0;
	vector <MST_vertex*> verticesToVisit;	
	bool ContainsVertex(int vertex_id);
	map <pair<int,int>,int> * allEdgeWeights;
	chrono::system_clock::time_point current_time;
	chrono::system_clock::time_point start_time;
	chrono::system_clock::time_point time_to_compute_MST;	
		
public:
	int maxDegree;	
	unsigned int num_duplicated_sequences = 0;
    int sequence_length;
	vector <int> siteWeights;
	string sequenceFileName;
	int v_ind;
	int numberOfLargeEdgesThreshold_input = 0;
	int numberOfLargeEdgesThreshold = 0;
	int numberOfNonZeroWeightEdges = 0;
    int numberOfInputSequences = 0;
	vector <int> idsOfExternalVertices;
	vector <MST_vertex *> leader_ptrs;  
	void SetLeaders();
	map <string, int> mapDNAtoInteger;	
	map <string, vector <string>> unique_seq_id_2_dupl_seq_ids;
	map <vector <int>, MST_vertex *> unique_seq_2_MST_vertex_ptr;
	MST_vertex * subtree_v_ptr;
	map <int, MST_vertex *> * vertexMap;
	map <pair <int, int>, int> edgeWeightsMap;
	string EncodeAsDNA(vector<int> sequence);
	vector<int> DecompressSequence(vector<int>* compressedSequence, vector<vector<int>>* sitePatternRepeats);
	void AddEdgeWeight(int u_id, int v_id, int edgeWeight);
	void RemoveEdgeWeight(int u_id, int v_id);
	void AddEdgeWeightToDistanceGraph(int u_id, int v_id, int edgeWeight);
	void RemoveWeightedEdgesIncidentToVertexInDistanceGraph(int u_id);
	void SetCompressedSequencesAndSiteWeightsForInputSequences();
	bool IsSequenceDuplicated(vector<int> sequence);
	void AddDuplicatedSequenceName(string name, vector<int> sequence);
	void SetNumberOfLargeEdgesThreshold(int numberOfLargeEdges);
	void SetEdgeWeightThreshold(int edgeWeight){edgeWeightThreshold = edgeWeight;}
	void AddVertex(string name, vector <int> sequence);
	void AddVertexWithId(int id, string name, vector <int> sequence);
	void RemoveVertex(int vertex_id);
	void AddEdge(int u_id, int v_id, int edgeWeight);
	void RemoveEdge(int u_id, int v_id);
	void ResetVertexAttributesForSubtreeSearch();
	void UpdateMSTWithMultipleExternalVertices(vector <int> idsOfVerticesToKeep, vector <int> idsOfVerticesToRemove, vector <tuple<int,string,vector<int>>> idAndNameAndSeqTupleForVerticesToAdd, vector <int> idsOfExternalVertices);
	void UpdateMaxDegree();
	void UpdateMSTWithOneExternalVertex(vector <int> idsOfVerticesToRemove, string nameOfSequenceToAdd, vector <int> sequenceToAdd);
	bool ContainsEdge(int u_id, int v_id);
	int GetEdgeWeight(int u_id, int v_id);
	int GetEdgeIndex (int vertexIndex1, int vertexIndex2, int numberOfVertices);
	int GetNumberOfVertices();
	void ReadFasta(string sequenceFileNameToSet);
    void ReadPhyx(string sequenceFileNameToSet);
	void ComputeMST();
	void CLGrouping();		
	void ResetSubtreeSizeThreshold();	
	void doubleSubtreeSizeThreshold();
	int ComputeHammingDistance(vector <int> recodedSeq1, vector <int> recodedSeq2);
	pair <vector <int>, vector <int>> GetIdsForSubtreeVerticesAndExternalVertices();
	pair <bool, MST_vertex *> GetPtrToVertexSubtendingSubtree();
	pair <vector <int>,vector <int>> GetSubtreeVerticesAndExternalVertices();
	tuple <vector <string>, vector <vector <int>>, vector <int>, vector <vector<int>>> GetCompressedSequencesSiteWeightsAndSiteRepeats(vector <int> vertexIdList);
	vector <int> GetIdsOfClosestUnvisitedVertices(MST_vertex* u_ptr);
	void SetIdsOfExternalVertices();
	bool ShouldIComputeALocalPhylogeneticTree();
	void WriteToFile(string fileName);
	int ConvertDNAToChar(char dna);	
	MST() {		
		this->v_ind = 0;		
		vector <int> emptySequence;
		this->allEdgeWeights = new map <pair<int,int>,int> ; 
		this->vertexMap = new map <int, MST_vertex *>;
		this->mapDNAtoInteger["A"] = 0;
		this->mapDNAtoInteger["C"] = 1;
		this->mapDNAtoInteger["G"] = 2;
		this->mapDNAtoInteger["T"] = 3;		
		this->mapDNAtoInteger["-"] = 4;
		this->mapDNAtoInteger["N"] = 4;
		this->mapDNAtoInteger["W"] = 4;
		this->mapDNAtoInteger["S"] = 4;
		this->mapDNAtoInteger["M"] = 4;
		this->mapDNAtoInteger["K"] = 4;
		this->mapDNAtoInteger["R"] = 4;
		this->mapDNAtoInteger["Y"] = 4;
		this->mapDNAtoInteger["B"] = 4;
		this->mapDNAtoInteger["D"] = 4;
		this->mapDNAtoInteger["H"] = 4;
		this->mapDNAtoInteger["V"] = 4;		
	
	}
	~MST() {		
		for (pair<int,MST_vertex*> VptrMap: *this->vertexMap){			
			delete VptrMap.second;
		}
		delete this->vertexMap;
		delete this->allEdgeWeights;
	}
};

void MST::SetLeaders() {
	leader_ptrs.clear();
	for (pair<int,MST_vertex*> VptrMap: *this->vertexMap) {					
		if (VptrMap.second->degree > 1) {
			leader_ptrs.push_back(VptrMap.second);
		}			
	}
}

bool MST::IsSequenceDuplicated(vector<int> query_seq) {
	if (this->unique_seq_2_MST_vertex_ptr.find(query_seq) != this->unique_seq_2_MST_vertex_ptr.end()) {
		return (true);
	} else {
		return (false);
	}
}

void MST::AddDuplicatedSequenceName(string dupl_seq_name, vector <int> sequence) {	
	MST_vertex * v = this->unique_seq_2_MST_vertex_ptr[sequence];
	this->unique_seq_id_2_dupl_seq_ids[v->name].push_back(dupl_seq_name);
	this->num_duplicated_sequences += 1;
}

void MST::AddEdgeWeightToDistanceGraph(int u_id, int v_id, int edgeWeight) {
	if (u_id < v_id) {
		this->allEdgeWeights->insert(make_pair(make_pair(u_id,v_id),edgeWeight));
	} else {
		this->allEdgeWeights->insert(make_pair(make_pair(v_id,u_id),edgeWeight));
	}
}

void MST::SetIdsOfExternalVertices() {
	this->idsOfExternalVertices.clear();
	this->idsOfExternalVertices = this->GetIdsOfClosestUnvisitedVertices(this->subtree_v_ptr);
}

int MST::GetEdgeIndex (int vertexIndex1, int vertexIndex2, int numberOfVertices) {
	int edgeIndex;
	edgeIndex = numberOfVertices * (numberOfVertices-1)/2;
	edgeIndex -= (numberOfVertices - vertexIndex1) * (numberOfVertices-vertexIndex1-1)/2;
	edgeIndex += vertexIndex2 - vertexIndex1 - 1;
	return edgeIndex;
}

void MST::SetNumberOfLargeEdgesThreshold(int numberOfLargeEdges_toSet) {
	this->numberOfLargeEdgesThreshold_input = numberOfLargeEdges_toSet;
	this->numberOfLargeEdgesThreshold = numberOfLargeEdges_toSet;
}

bool MST::ShouldIComputeALocalPhylogeneticTree() {
	bool valueToReturn;
	bool verbose = 1;
	bool subtreeExtractionPossible;
	int numberOfNonZeroWeightEdgesInVmWithoutVs;
	tie (subtreeExtractionPossible, this->subtree_v_ptr) = this->GetPtrToVertexSubtendingSubtree();	
	if (subtreeExtractionPossible) {
		numberOfNonZeroWeightEdgesInVmWithoutVs = this->numberOfNonZeroWeightEdges - this->subtree_v_ptr->numberOfLargeEdgesInSubtree;
		if (numberOfNonZeroWeightEdgesInVmWithoutVs > this->numberOfLargeEdgesThreshold) {
			valueToReturn = 1;			
		} else {
			// if (verbose) {
			// 	cout << "Case 1: subtree extraction possible but number of external vertices is too small" << endl;
			// }
			valueToReturn = 0;
			
		}
	} else {
		// if (verbose) {
		// 	cout << "Case 2: subtree extraction is not possible" << endl;
		// }
		valueToReturn = 0;
	}
	return (valueToReturn);
}



void MST::ResetSubtreeSizeThreshold() {
	this->numberOfLargeEdgesThreshold = this->numberOfLargeEdgesThreshold_input;
}

void MST::doubleSubtreeSizeThreshold() {
	this->numberOfLargeEdgesThreshold = this->numberOfLargeEdgesThreshold * 2;
}

int MST::GetNumberOfVertices() {
	return this->vertexMap->size();
};


bool MST::ContainsVertex(int vertex_id) {
	return this->vertexMap->find(vertex_id)!=vertexMap->end();
}

void MST::AddVertex(string name, vector <int> sequence) {
	MST_vertex * v = new MST_vertex(this->v_ind, name, sequence);
	this->vertexMap->insert(pair<int,MST_vertex*>(this->v_ind,v));
	this->unique_seq_2_MST_vertex_ptr[sequence] = v;	
	this->v_ind += 1;
}

void MST::AddVertexWithId(int id, string name, vector <int> sequence) {
	MST_vertex * v = new MST_vertex(id, name, sequence);
	this->vertexMap->insert(pair<int,MST_vertex*>(id,v));
}

void MST::RemoveVertex(int vertex_id) {
	MST_vertex* v = (*this->vertexMap)[vertex_id];
	for (MST_vertex* n: v->neighbors) {
		if (n->id < v->id) {
			this->edgeWeightsMap.erase(pair<int,int>(n->id,v->id));
		} else {
			this->edgeWeightsMap.erase(pair<int,int>(v->id,n->id));
		}
		n->neighbors.erase(remove(n->neighbors.begin(),n->neighbors.end(),v),n->neighbors.end());
		n->degree -= 1;
	}
	v->neighbors.clear();
	v->sequence.clear();
	v->idsOfVerticesInSubtree.clear();
	this->vertexMap->erase(vertex_id);
	delete v;
}

void MST::AddEdge(int u_id, int v_id, int edgeWeight) {
	MST_vertex* u_ptr = (*this->vertexMap)[u_id];
	MST_vertex* v_ptr = (*this->vertexMap)[v_id];
	u_ptr->AddNeighbor(v_ptr);
	v_ptr->AddNeighbor(u_ptr);
	this->AddEdgeWeight(u_id,v_id,edgeWeight);
	if (edgeWeight > 0) {
		this->numberOfNonZeroWeightEdges += 1;
	}
};

void MST::RemoveEdge(int u_id, int v_id) {	
	if (u_id < v_id) {
		this->edgeWeightsMap.erase(pair<int,int>(u_id, v_id));
	} else {
		this->edgeWeightsMap.erase(pair<int,int>(v_id, u_id));
	}
	MST_vertex * u = (*this->vertexMap)[u_id];
	MST_vertex * v = (*this->vertexMap)[v_id];
	u->neighbors.erase(remove(u->neighbors.begin(),u->neighbors.end(),v),u->neighbors.end());
	u->degree -= 1;
	v->neighbors.erase(remove(v->neighbors.begin(),v->neighbors.end(),u),v->neighbors.end());
	v->degree -= 1;
}

bool MST::ContainsEdge(int u_id, int v_id) {
	if (u_id < v_id) {
		return (this->edgeWeightsMap.find(pair<int,int>(u_id,v_id)) != this->edgeWeightsMap.end());
	} else {
		return (this->edgeWeightsMap.find(pair<int,int>(v_id,u_id)) != this->edgeWeightsMap.end());
	}	
}

int MST::GetEdgeWeight(int u_id, int v_id) {
	if (u_id < v_id) {
		return this->edgeWeightsMap[pair<int,int>(u_id,v_id)];
	} else {
		return this->edgeWeightsMap[pair<int,int>(v_id,u_id)];
	}
}

void MST::AddEdgeWeight(int u_id, int v_id, int edgeWeight) {
	pair<int,int> edge ;
	if (u_id < v_id){
		edge = make_pair(u_id,v_id);
	} else {
		edge = make_pair(v_id,u_id);
	}
	if (this->edgeWeightsMap.find(edge) != this->edgeWeightsMap.end()) {
		this->edgeWeightsMap[edge] = edgeWeight;
	} else {
		this->edgeWeightsMap.insert(make_pair(edge,edgeWeight));
	}	
}

void MST::RemoveEdgeWeight(int u_id, int v_id) {
	pair <int, int> edge;
	if (u_id < v_id){
		edge = make_pair(u_id,v_id);
	} else {
		edge = make_pair(v_id,u_id);
	}
	this->edgeWeightsMap.erase(edge);
}

vector<int> MST::GetIdsOfClosestUnvisitedVertices(MST_vertex* v_ptr) {
	int numberOfLargeEdgesEncountered = 0;
	vector <int> idsOfClosestUnvisitedVertices;
	vector <MST_vertex*> verticesInCurrentLevel;
	for (MST_vertex* n_ptr: v_ptr->neighbors) {		
		if (find(v_ptr->idsOfVerticesInSubtree.begin(),v_ptr->idsOfVerticesInSubtree.end(),n_ptr->id)==v_ptr->idsOfVerticesInSubtree.end()) {
			idsOfClosestUnvisitedVertices.push_back(n_ptr->id);
			if (this->GetEdgeWeight(v_ptr->id,n_ptr->id)  > edgeWeightThreshold) {
				numberOfLargeEdgesEncountered+=1;
			}
			if (numberOfLargeEdgesEncountered < numberOfLargeEdgesThreshold) {
				verticesInCurrentLevel.push_back(n_ptr);
			}
		}
	}
	vector <MST_vertex *> verticesInNextLevel;
	while (numberOfLargeEdgesEncountered < numberOfLargeEdgesThreshold and verticesInCurrentLevel.size() > 0) {
		for(MST_vertex * x_ptr:verticesInCurrentLevel) {
			for (MST_vertex * n_ptr : x_ptr->neighbors) {
				if (find(idsOfClosestUnvisitedVertices.begin(),idsOfClosestUnvisitedVertices.end(),n_ptr->id)==idsOfClosestUnvisitedVertices.end() and n_ptr->id!=v_ptr->id) {
					idsOfClosestUnvisitedVertices.push_back(n_ptr->id);
					if(this->GetEdgeWeight(x_ptr->id,n_ptr->id) > edgeWeightThreshold) {
						numberOfLargeEdgesEncountered+=1;
					}
					if (numberOfLargeEdgesEncountered < numberOfLargeEdgesThreshold) {
						verticesInNextLevel.push_back(n_ptr);	
					}
				}
			}
		}
		verticesInCurrentLevel = verticesInNextLevel;
		verticesInNextLevel.clear();
	}
	return idsOfClosestUnvisitedVertices;
}

pair <bool, MST_vertex *> MST::GetPtrToVertexSubtendingSubtree() {
	this->ResetVertexAttributesForSubtreeSearch();	
	bool subTreeFound = 0;
	verticesToVisit.clear();
	for (pair <int, MST_vertex *> VptrMap: * this->vertexMap) {
		if (VptrMap.second->degree == 1) {
			verticesToVisit.push_back(VptrMap.second);
		}
	}
	vector <MST_vertex *> verticesVisited;
	int vertex_ind = verticesToVisit.size() -1;	
	while(vertex_ind > -1 and !subTreeFound) {
		this->subtree_v_ptr = verticesToVisit[vertex_ind];
		verticesToVisit.pop_back();
		vertex_ind -= 1;
		this->subtree_v_ptr->timesVisited += 1;
		for (MST_vertex * neighbor_ptr : this->subtree_v_ptr->neighbors) {
			if (neighbor_ptr->timesVisited < neighbor_ptr->degree) {
				neighbor_ptr->timesVisited += 1;
				for (int n_id : this->subtree_v_ptr->idsOfVerticesInSubtree) {
					neighbor_ptr->idsOfVerticesInSubtree.push_back(n_id);
				}
				neighbor_ptr->numberOfLargeEdgesInSubtree += this->subtree_v_ptr->numberOfLargeEdgesInSubtree;
				if (GetEdgeWeight(this->subtree_v_ptr->id,neighbor_ptr->id) > edgeWeightThreshold) {
					neighbor_ptr->numberOfLargeEdgesInSubtree+=1;
				}
				if (neighbor_ptr->degree - neighbor_ptr->timesVisited == 1) {
					if (neighbor_ptr->numberOfLargeEdgesInSubtree > numberOfLargeEdgesThreshold) {
						subTreeFound = 1;
						// set id to external vertex
						for (MST_vertex * v : neighbor_ptr->neighbors) {
							if (v->timesVisited < v->degree) {
								neighbor_ptr->idOfExternalVertex = v->id;
							}							
						}						
						return pair <bool, MST_vertex *> (subTreeFound,neighbor_ptr);						
					}
					verticesToVisit.push_back(neighbor_ptr);
					vertex_ind+=1;
				}
			}
		}
	}	
	return pair <bool, MST_vertex *> (subTreeFound,this->subtree_v_ptr);
}

void MST::ResetVertexAttributesForSubtreeSearch() {
	for (pair <int, MST_vertex *> VIdAndPtr: *this->vertexMap) {
		VIdAndPtr.second->numberOfLargeEdgesInSubtree = 0;
		VIdAndPtr.second->idsOfVerticesInSubtree.clear();
		VIdAndPtr.second->idsOfVerticesInSubtree.push_back(VIdAndPtr.second->id);
		VIdAndPtr.second->timesVisited = 0;
	}
}


void MST::UpdateMSTWithOneExternalVertex(vector<int> idsOfVerticesToRemove, string nameOfSequenceToAdd, vector <int> sequenceToAdd) {
	//	Remove vertices		
	for (int v_id: idsOfVerticesToRemove) {	
		this->RemoveVertex(v_id);
	}
	//	Remove neighbors and reset vertex attributes
	for (pair<int,MST_vertex*> VIdAndPtr: *this->vertexMap) {
		VIdAndPtr.second->numberOfLargeEdgesInSubtree = 0;
		VIdAndPtr.second->idsOfVerticesInSubtree.clear();
		VIdAndPtr.second->idsOfVerticesInSubtree.push_back(VIdAndPtr.second->id);
		VIdAndPtr.second->timesVisited = 0;
		VIdAndPtr.second->neighbors.clear();
		VIdAndPtr.second->degree = 0;
	}
	this->numberOfNonZeroWeightEdges = 0;
	// Add vertex
	int indexOfVertexToAdd = this->v_ind;
	this->AddVertex(nameOfSequenceToAdd, sequenceToAdd);
	MST_vertex * v_add; MST_vertex * v_inMST;
	v_add = (*this->vertexMap)[indexOfVertexToAdd];
	int edgeWeight;
	for (pair <int, MST_vertex *> idPtrPair : *this->vertexMap) {
		if (idPtrPair.first != indexOfVertexToAdd) {
			v_inMST = idPtrPair.second;
			edgeWeight = ComputeHammingDistance(v_add->sequence, v_inMST->sequence);
			if (v_inMST->id < v_add->id){					
				this->edgeWeightsMap[pair<int,int>(v_inMST->id,v_add->id)] = edgeWeight;					
			} else {					
				this->edgeWeightsMap[pair<int,int>(v_add->id,v_inMST->id)] = edgeWeight;					
			}
		}
	}	
	int numberOfVertices = int(this->vertexMap->size());	
	const int numberOfEdges = int(this->edgeWeightsMap.size());
	
	double * weights;
	weights = new double [numberOfEdges];
	
	typedef pair <int,int > E;
	E * edges;
	edges = new E [numberOfEdges];
	
	int edgeIndex = 0;
    for (pair<pair<int,int>,int> edgeAndWeight : this->edgeWeightsMap) {
		edges[edgeIndex] = E(edgeAndWeight.first.first,edgeAndWeight.first.second);
		weights[edgeIndex] = edgeAndWeight.second;
		edgeIndex += 1;		
	}
	
	vector<int> p(numberOfVertices); 

	emtr::prim_graph p_graph(numberOfVertices, edges, weights, numberOfEdges);
	
	emtr::prim(p_graph, &p[0]);
	
	delete[] edges;		
	delete[] weights;		
	vector <pair<int,int>> edgeWeightsToKeep;	
	vector <pair<int,int>> edgeWeightsToRemove;	
	for (size_t u = 0; u != p.size(); ++u) {
		if (p[u] != u) {		
			this->AddEdge(u,p[u],this->GetEdgeWeight(u,p[u]));
			if (u < p[u]) {
				edgeWeightsToKeep.push_back(pair<int,int>(u,p[u]));
			} else {
				edgeWeightsToKeep.push_back(pair<int,int>(p[u],u));
			}
		}
	}
	for (pair<pair<int,int>,int> edgeAndWeight : this->edgeWeightsMap) {
		if (find(edgeWeightsToKeep.begin(),edgeWeightsToKeep.end(),edgeAndWeight.first) == edgeWeightsToKeep.end()){
			edgeWeightsToRemove.push_back(edgeAndWeight.first);
		}
	}
	for (pair<int,int> edge: edgeWeightsToRemove) {
		this->edgeWeightsMap.erase(edge);
	}	
}

void MST::UpdateMaxDegree() {
	this->maxDegree = 0;
	for (pair <int,MST_vertex*> VIdAndPtr: *this->vertexMap) {
		if (this->maxDegree	< VIdAndPtr.second->degree) {
			this->maxDegree	= VIdAndPtr.second->degree;
		}
	}
}

void MST::UpdateMSTWithMultipleExternalVertices(vector<int> idsOfVerticesToKeep, vector<int> idsOfVerticesToRemove, vector<tuple<int,string,vector<int>>> idAndNameAndSeqTupleForVerticesToAdd, vector<int> idsOfExternalVertices) {
	// Remove weights of edges incident to vertex
	//	Remove vertices		
	for (int v_id: idsOfVerticesToRemove) {
		this->RemoveVertex(v_id);
	}
	//	Remove all edges in MST and reset all attributes for each vertex
	for (pair <int,MST_vertex*> VIdAndPtr: *this->vertexMap) {
		VIdAndPtr.second->numberOfLargeEdgesInSubtree = 0;
		VIdAndPtr.second->idsOfVerticesInSubtree.clear();
		VIdAndPtr.second->idsOfVerticesInSubtree.push_back(VIdAndPtr.second->id);
		VIdAndPtr.second->timesVisited = 0;
		VIdAndPtr.second->neighbors.clear();
		VIdAndPtr.second->degree = 0;
	}
	this->numberOfNonZeroWeightEdges = 0;
	int u_id; int v_id; int edgeWeight;
	vector <int> seq_u; vector <int> seq_v;
	string u_name; string v_name;
	
	int numberOfVerticesToKeep = int(idsOfVerticesToKeep.size());
	int numberOfVerticesToAdd = int(idAndNameAndSeqTupleForVerticesToAdd.size());
	int numberOfExternalVertices = int(idsOfExternalVertices.size());
	
	if (numberOfVerticesToAdd > 1) {
		for (int u_ind = 0; u_ind < numberOfVerticesToAdd -1; u_ind++) {
			tie (u_id, u_name, seq_u) = idAndNameAndSeqTupleForVerticesToAdd[u_ind];
			for (int v_ind = u_ind + 1; v_ind < numberOfVerticesToAdd; v_ind++) {
				tie (v_id, v_name, seq_v) = idAndNameAndSeqTupleForVerticesToAdd[v_ind];
				edgeWeight = ComputeHammingDistance(seq_u, seq_v);
				this->AddEdgeWeight(u_id,v_id,edgeWeight);
			}
		}
	}
// Add newly introduced vertices
	for (int u_ind = 0; u_ind < numberOfVerticesToAdd; u_ind++) {
		tie (u_id, u_name, seq_u) = idAndNameAndSeqTupleForVerticesToAdd[u_ind];
		this->AddVertexWithId(u_id, u_name, seq_u);
	}
// Add edge weights for vertices to add to external vertices
	for (int u_ind = 0; u_ind < numberOfVerticesToAdd; u_ind++) {
		tie (u_id, u_name, seq_u) = idAndNameAndSeqTupleForVerticesToAdd[u_ind];
		for (int v_ind = 0; v_ind < numberOfExternalVertices; v_ind++) {
			v_id = idsOfExternalVertices[v_ind];
			seq_v = (*this->vertexMap)[v_id]->sequence;
			edgeWeight = ComputeHammingDistance(seq_u,seq_v);
			this->AddEdgeWeight(u_id,v_id,edgeWeight);
		}
	}
		
	for (int u_ind = 0; u_ind < numberOfVerticesToAdd; u_ind++) {
		tie (u_id, u_name, seq_u) = idAndNameAndSeqTupleForVerticesToAdd[u_ind];		
		for (int v_ind = 0; v_ind < numberOfVerticesToKeep; v_ind++) {
			v_id = idsOfVerticesToKeep[v_ind];
			seq_v = (*this->vertexMap)[v_id]->sequence;
			edgeWeight = ComputeHammingDistance(seq_u,seq_v);
			this->AddEdgeWeight(u_id,v_id,edgeWeight);
		}		
	}

	for (int u_ind = 0; u_ind < numberOfVerticesToKeep; u_ind++) {
			u_id = idsOfVerticesToKeep[u_ind];
			seq_u = (*this->vertexMap)[u_id]->sequence;
		for (int v_ind = 0; v_ind < numberOfExternalVertices; v_ind++) {
			v_id = idsOfExternalVertices[v_ind];
			seq_v = (*this->vertexMap)[v_id]->sequence;
			edgeWeight = ComputeHammingDistance(seq_u,seq_v);
			this->AddEdgeWeight(u_id,v_id,edgeWeight);
		}
	}
		
	vector <int> mstIds;
	map<int,int> mstId2PrimId;
	int primId = 0;
	int mstId;
	for (pair <int,MST_vertex*> idPtrPair : *this->vertexMap) {
		mstId = idPtrPair.first;
		mstIds.push_back(mstId);
		mstId2PrimId.insert(make_pair(mstId,primId));
		primId += 1;
	}
	int numberOfVertices = int(this->vertexMap->size());	

	const int numberOfEdges = int(this->edgeWeightsMap.size());
	double * weights;
	weights = new double [numberOfEdges];
	
	typedef pair <int,int > E;
	E * edges;
	edges = new E [numberOfEdges];
		
	
	int edgeIndex = 0;	
    for (pair<pair<int,int>,int> edgeAndWeight : this->edgeWeightsMap) {
		tie (u_id, v_id) = edgeAndWeight.first;
		edges[edgeIndex] = E(mstId2PrimId[u_id],mstId2PrimId[v_id]);
		weights[edgeIndex] = edgeAndWeight.second;
		edgeIndex += 1;		
	}
	
	vector<int> p(numberOfVertices); 

	emtr::prim_graph p_graph(numberOfVertices, edges, weights, numberOfEdges);
	
	emtr::prim(p_graph, &p[0]);
				
	vector <pair<int,int>> edgeWeightsToKeep;	
	vector <pair<int,int>> edgeWeightsToRemove;	

	for (size_t u = 0; u != p.size(); ++u) {
		if (p[u] != u){
			u_id = mstIds[p[u]];
			v_id = mstIds[u];
			this->AddEdge(u_id,v_id,this->GetEdgeWeight(u_id,v_id));
			if (u_id < v_id){
				edgeWeightsToKeep.push_back(pair<int,int>(u_id,v_id));
			} else {
				edgeWeightsToKeep.push_back(pair<int,int>(v_id,u_id));
			}
		}
	}
	this->UpdateMaxDegree();
	
	if (this->maxDegree == 0) {
		ofstream edgeListFile;
		edgeListFile.open(sequenceFileName+".debugEdgeList");
		for (pair<pair<int,int>,int> edgeAndWeight : this->edgeWeightsMap) {
			edgeListFile << edgeAndWeight.first.first << "\t";
			edgeListFile << edgeAndWeight.first.second << "\t";
			edgeListFile << edgeAndWeight.second << endl;
		}
		edgeListFile.close();
	}
	
	delete[] edges;
	delete[] weights;
	for (pair<pair<int,int>,int> edgeAndWeight : this->edgeWeightsMap) {
		if (find(edgeWeightsToKeep.begin(),edgeWeightsToKeep.end(),edgeAndWeight.first) == edgeWeightsToKeep.end()){
			edgeWeightsToRemove.push_back(edgeAndWeight.first);
		}
	}
	for (pair<int,int> edge: edgeWeightsToRemove) {
		this->edgeWeightsMap.erase(edge);
	}
}

tuple <vector<string>,vector<vector<int>>,vector<int>,vector<vector<int>>> MST::GetCompressedSequencesSiteWeightsAndSiteRepeats(vector<int> vertexIdList){	
	vector <string> names;
	vector <vector<int>> compressedSequences;
	vector <int> sitePatternWeights_ptr;
	vector <vector <int>> sitePatternRepeats_ptr;
	vector <vector<int>> distinctPatterns;
	map <vector<int>,vector<int>> distinctPatternsToSitesWherePatternRepeats;
	vector <MST_vertex*> vertexPtrList;
	for (unsigned int i = 0; i < vertexIdList.size(); i++) {		
		MST_vertex* v_ptr = (*this->vertexMap)[vertexIdList[i]];
		vertexPtrList.push_back(v_ptr);
		vector <int> compressedSequence;
		compressedSequences.push_back(compressedSequence);
		names.push_back(v_ptr->name);
	}
	int numberOfSites = vertexPtrList[0]->sequence.size();
	vector<int> sitePattern;
	for(int site = 0; site < numberOfSites; site++){
		sitePattern.clear();
		for (MST_vertex* v_ptr: vertexPtrList) {
			sitePattern.push_back(v_ptr->sequence[site]);}
		if (find(distinctPatterns.begin(),distinctPatterns.end(),sitePattern)!=distinctPatterns.end()){
			distinctPatternsToSitesWherePatternRepeats[sitePattern].push_back(site);
			
		} else {
			distinctPatterns.push_back(sitePattern);	
			vector<int> sitePatternRepeats;
			sitePatternRepeats.push_back(site);
			distinctPatternsToSitesWherePatternRepeats[sitePattern] = sitePatternRepeats;						
			for (unsigned int i = 0; i < sitePattern.size(); i++){
				compressedSequences[i].push_back(sitePattern[i]);
			}
		}
	}
	for (vector<int> sitePattern : distinctPatterns){
		int sitePatternWeight = distinctPatternsToSitesWherePatternRepeats[sitePattern].size();
		sitePatternWeights_ptr.push_back(sitePatternWeight);		
		sitePatternRepeats_ptr.push_back(distinctPatternsToSitesWherePatternRepeats[sitePattern]);
	}
	return make_tuple(names,compressedSequences,sitePatternWeights_ptr,sitePatternRepeats_ptr);
}

void MST::SetCompressedSequencesAndSiteWeightsForInputSequences() {				
	vector <vector<int>> distinctPatterns;
	map <vector<int>,vector<int>> distinctPatternsToSitesWherePatternRepeats;	
	int numberOfSites = (*this->vertexMap)[0]->sequence.size();
	int numberOfInputSequences = this->vertexMap->size();
	int sitePatternWeight; int v_id; int site;
	vector<int> sitePattern;
	for(site=0; site < numberOfSites; site++) {
		sitePattern.clear();
		for (v_id = 0; v_id < numberOfInputSequences; v_id ++) {
			sitePattern.push_back((*this->vertexMap)[v_id]->sequence[site]);
			}
		if (find(distinctPatterns.begin(),distinctPatterns.end(),sitePattern)!=distinctPatterns.end()) {
			distinctPatternsToSitesWherePatternRepeats[sitePattern].push_back(site);			
		} else {
			distinctPatterns.push_back(sitePattern);	
			vector<int> sitePatternRepeats;
			sitePatternRepeats.push_back(site);
			distinctPatternsToSitesWherePatternRepeats[sitePattern] = sitePatternRepeats;						
			for (v_id = 0; v_id < numberOfInputSequences; v_id ++) {				
				(*this->vertexMap)[v_id]->globallyCompressedSequence.push_back(sitePattern[v_id]);
			}
		}
	}
	for (vector<int> sitePattern: distinctPatterns) {
		sitePatternWeight = distinctPatternsToSitesWherePatternRepeats[sitePattern].size();		
		this->siteWeights.push_back(sitePatternWeight);		
	}
}

void MST::WriteToFile(string FileName) {
	ofstream mstFile;	
	mstFile.open(FileName);	
	MST_vertex * v;
	for (pair <int, MST_vertex *> vIdAndPtr: *this->vertexMap) {
		v = vIdAndPtr.second;
		for (MST_vertex * n: v->neighbors) {
			if (v->id < n->id) {
				mstFile << v->name << "\t" << n->name << "\t" << this->GetEdgeWeight(v->id, n->id) << endl; 
			}
		}
	}
	mstFile.close();
}

int MST::ConvertDNAToChar(char dna) {
	string dna_upper = string(1,toupper(dna));
	int dna_char = 4;
	if (this->mapDNAtoInteger.find(dna_upper) != this->mapDNAtoInteger.end()) {
		dna_char = this->mapDNAtoInteger[dna_upper];
	} else {		
		cout << "DNA character " << dna_upper << " not recognized" << endl;
		throw mt_error("DNA character not recognized");
	}	
	return (dna_char);
}

void MST::ReadPhyx(string sequenceFileNameToSet) {
    this->sequenceFileName = sequenceFileNameToSet;
    ifstream inputFile(this->sequenceFileName.c_str());    
    vector <int> recodedSequence;
    string seqName;
    string seq = "";
    string line; getline(inputFile, line);
    vector <string> splitLine = emtr::split_ws(line);
    this->numberOfInputSequences = stoi(splitLine[0]);
    this->sequence_length = stoi(splitLine[1]);
    // cout << "sequence file contains " << this->numberOfInputSequences << " sequences of length " << this->sequence_length << endl;
    while(getline(inputFile,line)) {
        vector <string> splitLine = emtr::split_ws(line);
        seqName = splitLine[0];
        seq = splitLine[1];
        for (char const dna: seq) {
            recodedSequence.push_back(this->ConvertDNAToChar(dna));
        }
        this->AddVertex(seqName,recodedSequence);
        recodedSequence.clear();
    }    
	inputFile.close();
    this->numberOfInputSequences = this->vertexMap->size();
}

void MST::ReadFasta(string sequenceFileNameToSet) {
	this->sequenceFileName = sequenceFileNameToSet;
	vector <int> recodedSequence;
	recodedSequence.clear();
	unsigned int site = 0;
    unsigned int seq_len = 0;
	int dna_char;
	int num_amb = 0;
	int num_non_amb = 0;
	ifstream inputFile(this->sequenceFileName.c_str());
	string seqName;
	string seq = "";
	for(string line; getline(inputFile, line );) {
		if (line[0]=='>') {
			if (seq != "") {
				for (char const dna: seq) {
					if (!isspace(dna)) {
						dna_char = this->ConvertDNAToChar(dna);
						if (dna_char > 3) {
							num_amb += 1;
						} else {
							num_non_amb += 1;
						}
						recodedSequence.push_back(dna_char);					
						site += 1;							
						}						
				}
				if (this->IsSequenceDuplicated(recodedSequence)) {
					this->AddDuplicatedSequenceName(seqName,recodedSequence);
				} else {										
					this->AddVertex(seqName,recodedSequence);
				}
				recodedSequence.clear();
			} 
			seqName = line.substr(1,line.length());
			seq = "";
			site = 0;			
		}
		else {
			seq += line ;
		}		
	}		
	for (char const dna: seq) {
		if (!isspace(dna)) {
			dna_char = this->ConvertDNAToChar(dna);
			if (dna_char > 3) { // FIX_AMB
				num_amb += 1;
			} else {
				num_non_amb += 1;
			}
			recodedSequence.push_back(dna_char);
			site += 1;
		}
	}
	if (this->IsSequenceDuplicated(recodedSequence)) {
		this->AddDuplicatedSequenceName(seqName,recodedSequence);
	} else {
		this->AddVertex(seqName,recodedSequence);
	}
    seq_len = recodedSequence.size();
	recodedSequence.clear();
	inputFile.close();
    cout << "Number of sequences in fasta file is " << this->vertexMap->size() << endl;
    cout << "Sequence length is " << seq_len << endl;	
}

vector<int> MST::DecompressSequence(vector<int>* compressedSequence, vector<vector<int>>* sitePatternRepeats){
	int totalSequenceLength = 0;
	for (vector<int> sitePatternRepeat: *sitePatternRepeats){
		totalSequenceLength += int(sitePatternRepeat.size());
	}
	vector <int> decompressedSequence;
	for (int v_ind = 0; v_ind < totalSequenceLength; v_ind++){
		decompressedSequence.push_back(char(0));
	}
	int dnaToAdd;
	for (int sitePatternIndex = 0; sitePatternIndex < int(compressedSequence->size()); sitePatternIndex++){
		dnaToAdd = (*compressedSequence)[sitePatternIndex];
		for (int pos: (*sitePatternRepeats)[sitePatternIndex]){
			decompressedSequence[pos] = dnaToAdd;
		}
	}
	return (decompressedSequence);	
}

string MST::EncodeAsDNA(vector<int> sequence){
	string allDNA = "AGTC";
	string dnaSequence = "";
	for (int s : sequence){
		dnaSequence += allDNA[s];
	}
	return dnaSequence;
}


void MST::ComputeMST() {

	int numberOfVertices = (this->v_ind);		
	const int numberOfEdges = numberOfVertices*(numberOfVertices-1)/2;		
	
	double * weights;
	weights = new double [numberOfEdges];
		
	int edgeIndex = 0;
	for (int i=0; i<numberOfVertices; i++) {
		for (int j=i+1; j<numberOfVertices; j++) {			
			weights[edgeIndex] = this->ComputeHammingDistance((*this->vertexMap)[i]->sequence,(*this->vertexMap)[j]->sequence);
			if (weights[edgeIndex] == 0) {
				cout << this->EncodeAsDNA((*this->vertexMap)[i]->sequence) << endl;
			}
			edgeIndex += 1;
		}
	}
	typedef pair <int,int > E;

	E * edges;
	edges = new E [numberOfEdges];
	edgeIndex = 0;
	for (int i=0; i<numberOfVertices; i++) {
		for (int j=i+1; j<numberOfVertices; j++) {
			edges[edgeIndex] = E(i,j);
			edgeIndex += 1;
		}
	}	
	
	vector<int> p(numberOfVertices); 

	emtr::prim_graph p_graph(numberOfVertices, edges, weights, numberOfEdges);
	
	emtr::prim(p_graph, &p[0]);
	delete[] edges;		
	int edgeCount = 0;
	for (size_t u = 0; u != p.size(); u++) {
		if (p[u] != u) {
			edgeCount += 1;
			if (u < p[u]) {
				edgeIndex = GetEdgeIndex(u,p[u],numberOfVertices);
			} else {
				edgeIndex = GetEdgeIndex(p[u],u,numberOfVertices);
			}
			this->AddEdge(u, p[u], weights[edgeIndex]);
		}
	}
	this->UpdateMaxDegree();
	delete[] weights;
}

///...///...///...///...///...///...///...///...///... structural expectation maximization ///...///...///...///...///...///...///...///...///

class SEM_vertex {	
public:
	int degree = 0;
	int timesVisited = 0;
	bool observed = 0;
	string completeSequence;
    vector <int> DNArecoded;
	vector <int> DNAcompressed;	
	vector <int> AArecoded;
	vector <int> AAcompressed;	
	vector <string> dupl_seq_names;
	int id = -40;
	int global_id = -40;
	string newickLabel = "";
	string name = "";
	double logScalingFactors = 0;
	double vertexLogLikelihood = 0;
	double sumOfEdgeLogLikelihoods = 0;
	int rateCategory = 0;
	int GCContent = 0;
	bool AA_cond_like_initialized = 0;
	bool DNA_cond_like_initialized = 0;
	vector <SEM_vertex *> neighbors;
	vector <SEM_vertex *> children;
	array <double, 4> root_prob_hss;
	SEM_vertex * parent = this;
	void AddNeighbor(SEM_vertex * v_ptr);
	void RemoveNeighbor(SEM_vertex * v_ptr);
	void AddParent(SEM_vertex * v_ptr);
	void RemoveParent();
	void AddChild(SEM_vertex * v_ptr);
	void RemoveChild(SEM_vertex * v_ptr);
	void SetVertexLogLikelihood(double vertexLogLikelihoodToSet);
	int inDegree = 0;
	int outDegree = 0;
	emtr::Md transitionMatrix;
	Eigen::Matrix <double,20,20> AAtransitionMatrix;
	emtr::Md transitionMatrix_stored;	
	array <double, 4> rootProbability;
	array <double, 4> posteriorProbability;
	
	SEM_vertex (int idToAdd, vector <int> compressedSequenceToAdd) {
		this->id = idToAdd;
		this->DNArecoded = compressedSequenceToAdd;
		this->transitionMatrix = emtr::Md{};
		this->transitionMatrix_stored = emtr::Md{};
		for (int dna = 0; dna < 4; dna ++) {
			this->transitionMatrix[dna][dna] = 1.0;
			this->transitionMatrix_stored[dna][dna] = 1.0;
		}
		for (int i = 0; i < 4; i ++) {
			this->rootProbability[i] = 0;
			this->posteriorProbability[i] = 0;
			this->root_prob_hss[i] = 0;
		}
	}	
	~SEM_vertex () {
		this->neighbors.clear();
	}
};

void SEM_vertex::SetVertexLogLikelihood(double vertexLogLikelihoodToSet) {
	this->vertexLogLikelihood = vertexLogLikelihoodToSet;
}

void SEM_vertex::AddParent(SEM_vertex * v) {
	this->parent = v;
	this->inDegree += 1;
}

void SEM_vertex::RemoveParent() {
	this->parent = this;
	this->inDegree -=1;
}

void SEM_vertex::AddChild(SEM_vertex * v) {
	this->children.push_back(v);
	this->outDegree += 1;
}

void SEM_vertex::RemoveChild(SEM_vertex * v) {
	int ind = find(this->children.begin(),this->children.end(),v) - this->children.begin();
	this->children.erase(this->children.begin()+ind);
	this->outDegree -=1;
}

void SEM_vertex::AddNeighbor(SEM_vertex * v) {
	this->degree += 1;
	this->neighbors.push_back(v);
}

void SEM_vertex::RemoveNeighbor(SEM_vertex * v) {
	this->degree -= 1;
	int ind = find(this->neighbors.begin(),this->neighbors.end(),v) - this->neighbors.begin();
	this->neighbors.erase(this->neighbors.begin()+ind);
}

///...///...///...///...///...///...///... clique ...///...///...///...///...///...///...///

class clique {
	public:	
	map <clique *, double> logScalingFactorForMessages;
	double logScalingFactorForClique;
	map <clique *, std::array <double, 4>> messagesFromNeighbors;
    vector <int> compressedSequence;
	string name;
	int id;
	int site;
	int inDegree = 0;
	int outDegree = 0;
	int timesVisited = 0;
	clique * parent = this;
	SEM_vertex * root_variable;
	vector <clique *> children;
	void AddParent(clique * C);
	void AddChild(clique * C);
	void ComputeBelief();
	void ScaleFactor();
	SEM_vertex * x;
	SEM_vertex * y;
	std::array <double, 4> MarginalizeOverVariable(SEM_vertex * v);
//	emtr::Md DivideBeliefByMessageMarginalizedOverVariable(SEM_vertex * v);	
	// Clique is defined over the vertex pair (X,Y)
	// No of variables is always 2 for bifurcating tree-structured DAGs
	double scalingFactor;
	emtr::Md factor;
	emtr::Md initialPotential;	
	emtr::Md belief;
	// P(X,Y)
	
	void SetInitialPotentialAndBelief(int site);
	
	// If the clique contains an observed variable then initializing
	// the potential is the same as restricting the corresponding
	// CPD to row corresponding to observed variable
	void AddNeighbor(clique * C);
	clique (SEM_vertex * x, SEM_vertex * y) {
		this->x = x;
		this->y = y;
		this->name = to_string(x->id) + "-" + to_string(y->id);
		this->logScalingFactorForClique = 0;
	}
	
	~clique () {
		
	}
};



std::array <double, 4> clique::MarginalizeOverVariable(SEM_vertex * v) {
	std::array <double, 4> message;	
	if (this->x == v) {
		for (int dna_y = 0; dna_y < 4; dna_y ++) {
			message[dna_y] = 0;
			for (int dna_x = 0; dna_x < 4; dna_x ++) {
				message[dna_y] += this->belief[dna_x][dna_y];
			}
		}
	} else if (this->y == v) {
		for (int dna_x = 0; dna_x < 4; dna_x ++) {
			message[dna_x] = 0;
			for (int dna_y = 0; dna_y < 4; dna_y ++) {
				message[dna_x] += this->belief[dna_x][dna_y];
			}
		}
	} else {
		cout << "Check marginalization over variable" << endl;
	}
	return (message);
}

void clique::ComputeBelief() {
	this->factor = this->initialPotential;
	vector <clique *> neighbors = this->children;
	std::array <double, 4> messageFromNeighbor;
	bool debug = 0;		
	// if (this->y->name == "Mouse" && this->root_variable->name == "h_4" && this->site == 257) {debug = 1; cout << "..........................................................................." << endl;}
	if (this->parent != this) {
		neighbors.push_back(this->parent);
	}
	// cout << "6a" << endl;
	for (clique * C_neighbor : neighbors) {		
		this->logScalingFactorForClique += this->logScalingFactorForMessages[C_neighbor];
		messageFromNeighbor = this->messagesFromNeighbors[C_neighbor];
		for (int i = 0; i < 4; i ++) if (isnan(messageFromNeighbor[i])) throw mt_error("message contain nan");
		double message_sum = 0;		
				
		//
		if (this->y == C_neighbor->x or this->y == C_neighbor->y) {
		    //	factor_row_i = factor_row_i (dot) message
			for (int dna_x = 0 ; dna_x < 4; dna_x ++) for (int dna_y = 0 ; dna_y < 4; dna_y ++) this->factor[dna_x][dna_y] *= messageFromNeighbor[dna_y];
			// scale factor
			this->scalingFactor = 0;
			for (int dna_x = 0 ; dna_x < 4; dna_x ++) for (int dna_y = 0 ; dna_y < 4; dna_y ++) this->scalingFactor += this->factor[dna_x][dna_y];			
			assert(scalingFactor > 0);
			for (int dna_x = 0 ; dna_x < 4; dna_x ++) for (int dna_y = 0 ; dna_y < 4; dna_y ++) this->factor[dna_x][dna_y] /= this->scalingFactor;
			this->logScalingFactorForClique += log(scalingFactor);
			if (debug) cout << "Performing row-wise multiplication" << endl;
		} else if (this->x == C_neighbor->x or this->x == C_neighbor->y) {
			for (int dna_y = 0 ; dna_y < 4; dna_y ++) for (int dna_x = 0 ; dna_x < 4; dna_x ++) this->factor[dna_x][dna_y] *= messageFromNeighbor[dna_x];
			// scale factor
			this->scalingFactor = 0;
			for (int dna_x = 0 ; dna_x < 4; dna_x ++) for (int dna_y = 0 ; dna_y < 4; dna_y ++) this->scalingFactor += this->factor[dna_x][dna_y];			
			assert(scalingFactor > 0);
			for (int dna_x = 0 ; dna_x < 4; dna_x ++) for (int dna_y = 0 ; dna_y < 4; dna_y ++) this->factor[dna_x][dna_y] /= this->scalingFactor;
			this->logScalingFactorForClique += log(scalingFactor);
			if (debug) cout << "Performing column-wise multiplication" << endl;
		} else {
			cout << "Check product step" << endl;
            throw mt_error("check product step");
		}
	}
	// cout << "6b" << endl;	
	this->scalingFactor = 0;
	for (int dna_x = 0 ; dna_x < 4; dna_x ++) for (int dna_y = 0 ; dna_y < 4; dna_y ++) this->scalingFactor += this->factor[dna_x][dna_y];
	assert(scalingFactor > 0);
	for (int dna_x = 0 ; dna_x < 4; dna_x ++) for (int dna_y = 0 ; dna_y < 4; dna_y ++) this->factor[dna_x][dna_y] /= this->scalingFactor;
	this->logScalingFactorForClique += log(this->scalingFactor);
	// cout << "6c" << endl;		
	// cout << "6d" << endl;
	for (int dna_x = 0 ; dna_x < 4; dna_x ++) for (int dna_y = 0 ; dna_y < 4; dna_y ++) this->belief[dna_x][dna_y] = this->factor[dna_x][dna_y];
	
	for (int dna_x = 0 ; dna_x < 4; dna_x ++) {
		for (int dna_y = 0 ; dna_y < 4; dna_y ++) {			
			if (isnan(this->belief[dna_x][dna_y])){
				throw mt_error("belief has nan");
			}
		}
	}
	// cout << "6e" << endl;
}

void clique::AddParent(clique * C) {
	this->parent = C;
	this->inDegree += 1; 
}

void clique::AddChild(clique * C) {
	this->children.push_back(C);
	this->outDegree += 1;
}

void clique::SetInitialPotentialAndBelief(int site) {
	this->site = site;
	// Initialize psi
	// V = (X,Y) X->Y (wlog), X is always an unobserved vertex
	int matchingCase;
	// Case 1. Y is an observed vertex
	// Product factor psi = P(Y|X) restricted to observed value xe of X
	// psi (y|x) is set to 0 if x != xe
	// of y->DNAcompressed[site] is -1 (gap) then transition matrix is not conditioned
	// if (y->observed && y->DNAcompressed[site] == -1) {
	// 	cout << "dealing with gap" << endl;
	// }

	if (y->observed) {
		this->initialPotential = y->transitionMatrix;
		int dna_y = y->DNAcompressed[site];
		if (dna_y == -1) {
			// cout << "case 1" << endl;
			matchingCase = 1; // gap
		} else {
			// cout << "case 2" << endl;
			matchingCase = 2;			
			for (int dna_p = 0; dna_p < 4; dna_p++) {
				for (int dna_c = 0; dna_c < 4; dna_c++) {
					if (dna_c != dna_y) {
						this->initialPotential[dna_p][dna_c] *= 0;
					} else {
						this->initialPotential[dna_p][dna_c] *= 1;
					}
				}
			}	
		}
	}
	
	// Case 2. X and Y are hidden and X is not the root
	// psi = P(Y|X)
	if (!y->observed) {
		// cout << "case 3" << endl;
		matchingCase = 3;
		this->initialPotential = y->transitionMatrix;		
	}	
	
	// Case 3. X and Y are hidden and X is the root and "this" is not the root clique
	// psi = P(Y|X)
	if (!y->observed and (x->parent == x) and (this->parent != this)) {
		// cout << "case 4" << endl;
		matchingCase = 4;
		this->initialPotential = y->transitionMatrix;
	}
	
	// Case 4. X and Y are hidden and X is the root and "this" is the root clique
	// psi = P(X) * P(Y|X) 
	if (!y->observed and (x->parent == x) and (this->parent == this)) {	
		// cout << "case 5" << endl;
		matchingCase = 5;
		this->initialPotential = y->transitionMatrix;
		for (int dna_p = 0; dna_p < 4; dna_p++) {
			for (int dna_c = 0; dna_c < 4; dna_c++) {
				this->initialPotential[dna_p][dna_c] *= x->rootProbability[dna_c];
			}
		}
	}

	bool all_zero_initial_potential = 1;
	bool all_zero_transition = 1;
	for(int i = 0; i < 4; i ++) {
		for (int j = 0; j < 4; j ++) {
			if(this->initialPotential[i][j] != 0) all_zero_initial_potential = 0;
			if(y->transitionMatrix[i][j] != 0) all_zero_transition = 0;
		}
	}

	if (all_zero_transition) {
		cout  << "all zero transition" << endl;
		exit(-1);
	}
	
	if (all_zero_initial_potential) {
		cout << y->name << endl;		
		int dna_y = y->DNAcompressed[site];		
		cout << dna_y << endl;
		cout  << "all_zero_initial_potential" << endl;
		for (int i = 0; i < 4; i ++) {
			for (int j = 0; j < 4; j ++) {
				cout << i << "\t" << j << "\t" << y->transitionMatrix[i][j] << endl;
			}
		}
		exit(-1);
	}

	// double maxValue = 0;
	for (int i = 0; i < 4; i ++) {
		for (int j = 0; j < 4; j ++) {			
			if(isnan(this->initialPotential[i][j])){
				cout << "initial potential has nan" << endl;
				throw mt_error("initial potential has nan");
			}
		}
	}
	
	this->belief = this->initialPotential;
	this->logScalingFactorForClique = 0;
	this->logScalingFactorForMessages.clear();
	this->messagesFromNeighbors.clear();
}

///...///...///...///...///...///...///...///... clique tree ...///...///...///...///...///...///...///...///

class cliqueTree {
public:
	vector < pair <clique *, clique *> > edgesForPreOrderTreeTraversal;
	vector < pair <clique *, clique *> > edgesForPostOrderTreeTraversal;
	vector < pair <clique *, clique *> > cliquePairsSortedWrtLengthOfShortestPath;
	map < pair <SEM_vertex *, SEM_vertex *>, emtr::Md> marginalizedProbabilitiesForVariablePair;
	map < pair <clique *, clique *>, pair <SEM_vertex *, SEM_vertex *>> cliquePairToVariablePair;
	int site;
	clique * root;
	bool rootSet;
	vector <clique *> leaves;
	vector <clique *> cliques;
	void CalibrateTree();
	void SetMarginalProbabilitesForEachEdgeAsBelief();
	void ComputeMarginalProbabilitesForEachVariablePair();
	void ComputePosteriorProbabilitesForVariable();
	void ConstructSortedListOfAllCliquePairs();
	clique * GetLCA (clique * C_1, clique * C_2);
	int GetDistance(clique * C_1, clique * C_2);
	int GetDistanceToAncestor(clique * C_d, clique * C_a);
	void SetLeaves();
	void SetRoot();
	void AddEdge(clique * C_1, clique * C_2);
	void SendMessage(clique * C_1, clique * C_2);
	void AddClique(clique * C);
	void SetSite(int site);
	void InitializePotentialAndBeliefs();
	void SetEdgesForTreeTraversalOperations();
	void WriteCliqueTreeAndPathLengthForCliquePairs(string fileName);
	emtr::Md GetP_XZ(SEM_vertex * X, SEM_vertex * Y, SEM_vertex * Z);	
	SEM_vertex * GetCommonVariable(clique * Ci, clique * Cj);
	tuple <SEM_vertex *,SEM_vertex *,SEM_vertex *> GetXYZ(clique * Ci, clique * Cj);
	cliqueTree () {
		rootSet = 0;
	}
	~cliqueTree () {
		for (clique * C: this->cliques) {
			delete C;
		}
		this->cliques.clear();
		this->leaves.clear();
	}
};


tuple <SEM_vertex *,SEM_vertex *,SEM_vertex *> cliqueTree::GetXYZ(clique * Ci, clique * Cj) {
	SEM_vertex * X; SEM_vertex * Y; SEM_vertex * Z;
	SEM_vertex * Y_temp;
	clique * Cl;
	pair <clique *, clique *> cliquePairToCheck;
	if (Ci->parent == Cj or Cj->parent == Ci) {
		// Case 1: Ci and Cj are neighbors
		Y = this->GetCommonVariable(Ci, Cj);
		if (Ci->y == Y){
		X = Ci->x;		
		} else {
			X = Ci->y;			
		}
		if (Cj->y == Y){			
			Z = Cj->x;
		} else {
			Z = Cj->y;
		}
		
	} else {
		// Case 2: Ci and Cj are not neighbors
		// Ci-...-Cl-Cj
		vector <clique *> neighbors;
		if (Cj->parent != Cj) {
			neighbors.push_back(Cj->parent);
		}
		for (clique * C: Cj->children) {
			neighbors.push_back(C);
		}
		
		Cl = Ci;
		
		for (clique * C: neighbors) {
			if (C->name < Ci->name) {
				cliquePairToCheck = pair <clique*, clique*>(C,Ci);
			} else {
				cliquePairToCheck = pair <clique*, clique*>(Ci,C);
			}
			if (this->cliquePairToVariablePair.find(cliquePairToCheck) != this->cliquePairToVariablePair.end()) {
				if (Ci == cliquePairToCheck.first) {
					Cl = cliquePairToCheck.second;
				} else {
					Cl = cliquePairToCheck.first;
				}
				break;
			}
		}		
		
		// Scope(Ci,Cl) = {X,Y}
		if (Ci->name < Cl->name) {
			tie(X,Y) = this->cliquePairToVariablePair[pair <clique*, clique*>(Ci,Cl)];
		} else {
			tie(Y,X) = this->cliquePairToVariablePair[pair <clique*, clique*>(Cl,Ci)];
		}			
				
				
		if (Cj->x == Y or Cj->y == Y){
			// Case 2a
			// Scope(Cj) = {Y,Z}
			if (Cj->x == Y) {
				Z = Cj->y;
			} else {
				Z = Cj->x;
			}
		} else {
			// Case 2b
			// Scope(Cl,Cj) = {Y,Z}
			if (Cl->name < Cj->name) {
				tie(Y_temp,Z) = this->cliquePairToVariablePair[pair <clique*, clique*>(Cl,Cj)];
			} else {
				tie(Z,Y_temp) = this->cliquePairToVariablePair[pair <clique*, clique*>(Cj,Cl)];
			}
			if (Y_temp != Y){
                throw mt_error("check Case 2b");
            }
		}
	}	
	
	return (tuple <SEM_vertex *,SEM_vertex *,SEM_vertex *>(X,Y,Z));
}

emtr::Md cliqueTree::GetP_XZ(SEM_vertex * X, SEM_vertex * Y, SEM_vertex * Z) {
	emtr::Md P_XY; emtr::Md P_YZ;
	emtr::Md P_ZGivenY; emtr::Md P_XZ;
	
    if (X->id < Y->id) {
		P_XY = this->marginalizedProbabilitiesForVariablePair[pair<SEM_vertex *, SEM_vertex *>(X,Y)];
	} else {		
		P_XY = emtr::MT(this->marginalizedProbabilitiesForVariablePair[pair<SEM_vertex *, SEM_vertex *>(Y,X)]);
	}
	if (Y->id < Z->id) {
		P_YZ = this->marginalizedProbabilitiesForVariablePair[pair<SEM_vertex *, SEM_vertex *>(Y,Z)];
	} else {		
		P_YZ = emtr::MT(this->marginalizedProbabilitiesForVariablePair[pair<SEM_vertex *, SEM_vertex *>(Z,Y)]);
	}
//	cout << "P_XY is " << endl << P_XY << endl;
//	cout << "P_YZ is " << endl << P_YZ << endl;
	P_ZGivenY = emtr::Md{};
	double rowSum;
	for (int row = 0; row < 4; row ++) {		
		rowSum = 0;
		for (int col = 0; col < 4; col ++) {
			rowSum += P_YZ[row][col];
		}
		for (int col = 0; col < 4; col ++) {
			if (rowSum != 0){
				P_ZGivenY[row][col] = P_YZ[row][col]/rowSum;
			}			
		}
	}
	
//	cout << "P_ZGivenY is " << endl << P_ZGivenY << endl;
	
	for (int row = 0; row < 4; row ++) {		
		for (int col = 0; col < 4; col ++) {
			P_XZ[row][col] = 0;
		}
	}
	
	for (int dna_y = 0; dna_y < 4; dna_y ++) {		
		for (int dna_x = 0; dna_x < 4; dna_x ++) {
			for (int dna_z = 0; dna_z < 4; dna_z ++) {					
				// Sum over Y
				P_XZ[dna_x][dna_z] += P_XY[dna_x][dna_y] * P_ZGivenY[dna_x][dna_z];
			}
		}
	}
	
	return (P_XZ);
}

SEM_vertex * cliqueTree::GetCommonVariable(clique * Ci, clique * Cj) {
	SEM_vertex * commonVariable;
	if (Ci->x == Cj->x or Ci->x == Cj->y) {
		commonVariable = Ci->x;
	} else {
		commonVariable = Ci->y;
	}
	return (commonVariable);
}


void cliqueTree::ConstructSortedListOfAllCliquePairs() {
	this->cliquePairsSortedWrtLengthOfShortestPath.clear();
	vector < tuple <int, clique*, clique*>> sortedPathLengthAndCliquePair;
	int pathLength;
	for (clique * Ci : this->cliques) {
		for (clique * Cj : this->cliques) {
			if (Ci->name < Cj->name) {
				if (Ci->outDegree > 0 or Cj->outDegree > 0) {
					pathLength = this->GetDistance(Ci, Cj);
					sortedPathLengthAndCliquePair.push_back(make_tuple(pathLength,Ci,Cj));
				}				
			}
		}
	}
	sort(sortedPathLengthAndCliquePair.begin(),sortedPathLengthAndCliquePair.end());
	clique * Ci; clique * Cj;
	for (tuple <int, clique*, clique*> pathLengthCliquePair : sortedPathLengthAndCliquePair) {
		Ci = get<1>(pathLengthCliquePair);
		Cj = get<2>(pathLengthCliquePair);
		this->cliquePairsSortedWrtLengthOfShortestPath.push_back(pair <clique *, clique *> (Ci,Cj));
	}
}

int cliqueTree::GetDistance(clique * C_1, clique * C_2) {
	clique * lca = this->GetLCA(C_1, C_2);
	int d;
	d = this->GetDistanceToAncestor(C_1,lca) + this->GetDistanceToAncestor(C_2,lca);
	return (d);
}

int cliqueTree::GetDistanceToAncestor(clique * C_d, clique* C_a) {
	int d = 0;
	clique * C_p;
	C_p = C_d;
	while (C_p != C_a) {
		C_p = C_p->parent;
		d += 1;
	}
	return (d);
}

clique * cliqueTree::GetLCA(clique * C_1, clique * C_2) {
	vector <clique *> pathToRootForC1;
	vector <clique *> pathToRootForC2;
	clique * C1_p;
	clique * C2_p;
	C1_p = C_1;
	C2_p = C_2;
	
	clique * C_r = this->edgesForPreOrderTreeTraversal[0].first;
	
	while (C1_p->parent != C1_p) {
		pathToRootForC1.push_back(C1_p);
		C1_p = C1_p->parent;
	}
	pathToRootForC1.push_back(C1_p);
	if (C1_p != C_r) {
		cout << "Check get LCA for C1" << endl;
	}
	
	while (C2_p->parent != C2_p) {
		pathToRootForC2.push_back(C2_p);
		C2_p = C2_p->parent;
	}
	pathToRootForC2.push_back(C2_p);
	if (C2_p != C_r) {
		cout << "Check get LCA for C2" << endl;
	}
	
	clique * lca;
	lca = C_1;
	
	for (clique * C : pathToRootForC1) {
		if (find(pathToRootForC2.begin(),pathToRootForC2.end(),C)!=pathToRootForC2.end()) {
			lca = C;
			break;
		}
	}		
	return (lca);	
}

void cliqueTree::SetMarginalProbabilitesForEachEdgeAsBelief() {
	this->marginalizedProbabilitiesForVariablePair.clear();
	//	Store P(X,Y) for each clique
	for (clique * C: this->cliques) {
		if (C->x->id < C->y->id) {			
			this->marginalizedProbabilitiesForVariablePair.insert(pair<pair<SEM_vertex *, SEM_vertex *>, emtr::Md>(pair<SEM_vertex *, SEM_vertex *>(C->x,C->y),C->belief));
		} else {			
			this->marginalizedProbabilitiesForVariablePair.insert(pair<pair<SEM_vertex *, SEM_vertex *>, emtr::Md>(pair<SEM_vertex *, SEM_vertex *>(C->y,C->x),emtr::MT(C->belief)));
		}
	}
}

void cliqueTree::ComputeMarginalProbabilitesForEachVariablePair() {	
	this->marginalizedProbabilitiesForVariablePair.clear();
	this->cliquePairToVariablePair.clear();
	// For each clique pair store variable pair 
	// Iterate over clique pairs in order of increasing distance in clique tree	
	
	clique * Ci; clique * Cj;
	
	SEM_vertex * X; SEM_vertex * Z;
	SEM_vertex * Y;

	emtr::Md P_XZ;
		
	//	Store P(X,Y) for each clique
	for (clique * C: this->cliques) {
		if (C->x->id < C->y->id) {			
			this->marginalizedProbabilitiesForVariablePair.insert(pair<pair<SEM_vertex *, SEM_vertex *>, emtr::Md>(pair<SEM_vertex *, SEM_vertex *>(C->x,C->y),C->belief));
		} else {			
			this->marginalizedProbabilitiesForVariablePair.insert(pair<pair<SEM_vertex *, SEM_vertex *>, emtr::Md>(pair<SEM_vertex *, SEM_vertex *>(C->y,C->x),emtr::MT(C->belief)));
			
		}
	}
	
	for (pair <clique *, clique *> cliquePair : this->cliquePairsSortedWrtLengthOfShortestPath) {
		tie (Ci, Cj) = cliquePair;
		tie (X, Y, Z) = this->GetXYZ(Ci, Cj);		
		this->cliquePairToVariablePair.insert(pair <pair <clique *, clique *>,pair <SEM_vertex *, SEM_vertex *>>(pair <clique *, clique *>(Ci,Cj), pair <SEM_vertex *, SEM_vertex *>(X,Z)));
		P_XZ = this->GetP_XZ(X, Y, Z);
		if (X->id < Z->id) {
			this->marginalizedProbabilitiesForVariablePair.insert(pair<pair<SEM_vertex *, SEM_vertex *>, emtr::Md>(pair<SEM_vertex *, SEM_vertex *>(X,Z),P_XZ));			
		} else {			
			this->marginalizedProbabilitiesForVariablePair.insert(pair<pair<SEM_vertex *, SEM_vertex *>, emtr::Md>(pair<SEM_vertex *, SEM_vertex *>(Z,X),emtr::MT(P_XZ)));
		}
	}
}


void cliqueTree::SetRoot() {
	for (clique * C: this->cliques) {
		if (C->inDegree == 0) {
			this->root = C;
		}
	}
}

void cliqueTree::InitializePotentialAndBeliefs() {
	for (clique * C: this->cliques) {		
		C->SetInitialPotentialAndBelief(this->site);
	}
}

void cliqueTree::SetSite(int site) {
	this->site = site;
}

void cliqueTree::AddClique(clique * C) {
	this->cliques.push_back(C);
}

void cliqueTree::AddEdge(clique * C_1, clique * C_2) {
	C_1->AddChild(C_2);
	C_2->AddParent(C_1);
}

void cliqueTree::SetEdgesForTreeTraversalOperations() {
	for (clique * C : this->cliques) {
		C->timesVisited = 0;
	}
	this->edgesForPostOrderTreeTraversal.clear();
	this->edgesForPreOrderTreeTraversal.clear();
	vector <clique *> verticesToVisit;
	verticesToVisit = this->leaves;
	clique * C_child; clique * C_parent;
	int numberOfVerticesToVisit = verticesToVisit.size();
	
	while (numberOfVerticesToVisit > 0) {
		C_child = verticesToVisit[numberOfVerticesToVisit - 1];
		verticesToVisit.pop_back();
		numberOfVerticesToVisit -= 1;
		C_parent = C_child->parent;
		if (C_child != C_parent) {
			C_parent->timesVisited += 1;
			this->edgesForPostOrderTreeTraversal.push_back(make_pair(C_parent, C_child));
			if (C_parent->timesVisited == C_parent->outDegree) {				
				verticesToVisit.push_back(C_parent);
				numberOfVerticesToVisit += 1;				
			}
		}
	}
	
	for (int edgeInd = this->edgesForPostOrderTreeTraversal.size() -1; edgeInd > -1; edgeInd --) {
		this->edgesForPreOrderTreeTraversal.push_back(this->edgesForPostOrderTreeTraversal[edgeInd]);
	}
}

void cliqueTree::SetLeaves() {
	this->leaves.clear();
	for (clique * C: this->cliques) {
		if (C->outDegree == 0) {
			this->leaves.push_back(C);
		}		
	}
}


void cliqueTree::SendMessage(clique * C_from, clique * C_to) {	
	double logScalingFactor;
	double largestElement;	
	array <double, 4> messageFromNeighbor;
	array <double, 4> messageToNeighbor;
	bool verbose = 0;
	if (verbose) {
		cout << "Preparing message to send from " << C_from->x->name << "," << C_from->y->name << " to " ;
		cout << C_to->x->name << "," << C_to->y->name << " is " << endl;
	}
	
	// Perform the three following actions
	
	// A) Compute product: Multiply the initial potential of C_from
	// with messages from all neighbors of C_from except C_to, and
	
	// B) Compute sum: Marginalize over the variable that
	// is in C_from but not in C_to
	
	// C) Transmit: sending the message to C_to
	
	// Select neighbors
	vector <clique *> neighbors;
	if (C_from->parent != C_from and C_from->parent != C_to) {
		neighbors.push_back(C_from->parent);
	}
	
	for (clique * C_child : C_from->children) {
		if (C_child != C_to) {
			neighbors.push_back(C_child);
		}
	}
	
	emtr::Md factor;
	factor = C_from->initialPotential;
	
	logScalingFactor = 0;
		// A. PRODUCT: Multiply messages from neighbors that are not C_to
	// cout << "A step" << endl;
	bool allZero = 1;
	for (int dna_x = 0 ; dna_x < 4; dna_x ++) {
		for (int dna_y = 0 ; dna_y < 4; dna_y ++) {
			if (factor[dna_x][dna_y] != 0) allZero = 0;
			// cout << "initial factor for " << dna_x << "," << dna_y << " " << factor[dna_x][dna_y] << endl;
			
		}
	}
	if (allZero) {
		cout << "initial potential is all zero";
		throw mt_error("all zero in factor");
	}

	for (clique * C_neighbor : neighbors) {
		messageFromNeighbor = C_from->messagesFromNeighbors[C_neighbor];
		// for (int i = 0; i < 4; i ++) isnan(messageFromNeighbor)
		if (C_from->y == C_neighbor->x or C_from->y == C_neighbor->y) {
		// factor_row_i = factor_row_i (dot) message
			for (int dna_x = 0 ; dna_x < 4; dna_x ++) {
				for (int dna_y = 0 ; dna_y < 4; dna_y ++) {
					factor[dna_x][dna_y] *= messageFromNeighbor[dna_y];	
				}
			}
			if (verbose) {cout << "Performing row-wise multiplication" << endl;}
		} else if (C_neighbor->x == C_from->x or C_neighbor->y == C_from->x) {
		// factor_col_i = factor_col_i (dot) message
			for (int dna_y = 0 ; dna_y < 4; dna_y ++) {
				for (int dna_x = 0 ; dna_x < 4; dna_x ++) {
					factor[dna_x][dna_y] *= messageFromNeighbor[dna_x];
				}
			}
			if (verbose) {cout << "Performing column-wise multiplication" << endl;}
		} else {
			cout << "Check product step" << endl;
			throw mt_error("check product step");
		}		
		// Check to see if each entry in the factor is zero
		allZero = 1;	
		for (int dna_x = 0 ; dna_x < 4; dna_x ++) {
			for (int dna_y = 0 ; dna_y < 4; dna_y ++) {
				if (factor[dna_x][dna_y] != 0) {
					allZero = 0;
				}
			}
		}
		if (allZero) {
			cout << "all zero in factor" << endl;
			throw mt_error("all zero in factor");}		
		// Rescale factor
		largestElement = 0;
		for (int dna_x = 0 ; dna_x < 4; dna_x ++) {
			for (int dna_y = 0 ; dna_y < 4; dna_y ++) {
				if (largestElement < factor[dna_x][dna_y]) {
					largestElement = factor[dna_x][dna_y];
				}
			}
		}
		for (int dna_x = 0 ; dna_x < 4; dna_x ++) {
			for (int dna_y = 0 ; dna_y < 4; dna_y ++) {
				factor[dna_x][dna_y] /= largestElement;
			}
		}
		logScalingFactor += log(largestElement);
		logScalingFactor += C_from->logScalingFactorForMessages[C_neighbor];
	}

	// cout << "B step" << endl;
	allZero = 1;
	for (int dna_x = 0 ; dna_x < 4; dna_x ++) {
		for (int dna_y = 0 ; dna_y < 4; dna_y ++) {
			if (factor[dna_x][dna_y] != 0) allZero = 0;
		}
	}

	if (allZero) {
			cout << "all zero in factor" << endl;
			throw mt_error("all zero in factor");};
	// B. SUM
		// Marginalize factor by summing over common variable
	largestElement = 0;
	if (C_from->y == C_to->x or C_from->y == C_to->y) {
		// Sum over C_from->x		
		for (int dna_y = 0 ; dna_y < 4; dna_y ++) {
			messageToNeighbor[dna_y] = 0;
			for (int dna_x = 0 ; dna_x < 4; dna_x ++) {
				messageToNeighbor[dna_y] += factor[dna_x][dna_y];
			}
		}
		if (verbose) {
			cout << "Performing column-wise summation" << endl;
		}							
	} else if (C_from->x == C_to->x or C_from->x == C_to->y) {
		// Sum over C_from->y		
		for (int dna_x = 0 ; dna_x < 4; dna_x ++) {
			messageToNeighbor[dna_x] = 0;
			for (int dna_y = 0 ; dna_y < 4; dna_y ++) {				
				messageToNeighbor[dna_x] += factor[dna_x][dna_y];
			}
		}
		if (verbose) {
			cout << "Performing row-wise summation" << endl;
		}							
	} else {		
		cout << "Check sum step" << endl;
		throw mt_error("Check sum step");
	}
	// Rescale message to neighbor
	largestElement = 0;
	for (int dna_x = 0 ; dna_x < 4; dna_x ++) {
		if (largestElement < messageToNeighbor[dna_x]) {
			largestElement = messageToNeighbor[dna_x];
		}
	}
	if (largestElement == 0) {
		cout << "Division by zero" << endl;
		throw mt_error("Division by zero");
	}
	for (int dna_x = 0 ; dna_x < 4; dna_x ++) {
		messageToNeighbor[dna_x] /= largestElement;
	}
	logScalingFactor += log(largestElement);
	// cout << "check point 1" << endl;
	if (verbose) {
		cout << "Sending the following message from " << C_from->name << " to " << C_to->name << endl;
		for (int dna = 0; dna < 4; dna ++) {
			cout << "dna " << dna << endl;
			cout << to_string(messageToNeighbor[dna]) << endl;			
			cout << isnan(messageToNeighbor[dna]) << endl;
			if(isnan(messageToNeighbor[dna])) {
				cout << "messageToNeighbor[" << dna << "] is not a number";
				throw mt_error("Check message to neighbor");
			} else {
				cout << messageToNeighbor[dna] << endl;
			}
		}
	}
	// cout << "check point 2" << endl;
	// cout << endl;
	
	// cout << "C step" << endl;
	// C. TRANSMIT
	if(isnan(logScalingFactor)) {
		cout << "logScalingFactor is nan" << endl;
		throw mt_error("logScalingFactor is nan");
	}
	C_to->logScalingFactorForMessages.insert(make_pair(C_from,logScalingFactor));
	C_to->messagesFromNeighbors.insert(make_pair(C_from,messageToNeighbor));
}

void cliqueTree::CalibrateTree() {
	clique * C_p; clique * C_c;
	
	//	Send messages from leaves to root
	// cout << "13a" << endl;	
	for (pair <clique *, clique *> cliquePair : this->edgesForPostOrderTreeTraversal) {
		tie (C_p, C_c) = cliquePair;
		this->SendMessage(C_c, C_p);
	}
	
	// cout << "13b" << endl;

	//	Send messages from root to leaves	
	for (pair <clique *, clique *> cliquePair : this->edgesForPreOrderTreeTraversal) {
		tie (C_p, C_c) = cliquePair;
		this->SendMessage(C_p, C_c);
	}	
	//	Compute beliefs
	// cout << "13c" << endl;
	for (clique * C: this->cliques) {
		C->ComputeBelief();		
	}

	// cout << "13d" << endl;
}

///...///...///...///...///...///...///...///... Params Struct ///...///...///...///...///...///...///...///...///

struct EM_struct {
	string method;
	// init method  - store before starting rep
	int rep;
	// rep          - store at start of each rep
	int iter;   
	map <int, double> ecd_ll_per_iter;
	// ecd for each iteration of EM
	double ll_init;
	// ll initial store at completion of EM run
	double ll_final;
	// ll final store at completion of EM run
	string root;
	// init root_prob - store for each rep when params are initialized
	array <double,4> root_prob_init;
	// final root_prob - store if log likelihood score is maximum
	array <double,4> root_prob_final;
	// init trans_prob - store transition probability for each child node for each rep when params are initialized
	map <string, Md> trans_prob_init;	
	// init trans_prob - store transition probability for each child node for each rep when params are initialized
	map <string, Md> trans_prob_final;
};

struct EMTrifle_struct {
	vector <int> iter_trifle;
	vector <map <int,double>> ecd_ll_per_iter_for_trifle;
	vector <double> ll_initial_trifle;
	vector <double> ll_final_trifle;
	vector <map <string, Md>> trans_prob_final_trifle;
	vector <array <double, 4>> root_prob_final_trifle;
	vector <map <string, Md>> trans_prob_initial_trifle;
	vector <array <double, 4>> root_prob_initial_trifle;	
	string method;
	
	int rep;
	string root;
	EMTrifle_struct()
	: ecd_ll_per_iter_for_trifle(3),
	  iter_trifle(3),
	  ll_initial_trifle(3, numeric_limits<double>::quiet_NaN()),
	  ll_final_trifle(3, numeric_limits<double>::quiet_NaN()),
	  trans_prob_final_trifle(3),
	  root_prob_final_trifle(3),
	  trans_prob_initial_trifle(3),
	  root_prob_initial_trifle(3)
	{}
};

struct AA_opt_struct {
	int rep;
	int iter;
	map <int, double> aa_ll_per_iter;
	double aa_ll_init;
	double aa_ll_final;
	string root;
	vector <double> root_prob_init;
	vector <double> root_prob_final;
	Eigen::Matrix <double,20,20> exchangeability_matrix_init;
	Eigen::Matrix <double,20,20> exchangeability_matrix_final;
};

///...///...///...///...///...///...///...///...///... SEM ...///...///...///...///...///...///...///...///...///

class SEM {
private:	

public:
	void SetParameters(string parameters_json);
	Eigen::Matrix <double,20,20> Q_D;
	Eigen::Matrix <double,20,20> Vc, Vc_inv;	
    Eigen::Matrix <complex<double>,20,1> lambda;
	Eigen::Matrix <double,20,20> S_D;
	array <int,3> num_patterns_for_layer;
	vector <int> DNASitesByWeightDesc;
	vector <int> DNACumulativePatternWeight;
	vector <size_t> DNAPatternOrderByWeightDesc;
	vector <vector<int>> DNAPatternOriginalSites;
	vector <int> OriginalToCompressed;
	int compressedToFirstOriginalSite(size_t k) const;
	const vector<int>& compressedToAllOriginalSites(size_t k) const;
	int originalToCompressed(size_t origSite) const;
	bool spectral_ready_ = false;
	vector<double> pi_D;
	Graph T;
	void AssignLSBranchLengthsToTree();
	string dayhoff_rate_matrix_file_name;
	array <double, 4> alpha_pi;
	array <double, 4> alpha_M_row;
	array <double, 4> sample_dirichlet(const array <double, 4>& alpha, mt19937_64& gen);
	int largestIdOfVertexInMST = 1;
	default_random_engine generator;
	bool setParameters;
	bool verbose = 0;
	string modelForRooting;
	map <string,int> mapDNAtoInteger;
	map <int, SEM_vertex*> * vertexMap;
	vector <SEM_vertex*> vertices;
	map <pair<SEM_vertex *,SEM_vertex *>,emtr::Md> * M_hss;
	vector <int> DNAPatternWeights;
	vector <double> CumPercDNAPatternWeights;
	vector <int> AAPatternWeights;
	vector <pair<int,double>> AAsitecumweight; 
	vector <int> gaplesscompressedDNAsites;
	vector <bool> gapLessDNAFlag;
	vector <bool> gapLessAAFlag;
	int rep;	
	void addToZeroEntriesOfTransitionMatrix(emtr::Md& P);
	void addToZeroEntriesOfRootProbability(array<double,4>& p);
	int num_aa_patterns;
	int num_repetitions;	
	double thrCoarse = 50.0;
	double thrMedium = 30.0;
	double thrFine = 20.0;	
	double cum_pattern_weight_coarse;
	double cum_pattern_weight_medium;
	double cum_pattern_weight_fine = 100.0;
	vector<pair<int,int>> AAsite_and_weight;
	vector <double> cum_fractional_weight_AA;
	vector <vector <int> > sitePatternRepetitions;
	vector <int> sortedDeltaGCThresholds;
	int numberOfInputSequences;
	int numberOfVerticesInSubtree;
	int numberOfObservedVertices;
	int numberOfExternalVertices = 0;	
	int num_dna_patterns;
	int maxIter;
	double patW_coarse;
	double patW_medium;
	double patW_fine;
	double frac_iter_coarse;
	double frac_iter_medium;
	double frac_iter_fine;
	double ecdllConvergenceThreshold = 0.1;	
	vector <double> ecdllConvergenceThresholdForTrifle;
	double sumOfExpectedLogLikelihoods = 0;
	double maxSumOfExpectedLogLikelihoods = 0;
	int node_ind = 1;
	chrono::system_clock::time_point t_start_time;
	chrono::system_clock::time_point t_end_time;
	// ofstream * logFile;
	SEM_vertex * root;
	vector < pair <SEM_vertex *, SEM_vertex *>> edgesForPostOrderTreeTraversal;
	vector < pair <SEM_vertex *, SEM_vertex *>> edgesForPreOrderTreeTraversal;	
	vector < pair <SEM_vertex *, SEM_vertex *>> edgesForChowLiuTree;
	vector < pair <SEM_vertex *, SEM_vertex *>> directedEdgeList;
	map <pair <SEM_vertex *, SEM_vertex *>, double> edgeLengths;
	vector < SEM_vertex *> leaves;
	vector < SEM_vertex *> preOrderVerticesWithoutLeaves;
	vector < SEM_vertex *> postOrderVertices;
	map < pair <SEM_vertex * , SEM_vertex *>, emtr::Md > expectedCountsForVertexPair;
	map < pair <SEM_vertex * , SEM_vertex *>, emtr::Md > posteriorProbabilityForVertexPair;
	map < SEM_vertex *, array <double,4>> expectedCountsForVertex; 
	map < SEM_vertex *, array <double,4>> posteriorProbabilityForVertex;
	map <int, emtr::Md> rateMatrixPerRateCategory;
	map <int, emtr::Md> rateMatrixPerRateCategory_stored;
	map <int, double> scalingFactorPerRateCategory;
	map <int, double> scalingFactorPerRateCategory_stored;
	int numberOfRateCategories = 0;
	double maximumLogLikelihood;
	double max_log_likelihood_diri;
	double max_log_likelihood_pars;
	double max_log_likelihood_rep;
	double max_log_likelihood_hss;
	double max_log_likelihood_best;
	emtr::Md I4by4;	
	cliqueTree * cliqueT;
	bool debug;
	bool finalIterationOfSEM;
	bool flag_logDet = 0;
	bool flag_Hamming = 0;
	bool flag_JC = 0;
	bool flag_added_duplicated_sequences = 0;
	map <string, int> nameToIdMap;
	string DNAsequenceFileName;
	string AAsequenceFileName;
	string phylip_file_name;
	string topologyFileName;
	string probabilityFileName;
	string probabilityFileName_best;
	string probabilityFileName_pars;
	string probabilityFileName_hss;
	string probabilityFileName_diri;
	string prefix_for_output_files;
	string ancestralSequencesString = "";
	void SetFJTree(const fj::Graph& tree);
	double sequenceLength;
	// Add vertices (and compressed sequence for leaves)
	array <double, 4> rootProbability;
	array <double, 4> rootProbability_stored;
	SEM_vertex * root_stored;
	vector <int> compressedSequenceToAddToMST;
	string nameOfSequenceToAddToMST;
	double AAlogLikelihood;
	double logLikelihood;
	double logLikelihood_exp_counts;
	double logLikelihood_current;
	// Used for updating MST
	vector <int> indsOfVerticesOfInterest;
	vector <int> indsOfVerticesToKeepInMST;
	vector <int> idsOfVerticesOfInterest;
	vector <int> idsOfObservedVertices;	
	vector <int> idsOfVerticesToRemove;
	vector <int> idsOfVerticesToKeepInMST;	
	vector <int> idsOfExternalVertices;	
	vector < tuple <int, string, vector <int> > > idAndNameAndSeqTuple;
	vector<tuple <string,string,int,int,double,double,double,double>> EMTR_results;
	vector<tuple <string,string,int,int,int,double,double,double,double>> EMTrifle_results;
	// Used for updating global phylogenetic tree
	vector < pair <int, int> > edgesOfInterest_ind;	
	vector < pair < vector <int>, vector <int> > > edgesOfInterest_seq;
	string weightedEdgeListString;
	map < string, vector <int>> sequencesToAddToGlobalPhylogeneticTree;
	vector < tuple <string, string, double>> weightedEdgesToAddToGlobalPhylogeneticTree;
	vector < tuple <string, string, double>> edgeLogLikelihoodsToAddToGlobalPhylogeneticTree;		
	map <string, double> vertexLogLikelihoodsMapToAddToGlobalPhylogeneticTree;
	map <pair<SEM_vertex *,SEM_vertex *>,double> edgeLogLikelihoodsMap;
	SEM_vertex * externalVertex;
	void AddArc(SEM_vertex * from, SEM_vertex * to);
	void RemoveArc(SEM_vertex * from, SEM_vertex * to);
	void SetStream(ofstream& stream_to_set);
	void ClearDirectedEdges();
	void ClearUndirectedEdges();
	void ClearAllEdges();
	void AddVertex(string name, vector <int> compressedSequenceToAdd);
	void SetAASequenceForVertex(string name,vector <int> AARecodedSequenceToAdd);
	void AddFJVertices();
	void AddFJWeightedEdges();
	void AddWeightedEdges(vector<tuple<string,string,double>> weightedEdgesToAdd);
	void AddEdgeLogLikelihoods(vector<tuple<string,string,double>> edgeLogLikelihoodsToAdd);
	void AddExpectedCountMatrices(map < pair <SEM_vertex * , SEM_vertex *>, emtr::Md > expectedCountsForVertexPair);
	void AddVertexLogLikelihoods(map <string,double> vertexLogLikelihoodsMapToAdd);
	void SetNumberOfVerticesInSubtree(int numberOfVertices);	
	void AddSitePatternWeights(vector <int> sitePatternWeightsToAdd);
	void AddSitePatternRepeats(vector <vector <int> > sitePatternRepetitionsToAdd);
	void AddSequences(vector <vector <int>> sequencesToAdd);	
	void OpenAncestralSequencesFile();
	void AddRootVertex();
	void SetAASiteCumWeight();
	void BuildCumFractionIndicesAA();
	void SetVertexVector();
	void OptimizeRateMatrixNelderMead(int max_iter = 4000, double tol = 1e-5, double step = 0.1, int i0 = 0, int j0 = 1, int p_ref = 19);	
	void BuildSIndex(vector<pair<int,int>>& idxS, int i0, int j0) const;
	void PackTheta(const vector<pair<int,int>>& idxS, int p_ref,vector<double>& theta) const;
	void UnpackTheta(const vector<double>& theta, const vector<pair<int,int>>& idxS,int i0, int j0, int p_ref);
	double NMObjectiveNegLL(const vector<double>& theta, const vector<pair<int,int>>& idxS, int i0, int j0, int p_ref);
	void UNRESTForEachBranch();
//	void AddCompressedSequencesAndNames(map<string,vector<int>> sequencesList, vector <vector <int>> sitePatternRepeats);	
	void AddAllSequences(string sequencesFileName);
	void AddNames(vector <string> namesToAdd);
	void AddGlobalIds(vector <int> idsToAdd);		
	double ComputeDistance(int v_i, int v_j);
	void RootedTreeAlongAnEdgeIncidentToCentralVertex();
	void RootTreeAlongAnEdgePickedAtRandom();
	void RootTreeAtAVertexPickedAtRandom();
	void SetParsimonySites();	
	void EMTrifle_started_with_parsimony_rooted_at(SEM_vertex *v, int layer);
	void EMTrifle_started_with_dirichlet_rooted_at(SEM_vertex *v, int layer);
	void EMTrifle_started_with_HSS_rooted_at(SEM_vertex *v, int layer);
	tuple<int,double,double,double,double> EM_started_with_parsimony_rooted_at(SEM_vertex *v);
	tuple<int,double,double,double,double> EM_started_with_dirichlet_rooted_at(SEM_vertex *v);	
	tuple<int,double,double,double,double> EM_started_with_HSS_parameters_rooted_at(SEM_vertex *v);
	void StoreParamsInEMCurrent(string init_or_final);
	// void StoreParamsInEMTrifleCurrent(int layer);
	void StoreInitialParamsInEMTrifleCurrent(int layer);
	void StoreFinalParamsInEMTrifleCurrent(int layer);
	void ComputeSumOfExpectedLogLikelihoods();
	void RootTreeAlongEdge(SEM_vertex * u, SEM_vertex * v);
	void SelectEdgeIncidentToVertexViaMLForBH(SEM_vertex * v);
	void InitializeTransitionMatricesAndRootProbability();
	void ComputeMPEstimateOfAncestralSequences();
	void ComputeMPEstimateOfAncestralSequencesForTrifle(int layer);
	void ComputeMAPEstimateOfAncestralSequences();
	void ComputeMAPEstimateOfAncestralSequencesUsingCliques();
	void SetEdgesForPreOrderTraversal();
	void SetEdgesForPostOrderTraversal();
	void SetEdgesForTreeTraversalOperations();
	void SetLeaves();
	void SetVerticesForPreOrderTraversalWithoutLeaves();
	void SetVerticesForPostOrderTraversal();
	void SetObservedUnobservedStatus();
	void OptimizeParametersUsingMAPEstimates();
	void ComputeMLEOfRootProbability();
	void ComputeMLEOfTransitionMatrices();	
	void ComputeMarginalProbabilitiesUsingExpectedCounts();
	void ComputePosteriorProbabilitiesUsingMAPEstimates();
	void SetInfoForVerticesToAddToMST();
	void SetIdsOfExternalVertices();
	void ResetAncestralSequences();
	void WriteParametersOfBH(string BHparametersFileName);
	void RemoveEdgeLength(SEM_vertex * u, SEM_vertex * v);
	void AddEdgeLength(SEM_vertex * u, SEM_vertex * v, double t);
	double GetEdgeLength(SEM_vertex * u, SEM_vertex * v);
	double ComputeEdgeLength(SEM_vertex * u, SEM_vertex * v);
	void SetEdgeLength(SEM_vertex * u, SEM_vertex * v, double t);
	void SetEdgesFromTupleVector(vector<tuple <string, string, double>> edge_vector);
	void SetEdgesFromTopologyFile();
	string EncodeAsDNA(vector<int> sequence);
	vector<int> DecompressSequence(vector<int> compressedSequence, vector<vector<int>> sitePatternRepeats);	
	void ComputeChowLiuTree();
	void AddSubforestOfInterest(SEM * localPhylogeneticTree);
	void ReadRootedTree(string treeFileName);
	void SetBHparameters();
	void ReparameterizeBH();
	void Set_pi_for_neighbors_of_root();
	array<double,4> get_pi_child();
	void ReadProbabilities();
	void WriteProbabilities(string fileName);
	void ReadTransitionProbabilities(string fileName);
	int GetVertexId(string v_name);	
	SEM_vertex * GetVertex(string v_name);
	bool ContainsVertex(string v_name);
	int GetEdgeIndex (int vertexIndex1, int vertexIndex2, int numberOfVertices);
	emtr::Md GetP_yGivenx(emtr::Md P_xy);
	emtr::Md ComputeTransitionMatrixUsingAncestralStates(SEM_vertex * p, SEM_vertex * c);
	emtr::Md ComputeTransitionMatrixUsingAncestralStatesForTrifle(SEM_vertex * p, SEM_vertex * c, int layer);
	array <double, 4> GetBaseComposition(SEM_vertex * v);
	array <double, 4> GetBaseCompositionForTrifle(SEM_vertex * v, int layer);
	array <double, 4> GetObservedCountsForVariable(SEM_vertex * v);
	string modelSelectionCriterion;
	string distance_measure_for_NJ = "log-det";
	void SetModelSelectionCriterion(string modelSelectionCriterionToSet);
	void RootTreeAtVertex(SEM_vertex * r);
	void StoreEdgeListForChowLiuTree();
	void RestoreEdgeListForChowLiuTree();
	void StoreDirectedEdgeList();
	void RestoreDirectedEdgeList();
	void StoreBestProbability();
	void StoreRootAndRootProbability();
	void RestoreRootAndRootProbability();
	void StoreTransitionMatrices();	
	void RestoreTransitionMatrices();	
	void RestoreBestProbability();
	void StoreRateMatricesAndScalingFactors();
	void RestoreRateMatricesAndScalingFactors();
	void ResetPointerToRoot();
	void ResetTimesVisited();
	void SetIdsForObservedVertices(vector <int> idsOfObservedVerticesToAdd);
	void SetNumberOfInputSequences(int numOfInputSeqsToSet);
	void ComputeMLRootedTreeForFullStructureSearch();
	void SetNeighborsBasedOnParentChildRelationships();
	void ComputeMLRootedTreeForRootSearchUnderBH();	
	void ComputeMLEstimateOfBHGivenExpectedDataCompletion();		
	void SetMinLengthOfEdges();
	void SetParametersForRateMatrixForNelderMead(double x[], int rateCat);
	void NelderMeadForOptimizingParametersForRateCat(int rateCat, int n, double start[], double xmin[], 
		 double *ynewlo, double reqmin, double step[], int konvge,
		 int kcount, int *icount, int *numres, int *ifault);	
	void ComputeInitialEstimateOfModelParameters();
	void ComputeInitialEstimateOfModelParametersForTrifle(int layer);
	void BumpZeroEntriesOfModelParameters();
	void SetInitialEstimateOfModelParametersUsingDirichlet();
	void SetInitialEstimateOfModelParametersUsingHSS();
	void TransformRootedTreeToBifurcatingTree();
	void SwapRoot();
	void SuppressRoot();
	bool IsTreeInCanonicalForm();
	bool root_search;
	string init_criterion;
	string parameter_file;
	Matrix20d ExpQt_fromSpectrum(double t) const;
	void PrecomputeSpectrum();	
	void ComputeAALogLikelihood();
	void ComputeAALogLikelihood_maps();
	void ComputeLogLikelihood();
	void ComputeTrifleLogLikelihood(int layer);
	int siteIdx;
	std::tuple<double,double,double>
    LL_Grad_Hess_Counts_Reversible_Tuple(const Eigen::Matrix<double,20,20>& N, const Eigen::VectorXd& pis, const Eigen::VectorXd& invpis, const Eigen::MatrixXd& U, const Eigen::VectorXd& lambdas, double t) const;
	tuple<Eigen::Matrix<double,20,20>, Eigen::Matrix<double,20,20>, Eigen::Matrix<double,20,20>> EvalP_QP_Q2P_Tuple(const Eigen::VectorXd& pis, const Eigen::VectorXd& invpis, const Eigen::MatrixXd& U, const Eigen::VectorXd& lambdas, double t) const;
	tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd,Eigen::MatrixXd, Eigen::VectorXd> ReversibleSpecTuple(const Eigen::Matrix<double,20,20>& Q, const Eigen::VectorXd& pi) const;
	double GetMLDistance_NR_Reversible_Tuples(SEM_vertex* u, SEM_vertex* v, const Eigen::Matrix<double,20,20>& Q, const Eigen::VectorXd& pi, int max_iter = 50, double grad_tol = 1e-8, double ll_tol = 1e-10);
	void ComputeLogLikelihoodUsingExpectedDataCompletion();
	pair <bool, SEM_vertex *> CheckAndRetrieveSingletonHiddenVertex();
	pair <bool, SEM_vertex *> CheckAndRetrieveHiddenVertexWithOutDegreeZeroAndInDegreeOne();
	pair <bool, SEM_vertex *> CheckAndRetrieveHiddenVertexWithOutDegreeOneAndInDegreeOne();
	pair <bool, SEM_vertex *> CheckAndRetrieveHiddenVertexWithOutDegreeOneAndInDegreeZero();
	pair <bool, SEM_vertex *> CheckAndRetrieveHiddenVertexWithOutDegreeGreaterThanTwo();
	pair <bool, SEM_vertex *> CheckAndRetrieveObservedVertexThatIsTheRoot();
	pair <bool, SEM_vertex *> CheckAndRetrieveObservedVertexThatIsNotALeafAndIsNotTheRoot();
	double GetExpectedMutualInformation(SEM_vertex * u, SEM_vertex * v);
	void ResetLogScalingFactors();
	// Mutual information I(X;Y) is computed using 
	// P(X,Y), P(X), and P(Y), which in turn are computed using
	// MAP estimates
	void InitializeExpectedCounts();
	void InitializeExpectedCountsForEachVariable();
	void InitializeExpectedCountsForEachVariablePair();
	void InitializeExpectedCountsForEachEdge();
	void ResetExpectedCounts();
	void ConstructCliqueTree();	
	void ComputeExpectedCounts();
	void ComputeExpectedCountsForTrifle(int layer);
	emtr::Md GetObservedCounts(SEM_vertex * u, SEM_vertex * v);
	void AddToExpectedCounts();
	void AddToExpectedCountsForEachVariable();
	void AddToExpectedCountsForEachVariablePair();
	emtr::Md GetExpectedCountsForVariablePair(SEM_vertex * u, SEM_vertex * v);
	emtr::Md GetPosteriorProbabilityForVariablePair(SEM_vertex * u, SEM_vertex * v);
	void AddToExpectedCountsForEachEdge();
	DistanceMap DM_NR;
	DistanceMap DM_LS;
	Spectrum spectrum;
	// Mutual information I(X;Y) is computing using 
	// P(X,Y), P(X), and P(Y), which in turn are computed using
	// A calibrated clique tree
	// using P(X,Y) = Sum_{H\{X,Y}}{P(X,Y|H\{X,Y},O)}
	// where H is the set of hidden variables and O is the set of observed variables	
	void RootTreeUsingEstimatedParametersViaML();
	void SetFlagForFinalIterationOfSEM();
	int ConvertDNAtoIndex(char dna);
	int ConvertAAtoIndex(char aa);
	char GetDNAfromIndex(int dna_index);			
	double BIC;
	double AIC;
	void ComputeBIC();
	void ComputeAIC();
	void StoreEdgeListAndSeqToAdd();
	void SelectIndsOfVerticesOfInterestAndEdgesOfInterest();
	void RenameHiddenVerticesInEdgesOfInterestAndSetIdsOfVerticesOfInterest();
	void SetAncestralSequencesString();
	void SetWeightedEdgesToAddToGlobalPhylogeneticTree();	
	void ComputeVertexLogLikelihood(SEM_vertex * v);
	void ComputeEdgeLogLikelihood(SEM_vertex * u, SEM_vertex * v);
	void ComputeLSDistanceMatrix();
	void ComputeNRDistanceMatrix();
	void SetDayhoffExchangeabilityAndRootProbabilityParameters();	
	void SetRateMatrixParameters();
	void SetPerturbedDayhoffRateMatrix();
	void SetDirichletPerturbedDayhoffRateMatrix();
	void SetDNASequencesFromFile(string sequenceFileName);
	void SetAASequencesFromFile(string sequenceFileName);
	void CompressDNASequences();
	void CompressDNASequencesBkp();
	void CompressAASequences();
	void SetEdgeAndVertexLogLikelihoods();
	bool IsNumberOfNonSingletonComponentsGreaterThanZero();
	// chat ConvertDNAtoChar(char dna);
	// void WriteTree();
	void WriteAncestralSequences();
	void SetPrefixForOutputFiles(string prefix_for_output_files_to_set);
	void WriteRootedTreeInNewickFormat(string newickFileName);
	void WriteUnrootedTreeInNewickFormat(string newickFileName);
	void WriteCliqueTreeToFile(string cliqueTreeFileName);
	void WriteRootedTreeAsEdgeList(string fileName);
	void WriteUnrootedTreeAsEdgeList(string fileName);
	void ResetData();	
	void AddDuplicatedSequencesToRootedTree(MST * M);
	void AddDuplicatedSequencesToUnrootedTree(MST * M);
	void ReadFasta();
	Eigen::Matrix <double,20,20> GetCountMatrix(SEM_vertex * u, SEM_vertex * v);
	vector <double> GetCountVector (SEM_vertex * u);
	double DyadicLogLikelihood(SEM_vertex * u, SEM_vertex * v, double t);
	double DyadicLogLikelihoodFirstDer(SEM_vertex * u, SEM_vertex * v, double t);
	// double DyadicLogLikelihoodSecondDer(SEM_vertex * u, SEM_vertex * v, double t);
	double GetNewtonRaphsonDistance(SEM_vertex * u, SEM_vertex * v);
	double GetJukesCantorDistance(SEM_vertex * u, SEM_vertex * v);
	void ComputeMLDistancesForAAsequences();
	void GetOLSBranchLengths();
	void OptimizeQ_D();
	void SetSitePatternsAndWeights();
	void SetParameterFile();
	void initialize_BH(string init_criterion);
	void EM_DNA_rooted_at_each_internal_vertex_started_with_dirichlet_store_results(int num_repetitions);
	void EM_DNA_rooted_at_each_internal_vertex_started_with_parsimony_store_results(int num_repetitions);
	void EM_DNA_rooted_at_each_internal_vertex_started_with_HSS_store_results(int num_repetitions);
	void EMTrifle_DNA_for_replicate(int replicate);
	void EMTrifle_DNA_rooted_at_each_internal_vertex_started_with_parsimony_store_results();
	void EMtrifle_DNA_rooted_at_each_internal_vertex_started_with_dirichlet_store_results();
	void EMtrifle_DNA_rooted_at_each_internal_vertex_started_with_HSS_store_results();
	void Optimize_AA_RateMatrix_rooted_at_each_internal_vertex(int num_repetitions);
	void EM_rooted_at_each_internal_vertex_started_with_dirichlet(int num_repetitions);
	void EM_rooted_at_each_internal_vertex_started_with_parsimony(int num_repetitions);
	void EM_rooted_at_each_internal_vertex_started_with_HSS_par(int num_repetitions);
	void EM_root_search_at_each_internal_vertex_started_with_dirichlet(int num_repetitions);
	void EM_root_search_at_each_internal_vertex_started_with_HSS_par(int num_repetitions);
	void set_alpha_PI(double a1, double a2, double a3, double a4);
	void set_alpha_M_row(double a1, double a2, double a3, double a4);
	map <pair<SEM_vertex*, SEM_vertex*>,double> ML_distances;
	array <double, 4> sample_pi();
    array <double, 4> sample_M_row();
	vector <EM_struct> EM_DNA_runs_pars;
	vector <EMTrifle_struct> EMTrifle_DNA_runs_pars;
	vector <EM_struct> EM_DNA_runs_diri;
	vector <EMTrifle_struct> EMTrifle_DNA_runs_diri;
	vector <EM_struct> EM_DNA_runs_hss;
	vector <EMTrifle_struct> EMTrifle_DNA_runs_hss;
	vector <AA_opt_struct> AA_opt_runs;
	EM_struct EM_current{};
	EMTrifle_struct EMTrifle_current{};
	AA_opt_struct AA_struct{};	
	string emtrifle_to_json(const EMTrifle_struct& em) const;	
	int first_index_gt(const std::vector<double>& cum, double thr);
	// Select vertex for rooting Chow-Liu tree and update edges in T
	// Modify T such that T is a bifurcating tree and likelihood of updated
	// tree is equivalent to the likelihood of T
	SEM (double loglikelihood_conv_thresh, int max_EM_iter, bool verbose_flag_to_set) {
		// this->SetDayhoffRateMatrix();		
		this->alpha_pi = {100,100,100,100}; // default value
		this->alpha_M_row = {100,2,2,2};
		this->root_search = false;
		this->ecdllConvergenceThreshold = loglikelihood_conv_thresh;
		this->maxIter = max_EM_iter;
		this->distance_measure_for_NJ = "Hamming";		
		this->flag_Hamming = 1; this->flag_logDet = 0;
		this->verbose = verbose_flag_to_set;		
		this->node_ind = 0;
		this->vertexMap = new map <int, SEM_vertex *> ;
		this->M_hss = new map <pair<SEM_vertex*,SEM_vertex*>,emtr::Md>;
		// this->vertexName2IdMap = new map <string, int> ;
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		this->generator = default_random_engine(seed);
		this->I4by4 = emtr::Md{};
		for (int i = 0; i < 4; i++) {
			this->I4by4[i][i] = 1.0;
		}
		this->cliqueT = new cliqueTree;		
		mapDNAtoInteger["A"] = 0;
		mapDNAtoInteger["C"] = 1;
		mapDNAtoInteger["G"] = 2;
		mapDNAtoInteger["T"] = 3;
		this->finalIterationOfSEM = 0;
	}
	
	~SEM () {
		for (pair <int, SEM_vertex * > idPtrPair : * this->vertexMap) {
			delete idPtrPair.second;
		}
		this->vertexMap->clear();
		delete this->vertexMap;
		// delete this->vertexName2IdMap;
		delete this->cliqueT;
		delete this->M_hss;
	}
};

int SEM::compressedToFirstOriginalSite(size_t k) const {
    if (k >= DNAPatternOriginalSites.size()) return -1;
    const auto& v = DNAPatternOriginalSites[k];
    return v.empty() ? -1 : v.front();
}

const vector<int>& SEM::compressedToAllOriginalSites(size_t k) const {
    static const vector<int> EMPTY;
    if (k >= DNAPatternOriginalSites.size()) return EMPTY;
    return DNAPatternOriginalSites[k];
}

int SEM::originalToCompressed(size_t origSite) const {
    if (origSite >= OriginalToCompressed.size()) return -1;
    return OriginalToCompressed[origSite];
}


void SEM::addToZeroEntriesOfRootProbability(array <double,4> & baseComposition) {
	double sum = 0.0;

    for (size_t i = 0; i < 4; ++i) {
        if (baseComposition[i] == 0.0) {
            baseComposition[i] = 0.1;
        }
        sum += baseComposition[i];
    }

	for (size_t i = 0; i < 4; ++i) {
		baseComposition[i] /= sum;
	}
    
}

void SEM::addToZeroEntriesOfTransitionMatrix(array<array<double,4>,4>& P) {
	
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			if (P[i][j] == 0.0) {
				P[i][j] = (i == j) ? 0.97 : 0.01;
			}
		}
	}
	
	for (int i = 0; i < 4; ++i) {
		double rowSum = 0.0;
		for (int j = 0; j < 4; ++j) rowSum += P[i][j];

		if (rowSum == 0.0 || !isfinite(rowSum)) {		
			for (int j = 0; j < 4; ++j) P[i][j] = 0.01;
			P[i][i] = 0.97;
		} else {
			for (int j = 0; j < 4; ++j) P[i][j] /= rowSum;
		}
	}
}

int SEM::first_index_gt(const std::vector<double>& cum, double thr) {    
    constexpr double EPS = 1e-12;
    auto it = std::upper_bound(cum.begin(), cum.end(), thr + EPS);
    return static_cast<int>(it - cum.begin());
}


void SEM::SetParameters(string parameters_json) {
	json ph = json::parse(parameters_json);
	map<string, string> settings_map;
	flatten_json(ph, settings_map);
	for (const auto& kv : settings_map) {
		const string& key = kv.first;
		const string& value = kv.second;

		// cout << key << " = " << value << "\n";

		if (key.find("dna.reps") != string::npos) {
			try {this->num_repetitions = stoi(value);
				cout << "num_repetitions " << this->num_repetitions << endl; 
			} catch (const invalid_argument& e) {
				throw mt_error("couldn't parse dna.reps");
			}			 
    	} else if (key.find("dna.D_M[0]") != string::npos) {
			this->alpha_M_row[0] = stod(value);
		} else if (key.find("dna.D_M[1]") != string::npos) {
			this->alpha_M_row[1] = stod(value);
		} else if (key.find("dna.D_M[2]") != string::npos) {
			this->alpha_M_row[2] = stod(value);
		} else if (key.find("dna.D_M[3]") != string::npos) {
			this->alpha_M_row[3] = stod(value);
		} else if (key.find("dna.D_pi[0]") != string::npos) {
			this->alpha_pi[0] = stod(value);
		} else if (key.find("dna.D_pi[1]") != string::npos) {
			this->alpha_pi[1] = stod(value);
		} else if (key.find("dna.D_pi[2]") != string::npos) {
			this->alpha_pi[2] = stod(value);
		} else if (key.find("dna.D_pi[3]") != string::npos) {
			this->alpha_pi[3] = stod(value);
		} else if (key.find("dna.thr") != string::npos) {
			this->ecdllConvergenceThreshold = stod(value);
		} else if (key.find("dna.maxIter") != string::npos) {
			this->maxIter = stoi(value);
		} else if (key.find("dna.thrPerRegime.coarse") != string::npos) {
			this->thrCoarse = stod(value);
		} else if (key.find("dna.thrPerRegime.medium") != string::npos) {
			this->thrMedium = stod(value);
		} else if (key.find("dna.thrPerRegime.fine") != string::npos) {
			this->thrFine = stod(value);
		} else if (key.find("dna.patternCounts.coarse") != string::npos) {
			this->cum_pattern_weight_coarse = stod(value);
		} else if (key.find("dna.patternCounts.medium") != string::npos) {
			this->cum_pattern_weight_medium = stod(value);
		} 
	}
}
// void SEM::BuildCumFractionIndicesAA(double start, double end, double step) {
//     this->cum_fraction_indices_AA.clear();

//     if (cum_fractional_weight_AA.empty() || step <= 0.0) return;

//     // Clamp start/end to [0,1], allow tiny epsilon beyond 1.0 to catch rounding.
//     start = std::max(0.0, start);
//     end   = std::min(1.0000001, end);

//     // Use integer loop count to avoid FP drift.
//     int n_steps = static_cast<int>(std::floor((end - start) / step + 0.5)) + 1;
//     cum_fraction_indices_AA.reserve(n_steps);

//     for (int k = 0; k < n_steps; ++k) {
//         double thresh = start + k * step;
//         if (thresh > 1.0) thresh = 1.0; // final clamp

//         auto it = std::lower_bound(cum_fractional_weight_AA.begin(),
//                                    cum_fractional_weight_AA.end(),
//                                    thresh);
//         if (it != cum_fractional_weight_AA.end()) {
//             cum_fraction_indices_AA.push_back(
//                 static_cast<int>(std::distance(cum_fractional_weight_AA.begin(), it))
//             );
//         } else {
//             // Threshold not reached (shouldn’t happen if end <= 1), fallback to last index
//             cum_fraction_indices_AA.push_back(
//                 static_cast<int>(cum_fractional_weight_AA.size() - 1)
//             );
//         }
//     }
// }


void SEM::SetAASiteCumWeight() {
	this->AAsitecumweight.clear();	
	for (int site = 0; site < this->num_aa_patterns; site++ ) {
		this->AAsitecumweight.emplace_back(site,this->AAPatternWeights[site]);
				
	}
	sort(this->AAsitecumweight.begin(), this->AAsitecumweight.end(), [](const pair<int,int>& a, const pair<int,int>& b) {
        return a.second > b.second; });

	int total_weight = 0;
    for (const auto& p : this->AAsitecumweight) total_weight += p.second;
	this->cum_fractional_weight_AA.clear();
    this->cum_fractional_weight_AA.reserve(this->AAsitecumweight.size());

	double cumulative = 0.0;
    for (const auto& p : this->AAsitecumweight) {
        cumulative += p.second;
        this->cum_fractional_weight_AA.push_back(cumulative / total_weight);
    }
}

Matrix20d SEM::ExpQt_fromSpectrum(double t) const {
    if (!spectrum.ok) {    
        return (t * Q_D).exp();
    }
    
    Vector20cd d  = (spectrum.lambda.array() * t).exp();
    Matrix20cd Pc = (spectrum.V * d.asDiagonal()) * spectrum.V_inv;
    Matrix20d P = Pc.real();
    
    for (int i = 0; i < 20; ++i) {
        double sum = 0.0;
        for (int j = 0; j < 20; ++j) {
            double x = P(i, j);
            if (x < 0.0 && x > -1e-12) x = 0.0;
            P(i, j) = x;
            sum += x;
        }
        if (sum > 0.0) P.row(i) /= sum;
    }
    return P;
}


void SEM::PrecomputeSpectrum() {
    this->spectrum = Spectrum{};
    Eigen::ComplexEigenSolver<Matrix20d> ces;
    ces.compute(Q_D);
    if (ces.info() == Eigen::Success) {
        spectrum.V      = ces.eigenvectors();
        spectrum.lambda = ces.eigenvalues();
        spectrum.V_inv  = spectrum.V.partialPivLu().inverse();
        spectrum.ok = true;
    } else {
        spectrum.ok = false;
    }
}

void SEM::AssignLSBranchLengthsToTree() {
    if (this->DM_LS.empty()) {
		cout << "AssignOLSBranchLengthsToTree: DM_NR is empty. Build the DM_NR distance map first." << endl;
        throw std::runtime_error("AssignLSBranchLengthsToTree: DM_NR is empty. Build the DM_NR distance map first.");
    }
    // Refit OLS lengths on the current topology using NR distances
    FamilyJoining::AssignOLSLengthsToTree(this->T, this->DM_LS);
    FamilyJoining::SetNegativeLengthsToZero(this->T);
}

void SEM::SetFJTree(const fj::Graph& tree) {
    this->T = tree;
    this->T.rebuild();
}


void SEM::ComputeNRDistanceMatrix() {
    // --- settings ---
    const double T_MIN   = 1e-9;     // strictly positive lower bound
    const double T_MAX   = 10.0;     // hard cap for branch length
    const int    MAX_NR  = 60;       // max Newton iterations per pair
    const int    MAX_BT  = 60;       // max backtracking steps
    const double G_TOL   = 1e-10;    // gradient tolerance for stop
    const double X_TOL   = 1e-8;     // parameter change tolerance
    const double ALPHA0  = 1.0;      // initial step scale
    const double BT_SHR  = 0.5;      // backtracking shrink factor
    const double EPS_LOG = 1e-300;   // floor for log()
    const double EPS_DEN = 1e-16;    // floor to avoid /0 in derivatives

    this->DM_NR.clear();

    const int n = this->numberOfInputSequences;
    if (n <= 1) return;

    // Stationary distribution as VectorXd
    Eigen::VectorXd pi_vec(20);
    for (int k = 0; k < 20; ++k) pi_vec[k] = this->pi_D[k];

    // Precompute reversible spectral factors once for this Q_D
    Eigen::VectorXd piN, pis, invpis, lambdas;
    Eigen::MatrixXd U;
    std::tie(piN, pis, invpis, U, lambdas) = ReversibleSpecTuple(this->Q_D, pi_vec);

    // Canonical key for unordered pair
    const auto keyOf = [](const std::string& a, const std::string& b) {
        return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
    };

    // Build P, P', P'' via reversible spectral form at time t
    auto P_and_derivs = [&](double t,
                            Eigen::Matrix<double,20,20>& P,
                            Eigen::Matrix<double,20,20>& Pd,
                            Eigen::Matrix<double,20,20>& Pdd) {
        Eigen::Matrix<double,20,20> S0 = Eigen::Matrix<double,20,20>::Zero();
        Eigen::Matrix<double,20,20> S1 = Eigen::Matrix<double,20,20>::Zero();
        Eigen::Matrix<double,20,20> S2 = Eigen::Matrix<double,20,20>::Zero();

        for (int k = 0; k < 20; ++k) {
            const double lam = lambdas[k];
            const double e   = std::exp(lam * t);
            const auto   uk  = U.col(k);
            S0.noalias() += e           * (uk * uk.transpose());
            S1.noalias() += (lam*e)     * (uk * uk.transpose());
            S2.noalias() += (lam*lam*e) * (uk * uk.transpose());
        }

        for (int i = 0; i < 20; ++i) {
            for (int j = 0; j < 20; ++j) {
                const double pre  = invpis[i];
                const double post = pis[j];
                P  (i,j) = pre * S0 (i,j) * post;
                Pd (i,j) = pre * S1 (i,j) * post;
                Pdd(i,j) = pre * S2 (i,j) * post;
            }
        }
        // Do NOT renormalize rows here—keeps derivatives exact.
    };

    // log-likelihood, gradient, Hessian at t for counts N
    auto ll_g_h = [&](const Eigen::Matrix<double,20,20>& N, double t) {
        Eigen::Matrix<double,20,20> P, Pd, Pdd;
        P_and_derivs(t, P, Pd, Pdd);

        double ll = 0.0, g = 0.0, h = 0.0;

        for (int a = 0; a < 20; ++a) {
            const double log_pi_a = (pi_vec[a] > 0.0) ? std::log(pi_vec[a]) : std::log(EPS_LOG);
            for (int b = 0; b < 20; ++b) {
                const double n_ab = N(a,b);
                if (n_ab <= 0.0) continue;

                const double Pab      = P(a,b);
                const double Pab_log  = std::log((Pab > 0.0) ? Pab : EPS_LOG);
                const double Pab_safe = (std::abs(Pab) > EPS_DEN) ? Pab : (Pab >= 0.0 ? EPS_DEN : -EPS_DEN);

                const double Pdab = Pd(a,b);
                const double Pddb = Pdd(a,b);

                ll += n_ab * (log_pi_a + Pab_log);

                const double r = Pdab / Pab_safe;
                const double q = Pddb / Pab_safe - r*r;

                g += n_ab * r;
                h += n_ab * q;
            }
        }
        return std::tuple<double,double,double>(ll, g, h);
    };

    
    for (int i = 0; i < n; ++i) {
        SEM_vertex* u = this->vertices[i];
        const std::string& ui = u->name;

        for (int j = i + 1; j < n; ++j) {
            SEM_vertex* v = this->vertices[j];
            const std::string& vj = v->name;

            Eigen::Matrix<double,20,20> N = GetCountMatrix(u, v);
            if (N.sum() <= 0.0) { this->DM_NR[keyOf(ui, vj)] = 0.0; continue; }

            // Seed with JC(20); clamp to bounds
            double t = std::max(T_MIN, std::min(T_MAX, std::max(0.0, GetJukesCantorDistance(u, v))));

            double ll, g, h;
            std::tie(ll, g, h) = ll_g_h(N, t);

            // Newton–Raphson with backtracking line search
            for (int it = 0; it < MAX_NR; ++it) {
                if (std::abs(g) < G_TOL) break;

                double step;
                if (h < 0.0 && std::isfinite(h)) step = g / h; // true Newton (maximize)
                else                             step = g / ( (h==0.0)? 1.0 : std::abs(h) + 1e-12 );

                double alpha = ALPHA0;
                double t_new = t - alpha * step;

                for (int bt = 0; bt < MAX_BT; ++bt) {
                    if (!(t_new >= T_MIN && t_new <= T_MAX) || !std::isfinite(t_new)) {
                        alpha *= BT_SHR;
                        t_new = t - alpha * step;
                        continue;
                    }
                    double ll_new, g_new, h_new;
                    std::tie(ll_new, g_new, h_new) = ll_g_h(N, t_new);
                    if (ll_new >= ll || alpha < 1e-6) {
                        t = t_new; ll = ll_new; g = g_new; h = h_new;
                        break;
                    }
                    alpha *= BT_SHR;
                    t_new = t - alpha * step;
                }

                if (!std::isfinite(t)) break;
                if (std::abs(alpha * step) < X_TOL) break;
            }

            if (!(t >= 0.0) || !std::isfinite(t)) {
                double t_jc = std::max(0.0, GetJukesCantorDistance(u, v));
                t = std::max(T_MIN, std::min(T_MAX, t_jc));
            }

            this->DM_NR[keyOf(ui, vj)] = t;
        }
    }
}


// Builds an n x n symmetric matrix of MLE branch lengths between all pairs in `verts`.
// using line search
void SEM::ComputeLSDistanceMatrix() {
    // --- settings ---
    const double T_EPS   = 1e-9;                 // evaluate slightly > 0
    const double T_MAX   = 10.0;                 // hard cap for branch length
    const double T_TOL   = 1e-8;                 // interval tolerance
    const int    MAX_IT  = 120;                  // golden-section iterations
    const double PHI     = 0.6180339887498949;   // golden ratio conjugate
    const double IMPROVE = 1e-9;                 // min improvement when expanding bracket
    const double EPS_LL  = 1e-300;               // floor for probabilities in log

    this->DM_LS.clear();

    const int n = this->numberOfInputSequences; 
    if (n <= 1) return;

    // Stationary distribution as VectorXd
    Eigen::VectorXd pi_vec(20);
    for (int k = 0; k < 20; ++k) pi_vec[k] = this->pi_D[k];

    // Precompute reversible spectral factors once
    Eigen::VectorXd piN, pis, invpis, lambdas;
    Eigen::MatrixXd U;
    std::tie(piN, pis, invpis, U, lambdas) = ReversibleSpecTuple(this->Q_D, pi_vec);

    // Canonical key: (a,b) sorted
    const auto keyOf = [](const std::string& a, const std::string& b) {
        return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
    };

    // LL-only evaluator using spectral form (no g/H)
    auto ll_of_counts = [&](const Eigen::Matrix<double,20,20>& N, double t) -> double {
        // UdexpUt = sum_k e^{λ_k t} u_k u_k^T
        Eigen::Matrix<double,20,20> UdexpUt = Eigen::Matrix<double,20,20>::Zero();
        for (int k = 0; k < 20; ++k) {
            double e = std::exp(lambdas[k] * t);
            Eigen::Matrix<double,20,1> uk = U.col(k);
            UdexpUt.noalias() += e * (uk * uk.transpose());
        }
        // P = Π^{-1/2} * UdexpUt * Π^{1/2}
        Eigen::Matrix<double,20,20> P;
        for (int i = 0; i < 20; ++i)
            for (int j = 0; j < 20; ++j)
                P(i,j) = invpis[i] * UdexpUt(i,j) * pis[j];

        // hygiene: clip tiny negatives, renormalize each row
        for (int i = 0; i < 20; ++i) {
            for (int j = 0; j < 20; ++j)
                if (P(i,j) < 0.0 && P(i,j) > -1e-14) P(i,j) = 0.0;
            double rsum = P.row(i).sum();
            if (std::abs(rsum - 1.0) > 1e-10) P.row(i) /= rsum;
        }

        // log-likelihood
        double ll = 0.0;
        for (int a = 0; a < 20; ++a)
            for (int b = 0; b < 20; ++b) {
                double n_ab = N(a,b);
                if (n_ab <= 0.0) continue;
                double pab = P(a,b);
                if (pab < EPS_LL) pab = EPS_LL;
                ll += n_ab * std::log(pab);
            }
        return ll;
    };

    // All unordered pairs
    for (int i = 0; i < n; ++i) {
        SEM_vertex* u = this->vertices[i];
        const std::string& ui = u->name;

        for (int j = i + 1; j < n; ++j) {
            SEM_vertex* v = this->vertices[j];
            const std::string& vj = v->name;

            Eigen::Matrix<double,20,20> N = GetCountMatrix(u, v);
            if (N.sum() <= 0.0) { this->DM_LS[keyOf(ui, vj)] = 0.0; continue; }

            // Seed with JC(20)
            double t_seed = std::max(0.0, GetJukesCantorDistance(u, v));

            // ---- bracket [a,b] by expanding b until no improvement or T_MAX ----
            double a = 0.0;
            double b = std::min(T_MAX, std::max(0.25, 2.0 * t_seed + 1.0));
            double fb = ll_of_counts(N, b);
            while (b < T_MAX) {
                double b2 = std::min(T_MAX, 2.0 * b);
                double fb2 = ll_of_counts(N, b2);
                if (fb2 > fb + IMPROVE) { b = b2; fb = fb2; }
                else break;
            }

            // ---- golden-section search on [left, right] to maximize LL ----
            double left  = std::max(T_EPS, a);
            double right = std::max(left + T_TOL, b);

            double c = right - PHI * (right - left);
            double d = left  + PHI * (right - left);
            double fc = ll_of_counts(N, c);
            double fd = ll_of_counts(N, d);

            for (int it = 0; it < MAX_IT; ++it) {
                if (fc < fd) {
                    left = c;
                    c = d;
                    fc = fd;
                    d = left + PHI * (right - left);
                    fd = ll_of_counts(N, d);
                } else {
                    right = d;
                    d = c;
                    fd = fc;
                    c = right - PHI * (right - left);
                    fc = ll_of_counts(N, c);
                }
                if (std::abs(right - left) <= T_TOL) break;
            }

            // best of c, d, midpoint
            double tm = 0.5 * (left + right);
            double fm = ll_of_counts(N, tm);
            double t_star = tm;
            double fbest  = fm;
            if (fc > fbest) { fbest = fc; t_star = c; }
            if (fd > fbest) { fbest = fd; t_star = d; }

            if (!(t_star >= 0.0) || !std::isfinite(t_star)) {
                double t_jc = std::max(0.0, GetJukesCantorDistance(u, v));
                t_star = std::isfinite(t_jc) ? t_jc : 0.0;
            }

            this->DM_LS[keyOf(ui, vj)] = t_star;
        }
    }
}

void SEM::BuildSIndex(std::vector<std::pair<int,int>>& idxS, int i0, int j0) const {
	idxS.clear();
    for (int i = 0; i < 20; ++i)
        for (int j = i + 1; j < 20; ++j)
            if (!(i == i0 && j == j0))
                idxS.emplace_back(i, j);
}

void SEM::PackTheta(const std::vector<std::pair<int,int>>& idxS, int p_ref,
                    std::vector<double>& theta) const {
    theta.clear();
    theta.reserve(idxS.size() + 19);

    // S (positivity via exp)
    for (std::size_t k = 0; k < idxS.size(); ++k) {
        int i = idxS[k].first, j = idxS[k].second;
        double sij = this->S_D(i, j);
        if (!(sij > 0.0)) sij = 1e-6;
        theta.push_back(std::log(sij));
    }

    // pi (19 logits relative to p_ref)
    double pref = this->pi_D[p_ref];
    if (!(pref > 0.0)) pref = 1e-6;
    for (int i = 0; i < 20; ++i) {
        if (i == p_ref) continue;
        double pii = this->pi_D[i];
        if (!(pii > 0.0)) pii = 1e-12;
        theta.push_back(std::log(pii / pref));
    }
}

void SEM::UnpackTheta(const std::vector<double>& theta,
                      const std::vector<std::pair<int,int>>& idxS,
                      int i0, int j0, int p_ref) {
    const std::size_t mS = idxS.size();
    // S: fill only free pairs; keep (i0,j0) as-is
    for (std::size_t k = 0; k < mS; ++k) {
        int i = idxS[k].first, j = idxS[k].second;
        double sij = std::exp(theta[k]);
        if (!(sij > 0.0)) sij = 1e-12;
        this->S_D(i, j) = this->S_D(j, i) = sij;
    }

    // pi: anchored softmax
    double y[20]; for (int i=0;i<20;++i) y[i] = 0.0;
    y[p_ref] = 1.0;
    double sumy = 1.0;
    std::size_t t = mS;
    for (int i = 0; i < 20; ++i) {
        if (i == p_ref) continue;
        double yi = std::exp(theta[t++]);
        if (!(yi > 0.0)) yi = 1e-12;
        y[i] = yi; sumy += yi;
    }
    for (int i = 0; i < 20; ++i) this->pi_D[i] = y[i] / sumy;
}

double SEM::NMObjectiveNegLL(const std::vector<double>& theta,
                             const std::vector<std::pair<int,int>>& idxS,
                             int i0, int j0, int p_ref) {
    // Save
    Eigen::Matrix<double,20,20> S_save = this->S_D;
    Eigen::Matrix<double,20,20> Q_save = this->Q_D;
    std::vector<double>          pi_save = this->pi_D;

    double result = std::numeric_limits<double>::infinity();
    try {
        // Apply params
        this->UnpackTheta(theta, idxS, i0, j0, p_ref);

        // Build calibrated mean-rate-1 Q
        this->SetRateMatrixParameters();

        // Compute LL (stored in member)
        this->ComputeAALogLikelihood();
        double ll = this->AAlogLikelihood;
        if (std::isfinite(ll)) result = -ll;
    } catch (...) {
        // leave result as +inf
    }

    // Restore
    this->S_D  = S_save;
    this->Q_D  = Q_save;
    this->pi_D = pi_save;

    return result;
}

void SEM::OptimizeRateMatrixNelderMead(int max_iter, double tol, double step, int i0, int j0, int p_ref) {
    // Free S indices
    std::vector<std::pair<int,int>> idxS;
    this->BuildSIndex(idxS, i0, j0);
    const std::size_t D = idxS.size() + 19; // 189 + 19 = 208

    // Initial theta from current S_D, pi_D
    std::vector<double> x0;
    this->PackTheta(idxS, p_ref, x0);

    // Initial simplex of size D+1
    std::vector<std::vector<double>> X(D + 1, x0);
    for (std::size_t d = 0; d < D; ++d) X[d + 1][d] += step;

    auto f = [&](const std::vector<double>& x)->double {
        return this->NMObjectiveNegLL(x, idxS, i0, j0, p_ref);
    };

    std::vector<double> F(D + 1);
    for (std::size_t r = 0; r < D + 1; ++r) F[r] = f(X[r]);

    const double alpha = 1.0, gamma = 2.0, rho = 0.5, sigma = 0.5;

    std::vector<double> xc(D), xr(D), xe(D), xoc(D), xic(D);
    for (int it = 0; it < max_iter; ++it) {
        // Order simplex by objective
        std::vector<std::size_t> ord(D + 1);
        for (std::size_t r = 0; r < D + 1; ++r) ord[r] = r;
        std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b){ return F[a] < F[b]; });

        // Reorder
        std::vector<std::vector<double>> Xn(D + 1);
        std::vector<double> Fn(D + 1);
        for (std::size_t r = 0; r < D + 1; ++r) { Xn[r] = X[ord[r]]; Fn[r] = F[ord[r]]; }
        X.swap(Xn); F.swap(Fn);

        // Convergence checks
        double fspan = F.back() - F.front();
        double xspan = 0.0;
        for (std::size_t d = 0; d < D; ++d) {
            double vmin = X[0][d], vmax = X[0][d];
            for (std::size_t r = 1; r < D + 1; ++r) { vmin = std::min(vmin, X[r][d]); vmax = std::max(vmax, X[r][d]); }
            xspan = std::max(xspan, std::abs(vmax - vmin));
        }
        if (fspan < tol && xspan < tol) break;

        // Centroid of best D vertices
        for (std::size_t d = 0; d < D; ++d) {
            double s = 0.0; for (std::size_t r = 0; r < D; ++r) s += X[r][d];
            xc[d] = s / (double)D;
        }

        // Reflection
        for (std::size_t d = 0; d < D; ++d) xr[d] = xc[d] + alpha * (xc[d] - X[D][d]);
        double fr = f(xr);

        if (fr < F[0]) {
            // Expansion
            for (std::size_t d = 0; d < D; ++d) xe[d] = xc[d] + gamma * (xr[d] - xc[d]);
            double fe = f(xe);
            if (fe < fr) { X[D] = xe; F[D] = fe; }
            else         { X[D] = xr; F[D] = fr; }
        } else if (fr < F[D - 1]) {
            // Accept reflection
            X[D] = xr; F[D] = fr;
        } else {
            if (fr < F[D]) {
                // Outside contraction
                for (std::size_t d = 0; d < D; ++d) xoc[d] = xc[d] + rho * (xr[d] - xc[d]);
                double foc = f(xoc);
                if (foc <= fr) { X[D] = xoc; F[D] = foc; }
                else {
                    // Shrink towards best
                    for (std::size_t r = 1; r < D + 1; ++r) {
                        for (std::size_t d = 0; d < D; ++d) X[r][d] = X[0][d] + sigma * (X[r][d] - X[0][d]);
                        F[r] = f(X[r]);
                    }
                }
            } else {
                // Inside contraction
                for (std::size_t d = 0; d < D; ++d) xic[d] = xc[d] - rho * (xc[d] - X[D][d]);
                double fic = f(xic);
                if (fic < F[D]) { X[D] = xic; F[D] = fic; }
                else {
                    // Shrink
                    for (std::size_t r = 1; r < D + 1; ++r) {
                        for (std::size_t d = 0; d < D; ++d) X[r][d] = X[0][d] + sigma * (X[r][d] - X[0][d]);
                        F[r] = f(X[r]);
                    }
                }
            }
        }
    }

    // Apply best parameters to this SEM instance and finalize Q
    const std::vector<double>& xbest = X[0];
    this->UnpackTheta(xbest, idxS, i0, j0, p_ref);
    this->SetRateMatrixParameters();
}

double SEM::GetMLDistance_NR_Reversible_Tuples(SEM_vertex* u, SEM_vertex* v, const Eigen::Matrix<double,20,20>& Q, const Eigen::VectorXd& pi, int max_iter, double grad_tol, double ll_tol) {
    Eigen::Matrix<double,20,20> N = GetCountMatrix(u, v);
    if (N.sum() <= 0.0) return 0.0;
	Eigen::VectorXd piN;
	Eigen::VectorXd pis;
	Eigen::VectorXd invpis;
	Eigen::MatrixXd U;
	Eigen::VectorXd lambdas;
    
    tie (piN, pis, invpis, U, lambdas) = ReversibleSpecTuple(Q, pi);

    
    double t = std::max(0.0, GetJukesCantorDistance(u, v));
	double ll, g, H;

    tie (ll, g, H) = LL_Grad_Hess_Counts_Reversible_Tuple(N, pis, invpis, U, lambdas, t);

    constexpr double ALPHA_MIN = 1e-6;

    for (int it = 0; it < max_iter; ++it) {
        if (std::abs(g) <= grad_tol) break;

        
        double step;
        if (H < -1e-12) {
            step = g / H;    
        } else {
            double denom = std::max(std::abs(H), 1.0);
            step = -g / denom;
        }

        // Backtracking with t >= 0
        double alpha = 1.0;
        double t_new = t, ll_new = ll, g_new = g, H_new = H;

        for (int bt = 0; bt < 60; ++bt) {
            t_new = t - alpha * step;
            if (t_new < 0.0) t_new = 0.0;
			
			double ll_bt, g_bt, H_bt;
            tie (ll_bt, g_bt, H_bt) = LL_Grad_Hess_Counts_Reversible_Tuple(N, pis, invpis, U, lambdas, t_new);
            ll_new = ll_bt; g_new = g_bt; H_new = H_bt;

            if (ll_new >= ll || alpha < ALPHA_MIN) break;
            alpha *= 0.5;
        }

        if (ll_new < ll && std::abs(g) <= 10*grad_tol) break;

        // Accept step
        t  = t_new;
        ll = ll_new;
        g  = g_new;
        H  = H_new;

        if (std::abs(g) <= grad_tol) break;
        if (std::abs(ll_new - ll) < ll_tol) break;
    }

    if (!(t >= 0.0)) t = 0.0;
    return t;
}

tuple<double,double,double> SEM::LL_Grad_Hess_Counts_Reversible_Tuple(const Eigen::Matrix<double,20,20>& N, const Eigen::VectorXd& pis, const Eigen::VectorXd& invpis, const Eigen::MatrixXd& U, const Eigen::VectorXd& lambdas, double t) const {
    constexpr double EPS = 1e-300;
	Eigen::Matrix<double,20,20> P, QP, Q2P;
    tie(P, QP, Q2P) = EvalP_QP_Q2P_Tuple(pis, invpis, U, lambdas, t);

    double ll = 0.0, g = 0.0, H = 0.0;
    for (int a = 0; a < 20; ++a) {
        for (int b = 0; b < 20; ++b) {
            double n = N(a,b);
            if (n <= 0.0) continue;
            double pab = std::max(P(a,b), EPS);
            double r   = QP(a,b)  / pab;
            double s   = Q2P(a,b) / pab;
            ll += n * std::log(pab);
            g  += n * r;
            H  += n * (s - r*r);
        }
    }
    return std::make_tuple(ll, g, H);
}

tuple<Eigen::Matrix<double,20,20>, Eigen::Matrix<double,20,20>, Eigen::Matrix<double,20,20>> SEM::EvalP_QP_Q2P_Tuple(const Eigen::VectorXd& pis, const Eigen::VectorXd& invpis, const Eigen::MatrixXd& U, const Eigen::VectorXd& lambdas, double t) const {
    Eigen::Matrix<double,20,20> UdexpUt   = Eigen::Matrix<double,20,20>::Zero();
    Eigen::Matrix<double,20,20> UdlexpUt  = Eigen::Matrix<double,20,20>::Zero();
    Eigen::Matrix<double,20,20> UdllExpUt = Eigen::Matrix<double,20,20>::Zero();

    for (int k = 0; k < 20; ++k) {
        double e   = std::exp(lambdas[k] * t);
        double le  = lambdas[k] * e;
        double l2e = lambdas[k] * lambdas[k] * e;
        Eigen::VectorXd uk = U.col(k);
        Eigen::MatrixXd  outer = uk * uk.transpose(); // rank-1
        UdexpUt   += e   * outer;
        UdlexpUt  += le  * outer;
        UdllExpUt += l2e * outer;
    }

    Eigen::Matrix<double,20,20> P   = UdexpUt;
    Eigen::Matrix<double,20,20> QP  = UdlexpUt;
    Eigen::Matrix<double,20,20> Q2P = UdllExpUt;

    for (int i = 0; i < 20; ++i)
        for (int j = 0; j < 20; ++j) {
            P(i,j)   = invpis[i] * P(i,j)   * pis[j];
            QP(i,j)  = invpis[i] * QP(i,j)  * pis[j];
            Q2P(i,j) = invpis[i] * Q2P(i,j) * pis[j];
        }

    // Numerical hygiene: clip tiny negatives and row-normalize P
    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 20; ++j)
            if (P(i,j) < 0.0 && P(i,j) > -1e-14) P(i,j) = 0.0;
        double rsum = P.row(i).sum();
        if (std::abs(rsum - 1.0) > 1e-10) P.row(i) /= rsum;
    }

    return std::make_tuple(P, QP, Q2P);
}

tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::MatrixXd, Eigen::VectorXd> SEM::ReversibleSpecTuple(const Eigen::Matrix<double,20,20>& Q, const Eigen::VectorXd& pi_in) const {
    Eigen::VectorXd pi = pi_in;
    if (pi.size() != 20) throw std::invalid_argument("pi must have size 20");
    if ((pi.array() <= 0.0).any()) throw std::invalid_argument("pi must be strictly positive");
    double s = pi.sum();
    if (std::abs(s - 1.0) > 1e-10) pi /= s;

    Eigen::VectorXd pis    = pi.array().sqrt().matrix();
    Eigen::VectorXd invpis = pis.array().inverse().matrix();

    Eigen::Matrix<double,20,20> S;
    for (int i = 0; i < 20; ++i)
        for (int j = 0; j < 20; ++j)
            S(i,j) = pis[i] * Q(i,j) * invpis[j];

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,20,20>> es(S);
    if (es.info() != Eigen::Success)
        throw std::runtime_error("Eigen decomposition failed for symmetric S = Π^{1/2}QΠ^{-1/2}");

    Eigen::MatrixXd U     = es.eigenvectors();    // orthonormal
    Eigen::VectorXd evals = es.eigenvalues();     // real

    return std::make_tuple(pi, pis, invpis, U, evals);
}

double SEM::DyadicLogLikelihood(SEM_vertex * p, SEM_vertex * c, double t) {
	double ll = 0.0;	
	Eigen::Matrix <double,20,20> P_t;
	P_t = Q_D*t;
	P_t = P_t.exp();
	int aa_p; int aa_c;	
	for (int i = 0; i < num_aa_patterns; i++) {		
		aa_p = p->AAcompressed[i];
		aa_c = c->AAcompressed[i];
		if (aa_p > -1 && aa_c > -1) {
			ll += log(this->pi_D[aa_p]) * this->AAPatternWeights[i];
			ll += log(P_t(aa_p,aa_c)) * this->AAPatternWeights[i];
		}
	}
	return (ll);
}

Eigen::Matrix <double,20,20> SEM::GetCountMatrix(SEM_vertex * p, SEM_vertex * c) {
	Eigen::Matrix <double,20,20> N;
	N.setZero();
	int aa_p; int aa_c;
	for (int i = 0; i < this->num_aa_patterns; i++) {
		aa_p = p->AAcompressed[i];
		aa_c = c->AAcompressed[i];
		if (aa_p > -1 && aa_c > -1) {
			N(aa_p,aa_c) += this->AAPatternWeights[i];
		}
	}
	return (N);
}

vector <double> SEM::GetCountVector(SEM_vertex * p) {
	vector <double> C(20, 0.0);
	int aa_p;
	for (int i = 0; i < this->num_aa_patterns; i++) {
		aa_p = p->AAcompressed[i];
		if (aa_p > -1) C[aa_p] += this->AAPatternWeights[i];
	}
	return (C);	
}

double SEM::GetJukesCantorDistance(SEM_vertex * u, SEM_vertex * v) {
	double t, p;
	p = 0;
	Eigen::Matrix <double,20,20> N = GetCountMatrix(u,v);
	double total_count = 0;
	for (int i = 0; i < 20; i++) {
		for (int j = 0; j < 20; j++){
			total_count += N(i,j);
			if (i != j) p += N(i,j);
		}
	}
	p/=total_count;
	if (p > 19.0/20.0) throw mt_error("check p distance"); 
	t = -19.0/20.0 * log(1-((20.0*p)/19.0));
	return (t);
}



double SEM::DyadicLogLikelihoodFirstDer(SEM_vertex * p, SEM_vertex * c, double t) {
	double ll_first_der = 0.0;
	Eigen::Matrix<double,20,20> P_t;
	Eigen::Matrix<double,20,20> QP_t;
	P_t = Q_D*t;
	P_t = P_t.exp();
	QP_t = Q_D*P_t;
	int aa_p; int aa_c;
	for (int i = 0; i < num_aa_patterns; i++) {
		aa_p = p->AAcompressed[i];
		aa_c = c->AAcompressed[i];
		if (aa_p > -1 && aa_c > -1) {
			ll_first_der += (QP_t(aa_p,aa_c)/P_t(aa_p,aa_c)) * this->AAPatternWeights[i];
		}
	}
	return (ll_first_der);
}

double SEM::GetNewtonRaphsonDistance(SEM_vertex * u, SEM_vertex * v) {		
	float p = 0;
	double t;
	int aa_u; int aa_v;
	float sites_included = 0;
	double ll_fd;
	double ll;
	double nr_convergence = pow(10,-2);
	int max_iter = 100;
	int iter = 2;
	// compute initial estimate of evolutionary distance using Jukes-Cantor estimate
	// cout << u->name << "\t" << v->name << endl;
	for (int site = 0; site < num_aa_patterns; site++) {
		aa_u = u->AAcompressed[site];
		aa_v = v->AAcompressed[site];
		if (aa_u > -1 && aa_v > -1) {
			sites_included += this->AAPatternWeights[site];
			if (aa_u != aa_v) p += this->AAPatternWeights[site];
		}		
	}
	p /= sites_included;
	t = -19.0/20.0 * log(1-((20.0*p)/19.0));
	cout << "initial distance is " << t << endl;	
	ll = this->DyadicLogLikelihood(u,v,t);
	ll_fd = this->DyadicLogLikelihoodFirstDer(u,v,t);
	t -= ll/ll_fd;
	cout << "Distance after first iteration of NR is " << t << endl;
	cout << "Log likelihood after first iteration of NR is " << ll << endl;
	cout << "First derivate of ll after first iteration of NR is " << ll_fd << endl;	
	while (abs(ll_fd) > nr_convergence && iter < max_iter) {
		ll = this->DyadicLogLikelihood(u,v,t);
		ll_fd = this->DyadicLogLikelihoodFirstDer(u,v,t);
		t -= ll/ll_fd;
		cout << "Distance after iteration " << iter << " of NR is " << t << endl;
		cout << "Log likelihood after iteration " << iter << " of NR is " << ll << endl;
		cout << "First derivate after iteration " << iter << " of NR is " << ll_fd << endl;
		iter++;
	}
	return (t);
}

void SEM::ComputeMLDistancesForAAsequences() {
	double ML_d;
	for (SEM_vertex * u: this->vertices) {
		for (SEM_vertex * v: this->vertices) {
			if (u < v) {
				ML_d = this->GetNewtonRaphsonDistance(u,v);
				this->ML_distances[make_pair(u,v)] = ML_d;
				break;
			}
		}
	}
}



void SEM::SetDayhoffExchangeabilityAndRootProbabilityParameters() {
	ifstream Dayhoff_file(this->dayhoff_rate_matrix_file_name);

	if (!Dayhoff_file.is_open()){
		throw mt_error("Could not open Dayhoff rate matrix file");
	}
	this->Q_D.setZero();
	this->S_D.setZero();	
	string line;
	vector <string> line_split;
	int line_num = 0;
	int row; int col; int K = 20;
	double smallest_entry_in_dayhoff = 0;
	double rate_matrix_entry;
    while (getline(Dayhoff_file, line)) {
		line_num++;
		if (line.empty()) continue;
		// cout << "line " << line_num << " is " << line << endl;
		line_split = split_ws(line);
		if (line_num < K + 1) {
			if (line_split.size() != line_num -1) throw mt_error("Check the format of Dayhoff matrix");
			row = line_num - 1;
			for (col = 0; col < row; col++) {
				rate_matrix_entry = stod(line_split[col]);				
				this->S_D(row,col) = rate_matrix_entry;
				this->S_D(col,row) = rate_matrix_entry;
			}
		} else {
			if (line_split.size() != K) throw mt_error("Check the format of Dayhoff matrix");			
			for (int i = 0; i < K; i++) {
				this->pi_D[i] = stod(line_split[i]);
			}			
		}
	}
	Dayhoff_file.close();
	
	double sum_pi = 0.0;
	for (int i = 0; i < 20; i++) {
		sum_pi += this->pi_D[i];
	}

	for (int i = 0; i < 20; i++) {
		this->pi_D[i] /= sum_pi;
	}

}

void SEM::SetRateMatrixParameters() {
	int K = 20;

	// set rate matrix entries for off-diagonal elements
	// cout << "here 1" << endl;
	for (int row = 0; row < K; row++) {
		for (int col = 0; col < row; col++) {
			this->Q_D(row,col) = S_D(row,col) * this->pi_D[col];
			this->Q_D(col,row) = S_D(col,row) * this->pi_D[row];
		}
	}	
	// set rate matrix entries for diagonal elements
	// cout << "here 2" << endl;

	for (int row = 0; row < K; row++) {
    	double rowsum = 0.0;
    	for (int col = 0; col < K; col++) {
        	if (col == row) continue;
        	rowsum += this->Q_D(row, col);
    	}
    	this->Q_D(row, row) = -rowsum;
	}	

	// cout << "here 3" << endl;

	double sum_pi = 0.0;
	for (int i = 0; i < 20; i++) {
		sum_pi += this->pi_D[i];
	}
	
	// cout << "sum_pi is " << sum_pi << endl;

	// cout << " abs(sum_pi - 1.0) is " << abs(sum_pi - 1.0) << endl;

    if (abs(sum_pi - 1.0) > 1e-12) {
        throw mt_error("pi is not normalized");
    }

	// cout << "here 4" << endl;

    // compute mean rate: -sum_i pi[i] * Q[i,i]
    double mu = 0.0;
    for (int i = 0; i < K; i++) {
        mu += -1 * this->pi_D[i] * this->Q_D(i, i);
    }
	// cout << "mu is " << mu << endl;	
    if (mu <= 0.0) {
        throw mt_error("Non-positive mean rate; check Q and pi");
    }	
    this->Q_D /= mu;
	// cout << "Dayhoff rate matrix set" << endl;

}


array <double, 4> SEM::sample_pi() {
	return sample_dirichlet(this->alpha_pi, rng());
}

array <double, 4> SEM::sample_M_row() {
	return sample_dirichlet(this->alpha_M_row, rng());
}

array <double, 4> SEM::sample_dirichlet(const array<double, 4>& alpha, mt19937_64& gen) {
        array<double, 4> x{};
        double sum = 0.0;
        for (size_t i = 0; i < 4; ++i) {
            gamma_distribution<double> gamma(alpha[i], 1.0);
            x[i] = gamma(gen);
            sum += x[i];
        }
        for (auto& v : x) v /= sum;
        return x;
    }

void SEM::set_alpha_PI(double a1, double a2, double a3, double a4) {
	this->alpha_pi[0] = a1;
	this->alpha_pi[1] = a2;
	this->alpha_pi[2] = a3;
	this->alpha_pi[3] = a4;	
}

void SEM::set_alpha_M_row(double a1, double a2, double a3, double a4) {
	this->alpha_M_row[0] = a1;
	this->alpha_M_row[1] = a2;
	this->alpha_M_row[2] = a3;
	this->alpha_M_row[3] = a4;
}


void SEM::AddDuplicatedSequencesToRootedTree(MST * M) {
	// Store dupl seq names in uniq seq vertex
	double t;
	string uniq_seq_name;
	vector <string> dupl_seq_name_vec;
	SEM_vertex * u;
	SEM_vertex * p;
	SEM_vertex * d;
	SEM_vertex * h;
	vector <int> emptySequence;	
	int v_id = this->vertexMap->size() - 1;
	for (pair <string, vector <string> > uniq_seq_name_2_dupl_seq_name_vec : M->unique_seq_id_2_dupl_seq_ids) {
		uniq_seq_name = uniq_seq_name_2_dupl_seq_name_vec.first;
		dupl_seq_name_vec = uniq_seq_name_2_dupl_seq_name_vec.second;
		u = (*this->vertexMap)[this->nameToIdMap[uniq_seq_name]];
		p = u->parent;
		t = this->edgeLengths[make_pair(p,u)];
		
		v_id += 1;
		h = new SEM_vertex(v_id,emptySequence);
		this->vertexMap->insert(pair<int,SEM_vertex*>(h->id,h));
		h->name = "h_" + to_string(this->node_ind);
		this->nameToIdMap.insert(make_pair(h->name,h->id));
		this->node_ind += 1;
		
		u->RemoveParent();
		p->RemoveChild(u);
		this->RemoveEdgeLength(p,u);
		
		h->AddParent(p);
		p->AddChild(h);
		this->AddEdgeLength(p,h,t);

		u->AddParent(h);
		h->AddChild(u);
		this->AddEdgeLength(h,u,0.0);

		for (string dupl_seq_name: dupl_seq_name_vec) {
			v_id += 1;
			d = new SEM_vertex(v_id,emptySequence);
			d->name = dupl_seq_name;
			this->vertexMap->insert(pair<int,SEM_vertex*>(d->id,d));
			this->nameToIdMap.insert(make_pair(d->name,d->id));			
			d->AddParent(h);	
			h->AddChild(d);
			this->AddEdgeLength(h,d,0.0);			
		}
	}
	this->flag_added_duplicated_sequences = 1;
}

void SEM::AddDuplicatedSequencesToUnrootedTree(MST * M) {
	// Store dupl seq names in uniq seq vertex
	double t;
	string uniq_seq_name;
	vector <string> dupl_seq_name_vec;
	SEM_vertex * l;
	SEM_vertex * n;
	SEM_vertex * d;
	SEM_vertex * h;
	vector <int> emptySequence;
	// vector <SEM_vertex *> uniq_vertex_ptr_vec;
	int v_id = this->vertexMap->size() - 1;
	for (pair <string, vector <string> > uniq_seq_name_2_dupl_seq_name_vec : M->unique_seq_id_2_dupl_seq_ids) {
		uniq_seq_name = uniq_seq_name_2_dupl_seq_name_vec.first;
		dupl_seq_name_vec = uniq_seq_name_2_dupl_seq_name_vec.second;
		l = (*this->vertexMap)[this->nameToIdMap[uniq_seq_name]];
		n = l->neighbors[0];
		t = this->edgeLengths[make_pair(n,l)];
		
		v_id += 1;
		h = new SEM_vertex(v_id,emptySequence);
		this->vertexMap->insert(pair<int,SEM_vertex*>(h->id,h));
		h->name = "h_" + to_string(this->node_ind);
		this->nameToIdMap.insert(make_pair(h->name,h->id));
		this->node_ind += 1;
		
		l->RemoveNeighbor(n);
		n->RemoveNeighbor(l);
		this->RemoveEdgeLength(n,l);
		
		h->AddNeighbor(n);
		n->AddNeighbor(h);
		this->AddEdgeLength(n,h,t);

		h->AddNeighbor(l);
		l->AddNeighbor(h);
		this->AddEdgeLength(h,l,0.0);

		for (string dupl_seq_name: dupl_seq_name_vec) {
			v_id += 1;
			d = new SEM_vertex(v_id,emptySequence);
			d->name = dupl_seq_name;
			this->vertexMap->insert(pair<int,SEM_vertex*>(d->id,d));
			this->nameToIdMap.insert(make_pair(d->name,d->id));			
			d->AddNeighbor(h);	
			h->AddNeighbor(d);
			this->AddEdgeLength(h,d,0.0);			
		}
	}
}


void SEM::SetEdgeAndVertexLogLikelihoods() {
	SEM_vertex * u;	SEM_vertex * v;
	int u_id; int v_id;	
	this->vertexLogLikelihoodsMapToAddToGlobalPhylogeneticTree.clear();
	this->edgeLogLikelihoodsToAddToGlobalPhylogeneticTree.clear();
	for (pair <int, int> edge_ind : this->edgesOfInterest_ind) {
		tie (u_id, v_id) = edge_ind;
		u = (*this->vertexMap)[u_id];
		v = (*this->vertexMap)[v_id];
		if (this->vertexLogLikelihoodsMapToAddToGlobalPhylogeneticTree.find(u->name) == this->vertexLogLikelihoodsMapToAddToGlobalPhylogeneticTree.end()) {
			this->ComputeVertexLogLikelihood(u);
			this->vertexLogLikelihoodsMapToAddToGlobalPhylogeneticTree.insert(pair<string,double>(u->name,u->vertexLogLikelihood));
		}	
		if (this->vertexLogLikelihoodsMapToAddToGlobalPhylogeneticTree.find(v->name) == this->vertexLogLikelihoodsMapToAddToGlobalPhylogeneticTree.end()) {
			this->ComputeVertexLogLikelihood(v);
			this->vertexLogLikelihoodsMapToAddToGlobalPhylogeneticTree.insert(pair<string,double>(v->name,v->vertexLogLikelihood));
		}
		this->ComputeEdgeLogLikelihood(u,v);
		this->ComputeEdgeLogLikelihood(v,u);
	}	
}

void SEM::ComputeVertexLogLikelihood(SEM_vertex * v) {
	array <double, 4> prob = this->posteriorProbabilityForVertex[v];
	array <double, 4> Counts = this->expectedCountsForVertex[v];
	v->vertexLogLikelihood = 0;
	for (int i = 0; i < 4; i ++) {
		if (prob[i] > 0) {
			v->vertexLogLikelihood += (log(prob[i]) * Counts[i]);
		}		
	}	
}


void SEM::ComputeEdgeLogLikelihood(SEM_vertex* x, SEM_vertex * y) {	
	emtr::Md P_xy = this->GetPosteriorProbabilityForVariablePair(x,y);
	emtr::Md P_yGivenx = this->GetP_yGivenx(P_xy);
	emtr::Md Counts = this->GetExpectedCountsForVariablePair(x,y);
	double edgeLogLikelihood = 0;
	for (int dna_x = 0; dna_x < 4; dna_x ++) {
		for (int dna_y = 0; dna_y < 4; dna_y ++) {
			if (P_yGivenx[dna_x][dna_y] > 0) {
				edgeLogLikelihood += (log(P_yGivenx[dna_x][dna_y]) * Counts[dna_x][dna_y]);
			}
		}
	}
	this->edgeLogLikelihoodsToAddToGlobalPhylogeneticTree.push_back(tuple <string, string, double>(x->name,y->name,edgeLogLikelihood));
}

void SEM::SetWeightedEdgesToAddToGlobalPhylogeneticTree() {
	this->weightedEdgesToAddToGlobalPhylogeneticTree.clear();
	int u_id; int v_id;
	SEM_vertex * u; SEM_vertex * v; 
	string u_name; string v_name;
	double t;	
//	cout << "Adding the following edges to the global phylogenetic tree" << endl;
	for (pair <int, int> edge_ind : this->edgesOfInterest_ind) {
		tie (u_id, v_id) = edge_ind;
		u = (*this->vertexMap)[u_id];
		v = (*this->vertexMap)[v_id];
		u_name = u->name;
		v_name = v->name;
		t = this->ComputeEdgeLength(u,v);
//		cout << u_name << "\t" << v_name << "\t" << t << endl;
		this->weightedEdgesToAddToGlobalPhylogeneticTree.push_back(make_tuple(u_name,v_name,t));
	}
}

void SEM::SetAncestralSequencesString() {
	vector <SEM_vertex *> verticesOfInterest;
	int u_id; int v_id;
	vector <int> fullSeq;
	string DNAString;
	SEM_vertex * u;	SEM_vertex * v;
	this->ancestralSequencesString = "";
	for (pair <int, int> edge_ind : this->edgesOfInterest_ind) {
		tie(u_id, v_id) = edge_ind;		
		u = (*this->vertexMap)[u_id];		
		v = (*this->vertexMap)[v_id];
		if (!u->observed and find(verticesOfInterest.begin(),verticesOfInterest.end(),u)==verticesOfInterest.end()) {
			fullSeq = this->DecompressSequence(u->DNArecoded,this->sitePatternRepetitions);
			DNAString = EncodeAsDNA(fullSeq);	
			this->ancestralSequencesString += ">"; 
			this->ancestralSequencesString += u->name + "\n";
			this->ancestralSequencesString += DNAString + "\n";
			verticesOfInterest.push_back(u);
			
		}		
		if (!v->observed and find(verticesOfInterest.begin(),verticesOfInterest.end(),v)==verticesOfInterest.end()) {
			fullSeq = this->DecompressSequence(v->DNArecoded,this->sitePatternRepetitions);
			DNAString = EncodeAsDNA(fullSeq);
			this->ancestralSequencesString += ">"; 
			this->ancestralSequencesString += v->name + "\n";
			this->ancestralSequencesString += DNAString + "\n";
			verticesOfInterest.push_back(v);
		}		
	}	
}


void SEM::SetNeighborsBasedOnParentChildRelationships() {
	this->ClearUndirectedEdges();
	SEM_vertex * p; SEM_vertex * c;
	for (pair<SEM_vertex*, SEM_vertex*> edge : this->edgesForPreOrderTreeTraversal) {
		tie (p, c) = edge;
		p->AddNeighbor(c);
		c->AddNeighbor(p);
	}
}

void SEM::SetIdsForObservedVertices(vector <int> idsOfObservedVerticesToAdd) {
	this->idsOfObservedVertices = idsOfObservedVerticesToAdd;
}

void SEM::SetNumberOfInputSequences(int numOfInputSeqsToSet) {
	this->numberOfInputSequences = numOfInputSeqsToSet;	
}
void SEM::ComputeBIC() {
	this->ComputeLogLikelihood();
	this->BIC = -2.0 * this->logLikelihood;
	double n = this->sequenceLength;
	double numberOfFreeParameters = this->edgesForPostOrderTreeTraversal.size();
	numberOfFreeParameters += 11.0 * (this->numberOfRateCategories);
	bool rootHasDistinctRateCat = 1;
	for (SEM_vertex * v : this->root->children) {
		if (this->root->rateCategory == v->rateCategory) {
			rootHasDistinctRateCat = 0;
		}
	}
	if (rootHasDistinctRateCat) {
		numberOfFreeParameters += 3.0;
	}
	this->BIC += log(n) * numberOfFreeParameters;
}

void SEM::SetModelSelectionCriterion(string modelSelectionCriterionToSet) {
	this->modelSelectionCriterion = modelSelectionCriterionToSet;
}

void SEM::SetFlagForFinalIterationOfSEM() {
	this->finalIterationOfSEM = 1;
}


void SEM::ResetData() {
	for (pair<int,SEM_vertex*> idPtrPair : *this->vertexMap) {
		if (idPtrPair.first != -1) {
			delete idPtrPair.second;
		}
	}
	if(this->vertexMap->size()!=1){
        throw mt_error("vertexMap not set correctly");
    }
	(*this->vertexMap)[-1]->DNArecoded.clear();	
}

SEM_vertex * SEM::GetVertex(string v_name){
	bool contains_v = this->ContainsVertex(v_name);
	SEM_vertex * node;
	int node_id;
	if(!contains_v) {
        throw mt_error("v_name not found");
    }
	
	node_id = this->GetVertexId(v_name);		
	node = (*this->vertexMap)[node_id];
	return node;

}

int SEM::GetVertexId(string v_name) {
	SEM_vertex * v;
	int idToReturn = -10;
	for (pair<int,SEM_vertex*> idPtrPair : *this->vertexMap) {
		v = idPtrPair.second;
		if (v->name == v_name){
			idToReturn = v->id;						
		}
	}
	if (idToReturn == -10){
		cout << "Unable to find id for:" << v_name << endl;
	}
	return (idToReturn);
}

void SEM::SuppressRoot() {
	SEM_vertex * c_l;
	SEM_vertex * c_r;
	bool proceed = this->root->outDegree == 2;		
	if (proceed) {		
		c_l = this->root->children[0];		
		c_r = this->root->children[1];		
		c_l->AddNeighbor(c_r);
		c_r->AddNeighbor(c_l);
		c_l->RemoveNeighbor(this->root);
		c_r->RemoveNeighbor(this->root);
		this->RemoveArc(this->root,c_l);
		this->RemoveArc(this->root,c_r);
	}
}

void SEM::SwapRoot() {
	SEM_vertex * root_current;
	SEM_vertex * vertexNamedHRoot;
	vector <SEM_vertex *> childrenOfCurrentRoot;
	vector <SEM_vertex *> childrenOfVertexNamedHRoot;
	int n = this->numberOfObservedVertices;	
	if (this->root->name != "h_root") {
//		this->SetEdgesForPostOrderTraversal();		
//		this->ComputeLogLikelihood();		
		root_current = this->root;
		childrenOfCurrentRoot = root_current->children;
		
		vertexNamedHRoot = (*this->vertexMap)[((2*n)-2)];		
		childrenOfVertexNamedHRoot = vertexNamedHRoot->children;
		
		// Swap children of root
		for (SEM_vertex * c: childrenOfVertexNamedHRoot) {
			this->RemoveArc(vertexNamedHRoot,c);
			this->AddArc(root_current,c);
		}
		
		for (SEM_vertex * c: childrenOfCurrentRoot) {			
			this->RemoveArc(root_current,c);
			this->AddArc(vertexNamedHRoot,c);
		}
		
		vertexNamedHRoot->rootProbability = root_current->rootProbability;
		root_current->transitionMatrix = vertexNamedHRoot->transitionMatrix;
		vertexNamedHRoot->transitionMatrix = this->I4by4;
		
		this->AddArc(vertexNamedHRoot->parent,root_current);
		this->RemoveArc(vertexNamedHRoot->parent,vertexNamedHRoot);
		this->root = vertexNamedHRoot;
		this->SetLeaves();	
		this->SetEdgesForPreOrderTraversal();
		this->SetVerticesForPreOrderTraversalWithoutLeaves();
		this->SetEdgesForPostOrderTraversal();
//		this->ComputeLogLikelihood();
	}	
}

int SEM::GetEdgeIndex (int vertexIndex1, int vertexIndex2, int numberOfVertices) {
	int edgeIndex;
	edgeIndex = numberOfVertices*(numberOfVertices-1)/2;
	edgeIndex -= (numberOfVertices-vertexIndex1)*(numberOfVertices-vertexIndex1-1)/2;
	edgeIndex += vertexIndex2 - vertexIndex1 - 1;
	return edgeIndex;
}

char SEM::GetDNAfromIndex(int dna_index){
	char bases[4] = {'A', 'C', 'G', 'T'};
	if (dna_index > -1 && dna_index <4){
		return bases[dna_index];
	} else {
		return 'Z';
	}
}

// Dayhoff order indices (K = 20):
// 0:A 1:R 2:N 3:D 4:C 5:Q 6:E 7:G 8:H 9:I
// 10:L 11:K 12:M 13:F 14:P 15:S 16:T 17:W 18:Y 19:V
int SEM::ConvertAAtoIndex(char aa) {
    switch (std::toupper(static_cast<unsigned char>(aa))) {
    case 'A': return 0;  // Ala
    case 'R': return 1;  // Arg
    case 'N': return 2;  // Asn
    case 'D': return 3;  // Asp
    case 'C': return 4;  // Cys
    case 'Q': return 5;  // Gln
    case 'E': return 6;  // Glu
    case 'G': return 7;  // Gly
    case 'H': return 8;  // His
    case 'I': return 9;  // Ile
    case 'L': return 10; // Leu
    case 'K': return 11; // Lys
    case 'M': return 12; // Met
    case 'F': return 13; // Phe
    case 'P': return 14; // Pro
    case 'S': return 15; // Ser
    case 'T': return 16; // Thr
    case 'W': return 17; // Trp
    case 'Y': return 18; // Tyr
    case 'V': return 19; // Val
    default:
        // Ambiguous/unknown/stop/gap: B, Z, J, X, U, O, '*', '-',
        // or any other symbol -> treat as missing.
        return -1;
    }
}



int SEM::ConvertDNAtoIndex(char dna){
	int value;
	switch (dna)
	{
	case 'A':
		value = 0;
		break;
	case 'C':
		value = 1;
		break;
	case 'G':
		value = 2;
		break;
	case 'T':
		value = 3;
		break;
	default:
		value = -1;        
		break;
	}	
	return (value);
}

void SEM::SetVertexVector(){
	this->vertices.clear();
	SEM_vertex * v;
	for (pair <int, SEM_vertex *> idPtrPair : * this->vertexMap) {
		v = idPtrPair.second;
		this->vertices.push_back(v);
	}
}


void SEM::SetBHparameters() {
	this->ReadProbabilities();	
}

void SEM::ReparameterizeBH() {
	
	// Compute pi, P(u,v) and P(v,u) for each vertex and edge using the method described in hss paper	
	
	// 1. Set pi for root	
	for (int p_id = 0; p_id < 4; ++p_id) {
		if(this->rootProbability[p_id] == 0) {cout << "invalid probability"; throw mt_error("invalid probability");}
		this->root->root_prob_hss[p_id] = this->rootProbability[p_id];
	}
	// Store transition matrices in transition_prob_hss map, store root prob in node as root_prob_hss
	SEM_vertex * p; SEM_vertex * c; array <double,4> pi_p; array <double,4> pi_c;
	M_hss->clear();
	for (pair<SEM_vertex*, SEM_vertex*> edge : this->edgesForPreOrderTreeTraversal) {
		p = edge.first;
		c = edge.second;
		
		emtr::Md M_pc = c->transitionMatrix; // transition matrix of orig BH parameters 
		emtr::Md M_cp; 					     // transition matrix of reparameterized BH
						
		(*this->M_hss)[{p,c}] = M_pc;

		// 2. Initialize pi_p and pi_c
		for (int x = 0; x < 4; x ++) pi_p[x] = p->root_prob_hss[x]; // root prob already computed
		for (int x = 0; x < 4; x ++) pi_c[x] = 0;					// root prob to be computed
		
		// 3. Compute pi_c
		for (int x = 0; x < 4 ; x ++){
			for (int y = 0; y < 4; y ++){
				pi_c[x] += pi_p[y] * M_pc[y][x];
			}
		}

		// 4. Store pi_c
		for (int x = 0; x < 4; x++) {
			c->root_prob_hss[x] = pi_c[x];
		}

		// 5. Compute M_cp for root at child		
		for (int x = 0; x < 4; x++) {
			for (int y = 0; y < 4; y++) {
				M_cp[y][x] = M_pc[x][y] * pi_p[x]/pi_c[y];			// Bayes rule as described in HSS paper
			}
		}
		
		// 6. Store M_cp
		(*this->M_hss)[{c,p}] = M_cp;
	}	
}

void SEM::ReadProbabilities() {
    std::ifstream inputFile(this->probabilityFileName_best);
    if (!inputFile) {
        throw mt_error("Failed to open probability file: " + this->probabilityFileName_best);
    }

    std::string line;
    // skip first two header lines
    if (!std::getline(inputFile, line) || !std::getline(inputFile, line)) {
        throw mt_error("Probability file too short (missing headers)");
    }

    std::string node_name;
	// string node_parent_name;
    double prob; int i, j;

    while (std::getline(inputFile, line)) {
        vector<string> splitLine = emtr::split_ws(line);
        const int num_words = static_cast<int>(splitLine.size());
		// cout << num_words << endl;
        switch (num_words) {
            case 8: { 
				// node_parent_name = splitLine[5];
                node_name = splitLine[7];
                break;
            }
            case 16: {                
                SEM_vertex* n = this->GetVertex(node_name);
                for (int p_id = 0; p_id < 16; ++p_id) {
                    i = p_id / 4;
                    j = p_id % 4;
                    try	{
						prob = stod(splitLine[p_id]); 
					}
					catch(const exception& e) {
						prob = 0;
						cout << "setting to 0 small prob value " << splitLine[p_id] << " not converted by stod" << endl;
						// (*this->logFile) << "setting to 0 small prob value " << splitLine[p_id] << " not converted by stod" << endl;
					}
					n->transitionMatrix[i][j] = prob;
                }
                break;
            }
            case 9: {                
                node_name = splitLine[3];
                SEM_vertex* n = this->GetVertex(node_name);
                this->RootTreeAtVertex(n);
                break;
            }
            case 4: {                
                for (int p_id = 0; p_id < 4; ++p_id) {
                    prob = std::stod(splitLine[p_id]);
                    this->rootProbability[p_id] = prob;
                }
                break;
            }
            default:
                std::cerr << "ReadProbabilities: unexpected token count (" << num_words
                          << ") on line: " << line << "\n";
                break;
        }
    }
}

void SEM::WriteProbabilities(string fileName) {
	ofstream probabilityFile;
	probabilityFile.open(fileName);
	SEM_vertex * v;
	probabilityFile << "transition matrix for edge from " << "parent_name" << " to " << "child_name" << endl;
	
	for (int row = 0; row < 4; row++) {
		for (int col = 0; col < 4; col++) {
			probabilityFile << "p(" << this->GetDNAfromIndex(col) << "|" << this->GetDNAfromIndex(row) << ")";
			if (row ==3 and col ==3) {
				continue;
			} else {
				probabilityFile << " ";
			}
			
		}		
	}
	probabilityFile << endl;

	for (pair <int, SEM_vertex *> idPtrPair : * this->vertexMap) {
		v = idPtrPair.second;
		if (v != v->parent) {
			probabilityFile << "transition matrix for edge from " << v->parent->name << " to " << v->name << endl;
			for (int row = 0; row < 4; row++) {
				for (int col = 0; col < 4; col++) {
					probabilityFile << v->transitionMatrix[row][col];
					if (row ==3 and col ==3){
						continue;
					} else {
						probabilityFile << " ";
					}
				}				
			}
			probabilityFile << endl;
		}
	}

	probabilityFile << "probability at root " << this->root->name << " is ";
	for (int row = 0; row < 3; row++) {
		probabilityFile << "p(" << this->GetDNAfromIndex(row) << ") ";		
	}
	probabilityFile << "p(" << this->GetDNAfromIndex(3) << ") " << endl;

	for (int row = 0; row < 3; row++) {
		probabilityFile << this->rootProbability[row] << " ";
	}
	probabilityFile << this->rootProbability[3] << endl;
	probabilityFile.close();
}

void SEM::ReadRootedTree(string treeFileName) {
	string u_name;
	string v_name;
	int u_id;
	int v_id;
	SEM_vertex * u;
	SEM_vertex * v;
	vector <string> splitLine;
	vector <string> leafNames;
	vector <string> ancestorNames;
	vector <string> nonRootVertexNames;	
	string rootName = "";
	vector <int> emptySequence;
	v_id = 0;
	ifstream edgeListFile(treeFileName.c_str());
	for (string line; getline(edgeListFile, line);) {
		vector<string> splitLine = emtr::split_ws(line);
		u_name = splitLine[0];		
		v_name = splitLine[1];
		if (find(nonRootVertexNames.begin(),nonRootVertexNames.end(),v_name) == nonRootVertexNames.end()) {
			nonRootVertexNames.push_back(v_name);
		}		
		if (find(ancestorNames.begin(),ancestorNames.end(),u_name)==ancestorNames.end()) {
			ancestorNames.push_back(u_name);
		}
		if (find(leafNames.begin(),leafNames.end(),v_name)==leafNames.end()) {
			if(!emtr::starts_with(v_name, "h_")) {
				leafNames.push_back(v_name);
			}
		}
	}
	for (string name: leafNames) {
		SEM_vertex * v = new SEM_vertex(v_id, emptySequence);
		v->name = name;
		v->observed = 1;
		this->vertexMap->insert(pair<int,SEM_vertex*>(v_id,v));
		v_id += 1;
	}
	// Remove root from ancestor names
	for (string name: ancestorNames) {
		if (find(nonRootVertexNames.begin(),nonRootVertexNames.end(),name)==nonRootVertexNames.end()){
			rootName = name;
		}
	}
	this->numberOfObservedVertices = leafNames.size();
	int n = this->numberOfObservedVertices;			
	// Change root name
	
	ancestorNames.erase(remove(ancestorNames.begin(), ancestorNames.end(), rootName), ancestorNames.end());
	for (string name: ancestorNames) {
		SEM_vertex * v = new SEM_vertex(v_id,emptySequence);
		v->name = name;
		this->vertexMap->insert(pair <int,SEM_vertex*> (v_id,v));
		v_id += 1;
	}
	
	this->root = new SEM_vertex (((2 * n) - 2), emptySequence);	
	this->root->name = rootName;
	this->vertexMap->insert(pair <int,SEM_vertex*> (((2 * n) - 2), this->root));
	edgeListFile.clear();
	edgeListFile.seekg(0, ios::beg);
	for (string line; getline(edgeListFile, line);) {
		vector<string> splitLine = emtr::split_ws(line);
		u_name = splitLine[0];
		v_name = splitLine[1];
		u_id = this->GetVertexId(u_name);
		v_id = this->GetVertexId(v_name);
		u = (*this->vertexMap)[u_id];
		v = (*this->vertexMap)[v_id];
		u->AddChild(v);
		v->AddParent(u);
	}
	edgeListFile.close();	
	this->SetLeaves();
	// cout << "Number of leaves is " << this->leaves.size() << endl;
	this->SetEdgesForPreOrderTraversal();
	// cout << "Number of edges for pre order traversal is " << this->edgesForPreOrderTreeTraversal.size() << endl;
	this->SetVerticesForPreOrderTraversalWithoutLeaves();
	// cout << "Number of vertices for pre order traversal is " << this->preOrderVerticesWithoutLeaves.size() << endl;
	this->SetEdgesForPostOrderTraversal();
	// cout << "Number of edges for post order traversal is " << this->edgesForPostOrderTreeTraversal.size() << endl;
}

bool SEM::IsTreeInCanonicalForm() {
	bool valueToReturn = 1;
	SEM_vertex * v;
	for (pair<int,SEM_vertex*> idPtrPair : *this->vertexMap) {
		v = idPtrPair.second;
		if ((!v->observed) and v->outDegree != 2) {
			valueToReturn = 0;
		}
		if (v->observed and v->outDegree != 0) {
			valueToReturn = 0;
		}
	}
	return (valueToReturn);
}


void SEM::AddAllSequences(string fileName) {
	vector <int> recodedSequence;
	ifstream inputFile(fileName.c_str());
	string v_name;
	string seq = "";	
	int v_id;
	vector <string> vertexNames;	
	vector <vector <int>> allSequences;	
	for (string line; getline(inputFile, line );) {
		if (line[0]=='>') {
			if (seq != "") {				
				for (char const dna: seq) {
					recodedSequence.push_back(mapDNAtoInteger[string(1,toupper(dna))]);					
				}				
				v_id = this->GetVertexId(v_name);				
				(*this->vertexMap)[v_id]->DNArecoded = recodedSequence;				
				recodedSequence.clear();
			}
			v_name = line.substr(1,line.length());
			seq = "";
		} else {
			seq += line ;
		}
	}
	inputFile.close();
	
	for (char const dna: seq) {
		recodedSequence.push_back(mapDNAtoInteger[string(1,toupper(dna))]);
	}
	
	v_id = this->GetVertexId(v_name);
	(*this->vertexMap)[v_id]->DNArecoded = recodedSequence;
	
	int numberOfSites = recodedSequence.size();	
	this->num_dna_patterns = numberOfSites;
	this->sequenceLength = numberOfSites;
	recodedSequence.clear();
	
	this->DNAPatternWeights.clear();
	
	for (int i = 0; i < numberOfSites; i++) {
		this->DNAPatternWeights.push_back(1);
	}	
}

void SEM::ResetAncestralSequences() {
	vector<int> gappedSequence ;
	for (int i = 0; i < this->num_dna_patterns; i++) gappedSequence.push_back(-1);
	for (pair <int,SEM_vertex*> idPtrPair : *this->vertexMap) {
		if (!idPtrPair.second->observed) {
			// cout << " resetting for " << idPtrPair.second->name << endl;
			idPtrPair.second->DNAcompressed = gappedSequence;
		}
	}
}

void SEM::RemoveEdgeLength(SEM_vertex * u, SEM_vertex * v) {
	pair <SEM_vertex *, SEM_vertex *> vertexPair;
	if (u->id < v->id) {
		vertexPair = make_pair(u,v);
	} else {
		vertexPair = make_pair(v,u);
	}
	this->edgeLengths.erase(vertexPair);
}

void SEM::AddEdgeLength(SEM_vertex * u, SEM_vertex * v, double t) {	
	pair <SEM_vertex *, SEM_vertex *> vertexPair;
	if (u->id < v->id) {
		vertexPair = make_pair(u,v);
	} else {
		vertexPair = make_pair(v,u);
	}
	this->edgeLengths[vertexPair] = t;
}

double SEM::GetEdgeLength(SEM_vertex * u, SEM_vertex * v) {
	double t;
	pair <SEM_vertex *, SEM_vertex *> vertexPair;
	if (u->id < v->id) {
		vertexPair = make_pair(u,v);
	} else {
		vertexPair = make_pair(v,u);
	}
	t = this->edgeLengths[vertexPair];
	return (t);
}

double SEM::ComputeEdgeLength(SEM_vertex * u, SEM_vertex * v) {
	double t = 0;
	int dna_u; int dna_v; 
	for (int site = 0; site < this->num_dna_patterns; site++) {
		dna_u = u->DNArecoded[site];
		dna_v = v->DNArecoded[site];
		if (dna_u != dna_v) {
			t += this->DNAPatternWeights[site];
		}
	}
	t /= this->sequenceLength;	
	return (t);
}

void SEM::SetEdgeLength(SEM_vertex * u, SEM_vertex * v, double t) {
	pair <SEM_vertex *, SEM_vertex *> vertexPair;
	if (u->id < v->id) {
		vertexPair = make_pair(u,v);
	} else {
		vertexPair = make_pair(v,u);
	}
	this->edgeLengths[vertexPair] = t;
}

void SEM::StoreEdgeListAndSeqToAdd() {
	this->weightedEdgeListString = "";		
	this->sequencesToAddToGlobalPhylogeneticTree.clear();
	this->weightedEdgesToAddToGlobalPhylogeneticTree.clear();
	SEM_vertex * u; SEM_vertex * v;	
	double t;	
	for (pair <SEM_vertex *, SEM_vertex *> vertexPair : this->edgesForPostOrderTreeTraversal) {
		tie (u, v) = vertexPair;		
		if (u->parent != u) {
			if (v == this->externalVertex and !this->finalIterationOfSEM) {
				this->compressedSequenceToAddToMST = u->DNArecoded;
				this->nameOfSequenceToAddToMST = u->name;				
			} else {
				t = this->ComputeEdgeLength(u,v);			
//				cout << "Adding edge 1 " << u->name << "\t" << v->name << endl;
				this->weightedEdgeListString += u->name + "\t" + v->name + "\t" + to_string(t) + "\n";
				this->weightedEdgesToAddToGlobalPhylogeneticTree.push_back(make_tuple(u->name,v->name,t));
			}
		}
	}	
	u = this->root->children[0];
	v = this->root->children[1];
	if ((v != this->externalVertex and u!= this->externalVertex) or this->finalIterationOfSEM) {
		t = this->ComputeEdgeLength(u,v);
//		cout << "Adding edge 2 " << u->name << "\t" << v->name << endl;
		this->weightedEdgeListString += u->name + "\t" + v->name + "\t" + to_string(t) + "\n";
		this->weightedEdgesToAddToGlobalPhylogeneticTree.push_back(make_tuple(u->name,v->name,t));
	} else if (u == this->externalVertex) {
		this->compressedSequenceToAddToMST = v->DNArecoded;
		this->nameOfSequenceToAddToMST = v->name;
	} else {
		if (v != this->externalVertex){
            throw mt_error("v should be equal to external vertex");
        }
		this->compressedSequenceToAddToMST = u->DNArecoded;
		this->nameOfSequenceToAddToMST = u->name;
	}
//	cout << "Name of external sequence is " << this->externalVertex->name << endl;
//	cout << "Name of sequence to add to MST is " << this->nameOfSequenceToAddToMST << endl;
	// Add sequences of all vertices except the following vertices
	// 1) root, 2) external vertex
	for (pair <int, SEM_vertex * > idPtrPair : * this->vertexMap) {
		u = idPtrPair.second;		
		if (u->parent != u){
			if (u != this->externalVertex) {
				this->sequencesToAddToGlobalPhylogeneticTree[u->name] = u->DNArecoded;
			} else if (this->finalIterationOfSEM) {
				this->sequencesToAddToGlobalPhylogeneticTree[u->name] = u->DNArecoded;
			}
		}		
	}	
}



emtr::Md SEM::ComputeTransitionMatrixUsingAncestralStatesForTrifle(SEM_vertex * p, SEM_vertex * c, int layer) {	
	emtr::Md P = emtr::Md{};
	// for (int i = 0; i < 4; i ++) {		
	// 	P[i][i] = 1.0;
	// }			
	int dna_p; int dna_c;
	
	for (int site = 0; site < this->num_patterns_for_layer[layer]; site ++) {
		dna_p = p->DNAcompressed[site];
		dna_c = c->DNAcompressed[site];
		if (dna_p > -1 && dna_c > -1) {
			P[dna_p][dna_c] += this->DNAPatternWeights[site];
		}
	}

	double rowSum;
	for (int i = 0; i < 4; i ++) {
		rowSum = 0;
		for (int j = 0; j < 4; j ++) {
			rowSum += P[i][j];
		}
		if (rowSum > 0) {
			for (int j = 0; j < 4; j ++) {
				 P[i][j] /= rowSum;
			}
		}		
	}

	// add 0.01 to off diagonal zero entries
	// add 0.97 to diagonal zero entries

	for (int i = 0; i < 4; i ++) {
		for (int j = 0; j < 4; j ++) {
			if (P[i][j] == 0) {
				if (i == j) {
					P[i][j] = 0.97;
				} else {
					P[i][j] = 0.01;
				}
			}
		}
	}

	for (int i = 0; i < 4; i ++) {
		for (int j = 0; j < 4; j ++) {
			rowSum = 0;
			for (int j = 0; j < 4; j ++) {
				rowSum += P[i][j];
			}
			for (int j = 0; j < 4; j ++) {
				P[i][j] /= rowSum;
			}
		}
	}

	return P;
}

emtr::Md SEM::ComputeTransitionMatrixUsingAncestralStates(SEM_vertex * p, SEM_vertex * c) {	
	emtr::Md P = emtr::Md{};			
	int dna_p; int dna_c;
	// cout << p->name << "\t" << c->name << endl;
	for (int site = 0; site < this->num_dna_patterns; site ++) {		
		dna_p = p->DNAcompressed[site];
		dna_c = c->DNAcompressed[site];
		if (dna_p > -1 and dna_c > -1) {
			P[dna_p][dna_c] += this->DNAPatternWeights[site];
		}
	}

	double rowSum;
	for (int i = 0; i < 4; i ++) {
		rowSum = 0;
		for (int j = 0; j < 4; j ++) {
			rowSum += P[i][j];
		}
		for (int j = 0; j < 4; j ++) {
			 P[i][j] /= rowSum;
		}
	}
	return P;
}

void SEM::SetPrefixForOutputFiles(string prefix_for_output_files_to_set){
	this->prefix_for_output_files = prefix_for_output_files_to_set;
}

void SEM::WriteRootedTreeInNewickFormat(string newickFileName) {
	vector <SEM_vertex *> verticesToVisit;
	SEM_vertex * c;
	SEM_vertex * p;	
	double edgeLength;	
	for (pair <int, SEM_vertex *> idAndVertex : * this->vertexMap) {
		idAndVertex.second->timesVisited = 0;
		if (idAndVertex.second->children.size() == 0) {
			idAndVertex.second->newickLabel = idAndVertex.second->name;
			verticesToVisit.push_back(idAndVertex.second);
		} else {
			idAndVertex.second->newickLabel = "";
		}
	}
	
	pair <SEM_vertex *, SEM_vertex * > vertexPair;
	unsigned int numberOfVerticesToVisit = verticesToVisit.size();
	while (numberOfVerticesToVisit > 0) {
		c = verticesToVisit[numberOfVerticesToVisit -1];
		verticesToVisit.pop_back();
		numberOfVerticesToVisit -= 1;
		if (c->parent != c) {
			p = c->parent;
			if (p->id < c->id) {
				vertexPair = make_pair(p,c);
			} else {
				vertexPair = make_pair(c,p);
			}
			p->timesVisited += 1;			
			if (this->edgeLengths.find(vertexPair) == this->edgeLengths.end()) {
				edgeLength = 0.1;
			} else {
				edgeLength = this->edgeLengths[vertexPair];
			}
			if (p->timesVisited == int(p->children.size())) {
				p->newickLabel += "," + c->newickLabel + ":" + to_string(edgeLength) + ")";
				verticesToVisit.push_back(p);
				numberOfVerticesToVisit += 1;
			} else if (p->timesVisited == 1) {
				p->newickLabel += "(" + c->newickLabel + ":" + to_string(edgeLength);
			} else {
				p->newickLabel += "," + c->newickLabel + ":" + to_string(edgeLength);
			}			
		}
	}
	ofstream newickFile;
	newickFile.open(newickFileName);
	newickFile << this->root->newickLabel << ";" << endl;
	newickFile.close();
}

void SEM::WriteCliqueTreeToFile(string cliqueTreeFileName) {
	ofstream cliqueTreeFile;
	cliqueTreeFile.open(cliqueTreeFileName);
	clique * parentClique;
	for (clique * childClique : this->cliqueT->cliques) {
		if (childClique->parent != childClique) {
			parentClique = childClique->parent;
			cliqueTreeFile << parentClique->x->name + "_" + parentClique->y->name +"\t";
			cliqueTreeFile << childClique->x->name + "_" + childClique->y->name << "\t";
			cliqueTreeFile << "0.01" << endl;
		}
	}
	cliqueTreeFile.close();
}

void SEM::WriteUnrootedTreeAsEdgeList(string fileName) {
	ofstream treeFile;
	treeFile.open(fileName);
	SEM_vertex * v;
	double t;
	for (pair <int, SEM_vertex *> idPtrPair : * this->vertexMap) {
		v = idPtrPair.second;
		for (SEM_vertex * n : v->neighbors) {
			if (v->id < n->id) {
				t = this->GetEdgeLength(v,n);
				treeFile << v->name << "\t" << n->name << "\t" << t << endl;
			}
		}
	}
	treeFile.close();
}


void SEM::WriteParametersOfBH(string fileName) {
	ofstream parameterFile;
	parameterFile.open(fileName);		
	
	parameterFile << "Root probability for vertex " << this->root->name << " is " << endl;
	for (int i = 0; i < 4; i++) {
		parameterFile << this->rootProbability[i] << "\t";
	}
	parameterFile << endl;


	SEM_vertex * c; SEM_vertex * p;
	for (pair <int, SEM_vertex *> idPtrPair : * this->vertexMap){
		c = idPtrPair.second;
		p = c->parent;
		if (p != c) {								
			parameterFile << "Transition matrix for " << p->name << " to " << c->name << " is " << endl;
			string trans_par_string = "";
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					trans_par_string.append(to_string(c->transitionMatrix[i][j]) + " ");					
				}
			}
			if (!trans_par_string.empty() && trans_par_string.back() == ' ') trans_par_string.pop_back();    		
			parameterFile << trans_par_string << endl;
		}
	}	
	parameterFile.close();
}

void SEM::WriteRootedTreeAsEdgeList(string fileName) {
	ofstream treeFile;
	treeFile.open(fileName);
	double t;
	SEM_vertex * v;
	for (pair <int, SEM_vertex *> idPtrPair : * this->vertexMap) {
		v = idPtrPair.second;
		if (v != v->parent) {
			t = this->GetEdgeLength(v,v->parent);
			treeFile << v->parent->name << "\t" << v->name << "\t" << t << endl;
		}
	}
	treeFile.close();
}

void SEM::RootTreeAtAVertexPickedAtRandom() {
	cout << "num of observed vertices is " << this->numberOfObservedVertices << endl;	
	int n = this->numberOfObservedVertices;
	uniform_int_distribution <int> distribution_v(n,(2*n-3));
	int v_ind = distribution_v(generator);
	cout << "index of vertex selected for rooting is " << v_ind << endl;
	cout << "number of vertices is " << this->vertexMap->size() << endl;
	SEM_vertex * v = (*this->vertexMap)[v_ind];
	cout << "vertex selected for rooting is " << v->id << endl;
	this->RootTreeAtVertex(v);
	
}

void SEM::RootTreeAlongAnEdgePickedAtRandom() {
	int n = this->numberOfObservedVertices;
//	int numOfVertices = this->vertexMap->size();
	uniform_int_distribution <int> distribution_v(0,(2*n-3));
	int v_ind = distribution_v(generator);
	SEM_vertex * v = (*this->vertexMap)[v_ind];	
	int numOfNeighbors = v->neighbors.size();
//	cout << "Number of neighbors of v are " << numOfNeighbors << endl;
	uniform_int_distribution <int> distribution_u(0,numOfNeighbors-1);	
	int u_ind_in_neighborList = distribution_u(generator);
	SEM_vertex * u = v->neighbors[u_ind_in_neighborList];
//	cout << "Rooting tree along edge ";
//	cout << u->name << "\t" << v->name << endl;
	this->RootTreeAlongEdge(u,v);
}

void SEM::UNRESTForEachBranch() {

}

void SEM::RootTreeAlongEdge(SEM_vertex * u, SEM_vertex * v) {
	// Remove lengths of edges incident to root if necessary
	if (this->root->children.size() == 2) {
		SEM_vertex * c_l = this->root->children[0];
		SEM_vertex * c_r = this->root->children[1];
		this->edgeLengths.erase(make_pair(this->root,c_l));
		this->edgeLengths.erase(make_pair(this->root,c_r));
	}
	this->ClearDirectedEdges();
	SEM_vertex * c;
	this->root->AddChild(u);
	this->root->AddChild(v);
	
	SEM_vertex * c_l = this->root->children[0];
	SEM_vertex * c_r = this->root->children[1];
	this->edgeLengths.insert(make_pair(make_pair(this->root,c_l),0.001));
	this->edgeLengths.insert(make_pair(make_pair(this->root,c_r),0.001));	
	u->AddParent(this->root);
	v->AddParent(this->root);
	vector <SEM_vertex *> verticesToVisit;
	vector <SEM_vertex *> verticesVisited;
	verticesToVisit.push_back(u);
	verticesToVisit.push_back(v);
	verticesVisited.push_back(u);
	verticesVisited.push_back(v);
	int numberOfVerticesToVisit = verticesToVisit.size();
	while (numberOfVerticesToVisit > 0) {		
		c = verticesToVisit[numberOfVerticesToVisit - 1];
		verticesToVisit.pop_back();
		verticesVisited.push_back(c);
		numberOfVerticesToVisit -= 1;
		for (SEM_vertex * n: c->neighbors) {
			if (find(verticesVisited.begin(),verticesVisited.end(),n)==verticesVisited.end()) {
				verticesToVisit.push_back(n);
				numberOfVerticesToVisit += 1;
				c->AddChild(n);
				n->AddParent(c);
			}
		}
	}	
	this->SetLeaves();
//	cout << "Number of leaves is " << this->leaves.size() << endl;
	this->SetEdgesForPreOrderTraversal();
//	cout << "Number of edges for pre order traversal is " << this->edgesForPreOrderTreeTraversal.size() << endl;
	this->SetVerticesForPreOrderTraversalWithoutLeaves();
//	cout << "Number of vertices for pre order traversal is " << this->preOrderVerticesWithoutLeaves.size() << endl;
	this->SetEdgesForPostOrderTraversal();
//	cout << "Number of edges for post order traversal is " << this->edgesForPostOrderTreeTraversal.size() << endl;
}

void SEM::InitializeTransitionMatricesAndRootProbability() {
	// ASR via MP
	// MLE of transition matrices and root probability
		
	vector <SEM_vertex *> verticesToVisit;

	SEM_vertex * p; int numberOfPossibleStates; int pos;
	map <SEM_vertex *, vector <int>> VU;
	map <SEM_vertex *, int> V;
	for (int site = 0; site < this->num_dna_patterns; site++) {
		VU.clear(); V.clear();
		// Set VU and V for leaves;
		for (SEM_vertex * v : this->leaves) {
			V[v] = v->DNArecoded[site];
			VU[v].push_back(v->DNArecoded[site]);
		}
		// Set VU for ancestors
		for (int v_ind = preOrderVerticesWithoutLeaves.size()-1; v_ind > -1; v_ind--) {
			p = this->preOrderVerticesWithoutLeaves[v_ind];
			map <int, int> dnaCount;
			for (int dna = 0; dna < 4; dna++) {
				dnaCount[dna] = 0;
			}
			for (SEM_vertex * c : p->children) {
				for (int dna: VU[c]) {
					dnaCount[dna] += 1;
				}
			}
			int maxCount = 0;
			for (pair<int, int> dnaCountPair: dnaCount) {
				if (dnaCountPair.second > maxCount) {
					maxCount = dnaCountPair.second;
				}
			}			
			for (pair<int, int> dnaCountPair: dnaCount) {
				if (dnaCountPair.second == maxCount) {
					VU[p].push_back(dnaCountPair.first);					
				}
			}			
		}
		// Set V for ancestors
		for (SEM_vertex * v : this->preOrderVerticesWithoutLeaves) {
			if (v->parent == v) {
			// Set V for root
				if (VU[v].size()==1) {
					V[v] = VU[v][0];
				} else {
					numberOfPossibleStates = VU[v].size();
					uniform_int_distribution <int> distribution(0,numberOfPossibleStates-1);
					pos = distribution(generator);
					V[v] = VU[v][pos];
				}				
			} else {
				p = v->parent;
				if (find(VU[v].begin(),VU[v].end(),V[p])==VU[v].end()){
					numberOfPossibleStates = VU[v].size();
					uniform_int_distribution <int> distribution(0,numberOfPossibleStates-1);
					pos = distribution(generator);
					V[v] = VU[v][pos];					
				} else {
					V[v] = V[p];
				}				
			}
			// Push states to compressedSequence
			v->DNArecoded.push_back(V[v]);
		}
	}
}


void SEM::AddToExpectedCountsForEachVariable() {
	SEM_vertex * v;
	double siteWeight = this->DNAPatternWeights[this->cliqueT->site];	
	// Add to counts for each unobserved vertex (C->x) where C is a clique
	array <double, 4> marginalizedProbability;
	vector <SEM_vertex *> vertexList;
	for (clique * C: this->cliqueT->cliques) {
		v = C->x;
		if(v->observed){
            throw mt_error("v should not be am observed vertex");
        }
		if (find(vertexList.begin(),vertexList.end(),v) == vertexList.end()) {
			vertexList.push_back(v);
			marginalizedProbability = C->MarginalizeOverVariable(C->y);
			for (int i = 0; i < 4; i++) {
				this->expectedCountsForVertex[v][i] += marginalizedProbability[i] * siteWeight;
			}
		}
	}
	vertexList.clear();
}

void SEM::AddToExpectedCountsForEachVariablePair() {
	SEM_vertex * u; SEM_vertex * v;
	double siteWeight = this->DNAPatternWeights[this->cliqueT->site];	
	pair <SEM_vertex *, SEM_vertex*> vertexPair;
	emtr::Md countMatrixPerSite;
	for (pair<int,SEM_vertex*> idPtrPair_1 : *this->vertexMap) {
		u = idPtrPair_1.second;
		for (pair<int,SEM_vertex*> idPtrPair_2 : *this->vertexMap) {
			v = idPtrPair_2.second;			
			if (u->id < v->id) {
				if (!u->observed or !v->observed) {
					vertexPair = pair <SEM_vertex *, SEM_vertex *>(u,v);
					countMatrixPerSite = this->cliqueT->marginalizedProbabilitiesForVariablePair[vertexPair];
					for (int dna_u = 0; dna_u < 4; dna_u ++) {
						for (int dna_v = 0; dna_v < 4; dna_v ++) {
							this->expectedCountsForVertexPair[vertexPair][dna_u][dna_v] += countMatrixPerSite[dna_u][dna_v] * siteWeight;
						}
					}
				}
			}
		}
	}
}

void SEM::AddExpectedCountMatrices(map < pair <SEM_vertex * , SEM_vertex *>, emtr::Md > expectedCountsForVertexPairToAdd) {
	string u_name;
	string v_name;
	SEM_vertex * u;
	SEM_vertex * v;
	pair <SEM_vertex *, SEM_vertex *> edge;
	emtr::Md CountMatrix;
	for (pair <pair <SEM_vertex * , SEM_vertex *>, emtr::Md> mapElem: expectedCountsForVertexPairToAdd) {
		u_name = mapElem.first.first->name;
		v_name = mapElem.first.second->name;
		CountMatrix = mapElem.second;
		u = (*this->vertexMap)[this->nameToIdMap[u_name]];
		v = (*this->vertexMap)[this->nameToIdMap[v_name]];
		if (u->id < v->id) {
			this->expectedCountsForVertexPair[pair<SEM_vertex*,SEM_vertex*>(u,v)] = CountMatrix;
		} else {			
			this->expectedCountsForVertexPair[pair<SEM_vertex*,SEM_vertex*>(v,u)] = emtr::MT(CountMatrix);

		}
	}	//this->expectedCountsForVertexPair
}


emtr::Md SEM::GetExpectedCountsForVariablePair(SEM_vertex * u, SEM_vertex * v) {
	emtr::Md C_pc;
	if (u->id < v->id) {
		C_pc = this->expectedCountsForVertexPair[pair<SEM_vertex*,SEM_vertex*>(u,v)];	
	} else {		
		C_pc = emtr::MT(this->expectedCountsForVertexPair[pair<SEM_vertex*,SEM_vertex*>(v,u)]);
	}
	return (C_pc);
}

emtr::Md SEM::GetPosteriorProbabilityForVariablePair(SEM_vertex * u, SEM_vertex * v) {
	emtr::Md P;
	if (u->id < v->id) {
		P = this->posteriorProbabilityForVertexPair[pair<SEM_vertex *, SEM_vertex *>(u,v)];
	} else {		
		P = emtr::MT(this->posteriorProbabilityForVertexPair[pair<SEM_vertex *, SEM_vertex *>(v,u)]);
	}
	return (P);
}

void SEM::AddToExpectedCountsForEachEdge() {
	double siteWeight = this->DNAPatternWeights[this->cliqueT->site];
	emtr::Md countMatrixPerSite;
	pair <SEM_vertex *, SEM_vertex *> vertexPair;
	SEM_vertex * u; SEM_vertex * v;
	for (pair <SEM_vertex *, SEM_vertex *> edge : this->edgesForPreOrderTreeTraversal) {
		tie (u,v) = edge;
		if (u->id < v->id) {
			vertexPair.first = u; vertexPair.second = v;
		} else {
			vertexPair.second = u; vertexPair.first = v;
		}
		countMatrixPerSite = this->cliqueT->marginalizedProbabilitiesForVariablePair[vertexPair];
		for (int dna_u = 0; dna_u < 4; dna_u ++) {
			for (int dna_v = 0; dna_v < 4; dna_v ++) {
				this->expectedCountsForVertexPair[vertexPair][dna_u][dna_v] += countMatrixPerSite[dna_u][dna_v] * siteWeight;
			}
		}
	}
}

void SEM::AddToExpectedCounts() {
	SEM_vertex * u; SEM_vertex * v;
	double siteWeight = this->DNAPatternWeights[this->cliqueT->site];	
	// Add to counts for each unobserved vertex (C->x) where C is a clique
	array <double, 4> marginalizedProbability;
	vector <SEM_vertex *> vertexList;
	for (clique * C: this->cliqueT->cliques) {
		v = C->x;
		if(v->observed){
            throw mt_error("v should not be observed");
        }
		if (find(vertexList.begin(),vertexList.end(),v) == vertexList.end()) {
			vertexList.push_back(v);    
			marginalizedProbability = C->MarginalizeOverVariable(C->y);
			for (int i = 0; i < 4; i++) {
				this->expectedCountsForVertex[v][i] += marginalizedProbability[i] * siteWeight;
			}
		}
	}
	vertexList.clear();
	// Add to counts for each vertex pair
	pair <SEM_vertex *, SEM_vertex*> vertexPair;
	emtr::Md countMatrixPerSite;
	for (pair<int,SEM_vertex*> idPtrPair_1 : *this->vertexMap) {
		u = idPtrPair_1.second;
		for (pair<int,SEM_vertex*> idPtrPair_2 : *this->vertexMap) {
			v = idPtrPair_2.second;			
			if (u->id < v->id) {
				if (!u->observed or !v->observed) {
					vertexPair = pair <SEM_vertex *, SEM_vertex *>(u,v);
					countMatrixPerSite = this->cliqueT->marginalizedProbabilitiesForVariablePair[vertexPair];					
//					if (u->name == "l_1" or v->name == "l_1") {
//						cout << "Count matrix for " << u->name << ", " << v->name << " for site " << this->cliqueT->site << " is" << endl;
// 						cout << countMatrixPerSite << endl;
//					}
					for (int dna_u = 0; dna_u < 4; dna_u ++) {
						for (int dna_v = 0; dna_v < 4; dna_v ++) {
							this->expectedCountsForVertexPair[vertexPair][dna_u][dna_v] += countMatrixPerSite[dna_u][dna_v] * siteWeight;
						}
					}
				}
			}
		}
	}
}

emtr::Md SEM::GetObservedCounts(SEM_vertex * u, SEM_vertex * v) {	
	emtr::Md countMatrix = emtr::Md{};
	int dna_u; int dna_v;
	for (int i = 0; i < this->sequenceLength; i++) {
		dna_u = u->DNArecoded[i];
		dna_v = v->DNArecoded[i];
		countMatrix[dna_u][dna_v] += this->DNAPatternWeights[i];
	}
	return (countMatrix);
}

void SEM::ComputeExpectedCountsForTrifle(int layer) {
    // cout << "Initializing expected counts" << endl;
	// cout << "11a" << endl;
	this->InitializeExpectedCountsForEachVariable();
	// cout << "11b" << endl;
	this->InitializeExpectedCountsForEachEdge();
	// cout << "11c" << endl;
    //	this->ResetExpectedCounts();
    //	SEM_vertex * x; SEM_vertex * y;
	emtr::Md P_XY;
	//	int dna_x; int dna_y;
	bool debug = 0;
	if (debug) {
		cout << "Debug computing expected counts" << endl;
	}
	// Iterate over sites	
    // cout << "number of sites for layer " << layer << " is " << this->num_patterns_for_layer[layer] << endl;
	for (int site = 0; site < this->num_patterns_for_layer[layer]; site++) {
		// cout << "12a site" << site << " layer " << layer << " root " << this->root->name << endl;
		// cout << "computing expected counts for site " << site << endl;
		this->cliqueT->SetSite(site);		
		// cout << "12b" << endl;
		this->cliqueT->InitializePotentialAndBeliefs();	// set min root prob to 0.2 and min trans prob to 0.01
		// cout << "12c" << endl;
		this->cliqueT->CalibrateTree();
		// cout << "12d" << endl;
		this->cliqueT->SetMarginalProbabilitesForEachEdgeAsBelief(); // set min root prob to 0.2 and min trans prob to 0.01
		// cout << "12e" << endl;
		this->AddToExpectedCountsForEachVariable();
		// cout << "12f" << endl;
		this->AddToExpectedCountsForEachEdge();		
		// cout << "12g" << endl;
	}
	// cout << "11d" << endl;
}

void SEM::ComputeExpectedCounts() {
    // cout << "Initializing expected counts" << endl;
	// cout << "11a" << endl;
	this->InitializeExpectedCountsForEachVariable();
	// cout << "11b" << endl;
	this->InitializeExpectedCountsForEachEdge();
	// cout << "11c" << endl;
    //	this->ResetExpectedCounts();
    //	SEM_vertex * x; SEM_vertex * y; 
	emtr::Md P_XY;
//	int dna_x; int dna_y;
	bool debug = 0;
	if (debug) {
		cout << "Debug computing expected counts" << endl;
	}
// Iterate over sites
	// parallelize here if needed
	for (int site = 0; site < this->num_dna_patterns; site++) {
		// cout << "12a" << endl;
		// cout << "computing expected counts for site " << site << endl;
		this->cliqueT->SetSite(site);		
		// cout << "12b" << endl;
		this->cliqueT->InitializePotentialAndBeliefs();
		// cout << "12c" << endl;
		this->cliqueT->CalibrateTree();
		// cout << "12d" << endl;
		this->cliqueT->SetMarginalProbabilitesForEachEdgeAsBelief();
		// cout << "12e" << endl;
		this->AddToExpectedCountsForEachVariable();
		// cout << "12f" << endl;
		this->AddToExpectedCountsForEachEdge();		
		// cout << "12g" << endl;
	}
	// cout << "11d" << endl;
}

void SEM::ComputeMAPEstimateOfAncestralSequencesUsingCliques() {
	this->logLikelihood = 0;
	this->ResetAncestralSequences();
	this->ConstructCliqueTree();
	clique * rootClique = this->cliqueT->root;
	SEM_vertex * v;
	map <SEM_vertex *, int> verticesVisitedMap;
	array <double, 4> posteriorProbability;
	int maxProbState;
	double maxProb;
	for (int site = 0; site < this->num_dna_patterns; site++) {		
		this->cliqueT->SetSite(site);		
		this->cliqueT->InitializePotentialAndBeliefs();		
		this->cliqueT->CalibrateTree();
		this->logLikelihood += rootClique->logScalingFactorForClique * this->DNAPatternWeights[site];
//		logLikelihood_c0 + = C_1->logScalingFactorForClique * this->sitePatternWeights[site];
//		for (int i = 0; i < 4; i ++) {
//			for (int j = 0; j < 4; j ++) {
//				
//			}
//		}
		verticesVisitedMap.clear();
		for (clique * C: this->cliqueT->cliques) {
			v = C->x;
			if (verticesVisitedMap.find(v) == verticesVisitedMap.end()) {
				posteriorProbability = C->MarginalizeOverVariable(C->y);
				maxProb = -1; maxProbState = -1;
				for (int i = 0; i < 4; i ++) {
					if (posteriorProbability[i] > maxProb) {
						maxProb = posteriorProbability[i];
						maxProbState = i;
					}
				}
				if(maxProbState == -1) {
                    throw mt_error("Check prob assignment");
                }
				v->DNArecoded.push_back(maxProbState);
				verticesVisitedMap.insert(make_pair(v,v->id));
			}
		}
	}	
}


void SEM::ComputeMarginalProbabilitiesUsingExpectedCounts() {	
	SEM_vertex * v;
	double sum;
	// Compute posterior probability for vertex
	this->posteriorProbabilityForVertex.clear();
	array <double, 4> P_X;	
	for (pair <SEM_vertex *, array <double, 4>> vertexAndCountArray: this->expectedCountsForVertex) {
		v = vertexAndCountArray.first;
		P_X = vertexAndCountArray.second;
		sum = 0;
		for (int i = 0; i < 4; i++) {
			sum += P_X[i];
		}
		if (sum > 0) {
			for (int i = 0; i < 4; i++) {
				P_X[i] /= sum;
			}
		}
		for (int i = 0; i < 4; i++) { 
			if(isnan(P_X[i])) {
				cout << "P_X[i] is nan" << endl;
				throw mt_error("P_X[i] is nan");
			}
		}
		this->posteriorProbabilityForVertex.insert(pair<SEM_vertex * , array <double, 4>>(v,P_X));
	}
	// Compute posterior probability for vertex pair
	this->posteriorProbabilityForVertexPair.clear();
	emtr::Md P_XY;
	pair <SEM_vertex *, SEM_vertex *> vertexPair;
	for (pair <pair<SEM_vertex *, SEM_vertex *>, emtr::Md> vertexPairAndCountMatrix: this->expectedCountsForVertexPair) {
		vertexPair = vertexPairAndCountMatrix.first;
		P_XY = vertexPairAndCountMatrix.second;
		sum = 0;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				sum += P_XY[i][j];
			}
		}
		if (sum > 0) {
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					P_XY[i][j] /= sum;
				}
			}
		}

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				if(isnan(P_XY[i][j])) {
					cout << "P_XY" << i << j << "is nan" << endl;
					throw mt_error("P_XY[i][j] is nan");
				}
			}
		}
		this->posteriorProbabilityForVertexPair.insert(pair<pair<SEM_vertex *, SEM_vertex *>,emtr::Md>(vertexPair,P_XY));
	}
}

void SEM::ConstructCliqueTree() {
	this->cliqueT->rootSet = 0;
	SEM_vertex * root;
	for (clique * C : this->cliqueT->cliques) {
		delete C;
	}
	this->cliqueT->cliques.clear();
	for (pair <SEM_vertex *, SEM_vertex *> edge : this->edgesForPreOrderTreeTraversal) {
		clique * C = new clique(edge.first, edge.second);		
		this->cliqueT->AddClique(C);
		if (C->x->parent == C->x and !this->cliqueT->rootSet) {
			this->cliqueT->root = C;
			this->cliqueT->rootSet = 1;
			root = C->x;
		}
	}
	clique * C_i; clique * C_j;
	// Iterate over clique pairs and identify cliques
	// that have one vertex in common
	for (unsigned int i = 0; i < this->cliqueT->cliques.size(); i ++) {
		C_i = this->cliqueT->cliques[i];
		// Set Ci as the root clique if Ci.x is the root vertex
		for (unsigned int j = i+1; j < this->cliqueT->cliques.size(); j ++) {
			C_j = this->cliqueT->cliques[j];
			// Add edge Ci -> Cj if Ci.y = Cj.x;
			if (C_i->y == C_j->x) {
//				cout << "Case 1" << endl;
//				cout << "C_i.x, C_i.y is " << C_i->x->id << ", " << C_i->y->id << endl;
//				cout << "C_j.x, C_j.y is " << C_j->x->id << ", " << C_j->y->id << endl;
				this->cliqueT->AddEdge(C_i, C_j);
				// Add edge Cj -> Ci if Cj.y = Ci.x;
			} else if (C_j->y == C_i->x) {
//				cout << "Case 2" << endl;
//				cout << "C_i.x, C_i.y is " << C_i->x->id << ", " << C_i->y->id << endl;
//				cout << "C_j.x, C_j.y is " << C_j->x->id << ", " << C_j->y->id << endl;
				this->cliqueT->AddEdge(C_j, C_i);
				// If Ci->x = Cj->x 
				// add edge Ci -> Cj
			} else if (C_i->x == C_j->x and C_i->parent == C_i) {
//				cout << "Case 3" << endl;
//				cout << "C_i.x, C_i.y is " << C_i->x->id << ", " << C_i->y->id << endl;
				this->cliqueT->AddEdge(C_i, C_j);
				// Check to see that Ci is the root clique				
				if (this->cliqueT->root != C_i) {
					cout << "Check root of clique tree" << endl;
                    throw mt_error("Check root of clique tree");
				}
			}
			// Note that Cj can never be the root clique
			// because Ci is visited before Cj
		}
	}	
	this->cliqueT->SetLeaves();
	this->cliqueT->SetEdgesForTreeTraversalOperations();
	for (clique * C: this->cliqueT->cliques) {
		C->root_variable = root;
	}
}

void SEM::ResetExpectedCounts() {
	SEM_vertex* u; SEM_vertex* v; 
	// Reset counts for each unobserved vertex
	for (pair <int, SEM_vertex *> idPtrPair : * this->vertexMap) {
		v = idPtrPair.second;
		if (!v->observed) {
			for (int i = 0; i < 4; i++) {
				this->expectedCountsForVertex[v][i] = 0;
			}
		}
	}
	// Reset counts for each vertex pair such that at least one vertex is not observed
	for (pair <int, SEM_vertex *> idPtrPair_1 : * this->vertexMap) {
		u = idPtrPair_1.second;
		for (pair<int,SEM_vertex *> idPtrPair_2 : * this->vertexMap) {
			v = idPtrPair_2.second;
			if (!u->observed or !v->observed) {
				if (u->id < v->id) {
					this->expectedCountsForVertexPair[pair <SEM_vertex *, SEM_vertex *>(u,v)] = emtr::Md{};
				}	
			}			
		}
	}
}

void SEM::InitializeExpectedCountsForEachVariable() {
	SEM_vertex * v;
	// Initialize expected counts for each vertex
	this->expectedCountsForVertex.clear();
	array <double, 4> observedCounts;
	for (pair<int, SEM_vertex *> idPtrPair : * this->vertexMap) {
		v = idPtrPair.second;
		for (int i = 0; i < 4; i++) {
			observedCounts[i] = 0;
		}
		if (v->observed) {			
			observedCounts = this->GetObservedCountsForVariable(v);
		}
		this->expectedCountsForVertex.insert(pair<SEM_vertex *, array<double,4>>(v,observedCounts));
	}	
}

void SEM::InitializeExpectedCountsForEachVariablePair() {
	SEM_vertex * u; SEM_vertex * v;
	// Initialize expected counts for each vertex pair
	this->expectedCountsForVertexPair.clear();	
	emtr::Md countMatrix;	
	int dna_u;
	int dna_v;
	for (pair<int,SEM_vertex *> idPtrPair_1 : * this->vertexMap) {
		u = idPtrPair_1.second;
		for (pair<int,SEM_vertex *> idPtrPair_2 : * this->vertexMap) {
			v = idPtrPair_2.second;
			if (u->id < v->id) {
				countMatrix = emtr::Md{};			
				if (u->observed and v->observed) {
					for (int site = 0; site < this->num_dna_patterns; site++) {
						dna_u = u->DNArecoded[site];
						dna_v = v->DNArecoded[site];
						if (dna_u < 4 && dna_v < 4) { // FIX_AMB
							countMatrix[dna_u][dna_v] += this->DNAPatternWeights[site];
						}						
					}
				}
				this->expectedCountsForVertexPair.insert(make_pair(pair <SEM_vertex *, SEM_vertex *>(u,v), countMatrix));
			}
		}
	}
}


void SEM::InitializeExpectedCountsForEachEdge() {
	// Initialize expected counts for each vertex pair
	SEM_vertex * u; SEM_vertex * v;
	pair <SEM_vertex *, SEM_vertex *> vertexPair;
	this->expectedCountsForVertexPair.clear();
	emtr::Md countMatrix;
	for (pair <SEM_vertex *, SEM_vertex *> edge : this->edgesForPreOrderTreeTraversal) {
		countMatrix = emtr::Md{};
		tie (u,v) = edge;
		if (u->id < v->id) {
			vertexPair.first = u; vertexPair.second = v;
		} else {
			vertexPair.first = v; vertexPair.second = u;
		}
		this->expectedCountsForVertexPair.insert(make_pair(vertexPair, countMatrix));
	}
}

void SEM::InitializeExpectedCounts() {
	SEM_vertex * u; SEM_vertex * v;
	// Initialize expected counts for each vertex
	this->expectedCountsForVertex.clear();
	array <double, 4> observedCounts;
	for (pair<int, SEM_vertex *> idPtrPair : * this->vertexMap) {
		v = idPtrPair.second;
		for (int i = 0; i < 4; i++) {
			observedCounts[i] = 0;
		}
		if (v->observed) {			
			observedCounts = this->GetObservedCountsForVariable(v);
		}
		this->expectedCountsForVertex.insert(pair<SEM_vertex *, array<double,4>>(v,observedCounts));
	}
	
	// Initialize expected counts for each vertex pair
	this->expectedCountsForVertexPair.clear();	
	emtr::Md countMatrix;	
	int dna_u;
	int dna_v;
	for (pair<int,SEM_vertex *> idPtrPair_1 : * this->vertexMap) {
		u = idPtrPair_1.second;
		for (pair<int,SEM_vertex *> idPtrPair_2 : * this->vertexMap) {
			v = idPtrPair_2.second;
			if (u->id < v->id) {
				countMatrix = emtr::Md{};			
				if (u->observed and v->observed) {
					for (int site = 0; site < this->num_dna_patterns; site++) {
						dna_u = u->DNArecoded[site];
						dna_v = v->DNArecoded[site];
						countMatrix[dna_u][dna_v] += this->DNAPatternWeights[site];
					}
				}
				this->expectedCountsForVertexPair.insert(make_pair(pair <SEM_vertex *, SEM_vertex *>(u,v), countMatrix));
			}
		}
	}
}

void SEM::ResetPointerToRoot() {
	//	Make sure that the pointer this->root stores the location
	//  of the vertex with in degree 0
	for (pair<int,SEM_vertex *> idPtrPair : *this->vertexMap) {
		if (idPtrPair.second->inDegree == 0){
			this->root = idPtrPair.second;
		}
	}
}

void SEM::AddArc(SEM_vertex * from, SEM_vertex * to) {
	to->AddParent(from);
	from->AddChild(to);
}

void SEM::RemoveArc(SEM_vertex * from, SEM_vertex * to) {	
	to->parent = to;
	to->inDegree -= 1;
	from->outDegree -= 1;	
	int ind = find(from->children.begin(),from->children.end(),to) - from->children.begin();
	from->children.erase(from->children.begin()+ind);	
}

pair<bool,SEM_vertex *> SEM::CheckAndRetrieveHiddenVertexWithOutDegreeOneAndInDegreeOne() {
	bool containsVertex = 0;
	SEM_vertex* vPtrToReturn = (*this->vertexMap)[this->vertexMap->size()-1];
	for (pair<int, SEM_vertex *> idPtrPair: *this->vertexMap) {
		if (idPtrPair.second->outDegree == 1 and idPtrPair.second->inDegree == 1) {
			if (idPtrPair.second->id > this->numberOfObservedVertices-1) {
				containsVertex = 1;
				vPtrToReturn = idPtrPair.second;		
				break;
			}
		}		
	}
	return (make_pair(containsVertex,vPtrToReturn));
}

pair <bool,SEM_vertex *> SEM::CheckAndRetrieveHiddenVertexWithOutDegreeOneAndInDegreeZero() {
	bool containsVertex = 0;
	SEM_vertex* vPtrToReturn = (*this->vertexMap)[this->vertexMap->size()-1];
	for (pair<int, SEM_vertex*> idPtrPair: *this->vertexMap) {
		if (idPtrPair.second->outDegree == 1 and idPtrPair.second->inDegree == 0) {
			if (idPtrPair.second->id > this->numberOfObservedVertices-1) {
				containsVertex = 1;
				vPtrToReturn = idPtrPair.second;		
				break;
			}
		}		
	}
	return (make_pair(containsVertex,vPtrToReturn));
}

pair <bool, SEM_vertex*> SEM::CheckAndRetrieveSingletonHiddenVertex() {
	bool containsVertex = 0;
	SEM_vertex* vPtrToReturn = (*this->vertexMap)[this->vertexMap->size()-1];
	for (pair<int, SEM_vertex*> idPtrPair: *this->vertexMap) {
		if (!idPtrPair.second->observed and idPtrPair.second->outDegree == 0 and idPtrPair.second->inDegree == 0) {			
			containsVertex = 1;
			vPtrToReturn = idPtrPair.second;		
			break;
		}
	}
	return (make_pair(containsVertex,vPtrToReturn));
}


pair <bool,SEM_vertex*> SEM::CheckAndRetrieveHiddenVertexWithOutDegreeZeroAndInDegreeOne() {
	bool containsVertex = 0;
	SEM_vertex* vPtrToReturn = (*this->vertexMap)[this->vertexMap->size()-1];
	for (pair <int, SEM_vertex*> idPtrPair: *this->vertexMap) {
		if (!idPtrPair.second->observed and idPtrPair.second->outDegree == 0 and idPtrPair.second->inDegree == 1) {
			containsVertex = 1;
			vPtrToReturn = idPtrPair.second;
			break;
		}
	}
	return (make_pair(containsVertex,vPtrToReturn));
}

pair <bool, SEM_vertex *> SEM::CheckAndRetrieveHiddenVertexWithOutDegreeGreaterThanTwo() {
	bool containsVertex = 0;
	SEM_vertex* vPtrToReturn = (*this->vertexMap)[this->vertexMap->size()-1];
	for (pair <int, SEM_vertex *> idPtrPair: *this->vertexMap) {
		if (!idPtrPair.second->observed and idPtrPair.second->outDegree > 2) {
			if (idPtrPair.second->id > this->numberOfObservedVertices-1) {
				containsVertex = 1;
				vPtrToReturn = idPtrPair.second;		
				break;
			}
		}		
	}
	return (make_pair(containsVertex,vPtrToReturn));
}

pair <bool,SEM_vertex *> SEM::CheckAndRetrieveObservedVertexThatIsTheRoot() {
	bool containsVertex = 0;
	SEM_vertex * vPtrToReturn = (*this->vertexMap)[this->vertexMap->size()-1];
	for (pair <int, SEM_vertex *> idPtrPair: *this->vertexMap) {
		if (idPtrPair.second->observed and idPtrPair.second->outDegree > 0) {
			containsVertex = 1;
			vPtrToReturn = idPtrPair.second;
			break;
		}
	}
	return (make_pair(containsVertex,vPtrToReturn));
}

pair <bool,SEM_vertex *> SEM::CheckAndRetrieveObservedVertexThatIsNotALeafAndIsNotTheRoot() {
	bool containsVertex = 0;
	SEM_vertex * vPtrToReturn = (*this->vertexMap)[this->vertexMap->size()-1];
	for (pair <int, SEM_vertex *> idPtrPair: *this->vertexMap) {
		if (idPtrPair.second->observed and idPtrPair.second->outDegree > 0) {
			containsVertex = 1;
			vPtrToReturn = idPtrPair.second;
			break;
		}
	}
	return (make_pair(containsVertex,vPtrToReturn));
}

void SEM::StoreBestProbability() {
	this->StoreRootAndRootProbability();
	this->StoreTransitionMatrices();
}

void SEM::StoreRootAndRootProbability() {
	this->root_stored = this->root;
	this->rootProbability_stored = this->rootProbability;
}

void SEM::StoreTransitionMatrices() {
	for (pair <int,SEM_vertex*> idPtrPair : * this->vertexMap) {
		idPtrPair.second->transitionMatrix_stored = idPtrPair.second->transitionMatrix;
	}
}

void SEM::RestoreBestProbability() {
	this->RestoreRootAndRootProbability();
	this->RestoreTransitionMatrices();
	this->RootTreeAtVertex(this->root);
}

void SEM::RestoreRootAndRootProbability() {
	this->rootProbability = this->rootProbability_stored;
	this->root = this->root_stored;
	this->root->rootProbability = this->rootProbability;	
}

void SEM::RestoreTransitionMatrices() {
	for (pair<int,SEM_vertex*> idPtrPair : * this->vertexMap) {
		idPtrPair.second->transitionMatrix = idPtrPair.second->transitionMatrix_stored;
	}
}

void SEM::StoreRateMatricesAndScalingFactors() {
	this->rateMatrixPerRateCategory_stored = this->rateMatrixPerRateCategory;
	this->scalingFactorPerRateCategory_stored = this->scalingFactorPerRateCategory;
}

void SEM::RestoreRateMatricesAndScalingFactors() {
	this->rateMatrixPerRateCategory = this->rateMatrixPerRateCategory_stored;
	this->scalingFactorPerRateCategory = this->scalingFactorPerRateCategory_stored;
}


void SEM::ComputeLogLikelihoodUsingExpectedDataCompletion() {
	this->logLikelihood = 0;
	array <double, 4> S_r = this->expectedCountsForVertex[this->root];
	for (int dna_r = 0; dna_r < 4; dna_r ++) {
		if (this->rootProbability[dna_r] > 0) {
			this->logLikelihood += S_r[dna_r] * log(this->rootProbability[dna_r]);
		}	
	}
	// Contribution of edges
	SEM_vertex * p; SEM_vertex * c;
	emtr::Md S_pc;
	emtr::Md P;
	for (pair <int, SEM_vertex *> idPtrPair : *this->vertexMap) {
		c = idPtrPair.second;
		p = c->parent;
		if (p != c) {
			P = c->transitionMatrix;
			if (p->id < c->id) {
				S_pc = this->expectedCountsForVertexPair[pair<SEM_vertex*,SEM_vertex*>(p,c)];	
			} else {				
				S_pc = emtr::MT(this->expectedCountsForVertexPair[pair<SEM_vertex*,SEM_vertex*>(c,p)]);
			}
			for (int dna_p = 0; dna_p < 4; dna_p ++) {
				for (int dna_c = 0; dna_c < 4; dna_c ++) {
					if(S_pc[dna_p][dna_c] < 0) {
						cout << "Expected counts for " << p->name << "\t" << c->name << " is " << endl;						
						throw mt_error("Check counts");
					}
					if (P[dna_p][dna_c] > 0) {
						this->logLikelihood += S_pc[dna_p][dna_c] * log(P[dna_p][dna_c]);
					}
				}
			}
		}
	}
}

void SEM::ComputeAALogLikelihood_maps() {
	this->AAlogLikelihood = 0;
	constexpr size_t K = 20;
	map <SEM_vertex*,array<double,K>> conditionalLikelihoodMap;
	std::array <double,K> conditionalLikelihood;
	double partialLikelihood;
	double siteLikelihood;	
	double largestConditionalLikelihood = 0;
	double currentProb;			
	vector <SEM_vertex *> verticesToVisit;	
	SEM_vertex * p;
	SEM_vertex * c;
	Eigen::Matrix <double,20,20> P;
	double t;
	// this->num_aa_patterns = 10;
	for (int site = 0; site < this->num_aa_patterns; site++) {
    conditionalLikelihoodMap.clear();
    this->ResetLogScalingFactors();

		for (auto& edge : this->edgesForPostOrderTreeTraversal) {
			tie(p, c) = edge;
			t = this->GetEdgeLength(p,c);
			P = this->Q_D * t;
			P = P.exp();
			c->AAtransitionMatrix = P;

			p->logScalingFactors += c->logScalingFactors;

			// Initialize leaf child
			if (c->outDegree == 0) {
				if (c->AAcompressed[site] > -1) {
					conditionalLikelihood.fill(0.0);
					conditionalLikelihood[c->AAcompressed[site]] = 1.0;
				} else {
					conditionalLikelihood.fill(1.0);
				}
				conditionalLikelihoodMap.insert({c, conditionalLikelihood});
			}

			// Initialize parent p if absent
			if (conditionalLikelihoodMap.find(p) == conditionalLikelihoodMap.end()) {
				if (p->id > this->numberOfObservedVertices - 1) {
					conditionalLikelihood.fill(1.0);      // latent
				} else {
					conditionalLikelihood.fill(0.0);      // observed
					conditionalLikelihood[p->AAcompressed[site]] = 1.0;
				}
				conditionalLikelihoodMap.insert({p, conditionalLikelihood});
			}

			// DP update using safer access for child
			double largestConditionalLikelihood = 0.0;
			const auto& childCL = conditionalLikelihoodMap.at(c);

			for (int dna_p = 0; dna_p < K; ++dna_p) {
				double partialLikelihood = 0.0;
				for (int dna_c = 0; dna_c < K; ++dna_c) {
					partialLikelihood += P(dna_p,dna_c) * childCL[dna_c];
				}
				conditionalLikelihoodMap[p][dna_p] *= partialLikelihood;
				largestConditionalLikelihood = std::max(largestConditionalLikelihood, conditionalLikelihoodMap[p][dna_p]);
			}

			if (largestConditionalLikelihood > 0.0) {
				for (int dna_p = 0; dna_p < K; ++dna_p) {
					conditionalLikelihoodMap[p][dna_p] /= largestConditionalLikelihood;
				}
				p->logScalingFactors += std::log(largestConditionalLikelihood);
			} else {
				throw mt_error("Largest conditional likelihood value is zero");
			}
		}
		double siteLikelihood = 0.0;
		const auto& rootCL = conditionalLikelihoodMap.at(this->root);
		for (int dna = 0; dna < K; ++dna) {
			siteLikelihood += this->pi_D[dna] * rootCL[dna];
		}
		if (siteLikelihood <= 0.0) throw mt_error("siteLikelihood <= 0");

		this->AAlogLikelihood += (this->root->logScalingFactors + std::log(siteLikelihood)) * this->AAPatternWeights[site];
	}
}

void SEM::ComputeAALogLikelihood() {
	int num_states = 20;
	this->AAlogLikelihood = 0;		
	vector <vector<double>> conditionalLikelihoodVec;
	int num_vertices = this->vertices.size();
	int num_obs_vertices = this->numberOfObservedVertices;
	vector <double> logScalingFactorVec(num_vertices, 0.0);
	std::array <double,20> conditionalLikelihood;
	double partialLikelihood;
	double siteLikelihood;	
	double largestConditionalLikelihood = 0;
	double currentProb;			
	vector <SEM_vertex *> verticesToVisit;	
	SEM_vertex * p;
	SEM_vertex * c;	
	int root_id = this->root->id;
	double blen;
	
	for (pair<SEM_vertex *,SEM_vertex *> edge : this->edgesForPostOrderTreeTraversal) {
		tie (p, c) = edge;
		blen = this->GetEdgeLength(p,c);		
		c->AAtransitionMatrix = ExpQt_fromSpectrum(blen);
		cout << " transition matrix entry is " << c->AAtransitionMatrix(3,3) << endl;
	}

	for (int site = 0; site < this->num_aa_patterns; site++) {				
		for (int i = 0; i < num_obs_vertices; i++) {
			c = this->vertices[i];
			if (c->AAcompressed[site] > -1) {
				conditionalLikelihoodVec[i].assign(20,0);
				conditionalLikelihoodVec[i][c->AAcompressed[site]] = 1;
			} else {				
				conditionalLikelihoodVec[i].assign(20,1);
			}			
		}
		for (int i = num_obs_vertices; i < num_vertices; i++) conditionalLikelihoodVec[i].assign(20,1);

		fill(logScalingFactorVec.begin(), logScalingFactorVec.end(), 0.0);				
		int p_id, c_id;
		for (pair<SEM_vertex *,SEM_vertex *> edge : this->edgesForPostOrderTreeTraversal) {
			tie (p, c) = edge;
			p_id = p->id;
			c_id = c->id;			
			logScalingFactorVec[p_id] += logScalingFactorVec[c_id];
			largestConditionalLikelihood = 0;
			for (int aa_p = 0; aa_p < 20; aa_p++) {
				partialLikelihood = 0;
				for (int aa_c = 0; aa_c < 20; aa_c++) {
					partialLikelihood += c->AAtransitionMatrix(aa_p,aa_c)*conditionalLikelihoodVec[c_id][aa_c];
				}
				conditionalLikelihoodVec[p_id][aa_p] *= partialLikelihood;
				if (conditionalLikelihoodVec[p_id][aa_p] > largestConditionalLikelihood) {
					largestConditionalLikelihood = conditionalLikelihoodVec[p_id][aa_p];
				}
			}
			if (largestConditionalLikelihood != 0){
				for (int aa_p = 0; aa_p < 20; aa_p++) {
					conditionalLikelihoodVec[p_id][aa_p] /= largestConditionalLikelihood;
				}				
				logScalingFactorVec[p_id] += log(largestConditionalLikelihood);
			} else {
				cout << site;
				cout << " parent name " << p->name << " ";
				cout << " child name " << c->name << " ";
				for (int i = 0; i < 20; i ++) {
					for (int j = 0; j < 20; j++) {
						cout << c->AAtransitionMatrix(i,j) << endl;
					}
				}				
				cout << "Largest conditional likelihood value is zero" << endl;				
				throw mt_error("Largest conditional likelihood value is zero");
			}					
		}
		siteLikelihood = 0; 							
		for (int aa = 0; aa < 20; aa++) {
			currentProb = this->pi_D[aa]*conditionalLikelihoodVec[root_id][aa];
			siteLikelihood += currentProb;
		}		
		this->AAlogLikelihood += (logScalingFactorVec[root_id] + log(siteLikelihood)) * this->AAPatternWeights[site];				
	}
}

// Case 1: Observed vertices may have out degree > 0
// Case 2: Root may have out degree = one
// Case 3: Directed tree (rooted) with vertices with outdegree 2 or 0.
// layers 0, 1 and 2
// coarse: 0, coarse: 1, coarse: 2 
void SEM::ComputeTrifleLogLikelihood(int layer) {
	this->logLikelihood = 0;
	map <SEM_vertex*,array<double,4>> conditionalLikelihoodMap;
	std::array <double,4> conditionalLikelihood;
	double partialLikelihood;
	double siteLikelihood;	
	double largestConditionalLikelihood = 0;
	double currentProb;			
	vector <SEM_vertex *> verticesToVisit;	
	SEM_vertex * p;
	SEM_vertex * c;
	emtr::Md P;
	// cout << "14a" << endl;
	for (int site = 0; site < this->num_patterns_for_layer[layer]; site++) {
		// cout << "14b" << endl;
    	conditionalLikelihoodMap.clear();
		// cout << "14c" << endl;
    	this->ResetLogScalingFactors();
		// cout << "14d" << endl;
		for (auto& edge : this->edgesForPostOrderTreeTraversal) {
			// cout << "15a" << endl;
			tie(p, c) = edge;
			P = c->transitionMatrix;			
			p->logScalingFactors += c->logScalingFactors;

			// Initialize leaf child
			if (c->outDegree == 0) {
				if (c->DNAcompressed[site] > -1) {
					conditionalLikelihood.fill(0.0);
					conditionalLikelihood[c->DNAcompressed[site]] = 1.0;
				} else {
					conditionalLikelihood.fill(1.0);
				}
				conditionalLikelihoodMap.insert({c, conditionalLikelihood});
			}
			// cout << "15b" << endl;

			// Initialize parent p if absent
			if (conditionalLikelihoodMap.find(p) == conditionalLikelihoodMap.end()) {
				if (p->id > this->numberOfObservedVertices - 1) {
					conditionalLikelihood.fill(1.0);      // latent
				} else {
					conditionalLikelihood.fill(0.0);      // observed
					conditionalLikelihood[p->DNAcompressed[site]] = 1.0;
				}
				conditionalLikelihoodMap.insert({p, conditionalLikelihood});
			}
			// cout << "15c" << endl;
			
			double largestConditionalLikelihood = 0.0;
			const auto& childCL = conditionalLikelihoodMap.at(c);
			// cout << "15d" << endl;
			for (int dna_p = 0; dna_p < 4; ++dna_p) {
				double partialLikelihood = 0.0;
				for (int dna_c = 0; dna_c < 4; ++dna_c) {
					partialLikelihood += P[dna_p][dna_c] * childCL[dna_c];
				}
				conditionalLikelihoodMap[p][dna_p] *= partialLikelihood;
				largestConditionalLikelihood = std::max(largestConditionalLikelihood, conditionalLikelihoodMap[p][dna_p]);
			}
			// cout << "15e" << endl;						
			if (largestConditionalLikelihood > 0.0) {
				for (int dna_p = 0; dna_p < 4; ++dna_p) {
					conditionalLikelihoodMap[p][dna_p] /= largestConditionalLikelihood;
				}
				p->logScalingFactors += log(largestConditionalLikelihood);
			} else {
				cout << "conditional likelihood is zero " << largestConditionalLikelihood << endl;
				cout << "p dna " << p->DNAcompressed[site] << endl;
				cout << "c dna " << c->DNAcompressed[site] << endl;
				for (int i = 0; i < 4; i ++) {
					for (int j = 0; j < 4; j ++) {
						cout << "c->transitionMatrix[" << i << "," << j << "] = " << c->transitionMatrix[i][j] << endl;
					}
				}
				throw mt_error("Largest conditional likelihood value is zero");
			}
			// cout << "15f" << endl;
		}
		// cout << "14e" << endl;
		double siteLikelihood = 0.0;
		const auto& rootCL = conditionalLikelihoodMap.at(this->root);
		for (int dna = 0; dna < 4; ++dna) {
			siteLikelihood += this->rootProbability[dna] * rootCL[dna];
		}
		// cout << "14f" << endl;
		if (siteLikelihood <= 0.0) throw mt_error("siteLikelihood <= 0");
		// cout << "14g" << endl;
		this->logLikelihood += (this->root->logScalingFactors + log(siteLikelihood)) * this->DNAPatternWeights[site];
		// cout << "14h" << endl;
	}
}

std::string SEM::emtrifle_to_json(const EMTrifle_struct& em) const {
    std::ostringstream os;

    // Top-level init values from layer 0 (when usable)
    const bool has_ll_init0 = (!em.ll_initial_trifle.empty() &&
                               std::isfinite(em.ll_initial_trifle[0]));
    const bool has_rpi0 = (!em.root_prob_initial_trifle.empty() &&
                           is_vec4_finite(em.root_prob_initial_trifle[0]));
    const bool has_tpi0 = (!em.trans_prob_initial_trifle.empty() &&
                           !em.trans_prob_initial_trifle[0].empty());

    os << '{'
       << "\"method\":" << estr(em.method) << ','
       << "\"rep\":"    << em.rep << ','
       << "\"root\":"   << estr(em.root) << ',';

    if (has_ll_init0) {
        os << "\"ll_init\":" << num3(em.ll_initial_trifle[0]) << ',';
    }

    os << "\"root_prob_init\":";
    if (has_rpi0) os << pvec4(em.root_prob_initial_trifle[0]); else os << "[]";
    os << ',';

    os << "\"trans_prob_init\":";
    if (has_tpi0) os << hmap_mat4(em.trans_prob_initial_trifle[0]); else os << "{}";
    os << ',';

    // Layers array
    os << "\"layers\":[";
    bool firstLayer = true;
    for (int layer = 0; layer < 3; ++layer) {
        const size_t L = static_cast<size_t>(layer);

        // Availability checks
        const bool has_iter        = (L < em.iter_trifle.size());
        const bool has_ll_initial  = (L < em.ll_initial_trifle.size() &&
                                      std::isfinite(em.ll_initial_trifle[layer]));
        const bool has_ll_final    = (L < em.ll_final_trifle.size() &&
                                      std::isfinite(em.ll_final_trifle[layer]));
        const bool has_rpf         = (L < em.root_prob_final_trifle.size() &&
                                      is_vec4_finite(em.root_prob_final_trifle[layer]));
        const bool has_tpf         = (L < em.trans_prob_final_trifle.size() &&
                                      !em.trans_prob_final_trifle[layer].empty());
        const bool has_rpi_layer   = (L < em.root_prob_initial_trifle.size() &&
                                      is_vec4_finite(em.root_prob_initial_trifle[layer]));
        const bool has_tpi_layer   = (L < em.trans_prob_initial_trifle.size() &&
                                      !em.trans_prob_initial_trifle[layer].empty());
        const bool has_ecd         = (L < em.ecd_ll_per_iter_for_trifle.size());

        // If absolutely nothing for this layer, skip
        if (!(has_iter || has_ll_initial || has_ll_final ||
              has_rpf || has_tpf || has_rpi_layer || has_tpi_layer || has_ecd)) continue;

        if (!firstLayer) os << ',';
        firstLayer = false;

        os << '{'
           << "\"layer\":" << layer;

        if (has_iter)       os << ",\"iter\":"       << em.iter_trifle[layer];
        if (has_ll_initial) os << ",\"ll_initial\":" << num3(em.ll_initial_trifle[layer]);
        if (has_ll_final)   os << ",\"ll_final\":"   << num3(em.ll_final_trifle[layer]);

        // Per-layer initial + final (containers never null)
        os << ",\"root_prob_initial\":";
        if (has_rpi_layer) os << pvec4(em.root_prob_initial_trifle[layer]); else os << "[]";

        os << ",\"trans_prob_initial\":";
        if (has_tpi_layer) os << hmap_mat4(em.trans_prob_initial_trifle[layer]); else os << "{}";

        os << ",\"root_prob_final\":";
        if (has_rpf) os << pvec4(em.root_prob_final_trifle[layer]); else os << "[]";

        os << ",\"trans_prob_final\":";
        if (has_tpf) os << hmap_mat4(em.trans_prob_final_trifle[layer]); else os << "{}";

        os << ",\"ecd_ll_per_iter\":";
        if (has_ecd) os << dump_iter_ll_obj(em.ecd_ll_per_iter_for_trifle[layer]); else os << "{}";

        // raw_json slice (optional mirror)
        os << ",\"raw_json\":{"
           << "\"layer\":" << layer;
        if (has_iter)       os << ",\"iter\":"       << em.iter_trifle[layer];
        if (has_ll_initial) os << ",\"ll_initial\":" << num3(em.ll_initial_trifle[layer]);
        if (has_ll_final)   os << ",\"ll_final\":"   << num3(em.ll_final_trifle[layer]);

        os << ",\"root_prob_initial\":";
        if (has_rpi_layer) os << pvec4(em.root_prob_initial_trifle[layer]); else os << "[]";

        os << ",\"trans_prob_initial\":";
        if (has_tpi_layer) os << hmap_mat4(em.trans_prob_initial_trifle[layer]); else os << "{}";

        os << ",\"root_prob_final\":";
        if (has_rpf) os << pvec4(em.root_prob_final_trifle[layer]); else os << "[]";

        os << ",\"trans_prob_final\":";
        if (has_tpf) os << hmap_mat4(em.trans_prob_final_trifle[layer]); else os << "{}";

        os << ",\"ecd_ll_per_iter\":";
        if (has_ecd) os << dump_iter_ll_obj(em.ecd_ll_per_iter_for_trifle[layer]); else os << "{}";

        os << '}'; // end raw_json
        os << '}'; // end layer obj
    }
    os << ']'; // end layers

    os << '}'; // end root object
    return os.str();
}



// Case 1: Observed vertices may have out degree > 0
// Case 2: Root may have out degree = one
// Case 3: Directed tree (rooted) with vertices with outdegree 2 or 0.
void SEM::ComputeLogLikelihood() {
	this->logLikelihood = 0;
	map <SEM_vertex*,array<double,4>> conditionalLikelihoodMap;
	std::array <double,4> conditionalLikelihood;
	double partialLikelihood;
	double siteLikelihood;	
	double largestConditionalLikelihood = 0;
	double currentProb;			
	vector <SEM_vertex *> verticesToVisit;	
	SEM_vertex * p;
	SEM_vertex * c;
	emtr::Md P;
	for (int site = 0; site < this->num_dna_patterns; site++) {
    conditionalLikelihoodMap.clear();
    this->ResetLogScalingFactors();
		for (auto& edge : this->edgesForPostOrderTreeTraversal) {
			tie(p, c) = edge;
			P = c->transitionMatrix;

			p->logScalingFactors += c->logScalingFactors;

			// Initialize leaf child
			if (c->outDegree == 0) {
				if (c->DNAcompressed[site] > -1) {
					conditionalLikelihood.fill(0.0);
					conditionalLikelihood[c->DNAcompressed[site]] = 1.0;
				} else {
					conditionalLikelihood.fill(1.0);
				}
				conditionalLikelihoodMap.insert({c, conditionalLikelihood});
			}

			// Initialize parent p if absent
			if (conditionalLikelihoodMap.find(p) == conditionalLikelihoodMap.end()) {
				if (p->id > this->numberOfObservedVertices - 1) {
					conditionalLikelihood.fill(1.0);      // latent
				} else {
					conditionalLikelihood.fill(0.0);      // observed
					conditionalLikelihood[p->DNAcompressed[site]] = 1.0;
				}
				conditionalLikelihoodMap.insert({p, conditionalLikelihood});
			}

			// DP update using safer access for child
			double largestConditionalLikelihood = 0.0;
			const auto& childCL = conditionalLikelihoodMap.at(c);

			for (int dna_p = 0; dna_p < 4; ++dna_p) {
				double partialLikelihood = 0.0;
				for (int dna_c = 0; dna_c < 4; ++dna_c) {
					partialLikelihood += P[dna_p][dna_c] * childCL[dna_c];
				}
				conditionalLikelihoodMap[p][dna_p] *= partialLikelihood;
				largestConditionalLikelihood = std::max(largestConditionalLikelihood, conditionalLikelihoodMap[p][dna_p]);
			}

			if (largestConditionalLikelihood > 0.0) {
				for (int dna_p = 0; dna_p < 4; ++dna_p) {
					conditionalLikelihoodMap[p][dna_p] /= largestConditionalLikelihood;
				}
				p->logScalingFactors += std::log(largestConditionalLikelihood);
			} else {
				throw mt_error("Largest conditional likelihood value is zero");
			}
		}
		double siteLikelihood = 0.0;
		const auto& rootCL = conditionalLikelihoodMap.at(this->root);
		for (int dna = 0; dna < 4; ++dna) {
			siteLikelihood += this->rootProbability[dna] * rootCL[dna];
		}
		if (siteLikelihood <= 0.0) throw mt_error("siteLikelihood <= 0");

		this->logLikelihood += (this->root->logScalingFactors + std::log(siteLikelihood)) * this->DNAPatternWeights[site];
	}
}

vector<int> SEM::DecompressSequence(vector<int> compressedSequence, vector<vector<int>> sitePatternRepeats){
	int totalSequenceLength = 0;
	for (vector<int> sitePatternRepeat: sitePatternRepeats){
		totalSequenceLength += int(sitePatternRepeat.size());
	}
	vector <int> decompressedSequence;
	for (int v_ind = 0; v_ind < totalSequenceLength; v_ind++){
		decompressedSequence.push_back(char(0));
	}
	int dnaToAdd;
	for (int sitePatternIndex = 0; sitePatternIndex < int(compressedSequence.size()); sitePatternIndex++){
		dnaToAdd = (compressedSequence)[sitePatternIndex];
		for (int pos: (sitePatternRepeats)[sitePatternIndex]){
			decompressedSequence[pos] = dnaToAdd;
		}
	}
	return (decompressedSequence);
}

string SEM::EncodeAsDNA(vector<int> sequence){
	string allDNA = "AGTC";
	string dnaSequence = "";
	for (int s : sequence){
		dnaSequence += allDNA[s];
	}
	return dnaSequence;
}

void SEM::Optimize_AA_RateMatrix_rooted_at_each_internal_vertex(int num_repetitions) {
	int n = this->numberOfObservedVertices;
	int num_vertices = this->vertexMap->size();		
	vector<int> vertex_indices_to_visit;
	vertex_indices_to_visit.reserve(std::max(0, num_vertices - n));
	
	for (int v_ind = n; v_ind < num_vertices; ++v_ind) vertex_indices_to_visit.push_back(v_ind);

	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 rng(seed);

    // Shuffle the vector
	// cout << "randomizing the order in which nodes will be visited" << endl;
    shuffle(vertex_indices_to_visit.begin(), vertex_indices_to_visit.end(), rng);
    int num_vertices_to_visit = vertex_indices_to_visit.size();
	int v_ind;
	SEM_vertex * v;
	double aa_ll_init;
	double aa_ll_final;
	int max_iter = 50;	
	for (int v_i = 0; v_i < num_vertices_to_visit; v_i++) {		
		v_ind = vertex_indices_to_visit[v_i];
		v = (*this->vertexMap)[v_ind];
		// cout << "node " << v_i+1 << " of " << num_vertices_to_visit << ":" << v->name << endl;
		this->RootTreeAtVertex(v);		
		double aa_ll_previous;
		double aa_ll_current;
		double conv_threshold = 0.01;
		int iter;
		for (int rep = 0; rep < num_repetitions; rep++) {
			this->SetDayhoffExchangeabilityAndRootProbabilityParameters();
			this->SetRateMatrixParameters();
			this->PrecomputeSpectrum();
			this->ComputeNRDistanceMatrix();
			this->AssignLSBranchLengthsToTree();
			this->ComputeAALogLikelihood();		
			this->AA_struct = AA_opt_struct{};
			this->AA_struct.root = v->name;			
			this->AA_struct.aa_ll_init = this->AAlogLikelihood;
			this->AA_struct.exchangeability_matrix_init = this->S_D;
			this->AA_struct.root_prob_init = this->pi_D;
			this->AA_struct.rep = rep;
			aa_ll_previous = this->AAlogLikelihood;
			// cout << "AA log likelihood for Dayhoff is " << this->AA_struct.aa_ll_init << endl;				
			for (iter = 0; iter < max_iter; iter++) {
				this->OptimizeRateMatrixNelderMead();
				this->ComputeNRDistanceMatrix();
				this->AssignLSBranchLengthsToTree();				
				this->ComputeAALogLikelihood();
				aa_ll_current = this->AAlogLikelihood;
				this->AA_struct.aa_ll_per_iter[iter+1] = aa_ll_current;
				double diff_ll = aa_ll_current - aa_ll_previous;
				// cout << "current aa_ll " << aa_ll_current << endl;
				// cout << "aa_ll_current - aa_ll_previous for iter " << iter + 1 << " is " << aa_ll_current - aa_ll_previous << endl;
				if (fabs(diff_ll) < conv_threshold) {
					break;
				} else {
					// cout << "AA log likelihood for iter " << iter+1 << " is " << setprecision(8) << aa_ll_current << endl;
					// cout << "Change in AA log likelihood for iter " << iter+1 << " is " << setprecision(8) << (diff_ll) << endl;
					aa_ll_previous = aa_ll_current;
				}				
			}
			this->AA_struct.iter = iter+1;
			this->AA_struct.aa_ll_init = this->AAlogLikelihood;
			this->AA_struct.exchangeability_matrix_final = this->S_D;
			this->AA_struct.root_prob_final = this->pi_D;
			this->AA_opt_runs.push_back(this->AA_struct);
		}
	}
}

void SEM::EM_DNA_rooted_at_each_internal_vertex_started_with_dirichlet_store_results(int num_repetitions) {
	int n = this->numberOfObservedVertices;
	int num_vertices = this->vertexMap->size();	
		
	this->max_log_likelihood_diri = -1 * pow(10,5);
	
	vector<int> vertex_indices_to_visit;
	vertex_indices_to_visit.reserve(std::max(0, num_vertices - n));
	
	for (int v_ind = n; v_ind < num_vertices; ++v_ind) vertex_indices_to_visit.push_back(v_ind);

	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 rng(seed);

    // Shuffle the vector
	// cout << "randomizing the order in which nodes will be visited" << endl;
    shuffle(vertex_indices_to_visit.begin(), vertex_indices_to_visit.end(), rng);

    int num_vertices_to_visit = vertex_indices_to_visit.size();
	for (int v_i = 0; v_i < num_vertices_to_visit; v_i++) {
        const int v_ind = vertex_indices_to_visit[v_i];
		SEM_vertex * v = (*this->vertexMap)[v_ind];
		// cout << "node " << v_i+1 << " of " << num_vertices_to_visit << ":" << v->name << endl;
		if(v->degree != 3) throw mt_error("Expect internal nodes to have degree three");
		
		for (int rep = 0; rep < num_repetitions; rep++) {
			this->EM_current = EM_struct{};
			this->EM_current.method = "dirichlet";
			this->EM_current.root = v->name;
			this->EM_current.rep = rep + 1;			
			auto tup = this->EM_started_with_dirichlet_rooted_at(v);
			EM_struct EM_diri{this->EM_current};
			this->EM_DNA_runs_diri.push_back(EM_diri);
			const int    iter                    = std::get<0>(tup);
            const double logLikelihood_diri      = std::get<1>(tup);
            const double loglikelihood_ecd_first = std::get<2>(tup);
            const double loglikelihood_ecd_final = std::get<3>(tup);
            const double logLikelihood_final     = std::get<4>(tup);

			emtr::push_result(
                this->EMTR_results,                 // vector<tuple<string,string,int,int,double,double,double,double>>
                "dirichlet",                        // init_method
                v->name,                            // root
                rep + 1,                            // repetition (1-based)
                iter,
                logLikelihood_diri,                 // ll_initial
                loglikelihood_ecd_first,
                loglikelihood_ecd_final,
                logLikelihood_final
            );

			if (this->max_log_likelihood_diri < logLikelihood_final) {				
				this->max_log_likelihood_diri = logLikelihood_final;				
				if (this->max_log_likelihood_best < this->max_log_likelihood_diri) {
					this->max_log_likelihood_best = this->max_log_likelihood_diri;
					this->StoreRootAndRootProbability();
					this->StoreTransitionMatrices();
				}
			}
		}
	}		
	cout << "max log likelihood obtained using Dirichlet parameters is " << setprecision(10) << this->max_log_likelihood_diri << endl;
}

void SEM::EMtrifle_DNA_rooted_at_each_internal_vertex_started_with_dirichlet_store_results() {
  int n = this->numberOfObservedVertices;
  int num_vertices = this->vertexMap->size();
  this->max_log_likelihood_diri = -1 * std::pow(10.0, 5.0);

  std::vector<int> vertex_indices_to_visit;
  vertex_indices_to_visit.reserve(std::max(0, num_vertices - n));
  for (int v_ind = n; v_ind < num_vertices; ++v_ind) vertex_indices_to_visit.push_back(v_ind);

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::mt19937 rng(seed);

  const int nodeTotal = static_cast<int>(vertex_indices_to_visit.size());

  auto emit_layer_done = [&](int rep1, int nodeIndex1, int layerIdx, const std::string& rootName, double ll_val) {
    nlohmann::json j{
      {"type", "progress"},
      {"kind", "layer:done"},
      {"method", "dirichlet"},
      {"rep", rep1},
      {"nodeIndex", nodeIndex1},
      {"nodeTotal", nodeTotal},
      {"layer", layerIdx},
      {"rootName", rootName},
      {"ll", ll_val}
    };
    std::cout << "[EM_Method] " << j.dump() << "\n";
    std::fflush(stdout);
  };

  int v_ind;
  SEM_vertex* v;

  for (int rep = 0; rep < num_repetitions; rep++) {
    std::cout << "rep " << rep + 1 << std::endl;

    for (int v_i = 0; v_i < static_cast<int>(vertex_indices_to_visit.size()); v_i++) {
      v_ind = vertex_indices_to_visit[v_i];
      v = (*this->vertexMap)[v_ind];
      std::cout << "node " << v_i + 1 << " of " << nodeTotal << ":" << v->name << std::endl;

      if (v->degree != 3) throw mt_error("Expect internal nodes to have degree three");

      // this->EMTrifle_current = EMTrifle_struct{};
      // this->EMTrifle_current.method = "dirichlet";
      // this->EMTrifle_current.root = v->name;
      // this->EMTrifle_current.rep = rep + 1;

      for (int layer = 0; layer < 3; layer++) {
        this->EMTrifle_started_with_dirichlet_rooted_at(v, layer);

        double ll_val = this->EMTrifle_current.ll_final_trifle[layer]; // if available
        emit_layer_done(/*rep1*/ rep + 1, /*nodeIndex1*/ v_i + 1, /*layerIdx*/ layer, /*rootName*/ v->name, ll_val);
      }
	  cout << "[EM_Trifle] " << this->emtrifle_to_json(this->EMTrifle_current) << endl;

      if (this->max_log_likelihood_diri < EMTrifle_current.ll_final_trifle[2]) {
        this->max_log_likelihood_diri = EMTrifle_current.ll_final_trifle[2];
        if (this->max_log_likelihood_best < this->max_log_likelihood_diri) {
          this->max_log_likelihood_best = this->max_log_likelihood_diri;
          this->StoreRootAndRootProbability();
          this->StoreTransitionMatrices();
        }
      }
      // this->EMTrifle_DNA_runs_diri.emplace_back(std::move(this->EMTrifle_current));
    }
  }
  std::cout << "max log likelihood obtained using Dirichlet parameters is "
            << std::setprecision(10) << this->max_log_likelihood_diri << std::endl;
}




void SEM::EM_rooted_at_each_internal_vertex_started_with_dirichlet(int num_repetitions) {
	int n = this->numberOfObservedVertices;
	int num_vertices = this->vertexMap->size();	
	SEM_vertex * v;
	vector <double> loglikelihoodscoresForEachRepetition;
	ofstream loglikelihood_node_rep_file;
	loglikelihood_node_rep_file.open(this->prefix_for_output_files + ".dirichlet_rooting_initial_final_rep_loglik");
    loglikelihood_node_rep_file << "root" << "\t"										
                                << "rep" << "\t"
                                << "iter" << "\t"
                                << "ll dirichlet" << "\t"
                                << "ecd-ll first" << "\t"
                                << "ecd-ll final" << "\t"
                                << "ll final" << endl;
	this->max_log_likelihood_diri = -1 * pow(10,5);
	double logLikelihood_pars;
	double loglikelihood_ecd_first;
	double loglikelihood_ecd_final;
	double logLikelihood_final;
	int iter;	
	tuple <int,double,double,double,double> iter_dirill_edllfirst_edllfinal_llfinal;
	vector<int> vertex_indices_to_visit;
	
	for (int v_ind = n; v_ind < num_vertices; v_ind++) {
		vertex_indices_to_visit.push_back(v_ind);
	}
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 rng(seed);

    // Shuffle the vector
	// cout << "randomizing the order in which nodes will be visited" << endl;
    shuffle(vertex_indices_to_visit.begin(), vertex_indices_to_visit.end(), rng);
    int v_ind;	
	for (int v_i = 0; v_i < vertex_indices_to_visit.size(); v_i++){
        v_ind = vertex_indices_to_visit[v_i];
		v = (*this->vertexMap)[v_ind];
		cout << "node " << v_i+1 << ":" << v->name << endl;
		if(v->degree != 3){
            throw mt_error("Expect internal nodes to have degree three");
        }
		loglikelihoodscoresForEachRepetition.clear();
		for (int rep = 0; rep < num_repetitions; rep++) {					
			iter_dirill_edllfirst_edllfinal_llfinal = this->EM_started_with_dirichlet_rooted_at(v);			
			iter = get<0>(iter_dirill_edllfirst_edllfinal_llfinal);
			logLikelihood_pars = get<1>(iter_dirill_edllfirst_edllfinal_llfinal);
			loglikelihood_ecd_first = get<2>(iter_dirill_edllfirst_edllfinal_llfinal);
			loglikelihood_ecd_final = get<3>(iter_dirill_edllfirst_edllfinal_llfinal);
			logLikelihood_final = get<4>(iter_dirill_edllfirst_edllfinal_llfinal);
			loglikelihood_node_rep_file << v->name << "\t"										
										<< rep +1 << "\t"
										<< iter << "\t"
										<< setprecision(ll_precision) << logLikelihood_pars << "\t"
										<< setprecision(ll_precision) << loglikelihood_ecd_first << "\t"
										<< setprecision(ll_precision) << loglikelihood_ecd_final << "\t"
										<< setprecision(ll_precision) << logLikelihood_final << endl;
			if (this->max_log_likelihood_diri < logLikelihood_final) {							
				this->max_log_likelihood_diri = logLikelihood_final;
				if (this->max_log_likelihood_best < this->max_log_likelihood_diri) {
					this->max_log_likelihood_best = this->max_log_likelihood_diri;
					this->StoreRootAndRootProbability();
					this->StoreTransitionMatrices();
				}				
			}
		}
	}
	
	loglikelihood_node_rep_file.close();
	cout << "max log likelihood, precision 10, obtained using Dirichlet parameters is " << setprecision(10) << this->max_log_likelihood_diri << endl;	
}

void SEM::EMTrifle_DNA_for_replicate(int rep) {
  this->max_log_likelihood_rep = -1 * std::pow(10.0, 5.0);
  int n = this->numberOfObservedVertices;
  int num_vertices = this->vertexMap->size();  

  std::vector<int> vertex_indices_to_visit;
  vertex_indices_to_visit.reserve(std::max(0, num_vertices - n));
  string node_name;
  SEM_vertex * v;
  int v_ind;
  
  for (int i = 1; i < n-1; i++) {
	node_name = "h_" + to_string(i);	
	cout << node_name << endl;
	vertex_indices_to_visit.emplace_back(this->GetVertexId(node_name));
  }
  
    
  const int nodeTotal = static_cast<int>(vertex_indices_to_visit.size());

  // --- Emitters (define BEFORE use) -----------------------------------------
  auto emit_method_start = [&](const char* method) {
    nlohmann::json j{
      {"type","progress"},
      {"kind","method:start"},
      {"method", method},
      {"rep", rep}
    };
    std::cout << "[EM_Method] " << j.dump() << "\n";
    std::fflush(stdout);
  };

  auto emit_layer_done = [&](const char* method,
                             int nodeIndex1,
                             int layerIdx,
                             const std::string& rootName,
                             double ll_val) {
    nlohmann::json j{
      {"type","progress"},
      {"kind","layer:done"},
      {"method", method},
      {"rep", rep},
      {"nodeIndex", nodeIndex1},
      {"nodeTotal", nodeTotal},
      {"layer", layerIdx},
      {"rootName", rootName},
      {"ll", ll_val}
    };
    std::cout << "[EM_Method] " << j.dump() << "\n";
    std::fflush(stdout);
  };
  // -------------------------------------------------------------------------
  
  emit_method_start("dirichlet");
  for (int v_i = 0; v_i < nodeTotal; ++v_i) {	
    v_ind = vertex_indices_to_visit[v_i];
    v = (*this->vertexMap)[v_ind];
	// ---- Dirichlet (layers 0..2) ----   
	this->EMTrifle_current = EMTrifle_struct{};
	this->EMTrifle_current.method = "dirichlet";
	this->EMTrifle_current.root = v->name;
	this->EMTrifle_current.rep = rep;    
	for (int layer = 0; layer < 3; ++layer) {
      this->EMTrifle_started_with_dirichlet_rooted_at(v, layer);      
	  double ll_val = this->EMTrifle_current.ll_final_trifle[layer];
      emit_layer_done("dirichlet", v_i + 1, layer, v->name, ll_val);
	  cout << "dirichlet" << "\t" << v_i + 1 << "\t" << layer << "\t" << v->name << "\t" << ll_val << endl;
    }
	// cout << this->emtrifle_to_json(this->EMTrifle_current) << endl;
	cout << "[EM_Trifle] " << this->emtrifle_to_json(this->EMTrifle_current) << endl;
    
	if (this->max_log_likelihood_rep < EMTrifle_current.ll_final_trifle[2]) {
      	this->max_log_likelihood_rep = EMTrifle_current.ll_final_trifle[2];
      	if (this->max_log_likelihood_best < this->max_log_likelihood_rep) {
				this->max_log_likelihood_best = this->max_log_likelihood_rep;
				this->StoreRootAndRootProbability();
				this->StoreTransitionMatrices();
      		}
    	}
	}

	emit_method_start("parsimony");
	for (int v_i = 0; v_i < nodeTotal; ++v_i) {
		v_ind = vertex_indices_to_visit[v_i];
    	v = (*this->vertexMap)[v_ind];
		// ---- Parsimony (layers 0..2) ----
		this->EMTrifle_current = EMTrifle_struct{};
		this->EMTrifle_current.method = "parsimony";
		this->EMTrifle_current.root = v->name;
		this->EMTrifle_current.rep = rep;		
		for (int layer = 0; layer < 3; ++layer) {
			this->EMTrifle_started_with_parsimony_rooted_at(v, layer);
			double ll_val = this->EMTrifle_current.ll_final_trifle[layer];
			emit_layer_done("parsimony", v_i + 1, layer, v->name, ll_val);
			cout << "parsimony" << "\t" << v_i + 1 << "\t" << layer << "\t" << v->name << "\t" << ll_val << endl;
			}		
		// cout << this->emtrifle_to_json(this->EMTrifle_current) << endl;
		cout << "[EM_Trifle] " << this->emtrifle_to_json(this->EMTrifle_current) << endl;

		if (this->max_log_likelihood_rep < EMTrifle_current.ll_final_trifle[2]) {
			this->max_log_likelihood_rep = EMTrifle_current.ll_final_trifle[2];
			if (this->max_log_likelihood_best < this->max_log_likelihood_rep) {
				this->max_log_likelihood_best = this->max_log_likelihood_rep;
				this->StoreRootAndRootProbability();
				this->StoreTransitionMatrices();
			}
		}
	}
	
	cout << "Restoring BH model with best probability" << endl;
    this->RestoreBestProbability();	
    this->ReparameterizeBH();
	cout << "BH reparameterized" << endl;
	
	emit_method_start("hss");
	for (int v_i = 0; v_i < nodeTotal; ++v_i) {
		v_ind = vertex_indices_to_visit[v_i];
    	v = (*this->vertexMap)[v_ind];
		// ---- HSS (layers 0..2) ----
		this->EMTrifle_current = EMTrifle_struct{};
		this->EMTrifle_current.method = "hss";
		this->EMTrifle_current.root = v->name;
		this->EMTrifle_current.rep = rep;		
		for (int layer = 0; layer < 3; ++layer) {
			this->EMTrifle_started_with_HSS_rooted_at(v, layer);
			double ll_val = this->EMTrifle_current.ll_final_trifle[layer];
			emit_layer_done("hss", v_i + 1, layer, v->name, ll_val);
			cout << "hss" << "\t" << v_i + 1 << "\t" << layer << "\t" << v->name << "\t" << ll_val << endl;
		}		
		// cout << this->emtrifle_to_json(this->EMTrifle_current) << endl;
		cout << "[EM_Trifle] " << this->emtrifle_to_json(this->EMTrifle_current) << endl;
	}
}


void SEM::EMTrifle_DNA_rooted_at_each_internal_vertex_started_with_parsimony_store_results() {
  int n = this->numberOfObservedVertices;
  int num_vertices = this->vertexMap->size();
  this->max_log_likelihood_pars = -1 * std::pow(10.0, 5.0);

  std::vector<int> vertex_indices_to_visit;
  vertex_indices_to_visit.reserve(std::max(0, num_vertices - n));
  for (int v_ind = n; v_ind < num_vertices; ++v_ind) vertex_indices_to_visit.push_back(v_ind);

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::mt19937 rng(seed);

  const int nodeTotal = static_cast<int>(vertex_indices_to_visit.size());

  
  auto emit_layer_done = [&](int rep, int nodeIndex1, int layerIdx, const std::string& rootName, double ll_val) {
    nlohmann::json j{
      {"type", "progress"},
      {"kind", "layer:done"},
      {"method", "parsimony"},
      {"rep", rep},
      {"nodeIndex", nodeIndex1},
      {"nodeTotal", nodeTotal},
      {"layer", layerIdx},
      {"rootName", rootName},
      {"ll", ll_val}
    };
    std::cout << "[EM_Method] " << j.dump() << "\n";
    std::fflush(stdout);
  };

  int v_ind;
  SEM_vertex* v;

  for (int rep = 0; rep < num_repetitions; rep++) {
    std::cout << "rep " << rep + 1 << std::endl;

    for (int v_i = 0; v_i < static_cast<int>(vertex_indices_to_visit.size()); v_i++) {
      v_ind = vertex_indices_to_visit[v_i];
      v = (*this->vertexMap)[v_ind];
      std::cout << "node " << v_i + 1 << " of " << nodeTotal << ":" << v->name << std::endl;

      if (v->degree != 3) throw mt_error("Expect internal nodes to have degree three");

      for (int layer = 0; layer < 3; layer++) {
        this->EMTrifle_started_with_parsimony_rooted_at(v, layer);

        
        double ll_val = this->EMTrifle_current.ll_final_trifle[layer];
        emit_layer_done(rep + 1, v_i + 1, layer, v->name, ll_val);
      }

      if (this->max_log_likelihood_pars < EMTrifle_current.ll_final_trifle[2]) {
        this->max_log_likelihood_pars = EMTrifle_current.ll_final_trifle[2];
        if (this->max_log_likelihood_best < this->max_log_likelihood_pars) {
          this->max_log_likelihood_best = this->max_log_likelihood_pars;
          this->StoreRootAndRootProbability();
          this->StoreTransitionMatrices();
        }
      }      
    }
  }
  std::cout << "max log likelihood obtained using Parsimony parameters is "
            << std::setprecision(10) << this->max_log_likelihood_pars << std::endl;
}



void SEM::EM_DNA_rooted_at_each_internal_vertex_started_with_parsimony_store_results(int num_repetitions) {
	int n = this->numberOfObservedVertices;	
	int num_vertices = this->vertexMap->size();	
	this->max_log_likelihood_pars = -1 * pow(10,5);
	
	vector<int> vertex_indices_to_visit;
	vertex_indices_to_visit.reserve(std::max(0, num_vertices - n));
	
	for (int v_ind = n; v_ind < num_vertices; ++v_ind) vertex_indices_to_visit.push_back(v_ind);

	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 rng(seed);

    // Shuffle the vector
	int num_vertices_to_visit = vertex_indices_to_visit.size();
	cout << "randomizing the order in which " << num_vertices_to_visit << " nodes will be visited" << endl;
    shuffle(vertex_indices_to_visit.begin(), vertex_indices_to_visit.end(), rng);
    
	for (int v_i = 0; v_i < vertex_indices_to_visit.size(); v_i++) {
		const int v_ind = vertex_indices_to_visit[v_i];
		SEM_vertex * v = (*this->vertexMap)[v_ind];
		cout << "node " << v_i+1 << " of " << num_vertices_to_visit << ":" << v->name << endl;				
		if (v->degree != 3) throw mt_error("Expect internal nodes to have degree three");
		// cout << "here1" << endl;
		for (int rep = 0; rep < num_repetitions; rep++) {
			// cout << " size of EM_DNA_runs_pars is " << this->EM_DNA_runs_pars.size() << endl;
			this->EM_current = EM_struct{};
			this->EM_current.method = "parsimony";
			this->EM_current.root = v->name;
			this->EM_current.rep = rep + 1;
			auto tup = this->EM_started_with_parsimony_rooted_at(v);			
			// cout << string_EM_pars << endl;
			EM_struct EM_pars = this->EM_current;
			this->EM_DNA_runs_pars.push_back(EM_pars);			
			const int    iter                    = std::get<0>(tup);
            const double logLikelihood_pars      = std::get<1>(tup);
            const double loglikelihood_ecd_first = std::get<2>(tup);
            const double loglikelihood_ecd_final = std::get<3>(tup);
            const double logLikelihood_final     = std::get<4>(tup);
			// cout << "initial " << logLikelihood_pars << "iter " << iter << endl;
			// cout << "final " << logLikelihood_final << "iter " << iter << endl;

			emtr::push_result(
                this->EMTR_results,                 // vector<tuple<string,string,int,int,double,double,double,double>>
                "parsimony",                        // init_method
                v->name,                            // root
                rep + 1,                            // repetition (1-based)
                iter,
                logLikelihood_pars,                 // ll_initial
                loglikelihood_ecd_first,
                loglikelihood_ecd_final,
                logLikelihood_final
            );

			if (this->max_log_likelihood_pars < logLikelihood_final) {				
				this->max_log_likelihood_pars = logLikelihood_final;				
				if (this->max_log_likelihood_best < this->max_log_likelihood_pars) {
					this->max_log_likelihood_best = this->max_log_likelihood_pars;
					this->StoreRootAndRootProbability();
					this->StoreTransitionMatrices();
				}
			}
		}
	}
		
	cout << "max log likelihood obtained using Parsimony parameters is " << setprecision(10) << this->max_log_likelihood_pars << endl;	
}


void SEM::EM_rooted_at_each_internal_vertex_started_with_parsimony(int num_repetitions) {
	int n = this->numberOfObservedVertices;
	int num_vertices = this->vertexMap->size();	
	SEM_vertex * v;
	vector <double> loglikelihoodscoresForEachRepetition;
	ofstream loglikelihood_node_rep_file;
	loglikelihood_node_rep_file.open(this->prefix_for_output_files + ".pars_rooting_initial_final_rep_loglik");
    loglikelihood_node_rep_file << "root" << "\t"										
                                << "rep" << "\t"
                                << "iter" << "\t"
                                << "ll pars" << "\t"
                                << "ecd-ll first" << "\t"
                                << "ecd-ll final" << "\t"
                                << "ll final" << endl;
	this->max_log_likelihood_pars = -1 * pow(10,5);
	double logLikelihood_pars;
	double loglikelihood_ecd_first;
	double loglikelihood_ecd_final;
	double logLikelihood_final;
	int iter;	
	tuple <int,double,double,double,double> iter_parsll_edllfirst_edllfinal_llfinal;
	vector<int> vertex_indices_to_visit;
	
	for (int v_ind = n; v_ind < num_vertices; v_ind++) {
		vertex_indices_to_visit.push_back(v_ind);
	}
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 rng(seed);

    // Shuffle the vector
	cout << "randomizing the order in which nodes will be visited" << endl;
    shuffle(vertex_indices_to_visit.begin(), vertex_indices_to_visit.end(), rng);
    int v_ind;
	for (int v_i = 0; v_i < vertex_indices_to_visit.size(); v_i++){
        v_ind = vertex_indices_to_visit[v_i];
		v = (*this->vertexMap)[v_ind];
		cout << "node " << v_i+1 << ":" << v->name << endl;
		if(v->degree != 3){
            throw mt_error("Expect internal nodes to have degree three");
        }
		loglikelihoodscoresForEachRepetition.clear();
		for (int rep = 0; rep < num_repetitions; rep++) {				
			iter_parsll_edllfirst_edllfinal_llfinal = this->EM_started_with_parsimony_rooted_at(v);			
			iter = get<0>(iter_parsll_edllfirst_edllfinal_llfinal);
			logLikelihood_pars = get<1>(iter_parsll_edllfirst_edllfinal_llfinal);
			loglikelihood_ecd_first = get<2>(iter_parsll_edllfirst_edllfinal_llfinal);
			loglikelihood_ecd_final = get<3>(iter_parsll_edllfirst_edllfinal_llfinal);
			logLikelihood_final = get<4>(iter_parsll_edllfirst_edllfinal_llfinal);
			loglikelihood_node_rep_file << v->name << "\t"										
										<< rep +1 << "\t"
										<< iter << "\t"
										<< setprecision(ll_precision) << logLikelihood_pars << "\t"
										<< setprecision(ll_precision) << loglikelihood_ecd_first << "\t"
										<< setprecision(ll_precision) << loglikelihood_ecd_final << "\t"
										<< setprecision(ll_precision) << logLikelihood_final << endl;
			if (this->max_log_likelihood_pars < logLikelihood_final) {				
				this->max_log_likelihood_pars = logLikelihood_final;
				if (this->max_log_likelihood_best < this->max_log_likelihood_pars) {
					this->max_log_likelihood_best = this->max_log_likelihood_pars;
					this->StoreRootAndRootProbability();
					this->StoreTransitionMatrices();
				}
			}
		}
	}
	
	loglikelihood_node_rep_file.close();	
	cout << "max log likelihood obtained using Parsimony parameters is " << setprecision(10) << this->max_log_likelihood_pars << endl;
	// (*this->logFile) << "max log likelihood obtained using Parsimony parameters is " << setprecision(ll_precision) << max_log_likelihood << endl;	
}


void SEM::EM_DNA_rooted_at_each_internal_vertex_started_with_HSS_store_results(int num_repetitions) {
	int n = this->numberOfObservedVertices;
	int num_vertices = this->vertexMap->size();	
		
	this->max_log_likelihood_hss = -1 * pow(10,5);
	
	vector<int> vertex_indices_to_visit;
	vertex_indices_to_visit.reserve(std::max(0, num_vertices - n));
	
	for (int v_ind = n; v_ind < num_vertices; ++v_ind) vertex_indices_to_visit.push_back(v_ind);

	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 rng(seed);

    // Shuffle the vector
	cout << "randomizing the order in which nodes will be visited" << endl;
    shuffle(vertex_indices_to_visit.begin(), vertex_indices_to_visit.end(), rng);

    int num_vertices_to_visit = vertex_indices_to_visit.size();
	for (int v_i = 0; v_i < vertex_indices_to_visit.size(); v_i++){
        const int v_ind = vertex_indices_to_visit[v_i];
		SEM_vertex * v = (*this->vertexMap)[v_ind];
		cout << "node " << v_i+1 << " of " << num_vertices_to_visit << ":" << v->name << endl;
		if(v->degree != 3) throw mt_error("Expect internal nodes to have degree three");
		
		for (int rep = 0; rep < num_repetitions; rep++) {
			this->EM_current = EM_struct{};
			this->EM_current.method = "hss";
			this->EM_current.root = v->name;
			this->EM_current.rep = rep + 1;
			auto tup = this->EM_started_with_HSS_parameters_rooted_at(v);
			EM_struct EM_hss{this->EM_current};
			this->EM_DNA_runs_hss.push_back(EM_hss);
			const int    iter                    = std::get<0>(tup);
            const double logLikelihood_hss       = std::get<1>(tup);
            const double loglikelihood_ecd_first = std::get<2>(tup);
            const double loglikelihood_ecd_final = std::get<3>(tup);
            const double logLikelihood_final     = std::get<4>(tup);

			emtr::push_result(
                this->EMTR_results,
                "hss",
                v->name,
                rep + 1,
                iter,
                logLikelihood_hss,
                loglikelihood_ecd_first,
                loglikelihood_ecd_final,
                logLikelihood_final
            );


			if (this->max_log_likelihood_hss < logLikelihood_final) {				
				this->max_log_likelihood_hss = logLikelihood_final;				
				if (this->max_log_likelihood_best < this->max_log_likelihood_hss) {
					this->max_log_likelihood_best = this->max_log_likelihood_hss;
					this->StoreRootAndRootProbability();
					this->StoreTransitionMatrices();
				}
			}
		}
	}
		
	cout << "max log likelihood obtained using HSS parameters is " << setprecision(10) << this->max_log_likelihood_hss << endl;	
}

void SEM::EMtrifle_DNA_rooted_at_each_internal_vertex_started_with_HSS_store_results() {
  int n = this->numberOfObservedVertices;
  int num_vertices = this->vertexMap->size();
  this->max_log_likelihood_hss = -1 * std::pow(10.0, 5.0);

  std::vector<int> vertex_indices_to_visit;
  vertex_indices_to_visit.reserve(std::max(0, num_vertices - n));
  for (int v_ind = n; v_ind < num_vertices; ++v_ind) vertex_indices_to_visit.push_back(v_ind);

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::mt19937 rng(seed);

  const int nodeTotal = static_cast<int>(vertex_indices_to_visit.size());

  auto emit_layer_done = [&](int rep1, int nodeIndex1, int layerIdx, const std::string& rootName, double ll_val) {
    nlohmann::json j{
      {"type", "progress"},
      {"kind", "layer:done"},
      {"method", "hss"},
      {"rep", rep1},
      {"nodeIndex", nodeIndex1},
      {"nodeTotal", nodeTotal},
      {"layer", layerIdx},
      {"rootName", rootName},
      {"ll", ll_val}
    };
    std::cout << "[EM_Method] " << j.dump() << "\n";
    std::fflush(stdout);
  };

  int v_ind;
  SEM_vertex* v;

  for (int rep = 0; rep < num_repetitions; rep++) {
    std::cout << "rep " << rep + 1 << std::endl;

    for (int v_i = 0; v_i < static_cast<int>(vertex_indices_to_visit.size()); v_i++) {
      v_ind = vertex_indices_to_visit[v_i];
      v = (*this->vertexMap)[v_ind];
      std::cout << "node " << v_i + 1 << " of " << nodeTotal << ":" << v->name << std::endl;

      if (v->degree != 3) throw mt_error("Expect internal nodes to have degree three");

      // this->EMTrifle_current = EMTrifle_struct{};
      // this->EMTrifle_current.method = "hss";
      // this->EMTrifle_current.root = v->name;
      // this->EMTrifle_current.rep = rep + 1;

      for (int layer = 0; layer < 3; layer++) {
        this->EMTrifle_started_with_HSS_rooted_at(v, layer);

        double ll_val = this->EMTrifle_current.ll_final_trifle[layer]; // if available
        emit_layer_done(/*rep1*/ rep + 1, /*nodeIndex1*/ v_i + 1, /*layerIdx*/ layer, /*rootName*/ v->name, ll_val);
      }

      if (this->max_log_likelihood_hss < this->EMTrifle_current.ll_final_trifle[2]) {
        this->max_log_likelihood_hss = this->EMTrifle_current.ll_final_trifle[2];
        if (this->max_log_likelihood_best < this->max_log_likelihood_hss) {
          this->max_log_likelihood_best = this->max_log_likelihood_hss;
          this->StoreRootAndRootProbability();
          this->StoreTransitionMatrices();
        }
      }
      // this->EMTrifle_DNA_runs_hss.emplace_back(std::move(this->EMTrifle_current));
    }
  }
  std::cout << "max log likelihood obtained using HSS is "
            << std::setprecision(10) << this->max_log_likelihood_hss << std::endl;
}


void SEM::EM_rooted_at_each_internal_vertex_started_with_HSS_par(int num_repetitions) {
	// cout << "convergence threshold for EM is " << this->logLikelihoodConvergenceThreshold << endl;
	// (* this->logFile) << "convergence threshold for EM is " << this->logLikelihoodConvergenceThreshold << endl;
	// cout << "maximum number of EM iterations allowed is " << this->maxIter << endl;
	// (* this->logFile) << "maximum number of EM iterations allowed is " << this->maxIter << endl;
	int n = this->numberOfObservedVertices;
	int num_vertices = this->vertexMap->size();	
	SEM_vertex * v;
	vector <double> loglikelihoodscoresForEachRepetition;
	ofstream loglikelihood_node_rep_file;
	loglikelihood_node_rep_file.open(this->prefix_for_output_files + ".HSS_rooting_initial_final_rep_loglik");
    loglikelihood_node_rep_file << "root" << "\t"
                                << "rep" << "\t"
                                << "iter" << "\t"
                                << "ll HSS" << "\t"
                                << "ecd ll first" << "\t"
                                << "ecd ll final" << "\t"
                                << "ll final" << endl;
	double max_log_likelihood = -1 * pow(10,5);	
	double logLikelihood_pars;
	double loglikelihood_ecd_first;
	double loglikelihood_ecd_final;
	double logLikelihood_final;
	int iter;	
	tuple <int,double,double,double,double> iter_parsll_edllfirst_edllfinal_llfinal;
	vector<int> vertex_indices_to_visit;
	
	for (int v_ind = n; v_ind < num_vertices; v_ind++) {
		vertex_indices_to_visit.push_back(v_ind);
	}
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 rng(seed);

    // Shuffle the vector
	cout << "randomizing the order in which nodes will be visited" << endl;
    shuffle(vertex_indices_to_visit.begin(), vertex_indices_to_visit.end(), rng);
	int v_ind;
	for (int v_i = 0; v_i < vertex_indices_to_visit.size(); v_i++){
        v_ind = vertex_indices_to_visit[v_i];
		v = (*this->vertexMap)[v_ind];
		cout << "node " << v_i+1 << ":" << v->name << endl;		
		if(v->degree != 3){
            throw mt_error("Expect internal nodes to have degree three");
        }
		loglikelihoodscoresForEachRepetition.clear();
		for (int rep = 0; rep < num_repetitions; rep++) {					
			iter_parsll_edllfirst_edllfinal_llfinal = this->EM_started_with_HSS_parameters_rooted_at(v);
			iter = get<0>(iter_parsll_edllfirst_edllfinal_llfinal);
			logLikelihood_pars = get<1>(iter_parsll_edllfirst_edllfinal_llfinal);
			loglikelihood_ecd_first = get<2>(iter_parsll_edllfirst_edllfinal_llfinal);
			loglikelihood_ecd_final = get<3>(iter_parsll_edllfirst_edllfinal_llfinal);
			logLikelihood_final = get<4>(iter_parsll_edllfirst_edllfinal_llfinal);
			loglikelihood_node_rep_file << v->name << "\t"										
										<< rep +1 << "\t"
										<< iter << "\t"
										<< setprecision(ll_precision) << logLikelihood_pars << "\t"
										<< setprecision(ll_precision) << loglikelihood_ecd_first << "\t"
										<< setprecision(ll_precision) << loglikelihood_ecd_final << "\t"
										<< setprecision(ll_precision) << logLikelihood_final << endl;
			if (max_log_likelihood < logLikelihood_final) {				
				// this->WriteProbabilities();
				max_log_likelihood = logLikelihood_final;
			}
		}
	}
	loglikelihood_node_rep_file.close();
	cout << "max log likelihood obtained using HSS parameters is " << setprecision(10) << max_log_likelihood << endl;
	// (*this->logFile) << "max log likelihood obtained using HSS parameters is " << setprecision(ll_precision) << max_log_likelihood << endl;	
}

void SEM::EMTrifle_started_with_dirichlet_rooted_at(SEM_vertex *v, int layer) {
	// iterate over each internal node		
	this->RootTreeAtVertex(v);
	if (layer == 0) {		
		this->SetInitialEstimateOfModelParametersUsingDirichlet();		
	}
	this->BumpZeroEntriesOfModelParameters();
	this->StoreInitialParamsInEMTrifleCurrent(layer);
	// cout << "10d" << endl;
	this->ComputeTrifleLogLikelihood(layer);
	// cout << "10e" << endl;		
	this->EMTrifle_current.ll_initial_trifle[layer] = this->logLikelihood;

	cout << "initial log likelihood is " << this->EMTrifle_current.ll_initial_trifle[layer] << endl;
	
	double logLikelihood_exp_data_previous;
	double logLikelihood_exp_data_first;
	double logLikelihood_exp_data_final;
	int iter = 0;	
	bool continueIterations = 1;
	this->debug = 0;
	bool verbose = 0;		
	// (*this->logFile) << "log-likelihood computed by marginalization using parsimony parameters is " << setprecision(ll_precision) << logLikelihood_pars << endl;
	logLikelihood_exp_data_previous = -1 * pow(10,4);
	this->ConstructCliqueTree();
	while (continueIterations) {
		// t_start_time = chrono::high_resolution_clock::now();
		iter += 1;
		// cout << "iteration " << iter << endl;			
		// cout << "here 2" << endl;
		this->ComputeExpectedCountsForTrifle(layer);
		// cout << "here 3" << endl;
		this->ComputeMarginalProbabilitiesUsingExpectedCounts();
		// cout << "here 4" << endl;
		this->ComputeMLEstimateOfBHGivenExpectedDataCompletion();
		// cout << "here 5" << endl;
		this->ComputeLogLikelihoodUsingExpectedDataCompletion();
		// cout << "here 6" << endl;
		// cout << "log-likelihood computed using expected counts after EM iteration " << iter << " is " << setprecision(ll_precision) << this->logLikelihood << endl;
		// (*this->logFile) << "log-likelihood computed using expected counts after EM iteration " << iter << " is " << setprecision(ll_precision) << this->logLikelihood << endl;		
		if (iter == 1) {			
			logLikelihood_exp_data_first = this->logLikelihood;
            logLikelihood_exp_data_previous = this->logLikelihood;
			this->EMTrifle_current.ecd_ll_per_iter_for_trifle[layer][iter] = this->logLikelihood;
		} else if ((this->logLikelihood > logLikelihood_exp_data_previous + this->ecdllConvergenceThresholdForTrifle[layer]) and (iter < this->maxIter)) {
			logLikelihood_exp_data_previous = this->logLikelihood;
			this->EMTrifle_current.ecd_ll_per_iter_for_trifle[layer][iter] = this->logLikelihood;
		} else {
			continueIterations = 0;
			logLikelihood_exp_data_final = logLikelihood_exp_data_previous;
		}
	}
	this->ComputeTrifleLogLikelihood(layer);
	this->StoreFinalParamsInEMTrifleCurrent(layer);	
	this->EMTrifle_current.iter_trifle[layer] = iter;
	this->EMTrifle_current.ll_final_trifle[layer] = this->logLikelihood;	
	cout << "final log likelihood is " << this->EMTrifle_current.ll_final_trifle[layer] << endl;
	return ;
	// return tuple<int,double,double,double,double>(iter,logLikelihood_pars,logLikelihood_exp_data_first,logLikelihood_exp_data_final,this->logLikelihood);
}


void SEM::EMTrifle_started_with_HSS_rooted_at(SEM_vertex *v, int layer) {	
	// iterate over each internal node	
	this->RootTreeAtVertex(v);
	if (layer == 0) {
		this->SetInitialEstimateOfModelParametersUsingHSS();		
	}
	this->BumpZeroEntriesOfModelParameters();
	this->StoreInitialParamsInEMTrifleCurrent(layer);
	// cout << "10d" << endl;
	this->ComputeTrifleLogLikelihood(layer);
	// cout << "10e" << endl;		
	this->EMTrifle_current.ll_initial_trifle[layer] = this->logLikelihood;
	cout << "initial log likelihood is " << this->EMTrifle_current.ll_initial_trifle[layer] << endl;
	
	// this->StoreParamsInEMCurrent("init");
	// cout << "10d" << endl;
	this->ResetAncestralSequences();
	
	double logLikelihood_exp_data_previous;
	double logLikelihood_exp_data_first;
	double logLikelihood_exp_data_final;
	int iter = 0;	
	bool continueIterations = 1;
	this->debug = 0;
	bool verbose = 0;
	
	// cout << "-    -     -     -     -     -     -     -     -     -     -     -     -     -" << endl;
	// cout << "log-likelihood computed by marginalization using parsimony parameters is " << setprecision(ll_precision) << logLikelihood_pars << endl;
	// (*this->logFile) << "-    -     -     -     -     -     -     -     -     -     -     -     -     -" << endl;
	// (*this->logFile) << "log-likelihood computed by marginalization using parsimony parameters is " << setprecision(ll_precision) << logLikelihood_pars << endl;
	logLikelihood_exp_data_previous = -1 * pow(10,4);
	this->ConstructCliqueTree();
	while (continueIterations) {
		// t_start_time = chrono::high_resolution_clock::now();
		iter += 1;
		// cout << "iteration " << iter << endl;			
		// cout << "here 2" << endl;
		this->ComputeExpectedCountsForTrifle(layer);
		// cout << "here 3" << endl;
		this->ComputeMarginalProbabilitiesUsingExpectedCounts();
		// cout << "here 4" << endl;
		this->ComputeMLEstimateOfBHGivenExpectedDataCompletion();
		// cout << "here 5" << endl;
		this->ComputeLogLikelihoodUsingExpectedDataCompletion();
		// cout << "here 6" << endl;
		// cout << "log-likelihood computed using expected counts after EM iteration " << iter << " is " << setprecision(ll_precision) << this->logLikelihood << endl;
		// (*this->logFile) << "log-likelihood computed using expected counts after EM iteration " << iter << " is " << setprecision(ll_precision) << this->logLikelihood << endl;		
		if (iter == 1) {			
			logLikelihood_exp_data_first = this->logLikelihood;
            logLikelihood_exp_data_previous = this->logLikelihood;
			this->EMTrifle_current.ecd_ll_per_iter_for_trifle[layer][iter] = this->logLikelihood;
		} else if ((this->logLikelihood > logLikelihood_exp_data_previous + this->ecdllConvergenceThresholdForTrifle[layer]) and (iter < this->maxIter)) {
			logLikelihood_exp_data_previous = this->logLikelihood;
			this->EMTrifle_current.ecd_ll_per_iter_for_trifle[layer][iter] = this->logLikelihood;
		} else {
			continueIterations = 0;
			logLikelihood_exp_data_final = logLikelihood_exp_data_previous;
		}
	}
	this->ComputeTrifleLogLikelihood(layer);
	this->StoreFinalParamsInEMTrifleCurrent(layer);	
	this->EMTrifle_current.iter_trifle[layer] = iter;
	this->EMTrifle_current.ll_final_trifle[layer] = this->logLikelihood;	
	cout << "final log likelihood is " << this->EMTrifle_current.ll_final_trifle[layer] << endl;
	return ;
}

void SEM::EMTrifle_started_with_parsimony_rooted_at(SEM_vertex *v, int layer) {
	// iterate over each internal node	
	this->RootTreeAtVertex(v);
	if (layer == 0) {		
		this->ComputeMPEstimateOfAncestralSequencesForTrifle(layer);
		this->ComputeInitialEstimateOfModelParametersForTrifle(layer);
	}
	this->BumpZeroEntriesOfModelParameters();
	this->StoreInitialParamsInEMTrifleCurrent(layer);	
	// cout << "10d" << endl;
	this->ComputeTrifleLogLikelihood(layer);
	// cout << "10e" << endl;		
	this->EMTrifle_current.ll_initial_trifle[layer] = this->logLikelihood;
	cout << "initial log likelihood is " << this->EMTrifle_current.ll_initial_trifle[layer] << endl;
	
	// this->StoreParamsInEMCurrent("init");
	// cout << "10d" << endl;
	this->ResetAncestralSequences();
	
	double logLikelihood_exp_data_previous;
	double logLikelihood_exp_data_first;
	double logLikelihood_exp_data_final;
	int iter = 0;	
	bool continueIterations = 1;
	this->debug = 0;
	bool verbose = 0;
	
	// cout << "-    -     -     -     -     -     -     -     -     -     -     -     -     -" << endl;
	// cout << "log-likelihood computed by marginalization using parsimony parameters is " << setprecision(ll_precision) << logLikelihood_pars << endl;
	// (*this->logFile) << "-    -     -     -     -     -     -     -     -     -     -     -     -     -" << endl;
	// (*this->logFile) << "log-likelihood computed by marginalization using parsimony parameters is " << setprecision(ll_precision) << logLikelihood_pars << endl;
	logLikelihood_exp_data_previous = -1 * pow(10,4);
	this->ConstructCliqueTree();
	while (continueIterations) {
		// t_start_time = chrono::high_resolution_clock::now();
		iter += 1;
		// cout << "iteration " << iter << endl;			
		// cout << "here 2" << endl;
		this->ComputeExpectedCountsForTrifle(layer);
		// cout << "here 3" << endl;
		this->ComputeMarginalProbabilitiesUsingExpectedCounts();
		// cout << "here 4" << endl;
		this->ComputeMLEstimateOfBHGivenExpectedDataCompletion();
		// cout << "here 5" << endl;
		this->ComputeLogLikelihoodUsingExpectedDataCompletion();
		// cout << "here 6" << endl;
		// cout << "log-likelihood computed using expected counts after EM iteration " << iter << " is " << setprecision(ll_precision) << this->logLikelihood << endl;
		// (*this->logFile) << "log-likelihood computed using expected counts after EM iteration " << iter << " is " << setprecision(ll_precision) << this->logLikelihood << endl;		
		if (iter == 1) {			
			logLikelihood_exp_data_first = this->logLikelihood;
            logLikelihood_exp_data_previous = this->logLikelihood;
			this->EMTrifle_current.ecd_ll_per_iter_for_trifle[layer][iter] = this->logLikelihood;
		} else if ((this->logLikelihood > logLikelihood_exp_data_previous + this->ecdllConvergenceThresholdForTrifle[layer]) and (iter < this->maxIter)) {
			logLikelihood_exp_data_previous = this->logLikelihood;
			this->EMTrifle_current.ecd_ll_per_iter_for_trifle[layer][iter] = this->logLikelihood;
		} else {
			continueIterations = 0;
			logLikelihood_exp_data_final = logLikelihood_exp_data_previous;
		}
	}
	this->ComputeTrifleLogLikelihood(layer);
	this->StoreFinalParamsInEMTrifleCurrent(layer);	
	this->EMTrifle_current.iter_trifle[layer] = iter;
	this->EMTrifle_current.ll_final_trifle[layer] = this->logLikelihood;	
	cout << "final log likelihood is " << this->EMTrifle_current.ll_final_trifle[layer] << endl;
	return ;
	// return tuple<int,double,double,double,double>(iter,logLikelihood_pars,logLikelihood_exp_data_first,logLikelihood_exp_data_final,this->logLikelihood);
}

tuple <int,double,double,double,double> SEM::EM_started_with_parsimony_rooted_at(SEM_vertex *v) {
	//	cout << "10a" << endl;
	// iterate over each internal node	
	this->RootTreeAtVertex(v);
	// cout << "10b" << endl;
	this->ComputeMPEstimateOfAncestralSequences();	
	// cout << "10c" << endl;
	this->ComputeInitialEstimateOfModelParameters();
	this->StoreParamsInEMCurrent("init");
	this->ComputeLogLikelihood();
	// cout << "Initial value of log-likelihood is " << setprecision(ll_precision) << this->logLikelihood << endl;
	// cout << "10d" << endl;
	this->ResetAncestralSequences();
	double logLikelihood_pars;
	double logLikelihood_exp_data_previous;
	double logLikelihood_exp_data_first;
	double logLikelihood_exp_data_final;
	int iter = 0;	
	bool continueIterations = 1;
	this->debug = 0;
	bool verbose = 0;
	logLikelihood_pars = this->logLikelihood;	
	this->EM_current.ll_init = logLikelihood_pars;
	// cout << "-    -     -     -     -     -     -     -     -     -     -     -     -     -" << endl;
	// cout << "log-likelihood computed by marginalization using parsimony parameters is " << setprecision(ll_precision) << logLikelihood_pars << endl;
	// (*this->logFile) << "-    -     -     -     -     -     -     -     -     -     -     -     -     -" << endl;
	// (*this->logFile) << "log-likelihood computed by marginalization using parsimony parameters is " << setprecision(ll_precision) << logLikelihood_pars << endl;
	logLikelihood_exp_data_previous = -1 * pow(10,4);
	while (continueIterations) {
		// t_start_time = chrono::high_resolution_clock::now();
		iter += 1;					
		this->ConstructCliqueTree();			
		// 2. Compute expected counts
		this->ComputeExpectedCounts();

		this->ComputeMarginalProbabilitiesUsingExpectedCounts();
		
		this->ComputeMLEstimateOfBHGivenExpectedDataCompletion();	
		
		this->ComputeLogLikelihoodUsingExpectedDataCompletion();		
		
		// cout << "log-likelihood computed using expected counts after EM iteration " << iter << " is " << setprecision(ll_precision) << this->logLikelihood << endl;
		// (*this->logFile) << "log-likelihood computed using expected counts after EM iteration " << iter << " is " << setprecision(ll_precision) << this->logLikelihood << endl;		
		if (iter == 1){			
			logLikelihood_exp_data_first = this->logLikelihood;
            logLikelihood_exp_data_previous = this->logLikelihood;
			this->EM_current.ecd_ll_per_iter[iter] = this->logLikelihood;
		} else if ((this->logLikelihood > logLikelihood_exp_data_previous + this->ecdllConvergenceThreshold) and (iter < this->maxIter)) {
			logLikelihood_exp_data_previous = this->logLikelihood;
			this->EM_current.ecd_ll_per_iter[iter] = this->logLikelihood;
		} else {
			continueIterations = 0;
			logLikelihood_exp_data_final = logLikelihood_exp_data_previous;
		}
	}
	this->ComputeLogLikelihood();
	this->StoreParamsInEMCurrent("final");
	this->EM_current.iter = iter;
	this->EM_current.ll_final = this->logLikelihood;
	
	return tuple<int,double,double,double,double>(iter,logLikelihood_pars,logLikelihood_exp_data_first,logLikelihood_exp_data_final,this->logLikelihood);
}

/*
1. Number of iterations until EM converges
2. Initial log-likelihood with HSS based parameters
3. Expected-data log-likelihood after one EM iteration
4. Expected-data log-likelihood after last EM iteration
5. Final log-likelihood using EM based parameters
*/

tuple <int,double,double,double,double> SEM::EM_started_with_HSS_parameters_rooted_at(SEM_vertex *v) {
	//	cout << "10a" << endl;
	// iterate over each internal node	
	this->RootTreeAtVertex(v);
	// cout << "10b" << endl;	
	this->SetInitialEstimateOfModelParametersUsingHSS();
	this->StoreParamsInEMCurrent("init");
	this->ComputeLogLikelihood();	
	// cout << "10d" << endl;
	this->ResetAncestralSequences();
	double logLikelihood_hss;
	double logLikelihood_exp_data_previous;
	double logLikelihood_exp_data_first;
	double logLikelihood_exp_data_final;
	int iter = 0;	
	bool continueIterations = 1;
	this->debug = 0;
	bool verbose = 0;
	logLikelihood_hss = this->logLikelihood;
	this->EM_current.ll_init = logLikelihood_hss;
	// cout << "-    -     -     -     -     -     -     -     -     -     -     -     -     -" << endl;
	// cout << "log-likelihood computed by marginalization using HSS parameters is " << setprecision(ll_precision) << logLikelihood_hss << endl;
	// (*this->logFile) << "-    -     -     -     -     -     -     -     -     -     -     -     -     -" << endl;
	// (*this->logFile) << "log-likelihood computed by marginalization using HSS parameters is " << setprecision(ll_precision) << logLikelihood_hss << endl;
	logLikelihood_exp_data_previous = -1 * pow(10,4);
	while (continueIterations) {
		// t_start_time = chrono::high_resolution_clock::now();
		iter += 1;		
		this->ConstructCliqueTree();			
		// 2. Compute expected counts
		this->ComputeExpectedCounts();

		this->ComputeMarginalProbabilitiesUsingExpectedCounts();

		this->ComputeMLEstimateOfBHGivenExpectedDataCompletion();
		
		this->ComputeLogLikelihoodUsingExpectedDataCompletion();
		
		// cout << "log-likelihood computed using expected counts after EM iteration " << iter << " is " << setprecision(ll_precision) << this->logLikelihood << endl;
		// (*this->logFile) << "log-likelihood computed using expected counts after EM iteration " << iter << " is " << setprecision(ll_precision) << this->logLikelihood << endl;		
		if (iter == 1) {
			logLikelihood_exp_data_previous = this->logLikelihood;
			logLikelihood_exp_data_first = this->logLikelihood;
			this->EM_current.ecd_ll_per_iter[iter] = this->logLikelihood;
		} else if ((this->logLikelihood > logLikelihood_exp_data_previous + this->ecdllConvergenceThreshold) and (iter < this->maxIter)) {
			logLikelihood_exp_data_previous = this->logLikelihood;
			this->EM_current.ecd_ll_per_iter[iter] = this->logLikelihood;
		} else {
			continueIterations = 0;
			logLikelihood_exp_data_final = logLikelihood_exp_data_previous;
		}
	}
	this->ComputeLogLikelihood();
	this->StoreParamsInEMCurrent("final");
	this->EM_current.iter = iter;
	this->EM_current.ll_final = this->logLikelihood;
	// cout << "log-likelihood computed by marginalization after EM iteration " << iter << " is " << setprecision(ll_precision) << this->logLikelihood << endl;
	// cout << "- + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + -" << endl;
	// (*this->logFile) << "log-likelihood computed by marginalization after EM iteration " << iter << " is " << setprecision(ll_precision) << this->logLikelihood << endl;
	// (*this->logFile) << "- + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + -" << endl;
	
	
	return tuple<int,double,double,double,double>(iter,logLikelihood_hss,logLikelihood_exp_data_first,logLikelihood_exp_data_final,this->logLikelihood);
}


/*
1. Number of iterations until EM converges
2. Initial log-likelihood with parsimony based parameters
3. Expected-data log-likelihood after one EM iteration
4. Expected-data log-likelihood after last EM iteration
5. Final log-likelihood using EM based parameters
*/

tuple <int,double,double,double,double> SEM::EM_started_with_dirichlet_rooted_at(SEM_vertex *v) {
	//	cout << "10a" << endl;
	// iterate over each internal node	
	this->RootTreeAtVertex(v); // take outside
	// cout << "10c" << endl;
	this->SetInitialEstimateOfModelParametersUsingDirichlet();
	this->StoreParamsInEMCurrent("init"); // take outside
	this->ComputeLogLikelihood();
	// cout << "Initial value of log-likelihood is " << setprecision(ll_precision) << this->logLikelihood << endl;
	// cout << "10d" << endl;
	this->ResetAncestralSequences();
	double logLikelihood_diri;
	double logLikelihood_exp_data_current;
	double logLikelihood_exp_data_first;
	double logLikelihood_exp_data_final;
	int iter = 0;	
	bool continueIterations = 1;
	this->debug = 0;
	bool verbose = 0;
	logLikelihood_diri = this->logLikelihood;
	this->EM_current.ll_init = logLikelihood_diri;	
	// (*this->logFile) << "-    -     -     -     -     -     -     -     -     -     -     -     -     -" << endl;
	// (*this->logFile) << "log-likelihood computed by marginalization using parsimony parameters is " << setprecision(ll_precision) << logLikelihood_diri << endl;
	logLikelihood_exp_data_current = -1 * pow(10,4);
	while (continueIterations) {
		// t_start_time = chrono::high_resolution_clock::now();
		iter += 1;		
		this->ConstructCliqueTree();			
		// 2. Compute expected counts
		this->ComputeExpectedCounts();

		this->ComputeMarginalProbabilitiesUsingExpectedCounts();
		
		this->ComputeMLEstimateOfBHGivenExpectedDataCompletion();
		
		this->ComputeLogLikelihoodUsingExpectedDataCompletion();
		
		// cout << "log-likelihood computed using expected counts after EM iteration " << iter << " is " << setprecision(ll_precision) << this->logLikelihood << endl;
		// (*this->logFile) << "log-likelihood computed using expected counts after EM iteration " << iter << " is " << setprecision(ll_precision) << this->logLikelihood << endl;		
		if (iter == 1){
			logLikelihood_exp_data_first = this->logLikelihood;
			logLikelihood_exp_data_current = this->logLikelihood;			
			this->EM_current.ecd_ll_per_iter[iter] = this->logLikelihood;
		} else if ((this->logLikelihood > logLikelihood_exp_data_current + this->ecdllConvergenceThreshold) and (iter < this->maxIter)) {
			logLikelihood_exp_data_current = this->logLikelihood;
			this->EM_current.ecd_ll_per_iter[iter] = this->logLikelihood;
		} else {
			continueIterations = 0;
			logLikelihood_exp_data_final = logLikelihood_exp_data_current;
		}
	}
	this->ComputeLogLikelihood();
	this->StoreParamsInEMCurrent("final"); // switch to three regime
	this->EM_current.iter = iter;
	this->EM_current.ll_final = this->logLikelihood;	
	// cout << "log-likelihood computed by marginalization using EM parameters " << iter << " is " << setprecision(ll_precision) << this->logLikelihood << endl;
	// cout << "- + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + -" << endl;
	// (*this->logFile) << "log-likelihood computed by marginalization using EM parameters " << iter << " is " << setprecision(ll_precision) << this->logLikelihood << endl;
	// (*this->logFile) << "- + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + -" << endl;
	
	return tuple<int,double,double,double,double>(iter,logLikelihood_diri,logLikelihood_exp_data_first,logLikelihood_exp_data_final,this->logLikelihood);
}

/*
1. Number of iterations until EM converges
2. Initial log-likelihood with parsimony based parameters
3. Expected-data log-likelihood after one EM iteration
4. Expected-data log-likelihood after last EM iteration
5. Final log-likelihood using EM based parameters
*/

// tuple <int,double,double,double,double> SEM::EM_started_with_parsimony_rooted_at(SEM_vertex *v) {
// 	// iterate over each internal node	
// 	this->RootTreeAtVertex(v);
// 	// cout << "10a" << endl;
// 	// select sites that don't have gap
// 	this->ComputeMPEstimateOfAncestralSequences(); // iterate over all sites including gaps	
// 	// cout << "10b" << endl;
// 	this->ComputeInitialEstimateOfModelParameters(); // skip sites that have gaps
// 	// cout << "10c" << endl;
// 	this->StoreParamsInEMCurrent("init");
// 	// cout << "10d" << endl;
// 	// cout << this->EM_current.root_prob_init[0] << "\t" << this->EM_current.root_prob_init[1] 
// 	// << "\t" << this->EM_current.root_prob_init[2] 
// 	// << "\t" << this->EM_current.root_prob_init[3] << endl;
// 	this->ComputeLogLikelihood(); // set conditional likelihood to 1 for all char having a gap	
// 	// cout << "Initial value of log-likelihood is " << setprecision(ll_precision) << this->logLikelihood << endl;
// 	// cout << "10e" << endl;
// 	this->ResetAncestralSequences();
// 	// cout << "10f" << endl;
// 	double logLikelihood_pars;
// 	double logLikelihood_exp_data_current;
// 	double logLikelihood_exp_data_first;
// 	double logLikelihood_exp_data_final;
// 	int iter = 0;	
// 	bool continueIterations = 1;
// 	this->debug = 0;
// 	bool verbose = 0;
// 	logLikelihood_pars = this->logLikelihood;
// 	this->EM_current.ll_init = logLikelihood_pars;
// 	// cout << "-    -     -     -     -     -     -     -     -     -     -     -     -     -" << endl;
// 	// cout << "log-likelihood computed by marginalization using parsimony parameters is " << setprecision(ll_precision) << logLikelihood_pars << endl;
// 	// (*this->logFile) << "-    -     -     -     -     -     -     -     -     -     -     -     -     -" << endl;
// 	// (*this->logFile) << "log-likelihood computed by marginalization using parsimony parameters is " << setprecision(ll_precision) << logLikelihood_pars << endl;
// 	logLikelihood_exp_data_current = -1 * pow(10,4);
// 	while (continueIterations) {
// 		// t_start_time = chrono::high_resolution_clock::now();
// 		iter += 1;
// 		// if (verbose) {
// 		// 	cout << "Iteration no. " << iter << endl;
// 		// 	(*this->logFile) << "Iteration no. " << iter << endl;
// 		// }
// 		// 1. Construct clique tree		
// 		// if (verbose) {
// 		// 	cout << "Construct clique tree" << endl;
// 		// 	(*this->logFile) << "Iteration no. " << iter << endl;
			
// 		// }		
		
// 		this->ConstructCliqueTree();
// 		// cout << "10g" << endl;
					
// 		// 2. Compute expected counts
// 		this->ComputeExpectedCounts();
// 		// cout << "10h" << endl;

// 		this->ComputePosteriorProbabilitiesUsingExpectedCounts();
// 		// cout << "10i" << endl;

		
// 		// 3. Optimize model parameters
// 		// if (verbose) {
// 		// 	cout << "Optimize model parameters given expected counts" << endl;
// 		// }
		
// 		// cout << "10j" << endl;
		
// 		this->ComputeLogLikelihoodUsingExpectedDataCompletion();
// 		// cout << "10k" << endl;
		
// 		// cout << "log-likelihood computed using expected counts after EM iteration " << iter << " is " << setprecision(ll_precision) << this->logLikelihood << endl;
// 		// (*this->logFile) << "log-likelihood computed using expected counts after EM iteration " << iter << " is " << setprecision(ll_precision) << this->logLikelihood << endl;
// 		this->EM_current.ecd_ll_per_iter[iter] = this->logLikelihood;
// 		if (iter == 1){
// 			// cout << "10l" << endl;
// 			logLikelihood_exp_data_first = this->logLikelihood;
// 			logLikelihood_exp_data_current = this->logLikelihood;			
// 		} else if ((this->logLikelihood > logLikelihood_exp_data_current + this->ecdllConvergenceThreshold) and (iter < this->maxIter)) {
// 			logLikelihood_exp_data_current = this->logLikelihood;
// 			// cout << "10m" << endl;
// 		} else {
// 			continueIterations = 0;
// 			logLikelihood_exp_data_final = logLikelihood_exp_data_current;
// 			// cout << "10n" << endl;
// 		}
// 	}
// 	this->ComputeLogLikelihood();
// 	// cout << "10o" << endl;
// 	this->StoreParamsInEMCurrent("final");
// 	// cout << "10p" << endl;
// 	this->EM_current.iter = iter;
// 	this->EM_current.ll_final = this->logLikelihood;	
// 	// cout << "log-likelihood computed by marginalization using EM parameters " << iter << " is " << setprecision(ll_precision) << this->logLikelihood << endl;
// 	// cout << "- + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + -" << endl;
// 	// (*this->logFile) << "log-likelihood computed by marginalization using EM parameters " << iter << " is " << setprecision(ll_precision) << this->logLikelihood << endl;
// 	// (*this->logFile) << "- + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + -" << endl;
// 	// cout << "10q" << endl;
// 	return tuple<int,double,double,double,double>(iter,logLikelihood_pars,logLikelihood_exp_data_first,logLikelihood_exp_data_final,this->logLikelihood);
// }


void SEM::ComputeSumOfExpectedLogLikelihoods() {
	this->sumOfExpectedLogLikelihoods = 0;
	this->sumOfExpectedLogLikelihoods += this->root->vertexLogLikelihood;
	double edgeLogLikelihood;	
	for (pair <SEM_vertex *, SEM_vertex *> edge : this->edgesForPostOrderTreeTraversal) {
		if (this->edgeLogLikelihoodsMap.find(edge) == this->edgeLogLikelihoodsMap.end()) {
//			cout << edge.first->name << "\t" << edge.second->name << endl;
		} else {
//			cout << edge.first->name << "\t" << edge.second->name << endl;
			edgeLogLikelihood = this->edgeLogLikelihoodsMap[edge];
			this->sumOfExpectedLogLikelihoods += edgeLogLikelihood;
		}				
	}
}

void SEM::ComputeMLRootedTreeForRootSearchUnderBH() {
	vector < SEM_vertex *> verticesToVisit = this->preOrderVerticesWithoutLeaves;	
	double logLikelihood_max = 0;
	int numberOfVerticesVisited = 0;	
	for (SEM_vertex * v : verticesToVisit) {
		numberOfVerticesVisited += 1;
		this->RootTreeAtVertex(v);		
		this->ComputeMLEstimateOfBHGivenExpectedDataCompletion();
		this->ComputeLogLikelihoodUsingExpectedDataCompletion();
		if ((numberOfVerticesVisited < 2) or (logLikelihood_max < this->logLikelihood)) {
			logLikelihood_max = this->logLikelihood;
			this->StoreRootAndRootProbability();
			this->StoreTransitionMatrices();			
			this->StoreDirectedEdgeList();
		}
	}
	this->RestoreRootAndRootProbability();
	this->RestoreTransitionMatrices();
	this->RestoreDirectedEdgeList();
	this->SetEdgesForTreeTraversalOperations();
	this->logLikelihood = logLikelihood_max;
}


void SEM::ComputeMLRootedTreeForFullStructureSearch() {
	this->StoreEdgeListForChowLiuTree();
	// For each vertex v of the Chow-Liu tree
	SEM_vertex * v;
	double logLikelihood_max = 0;
	int verticesTried = 0;
	bool debug = 0;
	bool useExpectedLogLikForSelectingRoot = 1;
	string nonCanonicalRootedTreeFileName = "/home/pk/Projects/EMTRasedForests/data/trees/nonCanonicalRootedTree_test_numberOfLeaves_16_replicate_1";
	string canonicalRootedTreeFileName = "/home/pk/Projects/EMTRasedForests/data/trees/canonicalRootedTree_test_numberOfLeaves_16_replicate_1";
	for (pair <int, SEM_vertex *> idPtrPair : * this->vertexMap) {
		verticesTried += 1;
		// 	Root tree at v
		v = idPtrPair.second;
		// if (debug) {
		// 	cout << "Rooting tree at vertex " << v->name << endl;
		// }
		this->RootTreeAtVertex(v);		
		// if (debug) {
		// 	this->WriteRootedTreeAsEdgeList(nonCanonicalRootedTreeFileName);
		// 	// cout << "Is v an observed variable?" << endl;
		// 	if (v->observed) {
		// 		// cout << "Yes" << endl;
		// 	} else {
		// 		cout << "No" << endl;
		// 	}
		// 	cout << "Root name is " << this->root->name << endl;
		// }
		// Compute MLE of model parameters
		// assuming that posterior probabilities
		// P(V) (for each vertex) and 
		// P(V1, V2) (for each vertex pair) are available
//		cout << "Computing MLE of model parameters" << endl;
		this->ComputeMLEstimateOfBHGivenExpectedDataCompletion();
		// Transform to bifurcating rooted tree
//		cout << "Transforming to bifurcating leaf-labeled tree" << endl;		
		if (useExpectedLogLikForSelectingRoot) {
			this->ComputeLogLikelihoodUsingExpectedDataCompletion();
		} else {
			this->TransformRootedTreeToBifurcatingTree();
			if (debug) {
				this->WriteRootedTreeAsEdgeList(canonicalRootedTreeFileName);
			}	
			// Compute loglikelihood
			this->SetLeaves();
			this->SetEdgesForPostOrderTraversal();
			this->ComputeLogLikelihood();
		}
		if (logLikelihood_max < this->logLikelihood or verticesTried < 2) {
			logLikelihood_max = this->logLikelihood;
		//	cout << "Current max loglikelihood is " << logLikelihood_max << endl;
			// Store root probability that maximizes loglikelihood
			this->StoreRootAndRootProbability();
			// Store transition matrices that maximize loglikelihood
			this->StoreTransitionMatrices();
			// Store directed edge list rooted tree which maximizes loglikelihood
			this->StoreDirectedEdgeList();
		}
	//		this->RestoreEdgeListForChowLiuTree();
	}
	// Select bifurcating rooted tree and parameters that maximize loglikelihood	
	this->RestoreRootAndRootProbability();
	this->RestoreTransitionMatrices();
	this->RestoreDirectedEdgeList();
	if (useExpectedLogLikForSelectingRoot) {
		this->TransformRootedTreeToBifurcatingTree();
	}
	this->SetLeaves();
	this->SetEdgesForPostOrderTraversal();
	this->SetEdgesForPreOrderTraversal();
	this->SetVerticesForPreOrderTraversalWithoutLeaves();
//	this->ComputeLogLikelihood();
//	Following step computes MLE of parameters of general Markov model	
//	this->logLikelihood = logLikelihood_max;
//	cout << "Current max expected logLikelihood is " << this->logLikelihood << endl;
//	this->ComputeMLEstimatesViaHardEM();
//	cout << "Current logLikelihood is " << this->logLikelihood << endl;
}

emtr::Md SEM::GetP_yGivenx(emtr::Md P_xy) {
	emtr::Md P_yGivenx = emtr::Md{};
	array <double, 4> P_x;
	for (int dna_x = 0; dna_x < 4; dna_x ++) {
		P_x[dna_x] = 0;
		for (int dna_y = 0; dna_y < 4; dna_y ++) {
			P_x[dna_x] += P_xy[dna_x][dna_y];
		}
	}
	for (int dna_x = 0; dna_x < 4; dna_x ++) {
		for (int dna_y = 0; dna_y < 4; dna_y ++) {
			// P_yGivenx[dna_x][dna_y] = P_xy[dna_x][dna_y] / P_x[dna_x];
			if (P_xy[dna_x][dna_y] == 0 && P_x[dna_x] == 0) {
				P_yGivenx[dna_x][dna_y] = 0;				
			} else {						
				P_yGivenx[dna_x][dna_y] = P_xy[dna_x][dna_y] / P_x[dna_x];
			}						
			if (isnan(P_yGivenx[dna_x][dna_y])) {
				cout << "P_yGivenx[" << dna_x << "][" << dna_y << "]" << " is nan " << endl;
				cout << "P_xy[" << dna_x << "][" << dna_y << "]" << " is " << P_xy[dna_x][dna_y]  << endl;
				cout << "P_x[" << dna_x << "] is " << P_x[dna_x] << endl;				
				throw mt_error ("nan found");
			}
		}
	}
	return (P_yGivenx);
}


void SEM::SetMinLengthOfEdges() {	
	for (pair<pair<SEM_vertex * , SEM_vertex * >,double> edgeAndLengthsPair: this->edgeLengths){		
		if (edgeAndLengthsPair.second < pow(10,-7)){
			this->edgeLengths[edgeAndLengthsPair.first]  = pow(10,-7);
		}
	}
}

 
void SEM::ComputeMLEstimateOfBHGivenExpectedDataCompletion() {	
	SEM_vertex * x; SEM_vertex * y;
	emtr::Md P_xy;
	for (pair <int, SEM_vertex*> idPtrPair : *this->vertexMap) {
		y = idPtrPair.second;
		x = y->parent;
		if (x != y) {
			if (x->id < y->id) {
				P_xy = this->posteriorProbabilityForVertexPair[make_pair(x,y)];
			} else {
				// Check following step				
				P_xy = emtr::MT(this->posteriorProbabilityForVertexPair[make_pair(y,x)]);
				
			}
			// MLE of transition matrices
			// cout << "eight eight" << endl;
			y->transitionMatrix = this->GetP_yGivenx(P_xy);
			// cout << "nine nine" << endl;
		} else {			
			// MLE of root probability
			this->rootProbability = this->posteriorProbabilityForVertex[y];
			y->rootProbability = this->rootProbability;
			y->transitionMatrix = emtr::Md{};
			for (int i = 0; i < 4; i ++) {
				y->transitionMatrix[i][i] = 1.0;
			}
			// cout << "ten ten" << endl;
		}
	}
}




void SEM::StoreEdgeListForChowLiuTree() {
	this->edgesForChowLiuTree.clear();	
	SEM_vertex * v;
	for(pair <int, SEM_vertex*> idPtrPair : *this->vertexMap) {
		v = idPtrPair.second;
		for (SEM_vertex * n : v->neighbors) {
			if (v->id < n->id) {
				this->edgesForChowLiuTree.push_back(make_pair(v,n));
			} else {
				this->edgesForChowLiuTree.push_back(make_pair(n,v));
			}
		}
	}
}

void SEM::RestoreEdgeListForChowLiuTree() {
	this->ClearDirectedEdges();
	SEM_vertex *u; SEM_vertex *v; 
	for(pair<SEM_vertex *,SEM_vertex*> edge : this->edgesForChowLiuTree){
		tie (u, v) = edge;
		u->AddNeighbor(v);
		v->AddNeighbor(u);
	}
}

void SEM::StoreDirectedEdgeList() {
	SEM_vertex * p;
	SEM_vertex * c;
	this->directedEdgeList.clear();
	for(pair <int, SEM_vertex*> idPtrPair : *this->vertexMap){
		c = idPtrPair.second;
		if (c->parent != c) {			
			p = c->parent;
			this->directedEdgeList.push_back(make_pair(p,c));
		}
	}
}

void SEM::RestoreDirectedEdgeList() {
	this->ClearDirectedEdges();
	SEM_vertex * p;
	SEM_vertex * c;
	for (pair <SEM_vertex *, SEM_vertex *> edge : this->directedEdgeList) {
		tie (p, c) = edge;
		c->AddParent(p);
		p->AddChild(c);		
	}
}

void SEM::RootTreeAtVertex(SEM_vertex* r) {	
	this->ClearDirectedEdges();
	vector <SEM_vertex*> verticesToVisit;
	vector <SEM_vertex*> verticesVisited;
	SEM_vertex * p;	
	verticesToVisit.push_back(r);	
	int numberOfVerticesToVisit = verticesToVisit.size();
	while (numberOfVerticesToVisit > 0) {
		p = verticesToVisit[numberOfVerticesToVisit-1];
		verticesToVisit.pop_back();
		verticesVisited.push_back(p);
		numberOfVerticesToVisit -= 1;
		for (SEM_vertex* c : p->neighbors) {
			if (find(verticesVisited.begin(),verticesVisited.end(),c)==verticesVisited.end()) {
				p->AddChild(c);
				c->AddParent(p);				
				verticesToVisit.push_back(c);
				numberOfVerticesToVisit += 1;
			}
		}
	}
	this->root = r;
	this->SetEdgesForTreeTraversalOperations();
}

void SEM::TransformRootedTreeToBifurcatingTree() {
	bool containsMatchingVertex;
	bool containsSingletonHiddenVertex;
	bool debug = 0;
	SEM_vertex * matchingVertex;
	SEM_vertex * p; SEM_vertex * c; SEM_vertex * h;
	SEM_vertex * o; SEM_vertex * h_s;
	SEM_vertex * c_1; SEM_vertex * c_2; 
	emtr::Md P;
	array <double, 4> rootProb_orig;
	array <double, 4> rootProb_new;
	vector <SEM_vertex *> childrenToSwap;
	bool checkForCanonicalForm;
	checkForCanonicalForm = IsTreeInCanonicalForm();
	string nonCanonicalRootedTreeFileName = "/home/pk/Projects/EMTRasedForests/data/trees/debugNonCanonicalRootedTree_before";
	this->WriteRootedTreeAsEdgeList(nonCanonicalRootedTreeFileName);
	int numberOfTransformations = 0;
	while (!checkForCanonicalForm) {
		numberOfTransformations += 1;
		if (numberOfTransformations > 10000) {
			cout << "Check Transformation of rooted tree to canonical form" << endl;
            throw mt_error("Check Transformation of rooted tree to canonical form");            
		}
		// Case 1. x->h
		// Check for hidden vertices that are leaves
		// Remove the arc x->h. This creates a singleton vertex h.
		tie (containsMatchingVertex, h) = this->CheckAndRetrieveHiddenVertexWithOutDegreeZeroAndInDegreeOne();		
		while (containsMatchingVertex) {
			// if (debug and containsMatchingVertex) {
			// 	cout << "Case 1. There is a hidden vertex that is a leaf" << endl;
			// }
			this->RemoveArc(h->parent,h);
			tie (containsMatchingVertex, h) = this->CheckAndRetrieveHiddenVertexWithOutDegreeZeroAndInDegreeOne();
		}
		// Case 2. p->h->c
		// Check for hidden vertices h with out degree 1 and in degree 1 
		// Remove p->h and h->c and add p->c
		tie (containsMatchingVertex, h) = this->CheckAndRetrieveHiddenVertexWithOutDegreeOneAndInDegreeOne();
		while (containsMatchingVertex) {
			// if (debug and containsMatchingVertex) {
			// 	cout << "Case 2. There is a non-root hidden vertex that has out degree 1" << endl;
			// }			
			p = h->parent;
			if (h->children.size() != 1){
                throw mt_error("Check case 2");
            }
			c = h->children[0];
			c->transitionMatrix = emtr::MM(h->transitionMatrix,c->transitionMatrix);
			h->transitionMatrix = this->I4by4;
//			cout << "Removing edge (p,h)" << endl;
			this->RemoveArc(p,h);
//			cout << "Removing edge (h,c)" << endl;
			this->RemoveArc(h,c);			
//			cout << "Adding edge (p,c)" << endl;			
			this->AddArc(p,c);
//			cout << "Edge (p,c) added" << endl;			
			tie (containsMatchingVertex, h) = this->CheckAndRetrieveHiddenVertexWithOutDegreeOneAndInDegreeOne();
		}
		// Case 3. The root is a hidden vertex with out degree 1
		if (!this->root->observed and (this->root->outDegree == 1)) {
			// if (debug) {
			// 	cout << "Case 3. The root is a hidden vertex with out degree 1" << endl;
			// }
			rootProb_orig = this->rootProbability;
			if (this->root->children.size()!=1){
                throw mt_error("Check case 3");
            }
			p = this->root;
			c = this->root->children[0];
			P = c->transitionMatrix;
			for (int y = 0; y < 4; y ++) {				
				rootProb_new[y] = 0;
				for (int x = 0; x < 4; x ++){
					rootProb_new[y] += rootProb_orig[x] * P[x][y];
				}
			}
			this->RemoveArc(p,c);
			this->root = c;
			this->rootProbability = rootProb_new;
			this->root->rootProbability = this->rootProbability;
			c->transitionMatrix = this->I4by4;
		}
		// Case 4. The root is an observed vertex
		if (this->root->observed) {
			// if (debug) {
			// 	cout << "Case 4. The root is an observed vertex" << endl;
			// }
			tie (containsSingletonHiddenVertex, h_s) = this->CheckAndRetrieveSingletonHiddenVertex();
			if(!containsSingletonHiddenVertex){
                throw mt_error("Check case 4");
            }
            
			if(h_s->children.size() != 0) {
                throw mt_error("Check case 4");
            }
			p = h_s;
			c = this->root;
			childrenToSwap = c->children;
			for (SEM_vertex * child : childrenToSwap) {
				this->RemoveArc(c, child);
				this->AddArc(p, child);
			}
			this->root = p;
			this->root->rootProbability = this->rootProbability;
			this->AddArc(p,c);
			c->transitionMatrix = this->I4by4;
		}			
		// Case 5. p->o, o->c1, ..., o->ck
		// Check for non-leaf non-root observed vertex
		tie (containsMatchingVertex, matchingVertex) = this->CheckAndRetrieveObservedVertexThatIsNotALeafAndIsNotTheRoot();
		tie (containsSingletonHiddenVertex, h_s) = this->CheckAndRetrieveSingletonHiddenVertex();
		while (containsMatchingVertex and containsSingletonHiddenVertex) {
			// if (debug) {
			// 	cout << "Case 5. There is a non-leaf non-root observed vertex" << endl;
			// }
			o = matchingVertex;
//			cout << o->name  << endl;
			// Swap children of o and h
			childrenToSwap = o->children;
			for (SEM_vertex * c: childrenToSwap){					
				this->RemoveArc(o,c);
				this->AddArc(h_s,c);
			}
			// Set parent of h to parent of o
			this->AddArc(o->parent,h_s);
			// Set P(h|p) to P(o|p)
			h_s->transitionMatrix = o->transitionMatrix;
			this->RemoveArc(o->parent,o);
			this->AddArc(h_s,o);
			// Set P(o|h) to Identity matrix I4by4
			o->transitionMatrix = this->I4by4;						
			tie (containsMatchingVertex, o) = this->CheckAndRetrieveObservedVertexThatIsNotALeafAndIsNotTheRoot();
			tie (containsSingletonHiddenVertex, h_s) = this->CheckAndRetrieveSingletonHiddenVertex();
		}
		// Case 6. p->h, h->c1, h->c2, ... h->ck
		// Check for a hidden vertex h with outdegree greater than two
		// and for a singleton hidden vertex h_s
		tie (containsMatchingVertex, h) = this->CheckAndRetrieveHiddenVertexWithOutDegreeGreaterThanTwo();
		tie (containsSingletonHiddenVertex, h_s) = this->CheckAndRetrieveSingletonHiddenVertex();			
		while (containsMatchingVertex and containsSingletonHiddenVertex) {			
			// if (debug) {
			// 	cout << "Case 6. There is a multifurcation" << endl;
			// }
			childrenToSwap = h->children;
			sort(childrenToSwap.begin(),childrenToSwap.end(), [](SEM_vertex * u, SEM_vertex * v) {
				return u->id < v->id;
			});
			// Select children c_1 and c_2 are by sorting the children of h in ascending order of id
			// and selecting the first two chilren in the sorted list
			c_1 = childrenToSwap[0];
			c_2 = childrenToSwap[1];
			// Remove children c_1 and c_2 of h and set them as children of h_s
			this->RemoveArc(h,c_1);
			this->RemoveArc(h,c_2);
			this->AddArc(h_s,c_1);
			this->AddArc(h_s,c_2);
			// Set h as the parent of h_s
			this->AddArc(h,h_s);
			h_s->transitionMatrix = this->I4by4;
			tie (containsMatchingVertex, h) = this->CheckAndRetrieveHiddenVertexWithOutDegreeGreaterThanTwo();
			tie (containsSingletonHiddenVertex, h_s) = this->CheckAndRetrieveSingletonHiddenVertex();								
		}
		checkForCanonicalForm = IsTreeInCanonicalForm();
	}
}

void SEM::ClearDirectedEdges() {
	// cout << "Resetting times visited " << endl;
	this->ResetTimesVisited();
	SEM_vertex * v;
	for (pair <int, SEM_vertex *> idPtrPair : * this->vertexMap) {
		v = idPtrPair.second;
		v->parent = v;
		v->children.clear();		
		v->inDegree = 0;
		v->outDegree = 0;		
	}
}

void SEM::ClearUndirectedEdges() {
	SEM_vertex * v;
	for (pair <int, SEM_vertex *> idPtrPair : * this->vertexMap) {
		v = idPtrPair.second;
		v->degree = 0;
		v->neighbors.clear();
	}
}

void SEM::ClearAllEdges() {
	this->ResetTimesVisited();
	SEM_vertex * v;
	for (pair <int, SEM_vertex *> idPtrPair : * this->vertexMap) {
		v = idPtrPair.second;
		v->parent = v;
		v->children.clear();
		v->neighbors.clear();
		v->degree = 0;
		v->inDegree = 0;
		v->outDegree = 0;	
	}
}

void SEM::ResetTimesVisited() {
	SEM_vertex * v;
	for (pair <int, SEM_vertex *> idPtrPair : * this->vertexMap) {
		v = idPtrPair.second;		
		v->timesVisited = 0;		
	}
}

array <double, 4> SEM::GetBaseCompositionForTrifle(SEM_vertex * v, int layer) {
	array <double, 4> baseCompositionArray;
	for (int dna = 0; dna < 4; dna ++) {
		baseCompositionArray[dna] = 0;
	}
	int dna_v;
	// cout << "here 12a" << endl;
	for (int site = 0; site < this->num_patterns_for_layer[layer]; site ++) {
		dna_v = v->DNAcompressed[site];
		if (dna_v > -1) baseCompositionArray[dna_v] += this->DNAPatternWeights[site];
	}
	// cout << "here 12b" << endl;
	int non_gap_sites = 0;
	for (int dna = 0; dna < 4; dna ++) {
		non_gap_sites += baseCompositionArray[dna];
	}
	// cout << "here 12c" << endl;
	for (int dna = 0; dna < 4; dna ++) {
		baseCompositionArray[dna] /= float(non_gap_sites);		
	}	

	return (baseCompositionArray);
}

array <double, 4> SEM::GetBaseComposition(SEM_vertex * v) {
	array <double, 4> baseCompositionArray;
	for (int dna = 0; dna < 4; dna ++) {
		baseCompositionArray[dna] = 0;
	}
	int dna_v;
	for (int site = 0; site < this->num_dna_patterns; site ++) {
		dna_v = v->DNAcompressed[site];
		if (dna_v > -1) baseCompositionArray[dna_v] += this->DNAPatternWeights[site];
	}
	int non_gap_sites = 0;
	for (int dna = 0; dna < 4; dna ++) {
		non_gap_sites += baseCompositionArray[dna_v];
	}
	for (int dna = 0; dna < 4; dna ++) {
		baseCompositionArray[dna] /= float(non_gap_sites);
	}
	return (baseCompositionArray);
}

array <double, 4> SEM::GetObservedCountsForVariable(SEM_vertex * v) {
	array <double, 4> observedCounts;
	for (int i = 0; i < 4; i++) {
		observedCounts[i] = 0;
 	}
	for (int site = 0; site < this->num_dna_patterns; site ++) {
		if (v->DNArecoded[site] < 4) { // FIX_AMB
			observedCounts[v->DNArecoded[site]] += this->DNAPatternWeights[site];
		}		
	}
	return (observedCounts);
}


void SEM::ComputeMLEOfRootProbability() {
	this->rootProbability = GetBaseComposition(this->root);
	this->root->rootProbability = this->rootProbability;
}

void SEM::ComputeMLEOfTransitionMatrices() {
	SEM_vertex * c; SEM_vertex * p;
	bool debug = 0;
	for (pair <int, SEM_vertex *> idPtrPair : * this->vertexMap){
		c = idPtrPair.second;
		p = c->parent;
		if (p != c) {
			c->transitionMatrix = this->ComputeTransitionMatrixUsingAncestralStates(p,c);
			// if (debug) {
			// 	cout << "Estimated transition matrix is" << endl;
			// 	cout << c->transitionMatrix << endl;
			// }
		}
	}
	
}

void SEM::SetInitialEstimateOfModelParametersUsingDirichlet() {

	array <double, 4> pi_diri = SEM::sample_pi();

	this->rootProbability = pi_diri;
	this->root->rootProbability = pi_diri;

	SEM_vertex * p; SEM_vertex * c; emtr::Md M_pc;
	for (pair<SEM_vertex*,SEM_vertex*> edge : this->edgesForPostOrderTreeTraversal) {		
		p = edge.first;
		c = edge.second;

		for (int row = 0; row < 4; row++) {
			array <double, 4> M_row_diri = SEM::sample_M_row();
			int diri_index = 1;
			for (int col = 0; col < 4; col++) {
				if (row == col){
					M_pc[row][col] = M_row_diri[0];
				} else {
					M_pc[row][col] = M_row_diri[diri_index++];
				}
			}
		}		
		c->transitionMatrix = M_pc;
	}
}

void SEM::StoreInitialParamsInEMTrifleCurrent(int layer) {	
	SEM_vertex * c;		
	for (int dna = 0; dna < 4; dna ++) this->EMTrifle_current.root_prob_initial_trifle[layer][dna] = this->rootProbability[dna];
	for (int dna = 0; dna < 4; dna ++) {
		if(this->rootProbability[dna] == 0) {
			cout << "root prob zero";
			exit(-1);
		}

	}
	for (pair <int, SEM_vertex *> idPtrPair : * this->vertexMap) {
		c = idPtrPair.second;
		this->EMTrifle_current.trans_prob_initial_trifle[layer][c->name] = c->transitionMatrix;
	}
}

void SEM::StoreFinalParamsInEMTrifleCurrent(int layer) {	
	SEM_vertex * c;		
	for (int dna = 0; dna < 4; dna ++) this->EMTrifle_current.root_prob_final_trifle[layer][dna] = this->rootProbability[dna];
	for (pair <int, SEM_vertex *> idPtrPair : * this->vertexMap) {
		c = idPtrPair.second;
		this->EMTrifle_current.trans_prob_final_trifle[layer][c->name] = c->transitionMatrix;
	}
}

// void SEM::StoreParamsInEMTrifleCurrent(int layer) {	
// 	SEM_vertex * c;	
// 	if (layer == -1) {
// 		for (int dna = 0; dna < 4; dna ++) this->EMTrifle_current.root_prob_init[dna] = this->rootProbability[dna];
// 		for (pair <int, SEM_vertex *> idPtrPair : * this->vertexMap) {
// 			c = idPtrPair.second;
// 			this->EMTrifle_current.trans_prob_init[c->name] =  c->transitionMatrix;
// 		}
// 	} else if (layer > -1 && layer < 3) {
// 		for (int dna = 0; dna < 4; dna ++) this->EMTrifle_current.root_prob_final_trifle[layer][dna] = this->rootProbability[dna];
// 		for (pair <int, SEM_vertex *> idPtrPair : * this->vertexMap) {
// 			c = idPtrPair.second;
// 			this->EMTrifle_current.trans_prob_final_trifle[layer][c->name] = c->transitionMatrix;
// 		}
// 	} else {
// 		throw mt_error("argument not recognized");
// 	}
// }

void SEM::StoreParamsInEMCurrent(string init_or_final) {
	for (int dna = 0; dna < 4; dna ++) {
		if (init_or_final == "init") {			
			this->EM_current.root_prob_init[dna] = this->rootProbability[dna];			
		} else if (init_or_final == "final") {
			this->EM_current.root_prob_final[dna] = this->rootProbability[dna];			
		} else {
			mt_error("argument not recognized");
		}		 	
	}

	SEM_vertex * c;
	for (pair <int, SEM_vertex *> idPtrPair : * this->vertexMap) {
		c = idPtrPair.second;		
		if (c->parent != c) {			
			if (init_or_final == "init") {
				this->EM_current.trans_prob_init[c->name] =  c->transitionMatrix;
			} else if (init_or_final == "final") {
				this->EM_current.trans_prob_final[c->name] = c->transitionMatrix;
			} else {
				mt_error("argument not recognized");
			}
		}
	}
}

void SEM::SetInitialEstimateOfModelParametersUsingHSS() {
	
	// cout << "root is set at " << this->root->name << endl;
	// set root probability 
	for (int x = 0; x < 4; x++) {
		this->rootProbability[x] = this->root->root_prob_hss[x];
		this->root->rootProbability[x] = this->root->root_prob_hss[x];
		// cout << "root probability for " << x << " is " << this->rootProbability[x] << endl;
	}


	SEM_vertex * p; SEM_vertex * c; emtr::Md M_pc;
	for (pair<SEM_vertex*,SEM_vertex*> edge : this->edgesForPostOrderTreeTraversal) {		
		p = edge.first;
		c = edge.second;		
		M_pc = (*this->M_hss)[{p,c}];	
		c->transitionMatrix = M_pc;
	}
}

void SEM::BumpZeroEntriesOfModelParameters() {
	// characters in next layer not present in previous layer may have
	// zero probabilities 
	bool debug = 0;
	// cout << "here 11a" << endl;
	this->addToZeroEntriesOfRootProbability(this->rootProbability);
	this->root->rootProbability = this->rootProbability;
	for (int i = 0; i < 4; i++) {
		if(isnan(this->rootProbability[i])) {
			cout << "root prob not valid" << endl;
			throw mt_error("root prob not valid");
		}
	}
	if (debug) {
		cout << "Root probability is " << endl;
		for (int i = 0; i < 4; i++) {
			// cout << this->rootProbability[i] << "\t";
		}
		cout << endl;
	}
	// cout << "here 11b" << endl;
	SEM_vertex * c; SEM_vertex * p;
	for (pair <int, SEM_vertex *> idPtrPair : * this->vertexMap) {
		c = idPtrPair.second;
		p = c->parent;
		if (p != c) {			
			this->addToZeroEntriesOfTransitionMatrix(c->transitionMatrix);			
			if (debug) {
				cout << "Transition matrix for " << p->name << " to " << c->name << " is " << endl;
				// cout << c->transitionMatrix[2][3] << endl;
			}			
		}
	}
}

void SEM::ComputeInitialEstimateOfModelParametersForTrifle(int layer) {
	bool debug = 0;
	// cout << "here 11a" << endl;
	this->rootProbability = GetBaseCompositionForTrifle(this->root, layer);
	this->root->rootProbability = this->rootProbability;
	for (int i = 0; i < 4; i++) {
		if(isnan(this->rootProbability[i])) {
			cout << "root prob not valid" << endl;
			throw mt_error("root prob not valid");
		}
	}
	if (debug) {
		cout << "Root probability is " << endl;
		for (int i = 0; i < 4; i++) {
			// cout << this->rootProbability[i] << "\t";
		}
		cout << endl;
	}
	// cout << "here 11b" << endl;
	SEM_vertex * c; SEM_vertex * p;
	for (pair <int, SEM_vertex *> idPtrPair : * this->vertexMap) {
		c = idPtrPair.second;
		p = c->parent;
		if (p != c) {
			c->transitionMatrix = this->ComputeTransitionMatrixUsingAncestralStatesForTrifle(p,c,layer);
			
			if (debug) {
				cout << "Transition matrix for " << p->name << " to " << c->name << " is " << endl;
				// cout << c->transitionMatrix[2][3] << endl;
			}			
		}
	}
	// cout << "here 11c" << endl;
	if (debug) {
		cout << "Transition matrices have been computed" << endl;
	}	
}

void SEM::ComputeInitialEstimateOfModelParameters() {
	bool debug = 0;
	this->rootProbability = GetBaseComposition(this->root);	
	this->root->rootProbability = this->rootProbability;
	if (debug) {
		// cout << "Root probability is " << endl;
		// for (int i = 0; i < 4; i++) {
		// 	cout << this->rootProbability[i] << "\t";
		// }
		// cout << endl;
	}

	SEM_vertex * c; SEM_vertex * p;
	for (pair <int, SEM_vertex *> idPtrPair : * this->vertexMap){
		c = idPtrPair.second;
		p = c->parent;
		if (p != c) {
			c->transitionMatrix = this->ComputeTransitionMatrixUsingAncestralStates(p,c);		
			if (debug) {
				cout << "Transition matrix for " << p->name << " to " << c->name << " is " << endl;
				// cout << c->transitionMatrix[2][3] << endl;
			}			
		}
	}
	if (debug) {
		cout << "Transition matrices have been computed" << endl;
	}	
}


void SEM::ResetLogScalingFactors() {
	for (pair <int, SEM_vertex * > idPtrPair : * this->vertexMap){
		idPtrPair.second->logScalingFactors = 0;
	}
}

void SEM::ComputeMPEstimateOfAncestralSequencesForTrifle(int layer) {
	SEM_vertex * p;	
	map <SEM_vertex * , int> V;
	map <SEM_vertex * , vector<int>> VU;					
	map <int, int> dnaCount;
	int pos;
	int maxCount; int numberOfPossibleStates;
	this->ResetAncestralSequences();
	if (this->preOrderVerticesWithoutLeaves.size() == 0) {
		// cout << "setting vertices for preorder traversal" << endl;
		this->SetVerticesForPreOrderTraversalWithoutLeaves();
	} else {
		// cout << "vertices for preorder traversal already set" << endl;
	}
	
	// cout << "p 1" << endl;

	// cout << "site patterns for layer is " << this->num_patterns_for_layer[layer] << endl;
	for (int site = 0; site < this->num_patterns_for_layer[layer]; site++) {		
		V.clear();
		VU.clear();		
	//	Compute V and VU for leaves	
		// if (site==0) cout << "p 2" << endl;
		for (SEM_vertex * c : this->leaves) {			
			V.insert(make_pair(c,c->DNAcompressed[site]));
			vector <int> vectorToAdd;
			vectorToAdd.push_back(c->DNAcompressed[site]);
			VU.insert(make_pair(c,vectorToAdd));	
		}		
	//	Set VU for ancestors
		for (SEM_vertex * c : this->preOrderVerticesWithoutLeaves) {
			vector <int> vectorToAdd;
			VU.insert(make_pair(c,vectorToAdd));
		}				
		for (int p_ind = this->preOrderVerticesWithoutLeaves.size()-1; p_ind > -1; p_ind--) {			
			p = preOrderVerticesWithoutLeaves[p_ind];
			// cout << p->name << " rf \t" << endl;			
			dnaCount.clear();
			for (int dna = -1; dna < 4; dna++) {
				dnaCount[dna] = 0;
			}
			for (SEM_vertex * c : p->children) {
				for (int dna: VU[c]) {
					dnaCount[dna] += 1;
				}
			}
			maxCount = 0;
			for (pair <int, int> dnaCountPair: dnaCount) {
				if (dnaCountPair.second > maxCount) {
					maxCount = dnaCountPair.second;
				}
			}			
			for (pair <int, int> dnaCountPair: dnaCount) {
				if (dnaCountPair.second == maxCount) {
					VU[p].push_back(dnaCountPair.first);
				}
			}			
		}			
		
	// Set V for ancestors
		for (SEM_vertex * c : preOrderVerticesWithoutLeaves) {
			if (c->parent == c) {			
			// Set V for root
				if (VU[c].size()==1) {
					V.insert(make_pair(c,VU[c][0]));
					// cout << c->name << "\t" << VU[c][0] << endl;
				} else {
					numberOfPossibleStates = VU[c].size();
					uniform_int_distribution <int> distribution(0,numberOfPossibleStates-1);
					pos = distribution(generator);
					V.insert(make_pair(c,VU[c][pos]));
				}				
			} else {
				p = c->parent;
				if (find(VU[c].begin(),VU[c].end(),V[p])==VU[c].end()) {
					numberOfPossibleStates = VU[c].size();
					uniform_int_distribution <int> distribution(0,numberOfPossibleStates-1);
					pos = distribution(generator);
					V.insert(make_pair(c,VU[c][pos]));
					// cout << c->name << " ancestor \t" << VU[c][pos] << endl;
				} else {
					V.insert(make_pair(c,V[p]));
					// cout << c->name << " ancestor \t" << V[p] << endl;
				}				
			}				
			c->DNAcompressed[site] = V[c];
		}				
	}	
}

void SEM::ComputeMPEstimateOfAncestralSequences() {
	SEM_vertex * p;	
	map <SEM_vertex * , int> V;
	map <SEM_vertex * , vector<int>> VU;					
	map <int, int> dnaCount;
	int pos;
	int maxCount; int numberOfPossibleStates;
	this->ResetAncestralSequences();
	if (this->preOrderVerticesWithoutLeaves.size() == 0) {
		this->SetVerticesForPreOrderTraversalWithoutLeaves();
	}
		
	for (int site = 0; site < this->num_dna_patterns; site++) {		
		V.clear();
		VU.clear();		
	//	Compute V and VU for leaves		
		for (SEM_vertex * c : this->leaves) {			
			V.insert(make_pair(c,c->DNAcompressed[site]));
			vector <int> vectorToAdd;
			vectorToAdd.push_back(c->DNAcompressed[site]);			
			VU.insert(make_pair(c,vectorToAdd));			
		}		
	//	Set VU for ancestors
		for (SEM_vertex* c : this->preOrderVerticesWithoutLeaves) {
			vector <int> vectorToAdd;
			VU.insert(make_pair(c,vectorToAdd));
		}				
		for (int p_ind = this->preOrderVerticesWithoutLeaves.size()-1; p_ind > -1; p_ind--) {
			p = preOrderVerticesWithoutLeaves[p_ind];			
			dnaCount.clear();
			for (int dna = -1; dna < 4; dna++) {
				dnaCount[dna] = 0;
			}
			for (SEM_vertex * c : p->children) {
				for (int dna: VU[c]) {
					dnaCount[dna] += 1;
				}
			}
			maxCount = 0;
			for (pair <int, int> dnaCountPair: dnaCount) {
				if (dnaCountPair.second > maxCount) {
					maxCount = dnaCountPair.second;
				}
			}			
			for (pair <int, int> dnaCountPair: dnaCount) { 
				if (dnaCountPair.second == maxCount) {
					VU[p].push_back(dnaCountPair.first);					
				}
			}			
		}			
		
	// Set V for ancestors
		for (SEM_vertex * c : preOrderVerticesWithoutLeaves) {
			if (c->parent == c) {			
			// Set V for root
				if (VU[c].size()==1) {
					V.insert(make_pair(c,VU[c][0]));
				} else {
					numberOfPossibleStates = VU[c].size();
					uniform_int_distribution <int> distribution(0,numberOfPossibleStates-1);
					pos = distribution(generator);
					V.insert(make_pair(c,VU[c][pos]));
				}				
			} else {
				p = c->parent;
				if (find(VU[c].begin(),VU[c].end(),V[p])==VU[c].end()) {
					numberOfPossibleStates = VU[c].size();
					uniform_int_distribution <int> distribution(0,numberOfPossibleStates-1);
					pos = distribution(generator);
					V.insert(make_pair(c,VU[c][pos]));
				} else {
					V.insert(make_pair(c,V[p]));
				}				
			}				
			c->DNAcompressed[site] = V[c];
		}				
	}	
}


void SEM::ComputeMAPEstimateOfAncestralSequences() {
	if (this->root->DNArecoded.size() > 0) {
		this->ResetAncestralSequences();
	}	
	this->logLikelihood = 0;
	double currentProbability;
	map <SEM_vertex*, array<double,4>> conditionalLikelihoodMap;
	array <double,4> conditionalLikelihood;
	double maxProbability;
	double stateWithMaxProbability;
	double partialLikelihood;
	double siteLikelihood;
	double largestConditionalLikelihood = 0;
	double currentProb;
	int dna_ind_p; int dna_ind_c;
	vector <SEM_vertex *> verticesToVisit;
	SEM_vertex * p;
	SEM_vertex * c;
	emtr::Md P;
	for (int site = 0 ; site < this->num_dna_patterns; site++){
		conditionalLikelihoodMap.clear();
		this->ResetLogScalingFactors();
		for (pair<SEM_vertex *,SEM_vertex *> edge : this->edgesForPostOrderTreeTraversal){
			tie (p, c) = edge;					
			P = c->transitionMatrix;
			p->logScalingFactors += c->logScalingFactors;				
			// Initialize conditional likelihood for leaves
			if (c->observed) {
				for (int dna_c = 0; dna_c < 4; dna_c ++){
					conditionalLikelihood[dna_c] = 0;
				}
				conditionalLikelihood[c->DNArecoded[site]] = 1;
				conditionalLikelihoodMap.insert(pair<SEM_vertex *,array<double,4>>(c,conditionalLikelihood));
			}
			// Initialize conditional likelihood for ancestors
			if (conditionalLikelihoodMap.find(p) == conditionalLikelihoodMap.end()){
				for (int dna_c = 0; dna_c < 4; dna_c++){
				conditionalLikelihood[dna_c] = 1;
				}
				conditionalLikelihoodMap.insert(pair<SEM_vertex *,array<double,4>>(p,conditionalLikelihood));
			}
			largestConditionalLikelihood = 0;
			for (int dna_p = 0; dna_p < 4; dna_p++) {
				partialLikelihood = 0;
				for (int dna_c = 0; dna_c < 4; dna_c++) {
					partialLikelihood += P[dna_p][dna_c]*conditionalLikelihoodMap[c][dna_c];
				}
				conditionalLikelihoodMap[p][dna_p] *= partialLikelihood;
				if (conditionalLikelihoodMap[p][dna_p] > largestConditionalLikelihood) {
					largestConditionalLikelihood = conditionalLikelihoodMap[p][dna_p];
				}
			}
			if (largestConditionalLikelihood != 0){
				for (int dna_p = 0; dna_p < 4; dna_p++) {
					conditionalLikelihoodMap[p][dna_p] /= largestConditionalLikelihood;
				}
				p->logScalingFactors += log(largestConditionalLikelihood);
			} else {
				cout << "Largest conditional likelihood value is zero" << endl;
                throw mt_error("Largest conditional likelihood value is zero");                
			}					
		}
		maxProbability = -1; stateWithMaxProbability = 10;	
		for (dna_ind_c = 0; dna_ind_c < 4; dna_ind_c ++) {
			currentProbability = this->rootProbability[dna_ind_c];
			currentProbability *= conditionalLikelihoodMap[this->root][dna_ind_c];
			if (currentProbability > maxProbability) {
				maxProbability = currentProbability;
				stateWithMaxProbability = dna_ind_c;
			}
		}
		if (stateWithMaxProbability > 3) {
			cout << maxProbability << "\tError in computing maximum a posterior estimate for ancestor vertex\n";
		} else {
			this->root->DNArecoded.push_back(stateWithMaxProbability);
		}
//		Compute MAP estimate for each ancestral sequence
		for (pair <SEM_vertex *,SEM_vertex *> edge : this->edgesForPreOrderTreeTraversal) {			
			tie (p, c) = edge;
			P = c->transitionMatrix;
			if (!c->observed) {
				maxProbability = -1; stateWithMaxProbability = 10;
				dna_ind_p = p->DNArecoded[site];
				for (dna_ind_c = 0; dna_ind_c < 4; dna_ind_c ++){ 
					currentProbability = P[dna_ind_p][dna_ind_c];
					currentProbability *= conditionalLikelihoodMap[c][dna_ind_c];
					if (currentProbability > maxProbability) {
						maxProbability = currentProbability;
						stateWithMaxProbability = dna_ind_c;
					}
				}
				if (stateWithMaxProbability > 3) {
//					cout << "Error in computing maximum a posterior estimate for ancestor vertex";
				} else {
					c->DNArecoded.push_back(stateWithMaxProbability);
				}
			}
		}		
		siteLikelihood = 0; 							
		for (int dna = 0; dna < 4; dna++) {
			currentProb = this->rootProbability[dna]*conditionalLikelihoodMap[this->root][dna];
			siteLikelihood += currentProb;					
		}
		this->logLikelihood += (this->root->logScalingFactors + log(siteLikelihood)) * this->DNAPatternWeights[site];				
	}
}

void SEM::ComputePosteriorProbabilitiesUsingMAPEstimates() {
	this->posteriorProbabilityForVertexPair.clear();
	emtr::Md P;
	SEM_vertex * u;
	SEM_vertex * v;
	double sum;
	int dna_u; int dna_v;	
	for (unsigned int u_id = 0; u_id < this->vertexMap->size()-1; u_id ++) {
		u = (*this->vertexMap)[u_id];		
		// Posterior probability for vertex u
		u->posteriorProbability = this->GetBaseComposition(u);
		// Posterior probabilies for vertex pair (u,v)
		for (unsigned int v_id = u_id + 1 ; v_id < this->vertexMap->size()-1; v_id ++) {			
			v = (*this->vertexMap)[v_id];
			P = emtr::Md{};
			for (int site = 0; site < this->num_dna_patterns; site++ ) {		
				dna_u = u->DNArecoded[site];
				dna_v = v->DNArecoded[site];
				P[dna_u][dna_v] += this->DNAPatternWeights[site];
			}
			sum = 0;
			for (dna_u = 0; dna_u < 4; dna_u ++) {				
				for (dna_v = 0; dna_v < 4; dna_v ++) {
					sum += P[dna_u][dna_v];
				}				
			}
			for (dna_u = 0; dna_u < 4; dna_u ++) {				
				for (dna_v = 0; dna_v < 4; dna_v ++) {
					P[dna_u][dna_v] /= sum;
				}				
			}			
			this->posteriorProbabilityForVertexPair.insert(make_pair(make_pair(u,v),P));
		}
		
	}	
}

double SEM::GetExpectedMutualInformation(SEM_vertex * x, SEM_vertex* y) {
	pair <SEM_vertex *, SEM_vertex *> vertexPair;
	if (x->id < y->id) {
		vertexPair = pair <SEM_vertex *, SEM_vertex *>(x,y);
	} else {
		vertexPair = pair <SEM_vertex *, SEM_vertex *>(y,x);
	}
	
	emtr::Md P = this->posteriorProbabilityForVertexPair[vertexPair];
//	cout << "Joint probability for vertex pair " << x->name << "\t" << y->name << " is " << endl;
//	cout << P << endl;
	std::array <double, 4> P_x;
	std::array <double, 4> P_y;
	for (int dna_x = 0; dna_x < 4; dna_x ++) {
		P_x[dna_x] = 0;
		P_y[dna_x] = 0;
		for (int dna_y = 0; dna_y < 4; dna_y ++) {
			P_x[dna_x] += P[dna_x][dna_y];
			P_y[dna_x] += P[dna_y][dna_x];
		}		
	}
//	cout << "P_x is ";
//	for (int i = 0; i < 4; i++) {
//		cout << P_x[i] << "\t";
//	}
//	cout << endl << "P_y is ";
//	for (int i = 0; i < 4; i++) {
//		cout << P_y[i] << "\t";
//	}
//	cout << endl;
	double mutualInformation = 0;
//	double inc;
	for (int dna_x = 0; dna_x < 4; dna_x ++) {
		for (int dna_y = 0; dna_y < 4; dna_y ++) {
			if (P[dna_x][dna_y] > 0) {
//				inc = P[dna_x][dna_y] * log(P[dna_x][dna_y]/(P_x[dna_x] * P_y[dna_y]));
//				cout << "Incrementing mutual information by " << inc << endl; 
				mutualInformation += P[dna_x][dna_y] * log(P[dna_x][dna_y]/(P_x[dna_x] * P_y[dna_y]));
			}
		}
	}
//	cout << "P_XY is " << endl << P << endl;
//	cout << "mutual information is " << mutualInformation << endl;
	return (mutualInformation);
}

void SEM::ComputeChowLiuTree() {
	this->ClearAllEdges();
	int numberOfVertices = this->vertexMap->size();
	for (int i = 0; i < numberOfVertices; i++) {
		if ((*this->vertexMap).find(i) == (*this->vertexMap).end()){
            throw mt_error("i should be present in vertex map");
        }
	}
	const int numberOfEdges = numberOfVertices * (numberOfVertices-1)/2;	
	double maxMutualInformation = 0;
	double mutualInformation;
	double * negMutualInformation;
	negMutualInformation = new double [numberOfEdges];		
	SEM_vertex * u; SEM_vertex * v;	
	int edgeIndex = 0;
	for (int i=0; i<numberOfVertices; i++) {
		u = (*this->vertexMap)[i];
		for (int j=i+1; j<numberOfVertices; j++) {
			v = (*this->vertexMap)[j];
			mutualInformation = this->GetExpectedMutualInformation(u,v);
			negMutualInformation[edgeIndex] = -1 * mutualInformation;
			if (mutualInformation > maxMutualInformation) {
				maxMutualInformation = mutualInformation;
			}
			edgeIndex += 1;
		}
	}
		
	
	typedef pair <int, int> E;

	E * edges;
	edges = new E [numberOfEdges];
	edgeIndex = 0;
	for (int i=0; i<numberOfVertices; i++) {
		for (int j=i+1; j<numberOfVertices; j++) {
			edges[edgeIndex] = E(i,j);
			edgeIndex += 1;
		}
	}

	vector<int> p(numberOfVertices); 

	emtr::prim_graph p_graph(numberOfVertices, edges, negMutualInformation, numberOfEdges);
	
	emtr::prim(p_graph, &p[0]);
	
	for (size_t i = 0; i != p.size(); i++) {
		if (p[i] != i) {
			u = (*this->vertexMap)[i];
			v = (*this->vertexMap)[p[i]];
			u->AddNeighbor(v);
			v->AddNeighbor(u);
			if (i < p[i]){
				edgeIndex = this->GetEdgeIndex(i,p[i],numberOfVertices);
			} else {
				edgeIndex = this->GetEdgeIndex(p[i],i,numberOfVertices);
			}
		}
	}	
	delete[] edges;
	delete[] negMutualInformation;
}

void SEM::SetEdgesForTreeTraversalOperations() {
	this->SetEdgesForPreOrderTraversal();
	this->SetEdgesForPostOrderTraversal();
	this->SetVerticesForPreOrderTraversalWithoutLeaves();	
	this->SetLeaves();
}

void SEM::SetEdgesForPreOrderTraversal() {
	this->edgesForPreOrderTreeTraversal.clear();
	vector <SEM_vertex*> verticesToVisit;	
	SEM_vertex * p;	
	verticesToVisit.push_back(this->root);
	int numberOfVerticesToVisit = verticesToVisit.size();
	while (numberOfVerticesToVisit > 0) {
		p = verticesToVisit[numberOfVerticesToVisit-1];
		verticesToVisit.pop_back();
		numberOfVerticesToVisit -= 1;
		for (SEM_vertex* c : p->children){
			this->edgesForPreOrderTreeTraversal.push_back(pair<SEM_vertex*, SEM_vertex*>(p,c));
			verticesToVisit.push_back(c);
			numberOfVerticesToVisit += 1;
		}
	}
}

void SEM::SetLeaves() {
	this->leaves.clear();
	for (pair<int, SEM_vertex*> idPtrPair : * this->vertexMap){
		if (idPtrPair.second->outDegree == 0) {
			this->leaves.push_back(idPtrPair.second);
		}
	}
}

void SEM::SetEdgesForPostOrderTraversal() {	
	vector <SEM_vertex*> verticesToVisit;	
	SEM_vertex* c;
	SEM_vertex* p;	
	for (pair<int,SEM_vertex*> idPtrPair : *this->vertexMap){
		idPtrPair.second->timesVisited = 0;		
	}	
	if (this->leaves.size()== 0) {
		this->SetLeaves();
	}
	this->edgesForPostOrderTreeTraversal.clear();	
	verticesToVisit = this->leaves;
	int numberOfVerticesToVisit = verticesToVisit.size();	
	while (numberOfVerticesToVisit > 0) {
		c = verticesToVisit[numberOfVerticesToVisit -1];
		verticesToVisit.pop_back();
		numberOfVerticesToVisit -= 1;
		if (c != c->parent) {
			p = c->parent;
			this->edgesForPostOrderTreeTraversal.push_back(pair<SEM_vertex*, SEM_vertex*>(p,c));
			p->timesVisited += 1;
			if (p->timesVisited == p->outDegree) {
				verticesToVisit.push_back(p);
				numberOfVerticesToVisit += 1;
			}
		}
	}
}

void SEM::SetVerticesForPostOrderTraversal() {

}

void SEM::SetVerticesForPreOrderTraversalWithoutLeaves() {
	this->preOrderVerticesWithoutLeaves.clear();
	for (pair <SEM_vertex*, SEM_vertex*> edge : this->edgesForPreOrderTreeTraversal) {
		if (find(this->preOrderVerticesWithoutLeaves.begin(),this->preOrderVerticesWithoutLeaves.end(),edge.first) == this->preOrderVerticesWithoutLeaves.end()){
			this->preOrderVerticesWithoutLeaves.push_back(edge.first);
		}
	}
	// for (SEM_vertex * v : this->preOrderVerticesWithoutLeaves) {
	// 	cout << v->name << endl;
	// }
}

void SEM::RootedTreeAlongAnEdgeIncidentToCentralVertex() {	
	// Identify a central vertex
	vector <SEM_vertex*> verticesToVisit;
	vector <SEM_vertex*> verticesVisited;
	SEM_vertex * u; SEM_vertex * v;
	int n_ind; int u_ind;
	for (pair <int, SEM_vertex *> idPtrPair : *this->vertexMap) {
		v = idPtrPair.second;
		v->timesVisited = 0;
		if (v->observed) {
			verticesToVisit.push_back(v);
		}
	}
	int numberOfVerticesToVisit = verticesToVisit.size();
	while (numberOfVerticesToVisit > 0) {
		v = verticesToVisit[numberOfVerticesToVisit-1];
		verticesToVisit.pop_back();
		verticesVisited.push_back(v);		
		numberOfVerticesToVisit -= 1;
		for (SEM_vertex* n: v->neighbors) {
			if (find(verticesVisited.begin(),verticesVisited.end(),n)==verticesVisited.end()) {
				n->timesVisited += 1;
				if ((n->degree - n->timesVisited) == 1) {
					verticesToVisit.push_back(n);
					numberOfVerticesToVisit += 1;
				}
			} else {				
				n_ind = find(verticesVisited.begin(),verticesVisited.end(),n) - verticesVisited.begin();
				verticesVisited.erase(verticesVisited.begin()+n_ind);
			}
		}
	}
	// v is a central vertex	
	// Root tree at a randomly selected neighbor u of v
	uniform_int_distribution <int> distribution(0,v->neighbors.size()-1);
	u_ind = distribution(generator);
	u = v->neighbors[u_ind];
	this->root->AddChild(u);
	this->root->AddChild(v);
	u->AddParent(this->root);
	v->AddParent(this->root);
	verticesToVisit.clear();
	verticesVisited.clear();
	verticesToVisit.push_back(u);
	verticesToVisit.push_back(v);
	verticesVisited.push_back(u);
	verticesVisited.push_back(v);
	numberOfVerticesToVisit = verticesToVisit.size() - 1;
	while (numberOfVerticesToVisit > 0) {
		v = verticesToVisit[numberOfVerticesToVisit-1];
		verticesToVisit.pop_back();
		verticesVisited.push_back(v);
		numberOfVerticesToVisit -= 1;
		for (SEM_vertex* n: v->neighbors) {
			if (find(verticesVisited.begin(),verticesVisited.end(),n)==verticesVisited.end()) {
				verticesToVisit.push_back(n);
				numberOfVerticesToVisit += 1;
				v->AddChild(n);
				n->AddParent(v);
			}
		}
	}
	this->SetLeaves();	
	this->SetEdgesForPreOrderTraversal();
	this->SetVerticesForPreOrderTraversalWithoutLeaves();
	this->SetEdgesForPostOrderTraversal();
}

bool SEM::IsNumberOfNonSingletonComponentsGreaterThanZero() {	
	bool valueToReturn;
	if (this->indsOfVerticesOfInterest.size() > 0) {
		valueToReturn = 1;
	} else {
		valueToReturn = 0;
	}
	return (valueToReturn);
}

void SEM::SelectIndsOfVerticesOfInterestAndEdgesOfInterest() {
	this->SuppressRoot();
	vector <SEM_vertex *> verticesToVisit;
	vector <SEM_vertex *> verticesVisited;
	int n_ind; SEM_vertex * v;
	pair <int, int> edgeToAdd;
	this->indsOfVerticesOfInterest.clear();
	this->indsOfVerticesToKeepInMST.clear();
	this->edgesOfInterest_ind.clear();
	bool vertex_n_NotVisited;
	for (pair <int, SEM_vertex *> idPtrPair : * this->vertexMap) {
		v = idPtrPair.second;
		v->timesVisited = 0;
		if (v->observed and v->id < this->numberOfVerticesInSubtree) {
			verticesToVisit.push_back(v);
		}
	}
	int numberOfVerticesToVisit = verticesToVisit.size();
	while (numberOfVerticesToVisit > 0) {
		v = verticesToVisit[numberOfVerticesToVisit-1];
		verticesToVisit.pop_back();
		verticesVisited.push_back(v);		
		numberOfVerticesToVisit -= 1;
		for (SEM_vertex * n: v->neighbors) {
			vertex_n_NotVisited = find(verticesVisited.begin(),verticesVisited.end(),n) == verticesVisited.end();
			if (vertex_n_NotVisited) {
				n->timesVisited += 1;
				if ((n->degree - n->timesVisited) == 1) {
					verticesToVisit.push_back(n);
					numberOfVerticesToVisit +=1;
				}
			} else {
				edgesOfInterest_ind.push_back(pair<int,int>(n->id,v->id));
				if (n->observed) {
					this->idsOfVerticesToRemove.push_back(n->global_id);
				}
				n_ind = find(verticesVisited.begin(),verticesVisited.end(),n) - verticesVisited.begin();
				verticesVisited.erase(verticesVisited.begin() + n_ind);					
			}
		}
	}
	for (SEM_vertex * v: verticesVisited) {	
		if (!v->observed) {
			this->indsOfVerticesOfInterest.push_back(v->id);
		} else {
			this->idsOfVerticesToKeepInMST.push_back(v->global_id);
		}
	}	
}

void SEM::RenameHiddenVerticesInEdgesOfInterestAndSetIdsOfVerticesOfInterest() {		
	vector <SEM_vertex *> vertices;	
	for (pair <int, int> edge: this->edgesOfInterest_ind) {		
		vertices.push_back((*this->vertexMap)[edge.first]);
		vertices.push_back((*this->vertexMap)[edge.second]);
		for (SEM_vertex * v : vertices) {
			if (v->global_id < 0) {				
				v->global_id = this->largestIdOfVertexInMST;
				this->largestIdOfVertexInMST += 1;
				v->name = "h_" + to_string(v->global_id - this->numberOfInputSequences +1);
			}
		}		
	}
	this->idsOfVerticesOfInterest.clear();
	SEM_vertex * v;	
	for (int v_id : this->indsOfVerticesOfInterest) {
		v = (*this->vertexMap)[v_id];
		this->idsOfVerticesOfInterest.push_back(v->global_id);
	}
}

void SEM::SetIdsOfExternalVertices() {	
	this->idsOfExternalVertices.clear();
	SEM_vertex * v;
	for (int i = this->numberOfVerticesInSubtree; i < this->numberOfObservedVertices; i++) {		
		v = (*this->vertexMap)[i];
		this->idsOfExternalVertices.push_back(v->global_id);
	}
}

void SEM::SetInfoForVerticesToAddToMST(){
	this->idAndNameAndSeqTuple.clear();	
	vector <int> fullSeq;
	SEM_vertex * v;	
	for (int i : this->indsOfVerticesOfInterest){
		v = (*this->vertexMap)[i];
		fullSeq = DecompressSequence(v->DNArecoded,this->sitePatternRepetitions);
		this->idAndNameAndSeqTuple.push_back(make_tuple(v->global_id,v->name,fullSeq));
	}	
}

void SEM::AddSitePatternWeights(vector <int> sitePatternWeightsToAdd) {
	this->DNAPatternWeights = sitePatternWeightsToAdd;
	this->num_dna_patterns = this->DNAPatternWeights.size();
	this->sequenceLength = 0;
	for (int sitePatternWeight : this->DNAPatternWeights) {
		this->sequenceLength += sitePatternWeight;
	}
}

void SEM::AddSitePatternRepeats(vector <vector <int> > sitePatternRepetitionsToAdd) {
	this->sitePatternRepetitions = sitePatternRepetitionsToAdd;
}

void SEM::SetNumberOfVerticesInSubtree(int numberOfVerticesToSet) {
	this->numberOfVerticesInSubtree = numberOfVerticesToSet;
}

void SEM::AddNames(vector <string> namesToAdd) {
	if (this->numberOfObservedVertices == 0) {
		this->numberOfObservedVertices = namesToAdd.size();
	}	
	for (int i = 0; i < this->numberOfObservedVertices; i++) {
		(*this->vertexMap)[i]->name = namesToAdd[i];
		this->nameToIdMap.insert(make_pair(namesToAdd[i],i));
	}
	this->externalVertex = (*this->vertexMap)[this->numberOfObservedVertices-1];
}

void SEM::AddGlobalIds(vector <int> idsToAdd) {
	if (this->numberOfObservedVertices == 0) {
		this->numberOfObservedVertices = idsToAdd.size();
	}
	SEM_vertex * v;
	for (int i = 0; i < this->numberOfObservedVertices; i++) {
		v = (*this->vertexMap)[i];
		v->global_id = idsToAdd[i];
	}
}

void SEM::AddSequences(vector <vector <int>> sequencesToAdd) {
	this->numberOfObservedVertices = sequencesToAdd.size();
	this->node_ind = this->numberOfObservedVertices;
	for (int i = 0 ; i < this->numberOfObservedVertices; i++) {
		SEM_vertex * v = new SEM_vertex(i,sequencesToAdd[i]);
		v->observed = 1;
		this->vertexMap->insert(make_pair(i,v));
	}	
}

void SEM::AddRootVertex() {
	int n = this->numberOfObservedVertices;
	vector <int> emptySequence;	
	this->root = new SEM_vertex (-1,emptySequence);
	this->root->name = "h_root";	
	this->root->id = ( 2 * n ) - 2;
	this->vertexMap->insert(pair<int,SEM_vertex*>((( 2 * n ) - 2 ),this->root));
	this->nameToIdMap.insert(make_pair(this->root->name,this->root->id));
}

void SEM::SetEdgesFromTopologyFile() {
	SEM_vertex * u; SEM_vertex * v;
	vector <string> splitLine;
	string u_name; string v_name; double t;
	t = 0.0;
	int num_edges = 0;
	vector <int> emptySequence;
	ifstream inputFile(this->topologyFileName.c_str());
	for(string line; getline(inputFile, line );) {
		num_edges++;
		vector<string> splitLine = emtr::split_ws(line);
		u_name = splitLine[0];
		v_name = splitLine[1];
		if (this->ContainsVertex(u_name)) {
			u = (*this->vertexMap)[this->nameToIdMap[u_name]];
		} else {
			u = new SEM_vertex(this->node_ind,emptySequence);
			u->name = u_name;
			u->id = this->node_ind;
			this->vertexMap->insert(pair<int,SEM_vertex*>(u->id,u));
			this->nameToIdMap.insert(make_pair(u->name,u->id));			
			this->node_ind += 1;
			if(!this->ContainsVertex(u_name)){
                throw mt_error("check why u is not in vertex map");
            }
		}

		if (this->ContainsVertex(v_name)) {
			v = (*this->vertexMap)[this->nameToIdMap[v_name]];
		} else {			
			v = new SEM_vertex(this->node_ind,emptySequence);
			v->name = v_name;
			v->id = this->node_ind;
			this->vertexMap->insert(pair<int,SEM_vertex*>(v->id,v));
			this->nameToIdMap.insert(make_pair(v->name,v->id));
			this->node_ind += 1;
			if(!this->ContainsVertex(v_name)){
                throw mt_error("Check why v is not in vertex map");
            }
		}
		u->AddNeighbor(v);
		v->AddNeighbor(u);		
		if (u->id < v->id) {
			this->edgeLengths.insert(make_pair(make_pair(u,v),t));			
		} else {
			this->edgeLengths.insert(make_pair(make_pair(v,u),t));
		}		
	}
	inputFile.close();
	// cout << "number of edges in topology file is " << num_edges << endl;
}

void SEM::CompressAASequences() {			
	this->AAPatternWeights.clear();	
	vector <vector<int>> uniquePatterns;
	this->gapLessAAFlag.clear();	
	map <vector<int>,vector<int>> uniquePatternsToSitesWherePatternRepeats;	
	for (unsigned int i = 0; i < this->leaves.size(); i++) {
		SEM_vertex * l_ptr = this->leaves[i];
		l_ptr->AAcompressed.clear();		
	}
	int numberOfSites = this->leaves[0]->AArecoded.size();
	vector<int> sitePattern;
	for(int site = 0; site < numberOfSites; site++) {
		sitePattern.clear();
		for (SEM_vertex* l_ptr: this->leaves) {
			sitePattern.push_back(l_ptr->AArecoded[site]);}
		if (find(uniquePatterns.begin(),uniquePatterns.end(),sitePattern)!=uniquePatterns.end()) {
			uniquePatternsToSitesWherePatternRepeats[sitePattern].push_back(site);
			
		} else {
			uniquePatterns.push_back(sitePattern);	
			vector<int> sitePatternRepeats;
			sitePatternRepeats.push_back(site);
			uniquePatternsToSitesWherePatternRepeats[sitePattern] = sitePatternRepeats;						
			for (unsigned int i = 0; i < sitePattern.size(); i++) {				
				SEM_vertex * l_ptr = this->leaves[i];
				l_ptr->AAcompressed.push_back(sitePattern[i]);
			}
		}
	}
	for (vector <int> sitePattern : uniquePatterns) {
		int sitePatternWeight = uniquePatternsToSitesWherePatternRepeats[sitePattern].size();
		this->AAPatternWeights.push_back(sitePatternWeight);
		bool gapLessPattern = 1;
		for (int character: sitePattern) if (character < 0) gapLessPattern = 0;
		this->gapLessAAFlag.push_back(gapLessPattern);			
	}	
	this->num_aa_patterns = this->AAPatternWeights.size();
}

#include <algorithm>
#include <numeric>
#include <map>
#include <vector>

void SEM::CompressDNASequences() {
    
    this->DNAPatternWeights.clear();
    this->gapLessDNAFlag.clear();

    vector<vector<int>> uniquePatterns;
    map<vector<int>, vector<int>> uniquePatternsToSitesWherePatternRepeats;

    // Clear per-leaf compressed storage
    for (unsigned int i = 0; i < this->leaves.size(); ++i) {
        SEM_vertex* l_ptr = this->leaves[i];
        l_ptr->DNAcompressed.clear();
    }

    // Build unique patterns and where they repeat
    const int numberOfSites = static_cast<int>(this->leaves[0]->DNArecoded.size());
    vector<int> sitePattern;
    sitePattern.reserve(this->leaves.size());

    for (int site = 0; site < numberOfSites; ++site) {
        sitePattern.clear();
        for (SEM_vertex* l_ptr : this->leaves) {
            sitePattern.push_back(l_ptr->DNArecoded[site]);
        }

        auto it = uniquePatternsToSitesWherePatternRepeats.find(sitePattern);
        if (it != uniquePatternsToSitesWherePatternRepeats.end()) {
            it->second.push_back(site);
        } else {
            uniquePatterns.push_back(sitePattern);
            uniquePatternsToSitesWherePatternRepeats[uniquePatterns.back()] = std::vector<int>(1, site);
        }
    }

    // Weights and gap-less flags (initial order)
    for (const std::vector<int>& pat : uniquePatterns) {
        int w = static_cast<int>(uniquePatternsToSitesWherePatternRepeats[pat].size());
        this->DNAPatternWeights.push_back(w);

        bool gapLess = true;
        for (int ch : pat) { if (ch < 0) { gapLess = false; break; } }
        this->gapLessDNAFlag.push_back(gapLess);
    }

    this->num_dna_patterns = static_cast<int>(this->DNAPatternWeights.size());
    if (this->num_dna_patterns == 0) {
        this->CumPercDNAPatternWeights.clear();
        return;
    }

    // Sort pattern indices by decreasing weight
    vector<size_t> order(this->DNAPatternWeights.size());
    iota(order.begin(), order.end(), 0);
    stable_sort(order.begin(), order.end(),
        [this](size_t a, size_t b) {
            return this->DNAPatternWeights[a] > this->DNAPatternWeights[b];
        }
    );

	this->DNAPatternOriginalSites.clear();
	this->DNAPatternOriginalSites.resize(order.size());

	this->OriginalToCompressed.assign(numberOfSites, -1);

	for (size_t k = 0; k < order.size(); ++k) {
    const size_t pidx = order[k];
    const std::vector<int>& sites =
        uniquePatternsToSitesWherePatternRepeats[ uniquePatterns[pidx] ];

    this->DNAPatternOriginalSites[k] = sites;

    for (int s : sites) {
        if (s >= 0 && s < numberOfSites) {
            this->OriginalToCompressed[s] = static_cast<int>(k);
        	}
    	}
	}

    // Rebuild DNAcompressed and metadata in sorted order
    for (SEM_vertex* l_ptr : this->leaves) {
        l_ptr->DNAcompressed.clear();
        l_ptr->DNAcompressed.reserve(order.size());
    }

    std::vector<int>  newWeights;  newWeights.reserve(order.size());
    std::vector<bool> newGapLess;  newGapLess.reserve(order.size());

    for (size_t k = 0; k < order.size(); ++k) {
        const size_t pidx = order[k];
        const std::vector<int>& pat = uniquePatterns[pidx];

        for (unsigned int i = 0; i < pat.size(); ++i) {
            SEM_vertex* l_ptr = this->leaves[i];
            l_ptr->DNAcompressed.push_back(pat[i]);
        }

        newWeights.push_back(this->DNAPatternWeights[pidx]);
        newGapLess.push_back(this->gapLessDNAFlag[pidx]);
    }

    this->DNAPatternWeights.swap(newWeights);
    this->gapLessDNAFlag.swap(newGapLess);
    
    this->CumPercDNAPatternWeights.clear();
    this->CumPercDNAPatternWeights.reserve(this->DNAPatternWeights.size());

    const double totalW = std::accumulate(
        this->DNAPatternWeights.begin(), this->DNAPatternWeights.end(), 0.0);

    double running = 0.0;
    if (totalW > 0.0) {
        for (size_t i = 0; i < this->DNAPatternWeights.size(); ++i) {
            running += static_cast<double>(this->DNAPatternWeights[i]);
            double pct = (running / totalW) * 100.0;
            if (pct > 100.0) pct = 100.0;
            this->CumPercDNAPatternWeights.push_back(pct);
        }        
        this->CumPercDNAPatternWeights.back() = 100.0;
    } else {
        // No weight -> all zeros
        for (size_t i = 0; i < this->DNAPatternWeights.size(); ++i) {
            this->CumPercDNAPatternWeights.push_back(0.0);
        }
    }
}

void SEM::CompressDNASequencesBkp() {			
	this->DNAPatternWeights.clear();	
	vector <vector<int>> uniquePatterns;
	this->gapLessDNAFlag.clear();
	map <vector<int>,vector<int>> uniquePatternsToSitesWherePatternRepeats;	
	for (unsigned int i = 0; i < this->leaves.size(); i++) {
		SEM_vertex * l_ptr = this->leaves[i];
		l_ptr->DNAcompressed.clear();		
	}
	int numberOfSites = this->leaves[0]->DNArecoded.size();
	vector<int> sitePattern;
	for(int site = 0; site < numberOfSites; site++) {
		sitePattern.clear();
		for (SEM_vertex* l_ptr: this->leaves) {
			sitePattern.push_back(l_ptr->DNArecoded[site]);}
		if (find(uniquePatterns.begin(),uniquePatterns.end(),sitePattern)!=uniquePatterns.end()) {
			uniquePatternsToSitesWherePatternRepeats[sitePattern].push_back(site);
			
		} else {
			uniquePatterns.push_back(sitePattern);	
			vector<int> sitePatternRepeats;
			sitePatternRepeats.push_back(site);
			uniquePatternsToSitesWherePatternRepeats[sitePattern] = sitePatternRepeats;						
			for (unsigned int i = 0; i < sitePattern.size(); i++) {				
				SEM_vertex * l_ptr = this->leaves[i];
				l_ptr->DNAcompressed.push_back(sitePattern[i]);
			}
		}
	}
	for (vector <int> sitePattern : uniquePatterns) {
		int sitePatternWeight = uniquePatternsToSitesWherePatternRepeats[sitePattern].size();
		this->DNAPatternWeights.push_back(sitePatternWeight);
		bool gapLessPattern = 1;
		for (int character: sitePattern) if (character < 0) gapLessPattern = 0;
		this->gapLessDNAFlag.push_back(gapLessPattern);
	}	
	this->num_dna_patterns = this->DNAPatternWeights.size();
	
	vector <size_t> order(this->DNAPatternWeights.size());
	iota(order.begin(), order.end(), 0);
	stable_sort(order.begin(), order.end(), [this](size_t a, size_t b) 
	{
		return this->DNAPatternWeights[a] > this->DNAPatternWeights[b];
	});
	
}

void SEM::SetDNASequencesFromFile(string sequenceFileName) {
	this->leaves.clear();
	this->DNAsequenceFileName = sequenceFileName;
	vector <int> recodedSequence;
	recodedSequence.clear();
	unsigned int site = 0;
    unsigned int seq_len = 0;
	int dna_index;	
	ifstream inputFile(sequenceFileName.c_str());
	string seqName;
	string seq = "";
	SEM_vertex * u;
	this->numberOfInputSequences = 0;
	this->numberOfObservedVertices = 0;
	for(string line; getline(inputFile, line );) {
		if (line[0]=='>') {
			if (seq != "") {
				for (char const dna: seq) {
					if (!isspace(dna)) {
						dna_index = this->ConvertDNAtoIndex(dna);						
						recodedSequence.push_back(dna_index);					
						site += 1;							
						}						
				}
				this->AddVertex(seqName,recodedSequence);
				u = this->GetVertex(seqName);
				this->leaves.push_back(u);
				u->observed = 1;
				this->numberOfInputSequences++;
				this->numberOfObservedVertices++;								
				recodedSequence.clear();
			} 
			seqName = line.substr(1,line.length());
			seq = "";
			site = 0;			
		}
		else {
			seq += line ;
		}		
	}		
	for (char const dna: seq) {
		if (!isspace(dna)) {
			dna_index = this->ConvertDNAtoIndex(dna);			
			recodedSequence.push_back(dna_index);
			site += 1;
		}
	}
	this->AddVertex(seqName,recodedSequence);
	u = this->GetVertex(seqName);
	this->leaves.push_back(u);
	u->observed = 1;
	this->numberOfInputSequences++;
	this->numberOfObservedVertices++;								
	recodedSequence.clear();
    seq_len = recodedSequence.size();
	recodedSequence.clear();
	inputFile.close(); 
	// cout << this->numberOfObservedVertices << endl;   
}

void SEM::SetAASequencesFromFile(string sequenceFileName) {
	this->AAsequenceFileName = sequenceFileName;
	vector <int> recodedSequence;
	recodedSequence.clear();
	unsigned int site = 0;
	int aa_index;
	ifstream inputFile(sequenceFileName.c_str());
	string seqName;
	string seq = "";
	// this->numberOfInputSequences = 0;
	// this->numberOfObservedVertices = 0;
	for(string line; getline(inputFile, line );) {
		if (line[0]=='>') {
			if (seq != "") {
				for (char const aa: seq) {
					if (!isspace(aa)) {
						aa_index = this->ConvertAAtoIndex(aa);
						recodedSequence.push_back(aa_index);
						site += 1;
						}
				}
				this->SetAASequenceForVertex(seqName,recodedSequence);
				recodedSequence.clear();
			}
			seqName = line.substr(1,line.length());
			seq = "";
			site = 0;
		}
		else {
			seq += line ;
		}
	}
	for (char const aa: seq) {
		if (!isspace(aa)) {
			aa_index = this->ConvertAAtoIndex(aa);
			recodedSequence.push_back(aa_index);
			site += 1;
		}
	}
	this->SetAASequenceForVertex(seqName,recodedSequence);
	// cout << "completed setting AA sequences" << endl;	
}

void SEM::AddWeightedEdges(vector <tuple <string,string,double> > weightedEdgesToAdd) {
	SEM_vertex * u; SEM_vertex * v;
	string u_name; string v_name; double t;
	vector <int> emptySequence;
	int edge_count_h1 = 0;
	for (tuple <string, string, double> weightedEdge : weightedEdgesToAdd) {
		tie (u_name, v_name, t) = weightedEdge;
		if (!this->ContainsVertex(u_name)) this->AddVertex(u_name, emptySequence);
		if (!this->ContainsVertex(v_name)) this->AddVertex(v_name, emptySequence);
		u = this->GetVertex(u_name); v = this->GetVertex(v_name);		
		
		u->AddNeighbor(v);
		v->AddNeighbor(u);
		if (u->id < v->id) {
			this->edgeLengths.insert(make_pair(make_pair(u,v),t));
		} else {
			this->edgeLengths.insert(make_pair(make_pair(v,u),t));
		}
	}
}

void SEM::AddEdgeLogLikelihoods(vector<tuple<string,string,double>> edgeLogLikelihoodsToAdd) {
	SEM_vertex * u; SEM_vertex * v; double edgeLogLikelihood;
	string u_name; string v_name;
	pair<SEM_vertex *, SEM_vertex *> vertexPair;
	for (tuple<string,string,double> edgeLogLikelihoodTuple : edgeLogLikelihoodsToAdd) {
		tie (u_name, v_name, edgeLogLikelihood) = edgeLogLikelihoodTuple;
		u = (*this->vertexMap)[this->nameToIdMap[u_name]];
		v = (*this->vertexMap)[this->nameToIdMap[v_name]];
		vertexPair = pair <SEM_vertex *, SEM_vertex *> (u,v);
		this->edgeLogLikelihoodsMap.insert(pair<pair <SEM_vertex *, SEM_vertex *>,double>(vertexPair, edgeLogLikelihood));
	}	
}

void SEM::AddVertexLogLikelihoods(map<string,double> vertexLogLikelihoodsMapToAdd) {
	string v_name; SEM_vertex * v; double vertexLogLikelihood;
	for (pair<string,double> vNameAndLogLik : vertexLogLikelihoodsMapToAdd) {
		tie (v_name, vertexLogLikelihood) = vNameAndLogLik;
		v = (*this->vertexMap)[this->nameToIdMap[v_name]];
		v->vertexLogLikelihood = vertexLogLikelihood;
	}
}

bool SEM::ContainsVertex(string v_name) {	
	if (this->nameToIdMap.find(v_name) == this->nameToIdMap.end()) {		
		return (0);
	} else {
		return (1);
	}
}

void SEM::SetAASequenceForVertex(string v_name, vector<int> aaRecoded){
	SEM_vertex * v = this->GetVertex(v_name);
	v->AArecoded = aaRecoded;
}

void SEM::AddVertex(string u_name, vector <int> emptySequence) {	
	SEM_vertex * u = new SEM_vertex(this->node_ind,emptySequence);
	u->name = u_name;
	u->id = this->node_ind;
	this->node_ind++;
	this->vertexMap->insert(pair<int,SEM_vertex*>(u->id,u));
	this->nameToIdMap.insert(make_pair(u->name,u->id));	
}

double SEM::ComputeDistance(int v_i, int v_j) {
	vector <int> seq_i; vector <int> seq_j;	
	seq_i = (*this->vertexMap)[v_i]->DNArecoded;
	seq_j = (*this->vertexMap)[v_j]->DNArecoded;
	double sequence_length = 0;
	for (int site = 0; site < num_dna_patterns; site++) {
		sequence_length += double(this->DNAPatternWeights[site]);
	}
	if(sequence_length <= 0){
        throw mt_error("Check pattern detection");
    }
	// logDet_distance
	double distance;	
	if (flag_Hamming) {
		distance = 0;
		for (int site = 0; site < num_dna_patterns; site++) {
			if (seq_i[site] != seq_j[site]) {
			distance += double(this->DNAPatternWeights[site])/sequence_length;
			}						
		}				
	} else {
        throw mt_error("Check distance flag");
    }
	return (distance);
}


///...///...///...///...///...///...///... Constructor and Destructor for EM manager ///...///...///...///...///...///...///...///...///

EMManager::EMManager(string DNAsequenceFileNameToSet,				 
					 string AAsequenceFileNameToSet,
					 const std::string& parameters_json) {
		this->parameters_json = parameters_json;		
		
		if (DNAsequenceFileNameToSet.empty()) {
        	DNAsequenceFileNameToSet = "/data/RAxML_DNA_test.fa";
		}
		if (AAsequenceFileNameToSet.empty()) {
			AAsequenceFileNameToSet = "/data/RAxML_AA_test.fa";
		}		

		this->prefix_for_output_files = "";
		this->supertree_method = "";        
		this->verbose = 0;
		this->distance_measure_for_NJ = "log-det";
		this->flag_topology = 1;		
		this->max_EM_iter = max_iter;
		this->numberOfLargeEdgesThreshold = 100;
		this->SetDNAMap();
		this->ancestralSequencesString = "";
		// make an alignment class	
		// parse fa, phy and nex

		this->P = new SEM(0.0,100,this->verbose);					
		this->F = new FamilyJoining(DNAsequenceFileNameToSet, 0.0);	
		this->P->SetFJTree(this->F->GetTree());			
		vector <tuple <string, string, double>> edge_vector = this->F->GetEdgeVector();
		{
			const auto& ev = this->F->GetEdgeVector();
			nlohmann::json j;

			// Build [[u,v,w], ...]
			std::vector<nlohmann::json> arr;
			arr.reserve(ev.size());
			for (const auto& e : ev) {
				const std::string& u = std::get<0>(e);
				const std::string& v = std::get<1>(e);
				const double       w = std::get<2>(e);
				arr.emplace_back(nlohmann::json::array({ u, v, w }));
			}

			j["edge_list"] = std::move(arr);
			j["source"]    = "familyjoining";
			j["method"]    = "epsilon_zero";
			j["label"]     = "final";
			j["n_edges"] = static_cast<int>(ev.size());			

			std::cout << "[FAMILYJOINING]\t" << j.dump() << "\n";
			fflush(stdout);
		}

		this->P->SetDNASequencesFromFile(DNAsequenceFileNameToSet);				
		this->P->AddWeightedEdges(edge_vector);				
		this->P->SetVertexVector();
		this->P->CompressDNASequences();		
		this->P->SetParameters(parameters_json);
		
		const auto& cum = this->P->CumPercDNAPatternWeights;				
		this->P->num_patterns_for_layer[0] = this->P->first_index_gt(cum, this->P->cum_pattern_weight_coarse);
		this->P->num_patterns_for_layer[1] = this->P->first_index_gt(cum, this->P->cum_pattern_weight_medium);
		this->P->num_patterns_for_layer[2] = this->P->first_index_gt(cum, 100);

		this->P->ecdllConvergenceThresholdForTrifle.resize(3);
		this->P->ecdllConvergenceThresholdForTrifle[0] = this->P->cum_pattern_weight_coarse/(100.0) * this->P->ecdllConvergenceThreshold;
		this->P->ecdllConvergenceThresholdForTrifle[1] = this->P->cum_pattern_weight_medium/(100.0) * this->P->ecdllConvergenceThreshold;
		this->P->ecdllConvergenceThresholdForTrifle[2] = this->P->ecdllConvergenceThreshold;
		}

		EMManager::~EMManager() {
			delete this->F;	
			delete this->P;		
		}

	void EMManager::EM_main() {
	this->num_repetitions = this->P->num_repetitions;
	this->P->max_log_likelihood_best = -1 * std::pow(10.0, 10.0);

	const int total_nodes = max(0, int(this->P->vertexMap->size()) - this->P->numberOfObservedVertices);

	auto emit_progress = [&](const char* kind,
							const char* method = nullptr,
							int rep1 = -1) {
		nlohmann::json j{
		{"type", "progress"},
		{"kind", kind}
		};
		if (method) j["method"] = method;
		if (rep1 >= 0) j["rep"] = rep1;

		
		if (std::strcmp(kind, "setup:reps") == 0) {
		j["totalReps"] = this->num_repetitions;
		} else if (std::strcmp(kind, "setup:edges") == 0) {
		j["nodeTotal"] = total_nodes;
		}

		std::cout << "[EM_Method] " << j.dump() << "\n";
		std::fflush(stdout);
	};

	
	emit_progress("setup:reps");
	emit_progress("setup:edges");

	for (int rep = 0; rep < this->num_repetitions; ++rep) {	
		this->P->EMTrifle_DNA_for_replicate(rep + 1);
		emit_progress("rep", nullptr, rep + 1);
	}
}

void EMManager::UNRESTforEachBranch(){
	// 

}


void EMManager::setDayhoffMatrixPath(const std::string& path) {
  this->rate_matrix_path_ = path;  
}

void EMManager::EMparsimony() {
    cout << "Starting EM with initial parameters set using parsimony" << endl;
	this->P->probabilityFileName_pars = this->prefix_for_output_files + ".pars_prob";
	this->probabilityFileName_pars = this->prefix_for_output_files + ".pars_prob";
    this->P->EM_rooted_at_each_internal_vertex_started_with_parsimony(this->num_repetitions);
}


void EMManager::EMdirichlet() {	
	cout << "Starting EM with initial parameters set using Dirichlet" << endl;
	this->P->probabilityFileName_diri = this->prefix_for_output_files + ".diri_prob";
	this->probabilityFileName_diri = this->prefix_for_output_files + ".diri_prob";
	this->P->EM_rooted_at_each_internal_vertex_started_with_dirichlet(this->num_repetitions);
}

void EMManager::SetprobFileforHSS() {
	if (this->max_log_lik_pars > this->max_log_lik_diri) {
		cout << "Initializing with Parsimony yielded higher log likelihood score" << endl;
		this->probabilityFileName_best = this->probabilityFileName_pars;
	} else {
		this->probabilityFileName_best = this->probabilityFileName_diri;
		cout << "Initializing with Dirichlet yielded higher log likelihood score" << endl;
	}
}

void EMManager::EMhss() {
	this->SetprobFileforHSS();
	cout << "Starting EM with initial parameters set using Bayes rule as described in HSS paper" << endl;
	this->P->probabilityFileName_best = this->probabilityFileName_best;	
    this->P->SetBHparameters();
	this->P->ReparameterizeBH();    
    this->P->EM_rooted_at_each_internal_vertex_started_with_HSS_par(this->num_repetitions);
}

void EMManager::SetDNAMap() {
	this->mapDNAtoInteger["A"] = 0;
	this->mapDNAtoInteger["C"] = 1;
	this->mapDNAtoInteger["G"] = 2;
	this->mapDNAtoInteger["T"] = 3;
}

void EMManager::EMgivenInputTopology() {
	// vector <string> names;
	// vector <vector <int> > sequences;
	// vector <int> sitePatternWeights;
	// vector <vector <int> > sitePatternRepetitions;
	// vector <int> idsOfVerticesForSEM;
	
	// int numberOfInputSequences = (int) this->M->vertexMap->size();
	// idsOfVerticesForSEM.clear();
	// for (pair <int, MST_vertex *> vIdAndPtr : * this->M->vertexMap) {
		// idsOfVerticesForSEM.push_back(vIdAndPtr.first);
	// }
	 // Transfer following feature to SEM
	// tie (names, sequences, sitePatternWeights, sitePatternRepetitions) = this->M->GetCompressedSequencesSiteWeightsAndSiteRepeats(idsOfVerticesForSEM);
	// cout << "setting sequence file name, topology file name, site pattern weights, number of input sequences" << endl;
    // cout << "number of site patterns is " << sitePatternWeights.size() << endl;	
	// this->P->sequenceFileName = this->fastaFileName;
	// this->P->topologyFileName = this->topologyFileName;
	// this->P->AddSequences(sequences);
	// this->P->AddNames(names);
	// this->P->AddSitePatternWeights(sitePatternWeights);
	// this->P->SetNumberOfInputSequences(numberOfInputSequences);
	// this->P->numberOfObservedVertices = numberOfInputSequences;
	// cout << "setting edges from topology file" << endl;	
	// this->P->SetEdgesFromTopologyFile();
}

int EMManager::GetEdgeIndex (int vertexIndex1, int vertexIndex2, int numberOfVertices){
	int edgeIndex;
	edgeIndex = numberOfVertices*(numberOfVertices-1)/2;
	edgeIndex -= (numberOfVertices-vertexIndex1)*(numberOfVertices-vertexIndex1-1)/2;
	edgeIndex += vertexIndex2 - vertexIndex1 - 1;
	return edgeIndex;
}

int EMManager::ComputeHammingDistance(string seq1, string seq2) {
	int hammingDistance = 0;
	for (unsigned int i=0;i<seq1.length();i++){
		if (seq1[i] != seq2[i]){
			hammingDistance+=1;
		}		
	}
	return (hammingDistance);
};
