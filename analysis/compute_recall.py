from typing import List, Tuple
import networkx as nx
from typing import Any, Set, FrozenSet
import networkx as nx
import matplotlib.pyplot as plt
import networkx as nx
from typing import Any, List
import numpy as np
from functools import lru_cache

def _parse_edge_line(line: str) -> Tuple[str, str]:
    """
    Parse a single line containing exactly two node labels separated by
    whitespace or a comma. Returns (u, v).
    """
    s = line.strip()
    if not s or s.startswith("#"):
        raise ValueError("blank/comment")  # handled by caller

    # Support either "u v" or "u,v"
    if "," in s:
        parts = [p.strip() for p in s.split(",")]
    else:
        parts = s.split()

    if len(parts) != 2:
        raise ValueError(f"Expected exactly two tokens, got: {s!r}")
    u, v = parts
    return u, v
import csv
import networkx as nx
from typing import List, Tuple

def read_directed_edges_from_csv(path: str) -> nx.DiGraph:
    """
    Read directed edges from a CSV file and return a NetworkX DiGraph.

    File format:
        - Two columns per row: parent,child
        - Delimiter can be comma or whitespace
        - Lines starting with '#' or blank lines are ignored

    Parameters
    ----------
    path : str
        Path to the edge list file.

    Returns
    -------
    nx.DiGraph
        Directed graph with the given edges.
    """
    edges: List[Tuple[str, str]] = []

    with open(path, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
            # Skip blank or malformed rows
            if not row or row[0].strip().startswith("#"):
                continue
            # Handle comma-delimited or space-delimited cases
            if len(row) == 1 and " " in row[0]:
                parts = row[0].split()
            else:
                parts = [c.strip() for c in row if c.strip()]
            if len(parts) != 2:
                raise ValueError(f"Bad edge line (need exactly 2 tokens): {row!r}")
            u, v = parts
            edges.append((u, v))

    G = nx.DiGraph()
    G.add_edges_from(edges)
    return G

def read_topology_to_graph(path: str, *, validate_tree: bool = False) -> nx.Graph:
    """
    Read an edge list from 'path' and return an undirected NetworkX Graph.

    The file should have two node labels per line, either whitespace- or
    comma-separated. Lines starting with '#' and blank lines are ignored.

    Parameters
    ----------
    path : str
        Path to the edge list file (e.g., 'ternary_topology.csv').
    validate_tree : bool
        If True, verify the result is a tree (connected & acyclic).

    Returns
    -------
    nx.Graph
        Undirected graph with the given edges.
    """
    edges: List[Tuple[str, str]] = []
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            try:
                u, v = _parse_edge_line(raw)
            except ValueError as e:
                # ignore blank/comment lines only
                if str(e) == "blank/comment":
                    continue
                else:
                    raise
            edges.append((u, v))

    G = nx.Graph()
    G.add_edges_from(edges)

    if validate_tree:
        if not nx.is_tree(G):
            if not nx.is_connected(G):
                comps = list(nx.connected_components(G))
                sizes = [len(c) for c in comps]
                raise ValueError(f"Graph is not connected: {len(comps)} components (sizes: {sizes})")
            else:
                # connected but has a cycle
                cycles = nx.cycle_basis(G)
                raise ValueError(f"Graph has cycles; example cycle: {cycles[0] if cycles else 'unknown'}")

    return G

def hierarchical_leaf_clusters_proper(T: nx.DiGraph) -> Set[FrozenSet[Any]]:
    """
    Given a rooted tree T (DiGraph, edges parent->child), return a set of
    frozensets, each the set of descendant leaf names under an internal node.
    
    Excludes:
      - singletons (size 1)
      - the full leaf set (entire tree)
    """
    if not isinstance(T, nx.DiGraph):
        raise TypeError("T must be a networkx.DiGraph")

    # infer the (unique) root
    roots = [n for n in T.nodes if T.in_degree(n) == 0]
    if len(roots) != 1:
        raise ValueError(f"Expected exactly one root, found {len(roots)}: {roots}")
    root = roots[0]

    # identify leaves
    leaves = {n for n in T.nodes if T.out_degree(n) == 0}
    if not leaves:
        raise ValueError("No leaves found in the tree.")

    full_leafset = frozenset(leaves)

    @lru_cache(maxsize=None)
    def leafset(u: Any) -> FrozenSet[Any]:
        if u in leaves:
            return frozenset([u])
        return frozenset().union(*(leafset(v) for v in T.successors(u)))

    clusters: Set[FrozenSet[Any]] = set()
    for n in T.nodes:
        if T.out_degree(n) > 0:  # internal node
            L = leafset(n)
            if len(L) > 1 and L != full_leafset:  # exclude singletons & full set
                clusters.add(L)

    return clusters

def _orient_undirected_tree(Gu: nx.Graph, root: Any) -> nx.DiGraph:
    """Orient an undirected tree Gu away from `root` and return a DiGraph."""
    if root not in Gu:
        raise ValueError(f"Root {root!r} not found in graph.")
    if not nx.is_tree(Gu):
        raise ValueError("Input graph must be a tree (connected, acyclic).")
    # Direct edges parent->child along a BFS from the root
    directed_edges: List[Tuple[Any, Any]] = list(nx.bfs_edges(Gu, root))
    Gd = nx.DiGraph()
    Gd.add_edges_from(directed_edges)
    return Gd

def clusters_when_rooted_at(
    T_unrooted: nx.Graph, root: Any
) -> Set[FrozenSet[Any]]:
    """
    Root the undirected tree T_unrooted at `root` and return the set of
    non-trivial hierarchical clusters (leaf-name frozensets), i.e.,
    excluding singletons and the full leaf set.
    """
    if T_unrooted.degree(root) <= 1:
        raise ValueError(f"Root {root!r} is not an internal node (degree={T_unrooted.degree(root)}).")
    T_rooted = _orient_undirected_tree(T_unrooted, root)
    return hierarchical_leaf_clusters_proper(T_rooted)

def _orient_undirected_tree(Gu: nx.Graph, root: Any) -> nx.DiGraph:
    """Orient an undirected tree Gu away from `root` and return a DiGraph."""
    if root not in Gu:
        raise ValueError(f"Root {root!r} not found in the undirected tree.")
    if not nx.is_tree(Gu):
        raise ValueError("T_unrooted must be a tree (connected & acyclic).")
    Gd = nx.DiGraph()
    Gd.add_edges_from(nx.bfs_edges(Gu, root))  # parent -> child
    return Gd

def cluster_recall_when_rooted_at(
    T_rooted: nx.DiGraph,
    T_unrooted: nx.Graph,
    root_candidate: Any,
    *,
    strict_leaf_match: bool = True,
) -> dict[str, Any]:
    """
    Compute the fraction (recall) of non-trivial clusters in T_rooted that
    are present in T_unrooted when rooted at `root_candidate`.

    Returns a dict with:
      - 'recall': float in [0,1]
      - 'n_common': int
      - 'n_ref': int  (non-trivial clusters in T_rooted)
      - 'n_test': int (non-trivial clusters in T_unrooted rooted at root_candidate)
      - 'missing_from_test': Set[FrozenSet[Any]] (clusters in ref not found in test)
      - 'common_clusters': Set[FrozenSet[Any]]   (intersection)
    """
    # Reference clusters from already-rooted tree
    ref_clusters: Set[FrozenSet[Any]] = hierarchical_leaf_clusters_proper(T_rooted)

    # Orient unrooted tree at the candidate and get its clusters
    T_test = _orient_undirected_tree(T_unrooted, root_candidate)
    test_clusters: Set[FrozenSet[Any]] = hierarchical_leaf_clusters_proper(T_test)

    # Optional: enforce identical leaf sets
    if strict_leaf_match:
        ref_leaves = {n for n in T_rooted.nodes if T_rooted.out_degree(n) == 0}
        test_leaves = {n for n in T_test.nodes if T_test.out_degree(n) == 0}
        if ref_leaves != test_leaves:
            raise ValueError(
                f"Leaf sets differ between T_rooted ({len(ref_leaves)}) and "
                f"T_unrooted@{root_candidate} ({len(test_leaves)}). "
                f"Diff: only-in-ref={sorted(ref_leaves - test_leaves)}, "
                f"only-in-test={sorted(test_leaves - ref_leaves)}"
            )

    common = ref_clusters & test_clusters
    missing = ref_clusters - test_clusters

    n_ref = len(ref_clusters)
    recall = (len(common) / n_ref) if n_ref > 0 else 0.0

    return {
        "recall": recall,
        "n_common": len(common),
        "n_ref": n_ref,
        "n_test": len(test_clusters),
        "missing_from_test": missing,
        "common_clusters": common,
    }

# --- Example usage ---
if __name__ == "__main__":
    T_unrooted = read_topology_to_graph("analysis/ternary_topology.csv", validate_tree=True)

    # print(f"Loaded graph with {T_unrooted.number_of_nodes()} nodes and {T_unrooted.number_of_edges()} edges.")
    # Example: list leaves (degree 1)
    # leaves = sorted([n for n in T_unrooted.nodes if T_unrooted.degree(n) == 1], key=str)
    # print(f"Leaves ({len(leaves)}): {leaves}")

    T_rooted = read_directed_edges_from_csv("analysis/Randall_real_rooted_tree.csv")
    # print(f"Loaded directed graph with {T_rooted.number_of_nodes()} nodes and {T_rooted.number_of_edges()} edges.")
    # roots = [n for n in T_rooted.nodes if T_rooted.in_degree(n) == 0]
    # # print(f"Root candidates: {roots}")
    
    
    # result = cluster_recall_when_rooted_at(T_rooted, T_unrooted, "h_30")

    # print(f"Recall@h_30 = {result['recall']:.3f} "
    #     f"({result['n_common']}/{result['n_ref']} clusters matched)")
    
    candidates = [f"h_{i}" for i in range(21, 38)]

    # Iterate and print recall values
    # for node in candidates:
    #     try:
    #         result = cluster_recall_when_rooted_at(T_rooted, T_unrooted, node)
    #         print(f"{node}: recall = {result['recall']:.3f} "
    #             f"({result['n_common']}/{result['n_ref']} clusters matched)")
    #     except Exception as e:
    #         print(f"{node}: ERROR -> {e}")


    # Generate node names h_21 ... h_37
    candidates = [f"h_{i}" for i in range(21, 38)]

    # Collect results
    recalls = []
    for node in candidates:
        try:
            result = cluster_recall_when_rooted_at(T_rooted, T_unrooted, node)
            recalls.append(result["recall"])
        except Exception as e:
            print(f"{node}: ERROR -> {e}")
            recalls.append(0.0)  # fallback so plot aligns

    # Plot
    # plt.figure(figsize=(10, 5))
    # plt.bar(candidates, recalls, color="skyblue", edgecolor="black")
    # plt.xticks(rotation=45, ha="right")
    # plt.ylabel("Recall")
    # plt.xlabel("Root candidate node")
    # plt.title("Cluster recall for different rootings of T_unrooted")
    # plt.ylim(0, 1)  # recall is always between 0 and 1
    # plt.tight_layout()
    # plt.show()

    # # Save to file
    # plt.savefig("recall_distribution.PNG", dpi=300)
    # plt.close()

def empirical_recall_pvalue(
    T_rooted: nx.DiGraph,
    T_unrooted: nx.Graph,
    observed_recall: float,
    n_samples: int | None = None,
    seed: int | None = None,
) -> float:
    """
    Compute an empirical non-parametric p-value for an observed recall value.

    Parameters
    ----------
    T_rooted : nx.DiGraph
        Reference rooted tree (edges parent->child).
    T_unrooted : nx.Graph
        Unrooted tree (undirected).
    observed_recall : float
        The recall value obtained by a specific rooting method.
    n_samples : int, optional
        If None, use all internal nodes of T_unrooted.
        Otherwise, sample this many internal nodes at random.
    seed : int, optional
        Random seed for reproducibility when sampling.

    Returns
    -------
    float
        Empirical p-value = fraction of rootings with recall >= observed_recall.
    """
    rng = np.random.default_rng(seed)

    # internal nodes = degree > 1
    internal_nodes = [n for n in T_unrooted.nodes if T_unrooted.degree(n) > 1]
    if n_samples is not None and n_samples < len(internal_nodes):
        internal_nodes = rng.choice(internal_nodes, size=n_samples, replace=False)

    recalls: List[float] = []
    for node in internal_nodes:
        try:
            result = cluster_recall_when_rooted_at(T_rooted, T_unrooted, node)
            recalls.append(result["recall"])
        except Exception:
            continue  # skip invalid rootings if any

    if not recalls:
        raise RuntimeError("No recall values could be computed.")

    # empirical p-value: P(R >= observed)
    more_extreme = sum(r >= observed_recall for r in recalls)
    pval = (more_extreme + 1) / (len(recalls) + 1)  # add-one correction

    return pval


print (empirical_recall_pvalue(T_rooted, T_unrooted,0.82,1000,33))

