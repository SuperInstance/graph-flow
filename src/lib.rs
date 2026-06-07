//! # graph-flow
//!
//! Network flow algorithms for Rust. Pure `std` — no external dependencies.
//!
//! ## Modules
//!
//! - [`max_flow`] — Ford-Fulkerson maximum flow
//! - [`edmonds_karp`] — Edmonds-Karp BFS-based max flow
//! - [`min_cost`] — Minimum-cost maximum flow
//! - [`bipartite`] — Bipartite matching via flow reduction
//! - [`circulation`] — Circulation with demands

pub mod max_flow;
pub mod edmonds_karp;
pub mod min_cost;
pub mod bipartite;
pub mod circulation;

/// A flow network represented as an adjacency list with capacity edges.
#[derive(Clone, Debug)]
pub struct FlowNetwork {
    n: usize,
    /// Each entry contains (to, capacity, rev_index) where rev_index is the
    /// index of the reverse edge in adj[to].
    adj: Vec<Vec<(usize, f64, usize)>>,
}

impl FlowNetwork {
    /// Create a new flow network with `n` vertices.
    pub fn new(n: usize) -> Self {
        Self {
            n,
            adj: vec![vec![]; n],
        }
    }

    /// Number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.n
    }

    /// Add a directed edge with given capacity (and a reverse edge with capacity 0).
    pub fn add_edge(&mut self, from: usize, to: usize, capacity: f64) {
        assert!(from < self.n && to < self.n, "Vertex index out of bounds");
        let rev_from = self.adj[to].len();
        let rev_to = self.adj[from].len();
        self.adj[from].push((to, capacity, rev_from));
        self.adj[to].push((from, 0.0, rev_to));
    }

    /// Get the adjacency list (for reading).
    pub fn adjacency(&self) -> &[Vec<(usize, f64, usize)>] {
        &self.adj
    }

    /// Get a mutable reference to the adjacency list.
    pub fn adjacency_mut(&mut self) -> &mut Vec<Vec<(usize, f64, usize)>> {
        &mut self.adj
    }
}
