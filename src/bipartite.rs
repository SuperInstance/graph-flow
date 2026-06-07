//! Bipartite matching via flow reduction.
//!
//! Maximum bipartite matching by reducing to a max-flow problem
//! with a super-source and super-sink.

use crate::{FlowNetwork, max_flow};

/// A bipartite graph with left and right vertex sets.
#[derive(Clone, Debug)]
pub struct BipartiteGraph {
    /// Number of left vertices.
    left_size: usize,
    /// Number of right vertices.
    right_size: usize,
    /// Edges from left to right.
    edges: Vec<(usize, usize)>,
}

impl BipartiteGraph {
    /// Create a new bipartite graph with `left_size` left vertices and `right_size` right vertices.
    pub fn new(left_size: usize, right_size: usize) -> Self {
        Self {
            left_size,
            right_size,
            edges: Vec::new(),
        }
    }

    /// Add an edge from left vertex `u` to right vertex `v`.
    pub fn add_edge(&mut self, u: usize, v: usize) {
        assert!(u < self.left_size, "Left vertex out of bounds");
        assert!(v < self.right_size, "Right vertex out of bounds");
        self.edges.push((u, v));
    }

    /// Get the number of left vertices.
    pub fn left_size(&self) -> usize {
        self.left_size
    }

    /// Get the number of right vertices.
    pub fn right_size(&self) -> usize {
        self.right_size
    }
}

/// Result of bipartite matching.
#[derive(Debug, Clone)]
pub struct MatchingResult {
    /// Number of matched pairs.
    pub matching_size: usize,
    /// Matched pairs (left, right).
    pub matches: Vec<(usize, usize)>,
}

/// Compute maximum bipartite matching via flow reduction.
///
/// Creates a flow network with super-source connected to all left vertices
/// (capacity 1), edges from left to right (capacity 1), and all right vertices
/// connected to super-sink (capacity 1).
pub fn maximum_matching(graph: &BipartiteGraph) -> MatchingResult {
    if graph.left_size == 0 || graph.right_size == 0 {
        return MatchingResult {
            matching_size: 0,
            matches: vec![],
        };
    }

    let n = graph.left_size + graph.right_size + 2;
    let source = n - 2;
    let sink = n - 1;

    let mut net = FlowNetwork::new(n);

    // Source to left vertices
    for u in 0..graph.left_size {
        net.add_edge(source, u, 1.0);
    }

    // Right vertices to sink
    for v in 0..graph.right_size {
        net.add_edge(graph.left_size + v, sink, 1.0);
    }

    // Left to right edges
    for &(u, v) in &graph.edges {
        net.add_edge(u, graph.left_size + v, 1.0);
    }

    let result = max_flow::ford_fulkerson(&net, source, sink);

    let matches: Vec<(usize, usize)> = result
        .flows
        .iter()
        .filter(|&&(u, v, f)| {
            u < graph.left_size && v >= graph.left_size && v < graph.left_size + graph.right_size && f > 0.5
        })
        .map(|&(u, v, _)| (u, v - graph.left_size))
        .collect();

    MatchingResult {
        matching_size: matches.len(),
        matches,
    }
}

/// Check if a bipartite graph has a perfect matching.
pub fn has_perfect_matching(graph: &BipartiteGraph) -> bool {
    let result = maximum_matching(graph);
    let min_side = graph.left_size.min(graph.right_size);
    result.matching_size == min_side
}

/// Compute a maximum independent set in a bipartite graph using König's theorem.
///
/// The maximum independent set size = n - maximum matching size.
pub fn maximum_independent_set_size(graph: &BipartiteGraph) -> usize {
    let result = maximum_matching(graph);
    graph.left_size + graph.right_size - result.matching_size
}

/// Compute the minimum vertex cover size (equals maximum matching by König's theorem).
pub fn minimum_vertex_cover_size(graph: &BipartiteGraph) -> usize {
    maximum_matching(graph).matching_size
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_perfect_matching() {
        let mut bg = BipartiteGraph::new(3, 3);
        bg.add_edge(0, 0);
        bg.add_edge(1, 1);
        bg.add_edge(2, 2);
        let result = maximum_matching(&bg);
        assert_eq!(result.matching_size, 3);
    }

    #[test]
    fn test_no_matching() {
        let bg = BipartiteGraph::new(2, 2);
        let result = maximum_matching(&bg);
        assert_eq!(result.matching_size, 0);
    }

    #[test]
    fn test_partial_matching() {
        let mut bg = BipartiteGraph::new(3, 3);
        bg.add_edge(0, 0);
        bg.add_edge(0, 1);
        bg.add_edge(1, 1);
        let result = maximum_matching(&bg);
        assert!(result.matching_size >= 2);
    }

    #[test]
    fn test_unequal_sides() {
        let mut bg = BipartiteGraph::new(2, 5);
        bg.add_edge(0, 0);
        bg.add_edge(1, 1);
        let result = maximum_matching(&bg);
        assert_eq!(result.matching_size, 2);
    }

    #[test]
    fn test_has_perfect_matching() {
        let mut bg = BipartiteGraph::new(2, 2);
        bg.add_edge(0, 0);
        bg.add_edge(1, 1);
        assert!(has_perfect_matching(&bg));
    }

    #[test]
    fn test_no_perfect_matching() {
        let mut bg = BipartiteGraph::new(2, 3);
        bg.add_edge(0, 0);
        assert!(!has_perfect_matching(&bg));
    }

    #[test]
    fn test_min_vertex_cover() {
        let mut bg = BipartiteGraph::new(3, 3);
        bg.add_edge(0, 0);
        bg.add_edge(1, 1);
        bg.add_edge(2, 2);
        assert_eq!(minimum_vertex_cover_size(&bg), 3);
    }

    #[test]
    fn test_max_independent_set() {
        let mut bg = BipartiteGraph::new(2, 2);
        bg.add_edge(0, 0);
        bg.add_edge(1, 1);
        let mis = maximum_independent_set_size(&bg);
        assert_eq!(mis, 2); // 4 - 2 = 2
    }
}
