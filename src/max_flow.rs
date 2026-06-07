//! Ford-Fulkerson maximum flow algorithm using DFS augmenting paths.
//!
//! Computes the maximum flow in a flow network using the classic
//! Ford-Fulkerson method with DFS-based path augmentation.

use crate::FlowNetwork;

/// Result of a max flow computation.
#[derive(Debug, Clone)]
pub struct MaxFlowResult {
    /// Maximum flow value.
    pub max_flow: f64,
    /// Flow on each edge (from, to, flow).
    pub flows: Vec<(usize, usize, f64)>,
}

/// Compute the maximum flow from `source` to `sink` using Ford-Fulkerson with DFS.
///
/// Returns the max flow value and the flow on each edge.
pub fn ford_fulkerson(network: &FlowNetwork, source: usize, sink: usize) -> MaxFlowResult {
    let n = network.vertex_count();
    if source == sink || n == 0 {
        return MaxFlowResult {
            max_flow: 0.0,
            flows: vec![],
        };
    }

    let mut adj = network.adjacency().to_vec();
    let mut total_flow = 0.0;
    let max_iter = n * n * 100; // Safety bound

    for _ in 0..max_iter {
        let mut visited = vec![false; n];
        let pushed = dfs_augment(&mut adj, source, sink, f64::INFINITY, &mut visited);
        if pushed < 1e-10 {
            break;
        }
        total_flow += pushed;
    }

    // Collect flow values
    let mut flows = Vec::new();
    #[allow(clippy::needless_range_loop)]
    for u in 0..n {
        for &(v, original_cap, _rev_idx) in &network.adjacency()[u] {
            let current_cap = adj[u].iter().find(|&&(vv, _, _)| vv == v).map(|&(_, c, _)| c).unwrap_or(0.0);
            let flow = original_cap - current_cap;
            if flow.abs() > 1e-10 {
                flows.push((u, v, flow));
            }
        }
    }

    MaxFlowResult {
        max_flow: total_flow,
        flows,
    }
}

fn dfs_augment(
    adj: &mut [Vec<(usize, f64, usize)>],
    u: usize,
    sink: usize,
    min_cap: f64,
    visited: &mut [bool],
) -> f64 {
    if u == sink {
        return min_cap;
    }
    visited[u] = true;

    for i in 0..adj[u].len() {
        let (v, cap, rev) = adj[u][i];
        if cap > 1e-10 && !visited[v] {
            let pushed = dfs_augment(adj, v, sink, min_cap.min(cap), visited);
            if pushed > 1e-10 {
                adj[u][i].1 -= pushed;
                adj[v][rev].1 += pushed;
                return pushed;
            }
        }
    }

    0.0
}

/// Compute the minimum s-t cut value (equals the max flow by the max-flow min-cut theorem).
pub fn min_cut_value(network: &FlowNetwork, source: usize, sink: usize) -> f64 {
    ford_fulkerson(network, source, sink).max_flow
}

/// Find the minimum s-t cut: returns (reachable vertices, unreachable vertices, cut capacity).
pub fn min_cut(network: &FlowNetwork, source: usize, sink: usize) -> (Vec<usize>, Vec<usize>, f64) {
    let n = network.vertex_count();
    if source == sink || n == 0 {
        return (vec![], (0..n).collect(), 0.0);
    }

    let result = ford_fulkerson(network, source, sink);

    // Build residual graph
    let mut adj = network.adjacency().to_vec();
    for &(u, v, f) in &result.flows {
        for entry in &mut adj[u] {
            if entry.0 == v {
                entry.1 -= f;
            }
        }
        for entry in &mut adj[v] {
            if entry.0 == u {
                entry.1 += f;
            }
        }
    }

    let mut visited = vec![false; n];
    let mut queue = std::collections::VecDeque::new();
    queue.push_back(source);
    visited[source] = true;

    while let Some(u) = queue.pop_front() {
        for &(v, cap, _) in &adj[u] {
            if cap > 1e-10 && !visited[v] {
                visited[v] = true;
                queue.push_back(v);
            }
        }
    }

    let reachable: Vec<usize> = (0..n).filter(|&i| visited[i]).collect();
    let unreachable: Vec<usize> = (0..n).filter(|&i| !visited[i]).collect();

    (reachable, unreachable, result.max_flow)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::FlowNetwork;

    #[test]
    fn test_simple_max_flow() {
        let mut net = FlowNetwork::new(2);
        net.add_edge(0, 1, 10.0);
        let result = ford_fulkerson(&net, 0, 1);
        assert!((result.max_flow - 10.0).abs() < 1e-6);
    }

    #[test]
    fn test_path_flow() {
        let mut net = FlowNetwork::new(4);
        net.add_edge(0, 1, 10.0);
        net.add_edge(1, 2, 5.0);
        net.add_edge(2, 3, 10.0);
        let result = ford_fulkerson(&net, 0, 3);
        assert!((result.max_flow - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_parallel_paths() {
        let mut net = FlowNetwork::new(4);
        net.add_edge(0, 1, 10.0);
        net.add_edge(0, 2, 10.0);
        net.add_edge(1, 3, 10.0);
        net.add_edge(2, 3, 10.0);
        let result = ford_fulkerson(&net, 0, 3);
        assert!((result.max_flow - 20.0).abs() < 1e-6);
    }

    #[test]
    fn test_bottleneck() {
        let mut net = FlowNetwork::new(4);
        net.add_edge(0, 1, 100.0);
        net.add_edge(1, 2, 1.0);
        net.add_edge(2, 3, 100.0);
        let result = ford_fulkerson(&net, 0, 3);
        assert!((result.max_flow - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_no_path() {
        let mut net = FlowNetwork::new(3);
        net.add_edge(0, 1, 10.0);
        let result = ford_fulkerson(&net, 0, 2);
        assert!((result.max_flow).abs() < 1e-6);
    }

    #[test]
    fn test_min_cut_value() {
        let mut net = FlowNetwork::new(4);
        net.add_edge(0, 1, 10.0);
        net.add_edge(0, 2, 10.0);
        net.add_edge(1, 3, 10.0);
        net.add_edge(2, 3, 10.0);
        let cut = min_cut_value(&net, 0, 3);
        assert!((cut - 20.0).abs() < 1e-6);
    }

    #[test]
    fn test_min_cut_partition() {
        let mut net = FlowNetwork::new(4);
        net.add_edge(0, 1, 10.0);
        net.add_edge(0, 2, 10.0);
        net.add_edge(1, 3, 10.0);
        net.add_edge(2, 3, 10.0);
        let (reachable, unreachable, cut_val) = min_cut(&net, 0, 3);
        assert!(reachable.contains(&0));
        assert!(unreachable.contains(&3));
        assert!((cut_val - 20.0).abs() < 1e-6);
    }

    #[test]
    fn test_source_equals_sink() {
        let net = FlowNetwork::new(2);
        let result = ford_fulkerson(&net, 0, 0);
        assert!((result.max_flow).abs() < 1e-6);
    }

    #[test]
    fn test_single_vertex() {
        let net = FlowNetwork::new(1);
        let result = ford_fulkerson(&net, 0, 0);
        assert!((result.max_flow).abs() < 1e-6);
    }
}
