//! Edmonds-Karp maximum flow algorithm using BFS.
//!
//! Computes maximum flow using shortest augmenting paths (BFS),
//! guaranteeing O(VE²) time complexity.

use crate::FlowNetwork;

/// Compute the maximum flow using the Edmonds-Karp algorithm.
///
/// Uses BFS to find shortest augmenting paths, providing polynomial time bounds.
pub fn edmonds_karp(network: &FlowNetwork, source: usize, sink: usize) -> f64 {
    let n = network.vertex_count();
    if source == sink || n == 0 {
        return 0.0;
    }

    let mut adj = network.adjacency().to_vec();
    let mut total_flow = 0.0;

    loop {
        // BFS to find augmenting path
        let mut parent = vec![None::<(usize, usize)>; n];
        let mut visited = vec![false; n];
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(source);
        visited[source] = true;

        while let Some(u) = queue.pop_front() {
            if u == sink {
                break;
            }
            #[allow(clippy::needless_range_loop)]
            for i in 0..adj[u].len() {
                let (v, cap, _) = adj[u][i];
                if cap > 1e-10 && !visited[v] {
                    visited[v] = true;
                    parent[v] = Some((u, i));
                    queue.push_back(v);
                }
            }
        }

        if !visited[sink] {
            break;
        }

        let mut bottleneck = f64::INFINITY;
        let mut v = sink;
        while let Some((u, i)) = parent[v] {
            bottleneck = bottleneck.min(adj[u][i].1);
            v = u;
        }

        if bottleneck < 1e-10 {
            break;
        }

        v = sink;
        while let Some((u, i)) = parent[v] {
            let rev = adj[u][i].2;
            adj[u][i].1 -= bottleneck;
            adj[v][rev].1 += bottleneck;
            v = u;
        }

        total_flow += bottleneck;
    }

    total_flow
}

/// Find the minimum cut using Edmonds-Karp.
///
/// Returns the reachable set from source in the residual graph after max flow.
pub fn min_cut_edmonds_karp(
    network: &FlowNetwork,
    source: usize,
    sink: usize,
) -> (Vec<usize>, f64) {
    let n = network.vertex_count();
    if source == sink || n == 0 {
        return (vec![], 0.0);
    }

    let mut adj = network.adjacency().to_vec();

    // Run EK to build residual graph
    let mut total_flow = 0.0;
    loop {
        let mut parent = vec![None::<(usize, usize)>; n];
        let mut visited = vec![false; n];
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(source);
        visited[source] = true;

        while let Some(u) = queue.pop_front() {
            if u == sink {
                break;
            }
            #[allow(clippy::needless_range_loop)]
            for i in 0..adj[u].len() {
                let (v, cap, _) = adj[u][i];
                if cap > 1e-10 && !visited[v] {
                    visited[v] = true;
                    parent[v] = Some((u, i));
                    queue.push_back(v);
                }
            }
        }

        if !visited[sink] {
            break;
        }

        let mut bottleneck = f64::INFINITY;
        let mut v = sink;
        while let Some((u, i)) = parent[v] {
            bottleneck = bottleneck.min(adj[u][i].1);
            v = u;
        }

        if bottleneck < 1e-10 {
            break;
        }

        v = sink;
        while let Some((u, i)) = parent[v] {
            let rev = adj[u][i].2;
            adj[u][i].1 -= bottleneck;
            adj[v][rev].1 += bottleneck;
            v = u;
        }

        total_flow += bottleneck;
    }

    // BFS in residual graph to find reachable
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
    (reachable, total_flow)
}

/// Compute the edge connectivity between source and sink.
pub fn edge_connectivity(network: &FlowNetwork, source: usize, sink: usize) -> usize {
    let flow = edmonds_karp(network, source, sink);
    flow.round() as usize
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::FlowNetwork;

    #[test]
    fn test_simple_flow() {
        let mut net = FlowNetwork::new(2);
        net.add_edge(0, 1, 10.0);
        let flow = edmonds_karp(&net, 0, 1);
        assert!((flow - 10.0).abs() < 1e-6);
    }

    #[test]
    fn test_path_bottleneck() {
        let mut net = FlowNetwork::new(4);
        net.add_edge(0, 1, 10.0);
        net.add_edge(1, 2, 5.0);
        net.add_edge(2, 3, 10.0);
        let flow = edmonds_karp(&net, 0, 3);
        assert!((flow - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_parallel_paths() {
        let mut net = FlowNetwork::new(4);
        net.add_edge(0, 1, 10.0);
        net.add_edge(0, 2, 10.0);
        net.add_edge(1, 3, 10.0);
        net.add_edge(2, 3, 10.0);
        let flow = edmonds_karp(&net, 0, 3);
        assert!((flow - 20.0).abs() < 1e-6);
    }

    #[test]
    fn test_no_path() {
        let mut net = FlowNetwork::new(3);
        net.add_edge(0, 1, 10.0);
        let flow = edmonds_karp(&net, 0, 2);
        assert!((flow).abs() < 1e-6);
    }

    #[test]
    fn test_min_cut() {
        let mut net = FlowNetwork::new(4);
        net.add_edge(0, 1, 10.0);
        net.add_edge(0, 2, 10.0);
        net.add_edge(1, 3, 10.0);
        net.add_edge(2, 3, 10.0);
        let (reachable, flow) = min_cut_edmonds_karp(&net, 0, 3);
        assert!(reachable.contains(&0));
        assert!((flow - 20.0).abs() < 1e-6);
    }

    #[test]
    fn test_edge_connectivity() {
        let mut net = FlowNetwork::new(4);
        net.add_edge(0, 1, 3.0);
        net.add_edge(0, 2, 2.0);
        net.add_edge(1, 3, 2.0);
        net.add_edge(2, 3, 3.0);
        let ec = edge_connectivity(&net, 0, 3);
        assert_eq!(ec, 4);
    }

    #[test]
    fn test_complex_network() {
        let mut net = FlowNetwork::new(6);
        net.add_edge(0, 1, 16.0);
        net.add_edge(0, 2, 13.0);
        net.add_edge(1, 2, 10.0);
        net.add_edge(1, 3, 12.0);
        net.add_edge(2, 1, 4.0);
        net.add_edge(2, 4, 14.0);
        net.add_edge(3, 2, 9.0);
        net.add_edge(3, 5, 20.0);
        net.add_edge(4, 3, 7.0);
        net.add_edge(4, 5, 4.0);
        let flow = edmonds_karp(&net, 0, 5);
        assert!((flow - 23.0).abs() < 1e-3);
    }

    #[test]
    fn test_source_equals_sink() {
        let net = FlowNetwork::new(2);
        let flow = edmonds_karp(&net, 0, 0);
        assert!((flow).abs() < 1e-6);
    }
}
