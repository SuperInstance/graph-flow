//! Minimum-cost maximum flow algorithm.
//!
//! Computes a flow that achieves maximum throughput at minimum total cost,
//! using successive shortest path augmentation with Bellman-Ford.


/// An edge with both capacity and cost.
#[derive(Clone, Debug)]
pub struct CostEdge {
    pub from: usize,
    pub to: usize,
    pub capacity: f64,
    pub cost: f64,
}

/// Result of min-cost max-flow computation.
#[derive(Debug, Clone)]
pub struct MinCostFlowResult {
    /// Total flow sent.
    pub total_flow: f64,
    /// Total cost of the flow.
    pub total_cost: f64,
}

/// A min-cost flow network.
#[derive(Clone, Debug)]
pub struct MinCostNetwork {
    n: usize,
    /// adj[u] = [(v, capacity, cost, rev_index)]
    adj: Vec<Vec<(usize, f64, f64, usize)>>,
}

impl MinCostNetwork {
    /// Create a new min-cost flow network with `n` vertices.
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

    /// Add a directed edge with capacity and cost per unit.
    pub fn add_edge(&mut self, from: usize, to: usize, capacity: f64, cost: f64) {
        assert!(from < self.n && to < self.n, "Vertex index out of bounds");
        let rev_from = self.adj[to].len();
        let rev_to = self.adj[from].len();
        self.adj[from].push((to, capacity, cost, rev_from));
        self.adj[to].push((from, 0.0, -cost, rev_to));
    }
}

/// Compute min-cost max-flow using successive shortest paths with Bellman-Ford.
pub fn min_cost_flow(
    network: &MinCostNetwork,
    source: usize,
    sink: usize,
    max_flow_target: f64,
) -> MinCostFlowResult {
    let n = network.vertex_count();
    let mut adj = network.adj.clone();
    let mut total_flow = 0.0;
    let mut total_cost = 0.0;

    while total_flow < max_flow_target - 1e-10 {
        // Bellman-Ford to find shortest path (minimum cost augmenting path)
        let mut dist = vec![f64::INFINITY; n];
        let mut parent = vec![None::<(usize, usize)>; n];
        dist[source] = 0.0;

        for _ in 0..n {
            let mut updated = false;
            for u in 0..n {
                if dist[u].is_infinite() {
                    continue;
                }
                #[allow(clippy::needless_range_loop)]
                for i in 0..adj[u].len() {
                    let (v, cap, cost, _) = adj[u][i];
                    if cap > 1e-10 && dist[u] + cost < dist[v] - 1e-10 {
                        dist[v] = dist[u] + cost;
                        parent[v] = Some((u, i));
                        updated = true;
                    }
                }
            }
            if !updated {
                break;
            }
        }

        if dist[sink].is_infinite() {
            break;
        }

        // Find bottleneck
        let mut bottleneck = max_flow_target - total_flow;
        let mut v = sink;
        while let Some((u, i)) = parent[v] {
            bottleneck = bottleneck.min(adj[u][i].1);
            v = u;
        }

        // Augment
        v = sink;
        while let Some((u, i)) = parent[v] {
            let rev = adj[u][i].3;
            adj[u][i].1 -= bottleneck;
            adj[v][rev].1 += bottleneck;
            total_cost += bottleneck * adj[u][i].2;
            v = u;
        }

        total_flow += bottleneck;
    }

    MinCostFlowResult {
        total_flow,
        total_cost,
    }
}

/// Compute the cost of a given flow (if feasible).
pub fn compute_flow_cost(network: &MinCostNetwork, source: usize, sink: usize) -> (f64, f64) {
    let result = min_cost_flow(network, source, sink, f64::INFINITY);
    (result.total_flow, result.total_cost)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_min_cost_flow() {
        let mut net = MinCostNetwork::new(2);
        net.add_edge(0, 1, 10.0, 1.0);
        let result = min_cost_flow(&net, 0, 1, 10.0);
        assert!((result.total_flow - 10.0).abs() < 1e-6);
        assert!((result.total_cost - 10.0).abs() < 1e-4);
    }

    #[test]
    fn test_two_paths_cheaper_first() {
        let mut net = MinCostNetwork::new(4);
        net.add_edge(0, 1, 5.0, 1.0); // cheap path
        net.add_edge(1, 3, 5.0, 1.0);
        net.add_edge(0, 2, 5.0, 3.0); // expensive path
        net.add_edge(2, 3, 5.0, 3.0);
        let result = min_cost_flow(&net, 0, 3, 10.0);
        assert!((result.total_flow - 10.0).abs() < 1e-4);
        // Cost: 5 * 2 (cheap) + 5 * 6 (expensive) = 10 + 30 = 40
        assert!((result.total_cost - 40.0).abs() < 1e-2);
    }

    #[test]
    fn test_limited_capacity() {
        let mut net = MinCostNetwork::new(3);
        net.add_edge(0, 1, 3.0, 1.0);
        net.add_edge(1, 2, 3.0, 1.0);
        let result = min_cost_flow(&net, 0, 2, 10.0);
        assert!((result.total_flow - 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_no_path() {
        let net = MinCostNetwork::new(3);
        let result = min_cost_flow(&net, 0, 2, 10.0);
        assert!((result.total_flow).abs() < 1e-6);
    }

    #[test]
    fn test_single_edge_zero_cost() {
        let mut net = MinCostNetwork::new(2);
        net.add_edge(0, 1, 5.0, 0.0);
        let result = min_cost_flow(&net, 0, 1, 5.0);
        assert!((result.total_flow - 5.0).abs() < 1e-6);
        assert!((result.total_cost).abs() < 1e-4);
    }

    #[test]
    fn test_compute_flow_cost() {
        let mut net = MinCostNetwork::new(3);
        net.add_edge(0, 1, 10.0, 2.0);
        net.add_edge(1, 2, 10.0, 3.0);
        let (flow, cost) = compute_flow_cost(&net, 0, 2);
        assert!((flow - 10.0).abs() < 1e-4);
        assert!((cost - 50.0).abs() < 1e-2);
    }
}
