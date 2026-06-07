//! Circulation with demands.
//!
//! Determines if a feasible circulation exists in a network
//! where each edge has a capacity and each vertex has a demand/supply.

use crate::{FlowNetwork, max_flow};

/// A circulation network with vertex demands.
#[derive(Clone, Debug)]
pub struct CirculationNetwork {
    n: usize,
    /// (from, to, lower_bound, upper_bound)
    edges: Vec<(usize, usize, f64, f64)>,
    /// Demand at each vertex. Positive = demand (needs flow in), negative = supply.
    demands: Vec<f64>,
}

impl CirculationNetwork {
    /// Create a new circulation network with `n` vertices.
    pub fn new(n: usize) -> Self {
        Self {
            n,
            edges: Vec::new(),
            demands: vec![0.0; n],
        }
    }

    /// Set the demand at vertex `v`. Positive = vertex needs flow, negative = vertex supplies flow.
    pub fn set_demand(&mut self, v: usize, demand: f64) {
        assert!(v < self.n, "Vertex index out of bounds");
        self.demands[v] = demand;
    }

    /// Add an edge with lower and upper capacity bounds.
    pub fn add_edge(&mut self, from: usize, to: usize, lower: f64, upper: f64) {
        assert!(from < self.n && to < self.n, "Vertex index out of bounds");
        assert!(lower <= upper, "Lower bound must not exceed upper bound");
        self.edges.push((from, to, lower, upper));
    }

    /// Number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.n
    }
}

/// Result of a circulation feasibility check.
#[derive(Debug, Clone)]
pub struct CirculationResult {
    /// Whether a feasible circulation exists.
    pub feasible: bool,
    /// Flow on each edge (satisfying lower/upper bounds and demands).
    pub flows: Vec<(usize, usize, f64)>,
}

/// Check if a feasible circulation exists.
///
/// Uses the standard reduction to max-flow with super-source and super-sink.
pub fn find_circulation(network: &CirculationNetwork) -> CirculationResult {
    let n = network.n;

    // Check demand balance
    let total_demand: f64 = network.demands.iter().sum();
    if total_demand.abs() > 1e-10 {
        return CirculationResult {
            feasible: false,
            flows: vec![],
        };
    }

    // Create flow network with super-source and super-sink
    let source = n;
    let sink = n + 1;
    let flow_n = n + 2;
    let mut net = FlowNetwork::new(flow_n);

    let mut excess = vec![0.0; n];
    let mut lower_flows = Vec::new();

    for &(from, to, lower, upper) in &network.edges {
        // Reduce to standard flow: capacity = upper - lower
        if upper - lower > 1e-10 {
            net.add_edge(from, to, upper - lower);
        }
        // Adjust demands
        excess[from] -= lower;
        excess[to] += lower;
        lower_flows.push((from, to, lower));
    }

    // Add demand edges to/from super-source/sink
    // Convention: demand[v] > 0 means vertex needs flow IN
    //             demand[v] < 0 means vertex has supply (flow OUT)
    // Standard reduction: supply vertices get flow from super-source,
    // demand vertices send flow to super-sink
    let mut required_flow = 0.0;
    #[allow(clippy::needless_range_loop)]
    for v in 0..n {
        let adjusted = network.demands[v] + excess[v];
        if adjusted < -1e-10 {
            // Vertex has net supply: super-source -> v
            net.add_edge(source, v, -adjusted);
            required_flow += -adjusted;
        } else if adjusted > 1e-10 {
            // Vertex has net demand: v -> super-sink
            net.add_edge(v, sink, adjusted);
        }
    }

    let result = max_flow::ford_fulkerson(&net, source, sink);

    let feasible = (result.max_flow - required_flow).abs() < 1e-6;

    // Build actual flows
    let mut flows = Vec::new();
    for (idx, &(from, to, _lower, _upper)) in network.edges.iter().enumerate() {
        let base_flow = lower_flows[idx].2;
        // Find the additional flow from the result
        let additional = result
            .flows
            .iter()
            .filter(|&&(u, v, f)| u == from && v == to && f > 0.0)
            .map(|&(_, _, f)| f)
            .sum::<f64>();
        flows.push((from, to, base_flow + additional));
    }

    CirculationResult { feasible, flows }
}

/// Check feasibility without computing the actual flow.
pub fn is_feasible(network: &CirculationNetwork) -> bool {
    find_circulation(network).feasible
}

/// Compute the minimum feasible flow value that satisfies all constraints.
pub fn total_flow(network: &CirculationNetwork) -> f64 {
    let result = find_circulation(network);
    if result.feasible {
        result.flows.iter().map(|(_, _, f)| f).sum()
    } else {
        -1.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_feasible_circulation() {
        let mut cn = CirculationNetwork::new(4);
        cn.add_edge(0, 1, 2.0, 5.0);
        cn.add_edge(1, 2, 1.0, 4.0);
        cn.add_edge(2, 3, 2.0, 6.0);
        cn.add_edge(3, 0, 1.0, 3.0);
        // No demands: just check if circulation exists
        let result = find_circulation(&cn);
        assert!(result.feasible);
    }

    #[test]
    fn test_infeasible_demands() {
        let mut cn = CirculationNetwork::new(2);
        cn.set_demand(0, 10.0);
        cn.set_demand(1, 5.0);
        // Demands don't balance
        assert!(!is_feasible(&cn));
    }

    #[test]
    fn test_balanced_demands() {
        let mut cn = CirculationNetwork::new(2);
        cn.set_demand(0, -5.0); // supply
        cn.set_demand(1, 5.0);  // demand
        cn.add_edge(0, 1, 0.0, 10.0);
        let result = find_circulation(&cn);
        assert!(result.feasible, "Should be feasible: {:?}", result);
    }

    #[test]
    fn test_capacity_constraint() {
        let mut cn = CirculationNetwork::new(2);
        cn.set_demand(0, -10.0);
        cn.set_demand(1, 10.0);
        cn.add_edge(0, 1, 0.0, 5.0); // capacity too low
        assert!(!is_feasible(&cn));
    }

    #[test]
    fn test_lower_bound() {
        let mut cn = CirculationNetwork::new(3);
        cn.add_edge(0, 1, 5.0, 10.0);
        cn.add_edge(1, 2, 5.0, 10.0);
        cn.add_edge(2, 0, 5.0, 10.0);
        let result = find_circulation(&cn);
        assert!(result.feasible);
        // All flows should be at least 5
        for (_, _, f) in &result.flows {
            assert!(*f >= 5.0 - 1e-6);
        }
    }

    #[test]
    fn test_empty_network() {
        let cn = CirculationNetwork::new(0);
        let result = find_circulation(&cn);
        assert!(result.feasible);
    }

    #[test]
    fn test_is_feasible() {
        let mut cn = CirculationNetwork::new(3);
        cn.add_edge(0, 1, 0.0, 10.0);
        cn.add_edge(1, 2, 0.0, 10.0);
        cn.add_edge(2, 0, 0.0, 10.0);
        assert!(is_feasible(&cn));
    }
}
