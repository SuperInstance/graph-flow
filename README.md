# graph-flow

Network flow algorithms for Rust. Pure `std` — no external dependencies.

## Features

- **Ford-Fulkerson** — Maximum flow via DFS augmenting paths
- **Edmonds-Karp** — BFS-based max flow with O(VE²) guarantee
- **Min-cost max-flow** — Successive shortest paths with Bellman-Ford
- **Bipartite matching** — Maximum matching via flow reduction, König's theorem
- **Circulation with demands** — Feasibility checking with lower/upper bounds

## Usage

```rust
use graph_flow::{FlowNetwork, max_flow, edmonds_karp};

let mut net = FlowNetwork::new(4);
net.add_edge(0, 1, 10.0);
net.add_edge(0, 2, 10.0);
net.add_edge(1, 3, 5.0);
net.add_edge(2, 3, 15.0);

let result = max_flow::ford_fulkerson(&net, 0, 3);
println!("Max flow: {}", result.max_flow);

let flow = edmonds_karp::edmonds_karp(&net, 0, 3);
println!("Edmonds-Karp flow: {flow}");
```

## License

MIT
