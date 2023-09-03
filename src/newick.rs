use petgraph::prelude::DiGraph;
use petgraph::Direction::{Incoming, Outgoing};
use petgraph::Graph;

#[derive(Debug)]
pub struct PhylogeneticTree<'a> {
    graph: Graph<&'a str, f64>,
}

impl PhylogeneticTree<'_> {
    pub fn into_newick(self) -> Option<String> {
        let root = self.graph.node_indices().next()?;
        let str = from_graph_to_newick(&self.graph, root);
        Some(format!("{};", str))
    }
    pub fn empty() -> Self {
        PhylogeneticTree {
            graph: Graph::new(),
        }
    }
}

fn from_graph_to_newick(
    graph: &DiGraph<&str, f64>,
    node: petgraph::stable_graph::NodeIndex,
) -> String {
    if graph.neighbors_directed(node, Outgoing).count() == 0 {
        let edge: Option<petgraph::graph::EdgeReference<'_, f64>> =
            graph.edges_directed(node, Incoming).next();
        if let Some(edge) = edge {
            return format!("{}:{:.1}", graph[node], edge.weight());
        }
    }
    let mut children = graph
        .neighbors_directed(node, Outgoing)
        .map(|child| from_graph_to_newick(graph, child))
        .collect::<Vec<String>>();
    children.reverse();
    let children = children.join(",");

    let edge = graph.edges_directed(node, Incoming).next();
    if let Some(edge) = edge {
        return format!("({}):{:.1}", children, edge.weight());
    } else {
        return format!("({})", children);
    }
}

mod tests {
    use super::*;
    #[test]
    fn test_happy_path() {
        let mut graph = Graph::new();
        let root = graph.add_node("");
        let node_a = graph.add_node("A");
        let node_b = graph.add_node("B");
        let internal_node = graph.add_node("");
        let node_c = graph.add_node("C");
        let node_d = graph.add_node("D");
        graph.extend_with_edges(&[
            (root, node_a, 0.1),
            (root, node_b, 0.2),
            (root, internal_node, 0.5),
            (internal_node, node_c, 0.3),
            (internal_node, node_d, 0.4),
        ]);
        let phylo = PhylogeneticTree { graph };
        assert_eq!(
            phylo.into_newick().unwrap_or_default(),
            "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
        )
    }
}
