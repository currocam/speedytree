use std::collections::HashMap;

use petgraph::{graph::UnGraph, stable_graph::NodeIndex};

#[derive(Debug, Clone)]
pub(crate) struct PhyloTree {
    pub tree: crate::Tree,
    pub nodes: HashMap<usize, NodeIndex>,
    n_nodes: usize,
}

impl PhyloTree {
    pub fn build(leafs: &[String]) -> PhyloTree {
        let mut tree: petgraph::Graph<String, f64, petgraph::Undirected> =
            UnGraph::new_undirected();
        let mut nodes = HashMap::new();
        for (index, leaf) in leafs.iter().enumerate() {
            let node: NodeIndex = tree.add_node(leaf.to_owned());
            nodes.insert(index, node);
        }
        let n_nodes = leafs.len();
        PhyloTree {
            tree,
            nodes,
            n_nodes,
        }
    }

    pub fn merge_neighbors(&mut self, a: usize, b: usize, dau: f64, dbu: f64) -> NodeIndex {
        // Get nodes to merge
        let a_node = self.nodes[&a];
        let b_node = self.nodes[&b];
        let u = self.tree.add_node("".to_string());

        self.nodes.insert(self.n_nodes, u);
        self.n_nodes += 1;

        // Add new edges
        self.tree.add_edge(u, a_node, dau);
        self.tree.add_edge(u, b_node, dbu);
        u
    }
}
