use std::collections::HashMap;

// Binary tree with edge lengths using petgraph
// graph from petagraph
use petgraph::{graph::UnGraph, stable_graph::NodeIndex};

#[derive(Debug, Clone)]
pub struct PhyloTree {
    pub tree: crate::Tree,
    n_unmerged_leaves: usize,
    pub nodes: HashMap<usize, NodeIndex>,
}

impl PhyloTree {
    pub fn new(tree: crate::Tree, nodes: HashMap<usize, NodeIndex>) -> PhyloTree {
        let n_unmerged_leaves = nodes.len();
        PhyloTree {
            tree,
            n_unmerged_leaves,
            nodes,
        }
    }
    pub fn build(leafs: &Vec<String>) -> PhyloTree {
        let mut tree: petgraph::Graph<String, f64, petgraph::Undirected> =
            UnGraph::new_undirected();
        let mut nodes = HashMap::new();
        for (i, leaf) in leafs.iter().enumerate() {
            let node: NodeIndex = tree.add_node(leaf.to_owned());
            //tree.add_edge(root, node, 0.0);
            nodes.insert(i, node);
        }
        let n_leaves: usize = leafs.len();
        let n_unmerged_leaves = n_leaves;
        PhyloTree {
            tree,
            nodes,
            n_unmerged_leaves,
        }
    }

    pub fn merge_neighbors(&mut self, a: usize, b: usize, dau: f64, dbu: f64) -> NodeIndex {
        // Get nodes to merge
        let n: &usize = &self.n_unmerged_leaves;
        let a_node = self.nodes.remove(&a).unwrap();
        let b_node = self.nodes.remove(&b).unwrap();
        // Swap nodes according to nj algorithm
        if b == n - 2 {
            let new_a = self.nodes.remove(&(n - 1)).unwrap();
            self.nodes.insert(a, new_a);
        } else {
            let new_a = self.nodes.remove(&(n - 2));
            if let Some(new_a) = new_a {
                self.nodes.insert(a, new_a);
            }
            let new_b = self.nodes.remove(&(n - 1));
            if let Some(new_b) = new_b {
                self.nodes.insert(b, new_b);
            }
        }
        let u = self.tree.add_node("".to_string());
        self.nodes.insert(self.n_unmerged_leaves - 2, u);
        self.n_unmerged_leaves -= 1;

        // Add new edges
        self.tree.add_edge(u, a_node, dau);
        self.tree.add_edge(u, b_node, dbu);
        u
    }
}

// Test
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_new() {
        let tree = PhyloTree::build(&vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
        ]);
        assert_eq!(tree.tree.node_count(), 4);
        assert_eq!(tree.tree.edge_count(), 0);
    }
}
