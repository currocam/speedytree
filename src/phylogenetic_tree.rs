use std::collections::HashMap;

// Binary tree with edge lengths using petgraph
// graph from petagraph
use petgraph::{graph::UnGraph, stable_graph::NodeIndex};
use rand::Rng;

#[derive(Debug, Clone)]
pub struct PhyloTree {
    pub tree: UnGraph<String, f64>,
    n_unmerged_leaves: usize,
    pub leaves: HashMap<usize, NodeIndex>,
    pub nodes: HashMap<usize, NodeIndex>,
}

impl PhyloTree {
    pub fn new(leafs: &Vec<String>) -> PhyloTree {
        let mut tree: petgraph::Graph<String, f64, petgraph::Undirected> =
            UnGraph::new_undirected();
        let mut nodes = HashMap::new();
        for (i, leaf) in leafs.iter().enumerate() {
            let node: NodeIndex = tree.add_node(leaf.to_owned());
            //tree.add_edge(root, node, 0.0);
            nodes.insert(i, node);
        }
        let n_leaves = leafs.len();
        let n_unmerged_leaves = n_leaves;
        let leaves = nodes.clone();
        PhyloTree {
            tree,
            leaves,
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
    pub fn merge_rapid_nj(&mut self, a: usize, b: usize, dau: f64, dbu: f64) -> NodeIndex {
        // Get nodes to merge
        let n = self.nodes.len();
        let u = self.tree.add_node("".to_string());
        self.nodes.insert(n, u);
        let a_node = self.nodes.get(&a).unwrap();
        let b_node = self.nodes.get(&b).unwrap();
        // Add new edges
        self.tree.add_edge(u, *a_node, dau);
        self.tree.add_edge(u, *b_node, dbu);
        u
    }

    pub fn random(n: usize) -> PhyloTree {
        let leafs = (1..n + 1)
            .map(|x| x.to_string())
            .rev()
            .collect::<Vec<String>>();
        // Shuffle leafs
        let mut tree = PhyloTree::new(&leafs);
        let n_nodes_end = 2 * tree.n_unmerged_leaves - 2;
        let mut previous_internal_node = None;
        let mut internal_node = None;
        for _ in 0..n_nodes_end - n {
            previous_internal_node = internal_node;
            internal_node = Some(tree.merge_neighbors(0, 1, 0.0, 0.0));
        }
        tree.tree
            .add_edge(previous_internal_node.unwrap(), internal_node.unwrap(), 0.0);
        // Iterate through all edges and assign random weights
        let mut rng = rand::thread_rng();
        for edge in tree.tree.edge_indices() {
            let weight = rng.gen_range(0.1..10.0);
            tree.tree[edge] = weight;
        }
        tree
    }
}

// Test
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_new() {
        let tree = PhyloTree::new(&vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
        ]);
        assert_eq!(tree.tree.node_count(), 4);
        assert_eq!(tree.tree.edge_count(), 0);
    }
}
