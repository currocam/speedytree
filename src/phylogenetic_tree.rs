use std::collections::HashMap;

// Binary tree with edge lengths using petgraph
// graph from petagraph
use petgraph::{graph::UnGraph, stable_graph::NodeIndex, visit::IntoNodeReferences};
use rand::Rng;

pub struct PhyloTree {
    pub tree: UnGraph<String, f64>,
    root: NodeIndex,
    n_leaves: usize,
    pub nodes: HashMap<usize, NodeIndex>,
}

impl PhyloTree {
    pub fn new(leafs: &Vec<String>) -> PhyloTree {
        let mut tree: petgraph::Graph<String, f64, petgraph::Undirected> =
            UnGraph::new_undirected();
        let root = tree.add_node("root".to_string());
        let mut nodes = HashMap::new();
        for (i, leaf) in leafs.iter().enumerate() {
            let node = tree.add_node(leaf.to_owned());
            tree.add_edge(root, node, 0.0);
            nodes.insert(i, node);
        }
        let n_leaves = leafs.len();
        PhyloTree {
            tree,
            root,
            nodes,
            n_leaves,
        }
    }
    pub fn merge_neighbors(&mut self, a: usize, b: usize) -> NodeIndex {
        // Get nodes to merge
        let n = &self.n_leaves;
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
        let u = self.merge_nodes(a_node, b_node);
        self.nodes.insert(self.n_leaves - 2, u);
        self.n_leaves -= 1;
        u
    }
    fn merge_nodes(&mut self, a: NodeIndex, b: NodeIndex) -> NodeIndex {
        // Create new node
        let new_node = self.tree.add_node("".to_string());
        let a_edge = self.tree.find_edge_undirected(a, self.root);
        let b_edge = self.tree.find_edge_undirected(b, self.root);

        let a_edge = a_edge.unwrap();
        let b_edge = b_edge.unwrap();
        // Get edge weights
        let a_weight = *self.tree.edge_weight(a_edge.0).unwrap();
        let b_weight = *self.tree.edge_weight(b_edge.0).unwrap();
        // Add new edges
        self.tree
            .add_edge(new_node, self.root, (a_weight + b_weight) / 2.0);
        self.tree.add_edge(new_node, a, a_weight / 2.0);
        self.tree.add_edge(new_node, b, b_weight / 2.0);
        // Remove old edges
        self.tree.remove_edge(a_edge.0);
        self.tree.remove_edge(b_edge.0);
        new_node
    }
    pub fn random(n: usize) -> PhyloTree {
        let leafs = (1..n).map(|x| x.to_string()).collect::<Vec<String>>();
        let mut tree = PhyloTree::new(&leafs);
        // Find a first node with more than 2 neighbors
        let mut node = tree.tree.neighbors(tree.root).next().unwrap();
        while tree.tree.neighbors(node).count() <= 2 {
            node = tree.tree.neighbors(node).next().unwrap();
        }
        // Merge nodes until only 2 neighbors are left
        while tree.tree.neighbors(node).count() > 2 {
            let a = tree.tree.neighbors(node).next().unwrap();
            let b = tree.tree.neighbors(node).nth(1).unwrap();
            tree.merge_nodes(a, b);
        }

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
        assert_eq!(tree.tree.node_count(), 5);
        assert_eq!(tree.tree.edge_count(), 4);
        assert_eq!(tree.tree.neighbors(tree.root).count(), 4);
    }
    #[test]
    fn test_merge_nodes() {
        let mut tree = PhyloTree::new(&vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
        ]);
        let a = tree.tree.neighbors(tree.root).next().unwrap();
        let b = tree.tree.neighbors(tree.root).nth(1).unwrap();
        tree.merge_nodes(a, b);
        assert_eq!(tree.tree.node_count(), 6);
        assert_eq!(tree.tree.edge_count(), 5);
        assert_eq!(tree.tree.neighbors(tree.root).count(), 3);
    }
}
