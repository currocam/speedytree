// Binary tree with edge lengths using petgraph
// graph from petagraph
use petgraph::{
    adj::EdgeReference,
    graph::Graph,
    stable_graph::NodeIndex,
    visit::{EdgeRef, IntoEdgesDirected},
    Direction::Incoming,
};
pub struct PhyloTree {
    pub tree: Graph<String, f64>,
    root: NodeIndex,
    pub nodes: Vec<NodeIndex>,
}

impl PhyloTree {
    pub fn new(leafs: &Vec<String>) -> PhyloTree {
        let mut tree = Graph::new();
        let root = tree.add_node("root".to_string());
        let mut nodes = Vec::new();
        for leaf in leafs {
            let node = tree.add_node(leaf.to_owned());
            tree.add_edge(root, node, 0.0);
            nodes.push(node);
        }
        PhyloTree { tree, root, nodes }
    }
    pub fn merge_nodes(&mut self, a: NodeIndex, b: NodeIndex) -> NodeIndex {
        // Create new node
        let new_node = self.tree.add_node("".to_string());
        // Get edge from a & b to root
        let a_edge = self.tree.find_edge_undirected(a, self.root).unwrap();
        let b_edge = self.tree.find_edge_undirected(b, self.root).unwrap();
        // Get edge weights
        let a_weight = *self.tree.edge_weight(a_edge.0).unwrap();
        let b_weight = *self.tree.edge_weight(b_edge.0).unwrap();
        // Add new edges
        self.tree
            .add_edge(self.root, new_node, (a_weight + b_weight) / 2.0);
        self.tree.add_edge(new_node, a, a_weight / 2.0);
        self.tree.add_edge(new_node, b, b_weight / 2.0);
        // Remove old edges
        self.tree.remove_edge(a_edge.0);
        self.tree.remove_edge(b_edge.0);
        new_node
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
    // Test  merge_nodes
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
        //assert_eq!(tree.tree.edge_count(), 5);
        assert_eq!(tree.tree.neighbors(tree.root).count(), 2);
    }
}
