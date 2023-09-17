// http://qspace.qu.edu.qa/bitstream/handle/10576/9747/A%20random%20binary%20trees%20generation%20method.pdf?sequence=10

use petgraph::graph::NodeIndex;
use petgraph::prelude::UnGraph;
use rand::Rng;

use crate::distances::DistanceMatrix;
use crate::Tree;

fn random_rooted_binary_tree(leaves: usize) -> UnGraph<String, f64> {
    let mut tree = UnGraph::new_undirected();
    let mut next_node = leaves;
    let mut nodes: Vec<NodeIndex> = (0..leaves).map(|_| tree.add_node("".to_string())).collect();
    while next_node <= 2 * leaves - 2 {
        let mut rng = rand::thread_rng();
        // Pop random node from range
        let a = nodes.swap_remove(rng.gen_range(0..nodes.len()));
        let b = nodes.swap_remove(rng.gen_range(0..nodes.len()));
        let u = tree.add_node("".to_string());
        nodes.push(u);
        tree.add_edge(u, a, rng.gen_range(0.1..100.0));
        tree.add_edge(u, b, rng.gen_range(0.1..100.0));
        next_node += 1;
    }
    for node in tree.node_indices() {
        // If has one edge, it's a leaf
        if tree.edges(node).count() == 1 {
            tree[node] = format!("L{}", node.index());
        }
    }
    tree
}

pub fn random_unrooted_binary_tree(n_leaves: usize) -> UnGraph<String, f64> {
    let mut t = random_rooted_binary_tree(n_leaves);
    // Remove root as node with degree 2
    let root = t
        .node_indices()
        .find(|node| t.edges(*node).count() == 2)
        .unwrap();
    let mut neighbors = t.neighbors_undirected(root);
    let a = neighbors.next().unwrap();
    let b = neighbors.next().unwrap();
    t.remove_node(root);
    // Add edge between a and b
    let mut rng = rand::thread_rng();
    t.add_edge(a, b, rng.gen_range(0.1..100.0));
    t
}

pub fn distance_matrix_from_tree(t: Tree) -> DistanceMatrix {
    // Find all leaves
    let leaves: Vec<NodeIndex> = t
        .node_indices()
        .filter(|node| t.edges(*node).count() == 1)
        .collect();
    let mut mat = vec![vec![0.0; leaves.len()]; leaves.len()];
    for (i, a) in leaves.iter().enumerate() {
        for (j, b) in leaves.iter().enumerate() {
            mat[i][j] =
                petgraph::algo::astar(&t, *a, |finish| finish == *b, |e| *e.weight(), |_| 0.0)
                    .unwrap()
                    .0;
        }
    }
    DistanceMatrix {
        names: leaves.iter().map(|node| t[*node].clone()).collect(),
        matrix: mat,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_random_rooted_binary_tree() {
        for _ in 0..50 {
            let mut rng = rand::thread_rng();
            let n_leaves = rng.gen_range(5..30);
            let tree = random_rooted_binary_tree(n_leaves);
            // Count number of non "" nodes
            let mut n = 0;
            for node in tree.node_indices() {
                if tree[node] != "" {
                    n += 1;
                }
            }
            assert_eq!(n, n_leaves);
            assert_eq!(tree.node_count(), 2 * n_leaves - 1);
        }
    }
    #[test]
    fn test_random_unrooted_binary_tree() {
        for _ in 0..50 {
            let mut rng = rand::thread_rng();
            let n_leaves = rng.gen_range(5..30);
            let tree = random_unrooted_binary_tree(n_leaves);
            // Count number of non "" nodes
            let mut n = 0;
            for node in tree.node_indices() {
                if tree[node] != "" {
                    n += 1;
                }
            }
            assert_eq!(n, n_leaves);
            assert_eq!(tree.node_count(), 2 * n_leaves - 2);
            // Assert any leave has degree 1 and any other node has degree 3
            for node in tree.node_indices() {
                if tree[node] != "" {
                    assert_eq!(tree.edges(node).count(), 1);
                } else {
                    assert_eq!(tree.edges(node).count(), 3);
                }
            }
        }
    }
}