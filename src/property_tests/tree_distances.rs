use std::collections::{HashMap, HashSet};

use bit_set::BitSet;
use bit_vec::BitVec;
use petgraph::stable_graph::EdgeIndex;

use crate::Tree;
/// Calculate the Branch-Score distance between two trees
pub fn branch_score(a: Tree, b: Tree, n_leaves: usize) -> f64 {
    let mut bits_a = HashMap::new();
    let mut bits_b = HashMap::new();
    a.edge_indices()
        .zip(a.edge_weights())
        .for_each(|(edge, w)| {
            bits_a.insert(collect_bit_vector(&a, edge, n_leaves), *w);
        });
    b.edge_indices()
        .zip(b.edge_weights())
        .for_each(|(edge, w)| {
            bits_b.insert(collect_bit_vector(&b, edge, n_leaves), *w);
        });

    // Get union of a and b keys
    let mut keys: HashSet<&BitVec> = HashSet::new();
    bits_a.keys().zip(bits_b.keys()).for_each(|(x, y)| {
        keys.insert(x);
        keys.insert(y);
    });

    let mut distance = 0.0;
    for key in keys {
        let a = bits_a.get(key).unwrap_or(&0.0);
        let b = bits_b.get(key).unwrap_or(&0.0);
        distance += (a - b).powi(2);
    }
    distance
}

/// Calculate the Robinson-Foulds distance between two trees
pub fn robinson_foulds(a: Tree, b: Tree, n_leaves: usize) -> usize {
    let bits_a: HashSet<BitVec> = HashSet::from_iter(
        a.edge_indices()
            .map(|edge| collect_bit_vector(&a, edge, n_leaves)),
    );
    let bits_b: HashSet<BitVec> = HashSet::from_iter(
        b.edge_indices()
            .map(|edge| collect_bit_vector(&b, edge, n_leaves)),
    );

    let mut distance = 0;
    for a in &bits_a {
        if !bits_b.contains(a) {
            distance += 1;
        }
    }
    for b in bits_b {
        if !bits_a.contains(&b) {
            distance += 1;
        }
    }
    distance
}

fn collect_bit_vector(
    tree: &petgraph::Graph<String, f64, petgraph::Undirected>,
    edge: EdgeIndex,
    n_leaves: usize,
) -> BitVec {
    // Collect all leaves
    let mut leaves = HashMap::new();
    tree.node_indices().enumerate().for_each(|(i, index)| {
        if !tree[index].is_empty() {
            leaves.insert(index, i);
        }
    });
    let mut left_nodes = HashSet::new();
    let (parent_left, parent_right) = tree.edge_endpoints(edge).expect("Valid edge");
    let mut left_queue = Vec::new();
    left_nodes.insert(parent_left);
    tree.neighbors(parent_left)
        .filter(|node| node != &parent_right)
        .for_each(|node| {
            left_nodes.insert(node);
            left_queue.push(node);
        });
    while let Some(node) = left_queue.pop() {
        for neighbor in tree.neighbors(node) {
            if !left_nodes.contains(&neighbor) {
                left_nodes.insert(neighbor);
                left_queue.push(neighbor);
            }
        }
    }
    let mut bit_vect = BitSet::with_capacity(n_leaves);
    for node in left_nodes {
        // If node is a leaf, then add it to the bit vector
        if let Some(i) = leaves.get(&node) {
            bit_vect.insert(*i);
        }
    }
    let mut bit_vect = bit_vect.into_bit_vec();
    if bit_vect[0] {
        bit_vect.negate();
    }
    bit_vect
}

#[cfg(test)]
mod tests {
    use rand::distributions::{Alphanumeric, DistString};

    use super::*;

    #[test]
    fn test_distance() {
        let names: Vec<String> = (0..4)
            .map(|_| Alphanumeric.sample_string(&mut rand::thread_rng(), 16))
            .collect();
        let mut t1 = Tree::new_undirected();
        {
            // Simulate random leaf node names

            let a = t1.add_node(names[0].to_owned());
            let b = t1.add_node(names[1].to_owned());
            let c = t1.add_node(names[2].to_owned());
            let d = t1.add_node(names[3].to_owned());
            let u = t1.add_node("".to_owned());
            let v = t1.add_node("".to_owned());
            t1.add_edge(a, v, 2.0);
            t1.add_edge(b, v, 1.0);
            t1.add_edge(c, u, 3.0);
            t1.add_edge(d, u, 10.0);
            t1.add_edge(v, u, 5.0);
        }
        let mut t2 = Tree::new_undirected();
        {
            let a = t2.add_node(names[0].to_owned());
            let b = t2.add_node(names[1].to_owned());
            let c = t2.add_node(names[2].to_owned());
            let d = t2.add_node(names[3].to_owned());
            let u = t2.add_node("".to_owned());
            let v = t2.add_node("".to_owned());
            t2.add_edge(a, v, 1.0);
            t2.add_edge(c, v, 1.0);
            t2.add_edge(b, u, 1.0);
            t2.add_edge(d, u, 1.0);
            t2.add_edge(v, u, 1.0);
        }
        assert_eq!(robinson_foulds(t1.clone(), t2.clone(), 4), 2);
        assert_eq!(robinson_foulds(t1.clone(), t1.clone(), 4), 0);
        assert_eq!(robinson_foulds(t2.clone(), t2.clone(), 4), 0);
        //
        assert_eq!(branch_score(t1.clone(), t2.clone(), 4), 112.0);
        assert_eq!(branch_score(t1.clone(), t1, 4), 0.0);
        assert_eq!(branch_score(t2.clone(), t2, 4), 0.0);
    }
}
