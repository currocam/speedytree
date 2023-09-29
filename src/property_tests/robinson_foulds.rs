use std::collections::HashMap;

use bit_set::BitSet;

use crate::Tree;

pub fn distance(a: Tree, b: Tree) -> usize {
    // Create a map from leaves to number according to their lexicographic order
    let mut leaves: Vec<String> = a
        .node_indices()
        .filter(|x| !a[*x].is_empty())
        .map(|x: petgraph::stable_graph::NodeIndex| a.node_weight(x).unwrap().clone())
        .collect();
    leaves.sort();
    // Create HashMap with leaf name as key and number as value
    let mut leaf_map: HashMap<String, usize> = HashMap::new();
    for (i, leaf) in leaves.iter().enumerate() {
        leaf_map.insert(leaf.clone(), i);
    }
    // Create closure that acts as a inner function to collec the bit vectors
    let collect_bit_vectors = |x: Tree| {
        // Let bit_vectors to be a set
        let mut bit_vct: Vec<BitSet> = Vec::new();
        for node in x.node_indices() {
            if x[node].is_empty() {
                let mut bit_vector = BitSet::with_capacity(leaves.len());
                // Set 1 to the neighbors leaves
                for neighbor in x.neighbors(node) {
                    let leaf = x.node_weight(neighbor).unwrap();
                    let leaf_number = leaf_map.get(leaf);
                    if let Some(leaf_number) = leaf_number {
                        bit_vector.insert(*leaf_number);
                    }
                }
                bit_vct.push(bit_vector);
            }
        }
        // bit_vct without empty sets
        bit_vct.retain(|x| !x.is_empty());
        bit_vct
    };
    // Create vector with collected bit vectors a and b, by concatenating them
    let vct_a = collect_bit_vectors(a);
    let vct_b = collect_bit_vectors(b);
    let mut bit_vct: Vec<BitSet> = vct_a.clone();
    bit_vct.append(&mut vct_b.clone());
    bit_vct.sort();
    // Count number of shared splits as number of doublets
    let mut shared_splits = 0;
    for i in 0..bit_vct.len() - 1 {
        if bit_vct[i] == bit_vct[i + 1] {
            shared_splits += 1;
        }
    }
    (vct_a.len() - shared_splits) + (vct_b.len() - shared_splits)
}

#[cfg(test)]
mod tests {
    use rand::distributions::{Alphanumeric, DistString};

    use super::*;

    #[test]
    fn test_distance() {
        let names: Vec<String> = (0..4)
            .into_iter()
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
            t1.add_edge(a, v, 1.0);
            t1.add_edge(b, v, 1.0);
            t1.add_edge(c, u, 1.0);
            t1.add_edge(d, u, 1.0);
            t1.add_edge(v, u, 1.0);
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
        assert_eq!(distance(t1.clone(), t2.clone()), 4);
        assert_eq!(distance(t1.clone(), t1), 0);
        assert_eq!(distance(t2.clone(), t2), 0);
    }
}
