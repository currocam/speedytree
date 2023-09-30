#[cfg(test)]
#[test]
fn test_random_additive_binary_trees_naive() {
    use crate::{
        naive_nj::naive_neighbor_joining,
        property_tests::random_additive_tree::{
            distance_matrix_from_tree, random_unrooted_binary_tree,
        },
    };
    for i in 4..20 {
        let original_tree = random_unrooted_binary_tree(i);
        let mut d = distance_matrix_from_tree(original_tree.clone());
        d.permutate();
        let tree = naive_neighbor_joining(d).unwrap();
        assert!(petgraph::algo::is_isomorphic(&original_tree, &tree));
    }
}
#[test]
fn test_random_additive_binary_trees_rapid() {
    use crate::property_tests::random_additive_tree::{
        distance_matrix_from_tree, random_unrooted_binary_tree,
    };
    use crate::rapid_nj::rapid_nj;
    for i in 4..20 {
        let original_tree = random_unrooted_binary_tree(i);
        let mut d = distance_matrix_from_tree(original_tree.clone());
        d.permutate();
        let tree = rapid_nj(d).unwrap();
        assert!(petgraph::algo::is_isomorphic(&original_tree, &tree));
    }
}

#[test]
fn test_random_additive_binary_trees_mix() {
    use crate::property_tests::random_additive_tree::{
        distance_matrix_from_tree, random_unrooted_binary_tree,
    };
    use crate::algo::neighbor_joining;
    let original_tree = random_unrooted_binary_tree(20);
    let mut d: crate::distances::DistanceMatrix = distance_matrix_from_tree(original_tree.clone());
    d.permutate();
    for i in 20..30 {
        let original_tree = random_unrooted_binary_tree(i);
        let d: crate::distances::DistanceMatrix = distance_matrix_from_tree(original_tree.clone());
        for _ in 0..5{
            let random = rand::random::<usize>() % (i + 1);
            dbg!(random, i);
            let tree = neighbor_joining(d.clone(), random).unwrap();
            assert!(petgraph::algo::is_isomorphic(&original_tree, &tree));
                
        }
    }
}

