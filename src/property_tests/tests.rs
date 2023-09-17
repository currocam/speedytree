#[cfg(test)]
#[test]
fn test_random_additive_binary_trees() {
    use crate::{
        naive_nj::naive_neighbor_joining,
        property_tests::random_additive_tree::{
            distance_matrix_from_tree, random_unrooted_binary_tree,
        },
    };
    for i in 4..20 {
        let original_tree = random_unrooted_binary_tree(i);
        let d = distance_matrix_from_tree(original_tree.clone());
        let tree = naive_neighbor_joining(d).unwrap();
        assert!(petgraph::algo::is_isomorphic(&original_tree, &tree));
    }
}
