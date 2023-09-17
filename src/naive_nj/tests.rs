#[cfg(test)]
mod tests {
    use crate::naive_nj::algorithm::naive_neighbor_joining;
    use crate::random_binary_trees::distance_matrix_from_tree;

    #[test]
    fn test_random_additive_binary_trees() {
        for i in 4..50 {
            let original_tree = crate::random_binary_trees::random_unrooted_binary_tree(i);
            let d = distance_matrix_from_tree(original_tree.clone());
            let tree = naive_neighbor_joining(d).unwrap();
            assert!(petgraph::algo::is_isomorphic(&original_tree, &tree));
        }
    }
}
