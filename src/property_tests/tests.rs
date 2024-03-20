#[cfg(test)]
fn assert_equal_tree(a: &crate::Tree, b: &crate::Tree, i: usize) {
    use crate::property_tests::tree_distances::{branch_score, robinson_foulds};
    assert_eq!(robinson_foulds(a.clone(), b.clone(), i), 0);
    assert!(petgraph::algo::is_isomorphic(&a.clone(), &b));
    assert!(branch_score(a.clone(), b.clone(), i) < f64::EPSILON);
}
#[test]
fn test_random_additive_binary_trees_naive() {
    use crate::{
        naive_nj::canonical_neighbor_joining,
        property_tests::random_additive_tree::{
            distance_matrix_from_tree, random_unrooted_binary_tree,
        },
    };
    for i in 4..20 {
        let original_tree = random_unrooted_binary_tree(i);
        let d = distance_matrix_from_tree(original_tree.clone());
        let tree = canonical_neighbor_joining(d).unwrap();
        assert_equal_tree(&original_tree, &tree, i)
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
        let d = distance_matrix_from_tree(original_tree.clone());
        let chunk_size = rand::random::<usize>() % (i + 1) + 1;
        let tree = rapid_nj(d, chunk_size).unwrap();
        assert_equal_tree(&original_tree, &tree, i)
    }
}

#[test]
fn test_random_additive_binary_trees_mix() {
    use crate::hybrid_nj::neighbor_joining;
    use crate::property_tests::random_additive_tree::{
        distance_matrix_from_tree, random_unrooted_binary_tree,
    };

    let original_tree = random_unrooted_binary_tree(20);
    let mut d: crate::distances::DistanceMatrix = distance_matrix_from_tree(original_tree.clone());
    d.permutate();
    for i in (25..100).step_by(25) {
        let original_tree = random_unrooted_binary_tree(i);
        let d: crate::distances::DistanceMatrix = distance_matrix_from_tree(original_tree.clone());
        for _ in 0..5 {
            let naive_steps = rand::random::<usize>() % (i + 1);
            let chunk_size = rand::random::<usize>() % (i + 1) + 1;
            let tree = neighbor_joining(d.clone(), naive_steps, chunk_size).unwrap();
            assert_equal_tree(&original_tree, &tree, i)
        }
    }
}
