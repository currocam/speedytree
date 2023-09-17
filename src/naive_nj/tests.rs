#[cfg(test)]
mod tests {
    use crate::random_binary_trees::distance_matrix_from_tree;
    use crate::{
        naive_nj::algorithm::naive_neighbor_joining, phylip_distance_matrix::DistanceMatrix,
    };

    #[test]
    fn test_example_wikipedia() {
        let mat: Vec<Vec<f64>> = vec![
            vec![0.0, 5.0, 9.0, 9.0, 8.0],
            vec![5.0, 0.0, 10.0, 10.0, 9.0],
            vec![9.0, 10.0, 0.0, 8.0, 7.0],
            vec![9.0, 10.0, 8.0, 0.0, 3.0],
            vec![8.0, 9.0, 7.0, 3.0, 0.0],
        ];
        let names = vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
            "E".to_string(),
        ];
        let phylo = naive_neighbor_joining(DistanceMatrix {
            matrix: mat,
            names: names,
        });
        assert!(phylo.is_ok());
        let phylo = phylo.unwrap();
        let a = phylo.leaves[&0];
        let b = phylo.leaves[&1];
        let c = phylo.leaves[&2];
        let d = phylo.leaves[&3];
        let e = phylo.leaves[&4];
        let u: petgraph::stable_graph::NodeIndex = phylo.tree.neighbors(a).next().unwrap();
        // Check b and u
        assert_eq!(
            phylo.tree.edge_weight(phylo.tree.find_edge(a, u).unwrap()),
            Some(&2.0)
        );
        assert_eq!(
            phylo.tree.edge_weight(phylo.tree.find_edge(b, u).unwrap()),
            Some(&3.0)
        );
        let w: petgraph::stable_graph::NodeIndex = phylo.tree.neighbors(e).next().unwrap();
        // Check a and b connected
        assert_eq!(
            phylo.tree.edge_weight(phylo.tree.find_edge(e, w).unwrap()),
            Some(&1.0)
        );
        assert_eq!(
            phylo.tree.edge_weight(phylo.tree.find_edge(d, w).unwrap()),
            Some(&2.0)
        );
        let v: petgraph::stable_graph::NodeIndex = phylo.tree.neighbors(c).next().unwrap();
        // Check a and b connected
        assert_eq!(
            phylo.tree.edge_weight(phylo.tree.find_edge(c, v).unwrap()),
            Some(&4.0)
        );
        assert_eq!(
            phylo.tree.edge_weight(phylo.tree.find_edge(v, u).unwrap()),
            Some(&3.0)
        );
        assert_eq!(
            phylo.tree.edge_weight(phylo.tree.find_edge(v, w).unwrap()),
            Some(&2.0)
        )
    }
    #[test]
    fn test_random_additive_binary_trees() {
        for i in 4..50 {
            let original_tree = crate::random_binary_trees::random_unrooted_binary_tree(i);
            let d = distance_matrix_from_tree(original_tree.clone());
            let tree = naive_neighbor_joining(d).unwrap().tree;
            assert!(petgraph::algo::is_isomorphic(&original_tree, &tree));
        }
    }
}
