use crate::{distances::DistanceMatrix, ResultBox, Tree};

use super::{phylo_tree::PhyloTree, qmatrix::QMatrix};

pub fn naive_neighbor_joining(dist: DistanceMatrix) -> ResultBox<Tree> {
    let mut t = PhyloTree::build(&dist.names);
    let mut q = QMatrix::build(dist);
    while q.n_leaves() > 3 {
        // Find the minimum element in the distance matrix
        let (i, j) = q.find_neighbors();
        let (dist_ui, dist_uj) = q.new_node_distances(i, j);
        t.merge_neighbors(i, j, dist_ui, dist_uj);
        q.update_distance_matrix(i, j);
    }

    Ok(terminate_nj(t, q))
}

pub fn terminate_nj(tree: PhyloTree, q: QMatrix) -> Tree {
    let (i, j, m) = (tree.nodes[&0], tree.nodes[&1], tree.nodes[&2]);
    let mut tree = tree.tree;

    let dvi = (q.distance(0, 1) + q.distance(0, 2) - q.distance(1, 2)) / 2.0;
    let dvj = (q.distance(0, 1) + q.distance(1, 2) - q.distance(0, 2)) / 2.0;
    let dvm = (q.distance(0, 2) + q.distance(1, 2) - q.distance(0, 1)) / 2.0;

    let v = tree.add_node("".to_owned());
    tree.add_edge(v, i, dvi);
    tree.add_edge(v, j, dvj);
    tree.add_edge(v, m, dvm);

    tree
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_example_wikipedia() {
        let d = DistanceMatrix {
            matrix: vec![
                vec![0.0, 5.0, 9.0, 9.0, 8.0],
                vec![5.0, 0.0, 10.0, 10.0, 9.0],
                vec![9.0, 10.0, 0.0, 8.0, 7.0],
                vec![9.0, 10.0, 8.0, 0.0, 3.0],
                vec![8.0, 9.0, 7.0, 3.0, 0.0],
            ],
            names: vec![
                "A".to_string(),
                "B".to_string(),
                "C".to_string(),
                "D".to_string(),
                "E".to_string(),
            ],
        };

        let phylo = naive_neighbor_joining(d);
        assert!(phylo.is_ok());

        let tree = phylo.unwrap();
        let mut node_indices = tree.node_indices();
        let (a, b, c, d, e) = (
            node_indices.next().unwrap(),
            node_indices.next().unwrap(),
            node_indices.next().unwrap(),
            node_indices.next().unwrap(),
            node_indices.next().unwrap(),
        );

        // Get internal nodes
        let u = tree.neighbors(a).next().unwrap();
        let w = tree.neighbors(e).next().unwrap();
        let v = tree.neighbors(c).next().unwrap();

        // Create array of node from, node to and distance
        let expected = [
            (a, u, 2.0),
            (b, u, 3.0),
            (e, w, 1.0),
            (d, w, 2.0),
            (c, v, 4.0),
            (v, u, 3.0),
            (v, w, 2.0),
        ];
        for (a, b, dist) in expected.into_iter() {
            assert_eq!(tree.edge_weight(tree.find_edge(a, b).unwrap()), Some(&dist));
        }
    }
}
