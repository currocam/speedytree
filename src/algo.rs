use std::collections::HashMap;
use crate::{distances::DistanceMatrix, ResultBox, Tree, naive_nj::DataNaiveNJ, rapid_nj::DataRapidNJ};

pub fn neighbor_joining(dist: DistanceMatrix, naive_iters: usize) -> ResultBox<Tree> {
    // If 0 < n < min_threshold, then we use naive neighbor joining
    if dist.size() < 4 || naive_iters >= dist.size(){
        return crate::naive_neighbor_joining(dist);
    }
    if naive_iters < 4  {
        return crate::rapid_nj(dist);
    }
    let mut q = crate::rapid_nj::QMatrix::from(&dist);
    let mut t = crate::rapid_nj::PhyloTree::build(&dist.names);
    while q.n_leaves() > naive_iters {
        // Find the minimum element in the distance matrix
        let (i, j) = q.find_neighbors();
        let (dist_ui, dist_uj) = q.new_node_distances(i, j);
        t.merge_neighbors(i, j, dist_ui, dist_uj);
        q.update(i, j);
    }
    let data = DataNaiveNJ::from(DataRapidNJ::new(q, t));
    let mut q = data.qmatrix;
    let mut t = data.phylo_tree;
    while q.n_leaves() > 3 {
        // Find the minimum element in the distance matrix
        let (i, j) = q.find_neighbors();
        let (dist_ui, dist_uj) = q.new_node_distances(i, j);
        t.merge_neighbors(i, j, dist_ui, dist_uj);
        q.update_distance_matrix(i, j);
    }
    Ok(crate::naive_nj::terminate_nj(t, q))
}

impl From<crate::rapid_nj::DataRapidNJ> for crate::naive_nj::DataNaiveNJ {
    fn from(data: crate::rapid_nj::DataRapidNJ) -> Self {
        let q = data.qmatrix;
        let tree = data.phylo_tree;
        let n = q.n_leaves();
        let mut nodes = HashMap::with_capacity(n);
        let mut sum_cols: Vec<f64> = Vec::with_capacity(n);
        let mut unmerged_index: Vec<usize> = Vec::with_capacity(n);
        for (index, elm) in q.sum_cols.iter().enumerate() {
            if let Some(elm) = elm {
                sum_cols.push(*elm);
                unmerged_index.push(index);
            }
        }
        debug_assert_eq!(sum_cols.len(), n);
        for (index, prev_index) in unmerged_index.iter().enumerate() {
            nodes.insert(index, tree.nodes[prev_index]);
        }
        // Create a vector of n vectors of n f64 with zeros
        let mut matrix: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
        for (i, prev_i) in unmerged_index.iter().enumerate() {
            for (j, prev_j) in unmerged_index.iter().enumerate() {
                matrix[i][j] = q.distance(*prev_i, *prev_j);
            }
        }
        let qmatrix = crate::naive_nj::QMatrix::new(matrix, sum_cols);
        let phylo_tree = crate::naive_nj::PhyloTree::new(tree.tree, nodes);
        DataNaiveNJ { qmatrix, phylo_tree }
    }
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

        let phylo = neighbor_joining(d, 4);
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
