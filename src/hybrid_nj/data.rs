use std::collections::HashMap;

use crate::naive_nj::DataNaiveNJ;

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
        DataNaiveNJ {
            qmatrix,
            phylo_tree,
        }
    }
}
