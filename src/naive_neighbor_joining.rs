use std::cmp::{max, min};

use petgraph::stable_graph::NodeIndex;

use crate::{phylip_distance_matrix::DistanceMatrix, phylogenetic_tree::PhyloTree, ResultBox};

#[derive(Debug)]
struct NjMatrix {
    matrix: Vec<Vec<f64>>,
    names: Vec<String>,
    sum_cols: Vec<f64>,
}

impl NjMatrix {
    fn new(d: DistanceMatrix) -> Self {
        let matrix = d.matrix;
        let names = d.names;
        let sum_cols = matrix
            .iter()
            .map(|row| row.iter().sum::<f64>())
            .collect::<Vec<f64>>();
        Self {
            matrix,
            names,
            sum_cols,
        }
    }
}

pub fn naive_neighbor_joining(dist: DistanceMatrix) -> ResultBox<PhyloTree> {
    let mut t = PhyloTree::new(&dist.names);
    dbg!(&t.nodes);
    let mut dist = NjMatrix::new(dist);
    // While n nodes is not 2n -2
    let n_nodes_end = 2 * t.n_leaves - 2;
    while t.tree.node_count() < n_nodes_end {
        // Find the minimum element in the distance matrix
        let (i, j) = find_neighbors(&mut dist);
        let _u: NodeIndex = t.merge_neighbors(i, j);
        update_distance_matrix(i, j, &mut dist);
    }
    t = terminate_nj(t, &mut dist);
    Ok(t)
}

fn terminate_nj(mut t: PhyloTree, dist: &mut NjMatrix) -> PhyloTree {
    t.merge_neighbors(0, 1);
    t
}

fn update_distance_matrix(i: usize, j: usize, d: &mut NjMatrix) {
    let matrix = &mut d.matrix;
    // Remove the ith and jth value to each row
    for k in 0..matrix.len() {
        d.sum_cols[k] -= matrix[i][k] + matrix[j][k];
    }

    let dij = matrix[i][j];
    let n = matrix.len();

    // Swap rows
    if j == n - 2 {
        matrix.swap(i, n - 1);
        d.sum_cols.swap(i, n - 1);
        for row in matrix.iter_mut() {
            row.swap(i, n - 1);
        }
    } else {
        matrix.swap(i, n - 2);
        matrix.swap(j, n - 1);
        d.sum_cols.swap(i, n - 2);
        d.sum_cols.swap(j, n - 1);
        for row in matrix.iter_mut() {
            row.swap(i, n - 2);
            row.swap(j, n - 1);
        }
    }

    // Update the row.len() - 2 row (aka u row)
    for k in 0..matrix.len() - 2 {
        matrix[n - 2][k] = (matrix[n - 2][k] + matrix[n - 1][k] - dij) / 2.0;
        matrix[k][n - 2] = matrix[n - 2][k];
    }

    // Remove the last row and every last column
    matrix.pop();
    d.sum_cols.pop();
    for row in matrix.iter_mut() {
        row.pop();
    }

    // Update the sum_cols with RS_i = RS'_i - x - y + z
    for index in 0..n - 2 {
        d.sum_cols[index] += matrix[n - 2][index];
    }

    // Compute the sum of the last row
    d.sum_cols[n - 2] = matrix[n - 2].iter().sum::<f64>();

    debug_assert_eq!(d.sum_cols.len(), d.matrix.len())
}

fn find_neighbors(d: &mut NjMatrix) -> (usize, usize) {
    let mut neighbors = (0, 0);
    let mut best_q = f64::INFINITY;
    let matrix = &d.matrix;
    let n = matrix.len();
    let sums = &d.sum_cols;

    for i in 0..n {
        for j in i + 1..n {
            let q = (matrix[i][j] * (n - 2) as f64) - sums[i] - sums[j];
            if q < best_q {
                best_q = q;
                neighbors = (i, j);
            }
        }
    }
    (min(neighbors.0, neighbors.1), max(neighbors.0, neighbors.1))
}

#[cfg(test)]
mod tests {
    use std::vec;

    use super::*;

    #[test]
    fn test_find_neighbors() {
        let mat = vec![
            vec![0.0, 5.0, 9.0, 9.0, 8.0],
            vec![5.0, 0.0, 10.0, 10.0, 9.0],
            vec![9.0, 10.0, 0.0, 8.0, 7.0],
            vec![9.0, 10.0, 8.0, 0.0, 3.0],
            vec![8.0, 9.0, 7.0, 3.0, 0.0],
        ];
        let mut mat = NjMatrix::new(DistanceMatrix {
            matrix: mat,
            names: vec![
                "A".to_string(),
                "B".to_string(),
                "C".to_string(),
                "D".to_string(),
                "E".to_string(),
            ],
        });
        let index = find_neighbors(&mut mat);
        assert_eq!(index, (0, 1));
    }
    // NJ should find a true binary tree if it exist
    #[test]
    fn test_naive_neighbor_joining() {
        for i in 4..8 {
            let original_tree = PhyloTree::random(i);
            let d = DistanceMatrix::from(original_tree.clone());
            let original_tree = original_tree.tree;
            let tree = naive_neighbor_joining(d).unwrap();
            assert!(petgraph::algo::is_isomorphic(&original_tree, &tree.tree));
        }
    }
}
