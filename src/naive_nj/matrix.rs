use std::cmp;

use crate::phylip_distance_matrix::DistanceMatrix;

#[derive(Debug)]
pub struct QMatrix {
    matrix: Vec<Vec<f64>>,
    sum_cols: Vec<f64>,
}

impl QMatrix {
    pub fn n_leaves(&self) -> usize {
        self.matrix.len()
    }
    pub fn distance(&self, i: usize, j: usize) -> f64 {
        self.matrix[i][j]
    }
    pub fn new_node_distances(&self, i: usize, j: usize) -> (f64, f64) {
        let s = (self.n_leaves() - 2) as f64;
        let dist_ui = self.distance(i, j) + self.sum_cols[i] / s - self.sum_cols[j] / s;
        (dist_ui / 2.0, self.distance(i, j) - dist_ui / 2.0)
    }

    pub fn new(d: DistanceMatrix) -> Self {
        let matrix = d.matrix;
        let sum_cols = matrix
            .iter()
            .map(|row| row.iter().sum::<f64>())
            .collect::<Vec<f64>>();
        Self { matrix, sum_cols }
    }
    pub fn find_neighbors(&self) -> (usize, usize) {
        let matrix = &self.matrix;
        let sums = &self.sum_cols;
        let n = matrix.len();
        let mut neighbors = (0, 0);
        let mut best_q = f64::INFINITY;

        for i in 0..n {
            for j in i + 1..n {
                let q = (matrix[i][j] * (n - 2) as f64) - sums[i] - sums[j];
                if q < best_q {
                    best_q = q;
                    neighbors = (i, j);
                }
            }
        }
        (
            cmp::min(neighbors.0, neighbors.1),
            cmp::max(neighbors.0, neighbors.1),
        )
    }

    pub fn update_distance_matrix(&mut self, i: usize, j: usize) {
        let matrix = &mut self.matrix;
        let sum_cols = &mut self.sum_cols;
        let dij = matrix[i][j];
        let n = matrix.len();
        // Remove the ith and jth value to each row
        for k in 0..matrix.len() {
            sum_cols[k] -= matrix[i][k] + matrix[j][k];
        }
        // Swap rows
        if j == n - 2 {
            matrix.swap(i, n - 1);
            sum_cols.swap(i, n - 1);
            for row in matrix.iter_mut() {
                row.swap(i, n - 1);
            }
        } else {
            matrix.swap(i, n - 2);
            matrix.swap(j, n - 1);
            sum_cols.swap(i, n - 2);
            sum_cols.swap(j, n - 1);
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
        sum_cols.pop();
        for row in matrix.iter_mut() {
            row.pop();
        }
        // Update the sum_cols with RS_i = RS'_i - x - y + z
        for index in 0..n - 2 {
            sum_cols[index] += matrix[n - 2][index];
        }
        // Compute the sum of the last row
        sum_cols[n - 2] = matrix[n - 2].iter().sum::<f64>();
    }
}
