use crate::distances::DistanceMatrix;
use crate::rapid_nj::node::Node;
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::sync::{Arc, Mutex};

pub struct QMatrix {
    distances: Vec<Option<Vec<f64>>>,
    sum_cols: Vec<Option<f64>>,
    trees: Vec<Option<BTreeSet<Node>>>,
    u_max: f64,
    n: usize,
    n_leaves: usize,
}

impl QMatrix {
    pub fn distance(&self, i: usize, j: usize) -> f64 {
        Self::distances_vec(&self.distances, i, j)
    }
    pub fn distances_vec(distances: &[Option<Vec<f64>>], i: usize, j: usize) -> f64 {
        match i.cmp(&j) {
            Ordering::Less => distances[i].as_ref().unwrap()[j - i - 1],
            Ordering::Greater => distances[j].as_ref().unwrap()[i - j - 1],
            Ordering::Equal => 0.0,
        }
    }

    pub fn find_neighbors(&self) -> (usize, usize) {
        let threads = rayon::current_num_threads();
        // Create slices of the self.trees
        let index = (0..self.n).collect::<Vec<usize>>();
        let qmin_shared = Arc::new(Mutex::new(f64::INFINITY)); // Shared q_min
        let min_index_shared = Arc::new(Mutex::new((0, 0))); // Shared min_index
        let minima: Vec<((usize, usize), f64)> = self
            .trees
            .par_chunks(threads)
            .zip(index.par_chunks(threads))
            .map(|(chunk, chunk_index)| {
                // Clone the current value of qmin and min_index
                let mut qmin = f64::INFINITY;
                let mut min_index = (0, 0);
                for (i, tree) in chunk_index.iter().zip(chunk.iter()) {
                    qmin = *qmin_shared.lock().unwrap();
                    min_index = *min_index_shared.lock().unwrap();
                    let i = *i;
                    if tree.is_none() {
                        continue;
                    }
                    let tree = tree.as_ref().unwrap();
                    for node in tree.iter() {
                        let j = node.index;
                        if (self.n_leaves as f64 - 2.0) * self.distance(i, j)
                            - self.sum_cols[i].unwrap()
                            - self.u_max
                            >= qmin
                        {
                            break;
                        }
                        let q = (self.n_leaves as f64 - 2.0) * self.distance(i, j)
                            - self.sum_cols[i].unwrap()
                            - self.sum_cols[j].unwrap();
                        if q < qmin {
                            qmin = q;
                            min_index = (i, j);
                            {
                                let mut qmin_shared = qmin_shared.lock().unwrap();
                                let mut min_index_shared = min_index_shared.lock().unwrap();
                                if qmin < *qmin_shared {
                                    *qmin_shared = qmin;
                                    *min_index_shared = min_index;
                                }
                            }
                        }
                    }
                }
                (min_index, qmin)
            })
            .collect();
        // Choose the minimum of the minima
        let mut qmin = f64::INFINITY;
        let mut min_index = (0, 0);
        for (index, q) in minima {
            if q < qmin {
                qmin = q;
                min_index = index;
            }
        }
        min_index
    }
    pub fn update(&mut self, i: usize, j: usize) {
        self.trees[i] = None;
        self.trees[j] = None;
        self.sum_cols[i] = None;
        self.sum_cols[j] = None;
        let distances = &mut self.distances;
        self.sum_cols.push(Some(0.0));
        for (m, row) in self.trees.iter_mut().enumerate() {
            if row.is_none() {
                continue;
            }
            let row = row.as_mut().unwrap();

            let dim = Self::distances_vec(distances, i, m);
            row.remove(&Node::new(i, dim));
            let djm = Self::distances_vec(distances, j, m);
            row.remove(&Node::new(j, djm));

            let new_distance = 0.5
                * (Self::distances_vec(distances, i, m) + Self::distances_vec(distances, j, m)
                    - Self::distances_vec(distances, i, j));
            row.insert(Node::new(self.n, new_distance));

            self.sum_cols[m] = Some(self.sum_cols[m].unwrap() - dim - djm + new_distance);
            self.sum_cols[self.n] = Some(self.sum_cols[self.n].unwrap() + new_distance);
            distances[m].as_mut().unwrap().push(new_distance);
        }
        self.n_leaves -= 1;
        self.n += 1;
        self.distances.push(Some(Vec::with_capacity(self.n_leaves))); // Maybe add with capacity
        self.distances[i] = None;
        self.distances[j] = None;
        self.trees.push(Some(BTreeSet::new_in(std::alloc::Global)));
    }

    pub fn n_leaves(&self) -> usize {
        self.n_leaves
    }
    pub fn new_node_distances(&self, i: usize, j: usize) -> (f64, f64) {
        let s = (self.n_leaves() - 2) as f64;
        let dist_ui =
            self.distance(i, j) + self.sum_cols[i].unwrap() / s - self.sum_cols[j].unwrap() / s;
        (dist_ui / 2.0, self.distance(i, j) - dist_ui / 2.0)
    }
    pub fn unmerged_nodes(&self) -> Vec<usize> {
        // Get index of all Some valyes in self.trees
        let mut unmerged = Vec::new();
        for (i, vct) in self.distances.iter().enumerate() {
            if vct.is_some() {
                unmerged.push(i);
            }
        }
        unmerged
    }
}

// Implement from DistanceMatrix

impl From<&DistanceMatrix> for QMatrix {
    fn from(d: &DistanceMatrix) -> Self {
        let n = d.size();
        let n_leaves = n;
        let matrix = &d.matrix;
        let sum_cols: Vec<Option<f64>> = matrix
            .par_iter()
            .map(|row| Some(row.iter().sum::<f64>()))
            .collect();
        let u_max = sum_cols
            .par_iter()
            .max_by(|a, b| a.unwrap().partial_cmp(&b.unwrap()).unwrap())
            .unwrap()
            .unwrap();
        let mut distances = Vec::with_capacity(n);
        for (i, whole_row) in matrix.iter().enumerate() {
            let mut row = Vec::with_capacity(n - i - 1);
            for distance in whole_row.iter().skip(i + 1) {
                row.push(*distance);
            }
            distances.push(Some(row));
        }
        let mut trees = Vec::with_capacity(n);
        for (row_index, row) in distances.iter().enumerate() {
            let mut tree = BTreeSet::new();
            for (col_index, value) in row.as_ref().unwrap().iter().enumerate() {
                tree.insert(Node::new(col_index + row_index + 1, *value));
            }
            trees.push(Some(tree));
        }
        // Sort according to sum_cols[i]
        QMatrix {
            distances,
            sum_cols,
            trees,
            u_max,
            n,
            n_leaves,
        }
    }
}

// Test QMatrix::from
#[cfg(test)]
mod tests {
    use super::QMatrix;
    use crate::{distances::DistanceMatrix, rapid_nj::qmatrix::Node};
    #[test]
    fn test_from_distance_matrix() {
        let d = wikipedia_distance_matrix();

        let q = QMatrix::from(&d);
        // Check column sums
        assert_eq!(
            q.sum_cols,
            vec![Some(31.0), Some(34.0), Some(34.0), Some(30.0), Some(27.0)]
        );
        // Check only the upper triangle is stored
        assert_eq!(
            &q.distances,
            &vec![
                Some(vec![5.0, 9.0, 9.0, 8.0]),
                Some(vec![10.0, 10.0, 9.0]),
                Some(vec![8.0, 7.0]),
                Some(vec![3.0]),
                Some(vec![])
            ]
        );
        for i in 0..5 {
            for j in 0..5 {
                assert_eq!(q.distance(i, j), d.matrix[i][j]);
            }
        }
        // Check tree one should be Node(1, 5.0), Node(2, 9.0), Node(3, 9.0)
        let expected_one = vec![
            Node::new(1, 5.0),
            Node::new(4, 8.0),
            Node::new(2, 9.0),
            Node::new(3, 9.0),
        ];
        for (node, expected) in q.trees[0].as_ref().unwrap().iter().zip(expected_one.iter()) {
            assert_eq!(node, expected);
        }
    }
    #[test]
    fn test_find_neighbors() {
        let d = wikipedia_distance_matrix();
        let q = QMatrix::from(&d);
        let neighbors = q.find_neighbors();
        assert_eq!(neighbors, (0, 1));
    }
    #[test]
    fn test_update() {
        let d = wikipedia_distance_matrix();
        let mut q = QMatrix::from(&d);
        q.update(0, 1);
        assert_eq!(
            &q.distances,
            &vec![
                None,
                None,
                Some(vec![8.0, 7.0, 7.0]),
                Some(vec![3.0, 7.0]),
                Some(vec![6.0]),
                Some(vec![])
            ]
        );
        assert_eq!(
            q.sum_cols,
            vec![None, None, Some(22.0), Some(18.0), Some(16.0), Some(20.0)]
        );
        assert_eq!(q.find_neighbors(), (3, 4));
        q.update(3, 4);
        assert_eq!(
            &q.distances,
            &vec![
                None,
                None,
                Some(vec![8.0, 7.0, 7.0, 6.0]),
                None,
                None,
                Some(vec![5.0]),
                Some(vec![])
            ]
        );
    }

    fn wikipedia_distance_matrix() -> DistanceMatrix {
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
        d
    }
}
