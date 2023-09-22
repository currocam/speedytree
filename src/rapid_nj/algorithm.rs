use std::collections::BTreeSet;

use crate::distances::DistanceMatrix;
#[derive(Debug)]
struct Node {
    index: usize,
    value: f64,
}

impl Node {
    fn new(index: usize, value: f64) -> Self {
        Self { index, value }
    }
}
// Define comparison for Node, so only the value is compared
impl PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}
impl Eq for Node {}
impl PartialOrd for Node {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.value.partial_cmp(&other.value)
    }
}
impl Ord for Node {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.value.partial_cmp(&other.value).unwrap()
    }
}

struct QMatrix {
    distances: Vec<Vec<f64>>,
    sum_cols: Vec<f64>,
    trees: Vec<BTreeSet<Node>>,
    u_max: f64,
}

impl QMatrix {
    pub fn distance(&self, i: usize, j: usize) -> f64 {
        if i == j {
            0.0
        } else if i < j {
            self.distances[i][j - i - 1]
        } else {
            self.distances[j][i - j - 1]
        }
    }
}

// Implement from DistanceMatrix

impl From<&DistanceMatrix> for QMatrix {
    fn from(d: &DistanceMatrix) -> Self {
        let n = d.size();
        let matrix = &d.matrix;
        let sum_cols: Vec<f64> = matrix.iter().map(|row| row.iter().sum::<f64>()).collect();
        let u_max = *sum_cols
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();
        let mut distances = Vec::with_capacity(n);
        for i in 0..n {
            let mut row = Vec::with_capacity(n - i - 1);
            for j in i + 1..n {
                row.push(matrix[i][j]);
            }
            distances.push(row);
        }
        let mut trees = Vec::with_capacity(n);
        for (row_index, row) in distances.iter().enumerate() {
            let mut tree = BTreeSet::new();
            for (col_index, value) in row.iter().enumerate() {
                tree.insert(Node {
                    index: col_index + row_index + 1,
                    value: *value,
                });
            }
            trees.push(tree);
        }
        QMatrix {
            distances,
            sum_cols,
            trees,
            u_max,
        }
    }
}

// Test QMatrix::from
#[cfg(test)]
mod tests {
    use super::QMatrix;
    use crate::{distances::DistanceMatrix, rapid_nj::algorithm::Node};
    #[test]
    fn test_from_distance_matrix() {
        let d = DistanceMatrix {
            matrix: vec![
                vec![0.0, 5.0, 9.0, 9.0],
                vec![5.0, 0.0, 10.0, 10.0],
                vec![9.0, 10.0, 0.0, 8.0],
                vec![9.0, 10.0, 8.0, 0.0],
            ],
            names: vec![
                "A".to_string(),
                "B".to_string(),
                "C".to_string(),
                "D".to_string(),
            ],
        };

        let q = QMatrix::from(&d);
        // Check column sums
        assert_eq!(q.sum_cols, vec![23.0, 25.0, 27.0, 27.0]);
        // Check only the upper triangle is stored
        assert_eq!(
            &q.distances,
            &vec![vec![5.0, 9.0, 9.0], vec![10.0, 10.0], vec![8.0], vec![]]
        );
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(q.distance(i, j), d.matrix[i][j]);
            }
        }
        // Check tree one should be Node(1, 5.0), Node(2, 9.0), Node(3, 9.0)
        let expected_one = vec![Node::new(1, 5.0), Node::new(2, 9.0), Node::new(3, 9.0)];
        for (node, expected) in q.trees[0].iter().zip(expected_one.iter()) {
            assert_eq!(node, expected);
        }
    }
}
