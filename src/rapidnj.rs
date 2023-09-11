use std::collections::BTreeSet;

use crate::phylip_distance_matrix::DistanceMatrix;
#[derive(Debug)]
struct Cell {
    index: usize,
    value: f64,
}

impl Ord for Cell {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.value.partial_cmp(&other.value) {
            Some(ordering) => ordering,
            None => std::cmp::Ordering::Equal, // Handle NaN values as equal
        }
    }
}
impl PartialEq for Cell {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}
impl Eq for Cell {}
impl PartialOrd for Cell {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

struct NjBtree {
    btrees: Vec<BTreeSet<Cell>>,
    sum_cols: Vec<f64>,
}

impl NjBtree {
    fn min(&self) -> (usize, usize) {
        let mut min = f64::INFINITY;
        let mut min_index = 0;
        let mut min_j = 0;
        for (i, bt) in self.btrees.iter().enumerate() {
            if let Some(cell) = bt.first() {
                if cell.value < min {
                    min = cell.value;
                    min_index = i;
                    min_j = cell.index;
                }
            }
        }
        (min_index, min_j)
    }
}

// Impl from DistanceMatrix

impl From<DistanceMatrix> for NjBtree {
    fn from(value: DistanceMatrix) -> Self {
        let matrix = value.matrix;

        let sum_cols = matrix
            .iter()
            .map(|row| row.iter().sum::<f64>())
            .collect::<Vec<f64>>();

        let n = matrix.len();

        let mut btrees = Vec::with_capacity(n);
        for i in 0..n {
            let mut bt = BTreeSet::new();
            for j in i..n {
                if i != j {
                    let q = (matrix[i][j] * (n - 2) as f64) - sum_cols[i] - sum_cols[j];
                    let cell = Cell { index: j, value: q };
                    bt.insert(cell);
                }
            }
            btrees.push(bt);
        }
        Self { btrees, sum_cols }
    }
}

mod tests {
    use crate::{
        naive_neighbor_joining::{self, NjMatrix},
        phylip_distance_matrix::DistanceMatrix,
        phylogenetic_tree::PhyloTree,
    };

    use super::*;

    #[test]
    fn test_find_neighbors() {
        let d = DistanceMatrix::from(PhyloTree::random(5));
        let mut naive = NjMatrix::new(d.clone());
        let expected = naive_neighbor_joining::find_neighbors(&mut naive);

        let rapid_nj = NjBtree::from(d);
        assert_eq!(rapid_nj.min(), expected)
    }
    #[test]
    fn test_nj_btree() {
        let d = DistanceMatrix::from(PhyloTree::random(5));
        let mut q = d.matrix.clone();
        let n = q.len();
        for i in 0..n {
            for j in 0..n {
                q[i][j] = (d.matrix[i][j] * (n - 2) as f64)
                    - d.matrix[i].iter().sum::<f64>()
                    - d.matrix[j].iter().sum::<f64>();
                if i >= j {
                    q[i][j] = f64::INFINITY;
                }
            }
        }
        let expected_argmin: Vec<usize> = q
            .iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .min_by(|(_, a), (_, b)| a.total_cmp(b))
                    .map(|(index, _)| index)
            })
            .map(|x| x.unwrap())
            .collect();
        let nj = NjBtree::from(d);
        assert_eq!(nj.btrees.len(), n);
        for index in 0..n {
            let argmin = nj.btrees[index].first();
            if argmin.is_none() {
                continue;
            }
            assert_eq!(argmin.unwrap().index, expected_argmin[index]);
        }
    }
}
