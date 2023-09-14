use std::cmp::{max, min};

use crate::{phylip_distance_matrix::DistanceMatrix, phylogenetic_tree::PhyloTree, ResultBox};

#[derive(Debug)]
struct NjMatrix {
    matrix: Vec<Vec<f64>>,
    sum_cols: Vec<f64>,
}

impl NjMatrix {
    fn new(d: DistanceMatrix) -> Self {
        let matrix = d.matrix;
        let sum_cols = matrix
            .iter()
            .map(|row| row.iter().sum::<f64>())
            .collect::<Vec<f64>>();
        Self { matrix, sum_cols }
    }
}

pub fn naive_neighbor_joining(dist: DistanceMatrix) -> ResultBox<PhyloTree> {
    let mut t = PhyloTree::new(&dist.names);
    //dbg!(&t.nodes);
    let mut dist = NjMatrix::new(dist);
    // While n nodes is not 2n -2
    while dist.matrix.len() > 3 {
        // Find the minimum element in the distance matrix
        let (i, j) = find_neighbors(&mut dist);
        let s = (dist.matrix.len() - 2) as f64;
        let dist_ui = (dist.matrix[i][j] + dist.sum_cols[i] / s - dist.sum_cols[j] / s) / 2.0;
        let dist_uj = dist.matrix[i][j] - dist_ui;
        t.merge_neighbors(i, j, dist_ui, dist_uj);
        update_distance_matrix(i, j, &mut dist);
    }
    t = terminate_nj(t, &mut dist);
    Ok(t)
}

fn terminate_nj(mut t: PhyloTree, d: &mut NjMatrix) -> PhyloTree {
    let (i, j, m) = (t.nodes[&0], t.nodes[&1], t.nodes[&2]);
    let dvi = (d.matrix[0][1] + d.matrix[0][2] - d.matrix[1][2]) / 2.0;
    let dvj = (d.matrix[0][1] + d.matrix[1][2] - d.matrix[0][2]) / 2.0;
    let dvm = (d.matrix[0][2] + d.matrix[1][2] - d.matrix[0][1]) / 2.0;
    let v = t.tree.add_node("".to_owned());
    t.tree.add_edge(v, i, dvi);
    t.tree.add_edge(v, j, dvj);
    t.tree.add_edge(v, m, dvm);
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

    use petgraph::visit::IntoNeighbors;

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
            let original_tree: PhyloTree = PhyloTree::random(i);
            let d = DistanceMatrix::from(original_tree.clone());
            let original_tree = original_tree.tree;
            let tree = naive_neighbor_joining(d).unwrap();
            assert!(petgraph::algo::is_isomorphic(&original_tree, &tree.tree));
        }
    }
}
