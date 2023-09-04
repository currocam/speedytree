use petgraph::stable_graph::NodeIndex;

use crate::{phylip_distance_matrix::DistanceMatrix, phylogenetic_tree::PhyloTree, ResultBox};

pub fn naive_neighbor_joining(dist: &mut DistanceMatrix) -> ResultBox<PhyloTree> {
    let mut t = PhyloTree::new(&dist.names);
    // Naive NJ
    while dist.matrix.len() > 2 {
        // Find the minimum element in the distance matrix
        let (i, j) = find_neighbors(&dist.matrix);
        // Merge the nodes
        let u = t.merge_nodes(t.nodes[i], t.nodes[j]);
        // Update the distance matrix
        update_distance_matrix(i, j, u, &mut t, &mut dist.matrix);
        dbg!(t.tree.node_count());
    }
    Ok(t)
}

fn update_distance_matrix(
    i: usize,
    j: usize,
    u: NodeIndex,
    t: &mut PhyloTree,
    matrix: &mut Vec<Vec<f64>>,
) {
    // Swap i and j so that i < j
    let (i, j) = if i < j { (i, j) } else { (j, i) };
    let dij = matrix[i][j];
    let n = matrix.len();
    // Swap i and j so those are the last two rows
    matrix.swap(i, n - 2);
    matrix.swap(j, n - 1);
    // Update nodes in tree
    t.nodes.swap(i, n - 2);
    t.nodes.swap(j, n - 1);
    t.nodes.pop();
    t.nodes.pop();
    t.nodes.push(u);

    // Swap i and j so those are the last two columns
    for row in matrix.iter_mut() {
        row.swap(i, n - 2);
        row.swap(j, n - 1);
    }
    // Update the row.len() - 2 row
    for k in 0..matrix.len() - 2 {
        matrix[n - 2][k] = (matrix[i][k] + matrix[j][k] - matrix[i][j]) / 2.0;
    }
    matrix[n - 2][n - 2] = dij;
    // Remove the last row and every last column
    matrix.pop();
    for row in matrix.iter_mut() {
        row.pop();
    }
}

fn find_neighbors(matrix: &Vec<Vec<f64>>) -> (usize, usize) {
    let mut neighbors = (0, 0);
    let mut best_q = f64::INFINITY;
    let n = matrix.len();
    let sums = matrix
        .iter()
        .map(|row| row.iter().sum::<f64>())
        .collect::<Vec<f64>>();

    for i in 0..n {
        for j in i + 1..n {
            let q = (matrix[i][j] * (n - 2) as f64) - sums[i] - sums[j];
            if q < best_q {
                best_q = q;
                neighbors = (i, j);
            }
        }
    }
    neighbors
}

#[cfg(test)]
mod tests {
    use std::vec;

    use ndarray::array;

    use super::*;

    #[test]
    fn test_find_neighbors() {
        let mut mat = vec![
            vec![0.0, 5.0, 9.0, 9.0, 8.0],
            vec![5.0, 0.0, 10.0, 10.0, 9.0],
            vec![9.0, 10.0, 0.0, 8.0, 7.0],
            vec![9.0, 10.0, 8.0, 0.0, 3.0],
            vec![8.0, 9.0, 7.0, 3.0, 0.0],
        ];
        let index = find_neighbors(&mat);
        assert_eq!(index, (0, 1));
    }
}
