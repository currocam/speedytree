mod btree;
mod cell;

use std::{
    cmp::{max, min},
    collections::BTreeSet,
    mem::swap,
};

use btree::NJRow;

use crate::{
    phylip_distance_matrix::DistanceMatrix, phylogenetic_tree::PhyloTree, rapidnj::cell::NJCell,
    ResultBox,
};

struct NjBtree {
    btrees: Vec<NJRow>,
    sum_cols: Vec<f64>,
    u_max: (f64, usize),
}

impl From<&DistanceMatrix> for NjBtree {
    fn from(value: &DistanceMatrix) -> Self {
        let matrix = &value.matrix;
        let sum_cols = matrix
            .iter()
            .map(|row| row.iter().sum::<f64>())
            .collect::<Vec<f64>>();
        let mut u_max = (f64::MIN, 0);
        for (i, val) in sum_cols.iter().enumerate() {
            if val > &u_max.0 {
                u_max = (*val, i);
            }
        }
        let btrees = matrix.iter().map(|x| NJRow::from(x.to_vec(), 0)).collect();
        Self {
            btrees,
            sum_cols,
            u_max,
        }
    }
}

fn find_neighbors(bt_tree: &NjBtree) -> (usize, usize, f64) {
    let n = bt_tree.btrees.len();
    let mut q_min = f64::INFINITY;
    let mut neighbors = (0, 0);
    let mut best_cell = NJCell::new(0, f64::NEG_INFINITY);

    for (i, row) in bt_tree.btrees.iter().enumerate() {
        for cell in row.val.iter() {
            let j = cell.index;
            let q = (cell.value * (n - 2) as f64) - bt_tree.sum_cols[i] - bt_tree.sum_cols[j];
            if i == j {
                continue;
            }
            if q < q_min {
                q_min = q;
                neighbors = (i, j);
                best_cell = cell.clone();
            } else {
                let best_q_min = q + bt_tree.sum_cols[j] - bt_tree.u_max.0;
                if best_q_min >= q_min {
                    break;
                }
            }
        }
    }
    (neighbors.0, neighbors.1, best_cell.value)
}

pub fn rapid_nj(mut dist: DistanceMatrix) -> ResultBox<PhyloTree> {
    let mut t = PhyloTree::new(&dist.names);
    //dbg!(&t.nodes);
    let mut bt = NjBtree::from(&dist);
    // While n nodes is not 2n -2
    let mut n = dist.matrix.len();
    while n > 3 {
        // Find the minimum element in the distance matrix
        let (i, j, dij) = find_neighbors(&bt);
        let s = (dist.matrix.len() - 2) as f64;
        let dist_ui = (dij + bt.sum_cols[i] / s - bt.sum_cols[j] / s) / 2.0;
        let dist_uj = dij - dist_ui;
        t.merge_rapid_nj(i, j, dist_ui, dist_uj);
        update_rapid_nj(i, j, dij, &mut bt);
        n -= 1;
    }
    let mut indices = vec![];
    for (i, row) in bt.btrees.iter().enumerate() {
        if !row.val.is_empty() {
            indices.push(i);
        }
    }
    let i = indices[0];
    let j = indices[1];
    let m = indices[2];

    let dij = bt.btrees[i]
        .val
        .iter()
        .find(|c| c.index == j)
        .unwrap()
        .value;
    let djm = bt.btrees[j]
        .val
        .iter()
        .find(|c| c.index == m)
        .unwrap()
        .value;
    let dim = bt.btrees[i]
        .val
        .iter()
        .find(|c| c.index == m)
        .unwrap()
        .value;
    let (i, j, m) = (t.nodes[&i], t.nodes[&j], t.nodes[&m]);
    let dvi = (dij + dim - djm) / 2.0;
    let dvj = (dij + djm - dim) / 2.0;
    let dvm = (dim + djm - dij) / 2.0;
    let v = t.tree.add_node("".to_owned());
    t.tree.add_edge(v, i, dvi);
    t.tree.add_edge(v, j, dvj);
    t.tree.add_edge(v, m, dvm);
    Ok(t)
}

fn update_rapid_nj(i: usize, j: usize, dij: f64, bts: &mut NjBtree) {
    let mut i_cells = NJRow::new();
    swap(&mut i_cells, &mut bts.btrees[i]);
    let mut j_cells = NJRow::new();
    swap(&mut j_cells, &mut bts.btrees[j]);
    // Update the col sums
    for cell in i_cells.val.iter() {
        bts.sum_cols[cell.index] -= cell.value;
    }
    for cell in j_cells.val.iter() {
        bts.sum_cols[cell.index] -= cell.value;
    }
    // Remove the ith and jth value to each row
    for bt in bts.btrees.iter_mut() {
        bt.val.retain(|c| c.index != i && c.index != j);
    }
    // Create new distances
    let mut new_distances: Vec<Option<f64>> = vec![None; bts.btrees.len()];
    for cell in i_cells.val {
        let prev = new_distances[cell.index].unwrap_or(0.0);
        new_distances[cell.index] = Some(prev + cell.value);
    }
    for cell in j_cells.val {
        let prev = new_distances[cell.index].unwrap_or(0.0);
        new_distances[cell.index] = Some(prev + cell.value);
    }
    new_distances[i] = None;
    new_distances[j] = None;
    bts.sum_cols[i] = f64::INFINITY;
    bts.sum_cols[j] = f64::INFINITY;

    let mut new_btree = BTreeSet::new();
    let n = bts.btrees.len();
    let mut new_col_sum = 0.0;
    for (index, val) in new_distances.iter().enumerate() {
        if let Some(dist) = val {
            let dist = (dist - dij) / 2.0;
            new_btree.insert(NJCell::new(index, dist));
            bts.btrees[index].val.insert(NJCell::new(n, dist));
            bts.sum_cols[index] += dist;
            new_col_sum += dist;
        }
    }
    bts.btrees.push(NJRow { val: new_btree });
    bts.sum_cols.push(new_col_sum);
}

fn terminate_nj(mut t: PhyloTree, dist: NjBtree) -> PhyloTree {
    let mut indices = vec![];
    for (i, row) in dist.btrees.iter().enumerate() {
        if !row.val.is_empty() {
            indices.push(i);
        }
    }
    let i = indices[0];
    let j = indices[1];
    let m = indices[2];
    t
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::vec;
    #[test]
    fn test_find_neighbors() {
        let mat = vec![
            vec![0.0, 5.0, 9.0, 9.0, 8.0],
            vec![5.0, 0.0, 10.0, 10.0, 9.0],
            vec![9.0, 10.0, 0.0, 8.0, 7.0],
            vec![9.0, 10.0, 8.0, 0.0, 3.0],
            vec![8.0, 9.0, 7.0, 3.0, 0.0],
        ];
        let bt_tree = NjBtree::from(&DistanceMatrix {
            matrix: mat,
            names: vec![
                "A".to_string(),
                "B".to_string(),
                "C".to_string(),
                "D".to_string(),
                "E".to_string(),
            ],
        });
        let index = find_neighbors(&bt_tree);
        assert_eq!(index, (0, 1, 5.0));
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
        let phylo = rapid_nj(DistanceMatrix {
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
}
