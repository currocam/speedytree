use crate::{phylip_distance_matrix::DistanceMatrix, phylogenetic_tree::PhyloTree, ResultBox};

use super::matrix::QMatrix;

pub fn naive_neighbor_joining(dist: DistanceMatrix) -> ResultBox<PhyloTree> {
    let mut t = PhyloTree::new(&dist.names);
    let mut q = QMatrix::new(dist);

    while q.n_leaves() > 3 {
        // Find the minimum element in the distance matrix
        let (i, j) = q.find_neighbors();
        let (dist_ui, dist_uj) = q.new_node_distances(i, j);
        t.merge_neighbors(i, j, dist_ui, dist_uj);
        q.update_distance_matrix(i, j);
    }

    t = terminate_nj(t, q);
    Ok(t)
}

fn terminate_nj(mut t: PhyloTree, q: QMatrix) -> PhyloTree {
    let (i, j, m) = (t.nodes[&0], t.nodes[&1], t.nodes[&2]);
    let dvi = (q.distance(0, 1) + q.distance(0, 2) - q.distance(1, 2)) / 2.0;
    let dvj = (q.distance(0, 1) + q.distance(1, 2) - q.distance(0, 2)) / 2.0;
    let dvm = (q.distance(0, 2) + q.distance(1, 2) - q.distance(0, 1)) / 2.0;
    let v = t.tree.add_node("".to_owned());
    t.tree.add_edge(v, i, dvi);
    t.tree.add_edge(v, j, dvj);
    t.tree.add_edge(v, m, dvm);
    t
}
