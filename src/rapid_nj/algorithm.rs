use crate::{distances::DistanceMatrix, ResultBox, Tree};

use super::{phylo_tree::PhyloTree, qmatrix::QMatrix};

pub fn rapid_nj(dist: DistanceMatrix, chunk_size: usize) -> ResultBox<Tree> {
    let mut q = QMatrix::from(&dist);
    q.set_chunk_size(chunk_size);
    let mut t = PhyloTree::build(&dist.names);
    while q.n_leaves() > 3 {
        // Find the minimum element in the distance matrix
        let (i, j) = q.find_neighbors();
        let (dist_ui, dist_uj) = q.new_node_distances(i, j);
        t.merge_neighbors(i, j, dist_ui, dist_uj);
        q.update(i, j);
    }
    Ok(terminate_nj(t, q))
}

fn terminate_nj(tree: PhyloTree, q: QMatrix) -> Tree {
    // Unmerged nodes are those that are not None in q.trees
    let unmerged = q.unmerged_nodes();
    let (i, j, m) = (unmerged[0], unmerged[1], unmerged[2]);

    let dvi = (q.distance(i, j) + q.distance(i, m) - q.distance(j, m)) / 2.0;
    let dvj = (q.distance(i, j) + q.distance(j, m) - q.distance(i, m)) / 2.0;
    let dvm = (q.distance(i, m) + q.distance(j, m) - q.distance(i, j)) / 2.0;

    let (i, j, m) = (tree.nodes[&i], tree.nodes[&j], tree.nodes[&m]);
    let mut tree = tree.tree;
    let v = tree.add_node("".to_owned());
    tree.add_edge(v, i, dvi);
    tree.add_edge(v, j, dvj);
    tree.add_edge(v, m, dvm);

    tree
}
