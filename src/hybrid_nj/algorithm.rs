use crate::{
    distances::DistanceMatrix, naive_nj::DataNaiveNJ, rapid_nj::DataRapidNJ, ResultBox, Tree,
};

/// Main function of the crate
/// This approach is a hybrid between the naive neighbor joining and the rapid neighbor joining.
/// If `naive_iters` is greater than n, then this function calls `naive_neighbor_joining` instead.
/// If `naive_iters` is less than 4, then this function calls `rapid_nj` instead.
/// Arguments:
/// * `dist` - Distance matrix
/// * `naive_iters` - Number of iterations to use the naive neighbor joining algorithm
///
/// Returns:
/// * `Ok(Tree)` - A phylogenetic tree
/// * `Err(Box<dyn error::Error>)` - An error
pub fn neighbor_joining(dist: DistanceMatrix, naive_iters: usize) -> ResultBox<Tree> {
    if dist.size() < 4 || naive_iters >= dist.size() {
        return crate::naive_neighbor_joining(dist);
    }
    if naive_iters < 4 {
        return crate::rapid_nj(dist);
    }
    let mut q = crate::rapid_nj::QMatrix::from(&dist);
    let mut t = crate::rapid_nj::PhyloTree::build(&dist.names);
    while q.n_leaves() > naive_iters {
        let (i, j) = q.find_neighbors();
        let (dist_ui, dist_uj) = q.new_node_distances(i, j);
        t.merge_neighbors(i, j, dist_ui, dist_uj);
        q.update(i, j);
    }
    // Convert to the inner data structure of the naive neighbor joining
    let data = DataNaiveNJ::from(DataRapidNJ::new(q, t));
    let mut q = data.qmatrix;
    let mut t = data.phylo_tree;
    while q.n_leaves() > 3 {
        let (i, j) = q.find_neighbors();
        let (dist_ui, dist_uj) = q.new_node_distances(i, j);
        t.merge_neighbors(i, j, dist_ui, dist_uj);
        q.update_distance_matrix(i, j);
    }
    Ok(crate::naive_nj::terminate_nj(t, q))
}
