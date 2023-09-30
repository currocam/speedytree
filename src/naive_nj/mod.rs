// Skeleton code for the Naive Neighbor Joining algorithm.
mod algorithm;
// QMatrix is a helper struct for the Naive Neighbor Joining algorithm.
mod matrix;
// PhyloTree is a helper struct for the Naive Neighbor Joining algorithm.
mod phylo_tree;
// Export the public interface of the Naive Neighbor Joining algorithm.
pub use algorithm::naive_neighbor_joining;
pub use matrix::QMatrix;
pub use phylo_tree::PhyloTree;
pub use algorithm::terminate_nj;
pub struct DataNaiveNJ {
    pub qmatrix : matrix::QMatrix,
    pub phylo_tree : phylo_tree::PhyloTree,
}
