// Skeleton code for the Naive Neighbor Joining algorithm.
mod algorithm;
// QMatrix is a helper struct for the Naive Neighbor Joining algorithm.
mod qmatrix;
// PhyloTree is a helper struct for the Naive Neighbor Joining algorithm.
mod phylo_tree;
// Export the public interface of the Naive Neighbor Joining algorithm.
pub use algorithm::canonical_neighbor_joining;
pub(crate) use algorithm::terminate_nj;
pub(crate) use phylo_tree::PhyloTree;
pub(crate) use qmatrix::QMatrix;
pub(crate) struct DataNaiveNJ {
    pub qmatrix: qmatrix::QMatrix,
    pub phylo_tree: phylo_tree::PhyloTree,
}
