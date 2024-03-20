mod algorithm;
mod node;
mod phylo_tree;
mod qmatrix;
pub use algorithm::rapid_nj;
pub(crate) use phylo_tree::PhyloTree;
pub(crate) use qmatrix::QMatrix;

pub(crate) struct DataRapidNJ {
    pub qmatrix: QMatrix,
    pub phylo_tree: phylo_tree::PhyloTree,
}

impl DataRapidNJ {
    pub fn new(qmatrix: QMatrix, phylo_tree: phylo_tree::PhyloTree) -> Self {
        Self {
            qmatrix,
            phylo_tree,
        }
    }
}
