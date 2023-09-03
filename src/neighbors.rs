use crate::{distances::DistanceMatrix, newick::PhylogeneticTree, ResultBox};

pub fn neighbor_joining(_d: &DistanceMatrix) -> ResultBox<PhylogeneticTree> {
    Ok(PhylogeneticTree::empty())
}
