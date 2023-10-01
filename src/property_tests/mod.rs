/// This module contains the property tests for the crate.¨
/// Random additive binary trees are generated and then converted to distance matrices.
/// The distance matrices are then permuted and the neighbor joining algorithm is run on them.
pub mod random_additive_tree;
/// Robinson-Foulds distance
pub mod robinson_foulds;
mod tests;
