/// This module contains the property tests for the crate.¨
/// Random additive binary trees are generated and then converted to distance matrices.
/// The distance matrices are then permuted and the neighbor joining algorithm is run on them.
pub mod random_additive_tree;
/// This module contains the property tests for the crate.¨
mod tests;
/// This module contains distance metrics for trees. Branch score-distance and Robinson Foulds are implemented.
pub mod tree_distances;
