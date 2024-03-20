//! speedytree: Command line tool for Neighbor Joining of biological sequences
//!
//! It implements different heuristics for fast Neighbor-Joining.
//!
//!   1. Naive Neighbor-Joining
//!   2. RapidNJ
//!   3. Hybrid
//!

#![warn(missing_docs)]
/// Hybrid neighbor joining algorithm
/// This approach is a hybrid between the naive neighbor joining and the rapid neighbor joining.
/// The idea is to use the rapidnj heuristic first to potentially stop the algorithm early,
/// and then use the naive neighbor joining to finish the algorithm, which is faster
/// in practice, but performs more comparisons in theory.
/// However, both algorithms are O(n^3), so the difference is not that big.
pub mod hybrid_nj;
/// Property tests for neighbor joining algorithm
pub mod property_tests;

/// Configuration of the program
pub mod configuration;
pub mod distances;
pub mod naive_nj;
pub mod newick;
pub mod rapid_nj;
use hybrid_nj::neighbor_joining;

use crate::distances::DistanceMatrix;
use crate::naive_nj::naive_neighbor_joining;
use crate::newick::to_newick;
use crate::rapid_nj::rapid_nj;
use std::{
    error,
    io::{self, Write},
    process,
};

type ResultBox<T> = std::result::Result<T, Box<dyn error::Error>>;
type Tree = petgraph::graph::UnGraph<String, f64>;

