//! speedytree: Command line tool for Neighbor Joining of biological sequences
//!
//! It implements different heuristics for fast Neighbor-Joining.
//!
//!   1. Naive Neighbor-Joining
//!   2. RapidNJ
//!   3. Hybrid
//!

#![feature(allocator_api)]
#![feature(btreemap_alloc)]
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
mod distances;
mod naive_nj;
mod newick;
mod rapid_nj;
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

/// Available algorithms in the program
#[derive(Debug, Clone)]
pub enum Algorithm {
    /// Naive neighbor joining
    Naive,
    /// Rapid neighbor joining
    RapidNJ,
    /// Hybrid neighbor joining
    Hybrid,
}
/// Main function of the crate
pub fn run(config: configuration::Config) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()
        .unwrap();

    let reader = io::stdin().lock();
    let d = DistanceMatrix::build_from_phylip(reader).unwrap_or_else(|err| {
        eprintln!("{err}");
        process::exit(1);
    });

    let d = match config.algo {
        Algorithm::Naive => naive_neighbor_joining(d),
        Algorithm::RapidNJ => rapid_nj(d, config.chunk_size),
        Algorithm::Hybrid => {
            let naive_steps = d.size() * config.naive_percentage / 100;
            dbg!(naive_steps);
            neighbor_joining(d, naive_steps, config.chunk_size)
        }
    };
    let graph = d.unwrap_or_else(|err| {
        eprintln!("{err}");
        process::exit(1);
    });
    io::stdout()
        .write_all(to_newick(&graph).as_bytes())
        .unwrap_or_else(|err| {
            eprintln!("{err}");
            process::exit(1);
        });
}
