//! Canonical and RapidNJ implementations of Neighbor-joining in Rust
//!
//! Provides Rust implementation of the Canonical algorithm and something in the spirit of RapidNJ but with B-trees. Some helper functions are also provided. 
//! 
//! This repository provides (a) a command line application that reads files in Phylip format and (b) a small library you can easily integrate in your projects. It relies on rayon for parallelization.
//! ## Example
//! A minimal example of the library is provided here. You can read more about the command line app by running speedytree -h
//! ```
//! use speedytree::DistanceMatrix;
//! use speedytree::robinson_foulds;
//! use speedytree::{NeighborJoiningSolver, Canonical, RapidBtrees, Hybrid};
//! use rayon;
//!// Raw Phylip format
//!let input = "5
//!    a    0    5    9    9    8
//!    b    5    0    10    10    9
//!    c    9    10    0    8    7
//!    d    9    10    8    0    3
//!    e    8    9    7    3    0
//!".as_bytes();
//!let d = DistanceMatrix::read_from_phylip(input).unwrap();
//! // Canonical Neighbor Joining. Optimal for small problems.
//! let tree1 = NeighborJoiningSolver::<Canonical>::default(d.clone())
//!   .solve()
//!   .unwrap();
//! // An algorithm in the spirit of RapidNJ using B-trees. Optimal for big problems.
//! let tree2 = NeighborJoiningSolver::<RapidBtrees>::default(d.clone())
//!   .solve()
//!   .unwrap();
//! assert_eq!(robinson_foulds(&tree1, &tree2), 0);
//!
//! // You can improve the speed a lot by using multiple threads (see rayon::ThreadPoolBuilder::new())
//! // If so, you may want to tune the chunk size every worker uses
//!   let tree3 = NeighborJoiningSolver::<RapidBtrees>::default(d.clone())
//!   .set_chunk_size(4)
//!   .solve()
//!   .unwrap();
//! // An algorithm that starts using RapidBtrees, later uses Canonical
//! // Optimal for intermediate problems (only if you tune the number of "canonical steps")
//! let tree4 = NeighborJoiningSolver::<Hybrid>::default(d.clone())
//!   .set_canonical_steps(2)
//!   .solve()
//!   .unwrap();
//! assert_eq!(robinson_foulds(&tree3, &tree4), 0);
//! ```

mod distances;
mod hybrid_nj;
mod naive_nj;
mod newick;
/// Property tests for neighbor joining algorithm
mod property_tests;
mod rapid_nj;
pub use distances::DistanceMatrix;
pub use newick::to_newick;
pub use property_tests::tree_distances::{branch_score, robinson_foulds};

use std::error;
type ResultBox<T> = std::result::Result<T, Box<dyn error::Error>>;
/// An undirected network built in top of [Petgraph](https://github.com/petgraph/petgraph). Internal nodes have empty names.
pub type Tree = petgraph::graph::UnGraph<String, f64>;

/// Generic solver
pub struct NeighborJoiningSolver<U> {
    algo: U,
    dist: DistanceMatrix,
}
/// Canonical Neighbor-Joining, similar to [QuickTree](https://github.com/khowe/quicktree). It runs on cubic time (worst and best case). It uses quadratic memory.  
pub struct Canonical {}
impl NeighborJoiningSolver<Canonical> {
    /// Construct solver from parameters
    pub fn build(dist: DistanceMatrix) -> Self {
        NeighborJoiningSolver {
            algo: Canonical {},
            dist,
        }
    }
    /// Default solver
    pub fn default(dist: DistanceMatrix) -> Self {
        Self::build(dist)
    }
    /// Solve the Neighbor-Joining problem
    pub fn solve(self) -> ResultBox<Tree> {
        naive_nj::canonical_neighbor_joining(self.dist)
    }
}
/// In the spirit of [RapidNJ](https://birc.au.dk/software/rapidnj/), but with B-trees. It runs on n^2 log(n) time best case and cubic time worst case.  It uses quadratic memory (with a higher constant).
pub struct RapidBtrees {
    chunk_size: usize,
}
impl NeighborJoiningSolver<RapidBtrees> {
    /// Construct solver from parameters
    pub fn build(dist: DistanceMatrix, chunk_size: usize) -> Self {
        NeighborJoiningSolver {
            algo: RapidBtrees { chunk_size },
            dist,
        }
    }
    /// Default solver (based on available rayon threads)
    pub fn default(dist: DistanceMatrix) -> Self {
        let n = dist.size();
        let threads = rayon::current_num_threads();
        let chunk_size = std::cmp::max(n / threads, 1);
        if chunk_size == 0 {
            panic!("Chunk size cannot be zero.");
        }
        NeighborJoiningSolver {
            algo: RapidBtrees { chunk_size },
            dist,
        }
    }
    /// Set chunk size (for every worker)
    pub fn set_chunk_size(self, chunk_size: usize) -> Self {
        if chunk_size < 1 {
            panic!("Chunk size  must be > 0.");
        }
        Self::build(self.dist, chunk_size)
    }
    /// Solve the Neighbor-Joining problem
    pub fn solve(self) -> ResultBox<Tree> {
        rapid_nj::rapid_nj(self.dist, self.algo.chunk_size)
    }
}

/// A mix of the Canonical and RapidBtrees. First, it starts with RapidBtrees (less lookups, but with an overhead), and then it changes the strategy. 
pub struct Hybrid {
    chunk_size: usize,
    canonical_iters: usize,
}
impl NeighborJoiningSolver<Hybrid> {
    /// Construct solver from parameters
    pub fn build(dist: DistanceMatrix, chunk_size: usize, canonical_iters: usize) -> Self {
        NeighborJoiningSolver {
            algo: Hybrid {
                chunk_size,
                canonical_iters,
            },
            dist,
        }
    }
    /// Default solver (based on available rayon threads and problem size)
    pub fn default(dist: DistanceMatrix) -> Self {
        let n = dist.size();
        let threads = rayon::current_num_threads();
        let chunk_size = std::cmp::max(n / threads, 1);
        let canonical_iters = std::cmp::max(n / 2, 1);
        NeighborJoiningSolver {
            algo: Hybrid {
                chunk_size,
                canonical_iters,
            },
            dist,
        }
    }
    /// Solve the Neighbor-Joining problem
    pub fn solve(self) -> ResultBox<Tree> {
        hybrid_nj::neighbor_joining(self.dist, self.algo.canonical_iters, self.algo.chunk_size)
    }
    /// Set chunk size (for every worker)
    pub fn set_chunk_size(self, chunk_size: usize) -> Self {
        if chunk_size < 1 {
            panic!("Chunk size  must be > 0.");
        }
        Self::build(self.dist, chunk_size, self.algo.canonical_iters)
    }
    /// Set number of canonical iterations will be done
    pub fn set_canonical_steps(self, n: usize) -> Self {
        if n < 1 {
            panic!("n must be > 0.");
        }
        Self::build(self.dist, self.algo.chunk_size, n)
    }
    /// Set fraction of canonical iterations will be done
    pub fn set_canonical_percentage(self, prop: f64) -> Self {
        if prop <= 0.0 || prop >= 1.0 {
            panic!("Proportion must be between 0 and 1.");
        }
        let n = self.dist.size() as f64 * prop / 100.0;
        Self::build(self.dist, self.algo.chunk_size, n as usize)
    }
}
