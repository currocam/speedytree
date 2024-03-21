//! Canonical and RapidNJ implementations of Neighbor-joining in Rust
//!
//! Provides Rust implementation of the Canonical algorithm and something in the spirit of RapidNJ but with B-trees. Work in progress.
//! A minimal example is provided here. 
//! ```
//! use speedytree::DistanceMatrix;
//! use speedytree::robinson_foulds;
//! use speedytree::{NeighborJoiningSolver, Canonical, RapidBtrees, Hybrid};
//! use rayon;
//!// Raw Phylip format
//!let input = "5
//!    a	0	5	9	9	8
//!    b	5	0	10	10	9
//!    c	9	10	0	8	7
//!    d	9	10	8	0	3
//!    e	8	9	7	3	0
//!".as_bytes();;
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

pub mod hybrid_nj;
/// Property tests for neighbor joining algorithm
pub mod property_tests;
pub mod distances;
pub mod naive_nj;
pub mod newick;
pub mod rapid_nj;
use std::error;
pub use distances::DistanceMatrix;
pub use property_tests::tree_distances::{branch_score, robinson_foulds};
type ResultBox<T> = std::result::Result<T, Box<dyn error::Error>>;
pub type Tree = petgraph::graph::UnGraph<String, f64>;

pub struct NeighborJoiningSolver<U> {algo: U, dist: DistanceMatrix}
pub struct Canonical{}
impl NeighborJoiningSolver<Canonical> {
  pub fn build(dist: DistanceMatrix) -> Self{
    NeighborJoiningSolver { algo: Canonical {}, dist: dist}
  }
  pub fn default(dist: DistanceMatrix) -> Self{
    Self::build(dist)
  }
  pub fn solve(self) -> ResultBox<Tree> {
    naive_nj::canonical_neighbor_joining(self.dist)
  }
}
pub struct RapidBtrees{chunk_size: usize}
impl NeighborJoiningSolver<RapidBtrees> {
  pub fn build(dist: DistanceMatrix, chunk_size: usize) -> Self{
    NeighborJoiningSolver { algo: RapidBtrees {chunk_size}, dist: dist}
  }
  pub fn default(dist: DistanceMatrix) -> Self {
    let n = dist.size();
    let threads = rayon::current_num_threads();
    let chunk_size = std::cmp::max(n/threads, 1);
    if chunk_size == 0 {
        panic!("Chunk size cannot be zero.");
    }
    NeighborJoiningSolver { algo: RapidBtrees {chunk_size}, dist: dist}
  }
  pub fn set_chunk_size(self, chunk_size: usize) -> Self {
    if chunk_size < 1 {
        panic!("Chunk size  must be > 0.");
    }
    Self::build(self.dist, chunk_size)
  }
  pub fn solve(self) -> ResultBox<Tree> {
    rapid_nj::rapid_nj(self.dist, self.algo.chunk_size)
  }
}

pub struct Hybrid{chunk_size: usize, canonical_iters: usize}
impl NeighborJoiningSolver<Hybrid> {
  pub fn build(dist: DistanceMatrix, chunk_size: usize, canonical_iters: usize) -> Self{
    NeighborJoiningSolver { algo: Hybrid {chunk_size, canonical_iters}, dist: dist}
  }
  pub fn default(dist: DistanceMatrix) -> Self {
    let n = dist.size();
    let threads = rayon::current_num_threads();
    let chunk_size = std::cmp::max(n/threads, 1);
    let canonical_iters = std::cmp::max(n/2, 1);
    NeighborJoiningSolver { algo: Hybrid {chunk_size, canonical_iters}, dist: dist}
  }
  pub fn solve(self) -> ResultBox<Tree> {
    hybrid_nj::neighbor_joining(self.dist, self.algo.canonical_iters, self.algo.chunk_size)
  }
  pub fn set_chunk_size(self, chunk_size: usize) -> Self {
    if chunk_size < 1 {
        panic!("Chunk size  must be > 0.");
    }
    Self::build(self.dist, chunk_size, self.algo.canonical_iters)
  }
  pub fn set_canonical_steps(self, n: usize) -> Self {
    if n < 1 {
        panic!("n must be > 0.");
    }
    Self::build(self.dist, self.algo.chunk_size, n)
  }
}
