//! Canonical and RapidNJ implementations of Neighbor-joining in Rust
//!
//! Provides Rust implementation of the Canonical algorithm and something in the spirit of RapidNJ but with B-trees. Work in progress.
//! A minimal example is provided here. 
//! ```
//! use speedytree::distances::DistanceMatrix;
//! use speedytree::canonical_neighbor_joining;
//! use speedytree::rapid_nj_neighbor_joining;
//! use speedytree::robinson_foulds;
//!// Raw Phylip format
//!let input = "5
//!    a	0	5	9	9	8
//!    b	5	0	10	10	9
//!    c	9	10	0	8	7
//!    d	9	10	8	0	3
//!    e	8	9	7	3	0
//!".as_bytes();;
//!let distances = DistanceMatrix::read_from_phylip(input).unwrap();
//! // Use canonical algorithm
//!let tree1 = canonical_neighbor_joining(distances.clone()).unwrap();
//! // Use RapidNJ with b-trees-
//!let tree2 = rapid_nj_neighbor_joining(distances.clone(), 2).unwrap();
//! assert_eq!(robinson_foulds(tree1, tree2, 5), 0);
//! ```

pub mod hybrid_nj;
/// Property tests for neighbor joining algorithm
pub mod property_tests;
pub mod distances;
pub mod naive_nj;
pub mod newick;
pub mod rapid_nj;
use std::error;

type ResultBox<T> = std::result::Result<T, Box<dyn error::Error>>;
type Tree = petgraph::graph::UnGraph<String, f64>;

pub use naive_nj::canonical_neighbor_joining as canonical_neighbor_joining;
pub use rapid_nj::rapid_nj as rapid_nj_neighbor_joining;
pub use property_tests::tree_distances::robinson_foulds as robinson_foulds;
pub use property_tests::tree_distances::branch_score as branch_score;