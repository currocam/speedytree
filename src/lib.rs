mod distances;
mod naive_nj;
mod newick;
mod random_binary_trees;

use crate::distances::DistanceMatrix;
use crate::naive_nj::naive_neighbor_joining;
use crate::newick::to_newick;
use std::{error, io, process};

type ResultBox<T> = std::result::Result<T, Box<dyn error::Error>>;
type Tree = petgraph::graph::UnGraph<String, f64>;

#[derive(Debug)]
pub struct Config;

impl Config {
    pub fn build() -> ResultBox<Config> {
        Ok(Config)
    }
}

pub fn run(_config: Config) {
    //dbg!(&config);
    let reader = io::stdin().lock();
    let d = DistanceMatrix::build_from_phylip(reader).unwrap_or_else(|err| {
        eprintln!("{err}");
        process::exit(1);
    });
    //dbg!(&distance_mat);
    let tree = naive_neighbor_joining(d).unwrap_or_else(|err| {
        eprintln!("{err}");
        process::exit(1);
    });
    //dbg!(&tree);
    let graph = &tree;
    println!("{}", to_newick(graph));
}
