pub mod naive_neighbor_joining;
pub mod phylip_distance_matrix;
pub mod phylogenetic_tree;
use petgraph::dot::{self, Dot};

use std::{error, io, process};

use crate::phylogenetic_tree::PhyloTree;

type ResultBox<T> = std::result::Result<T, Box<dyn error::Error>>;

#[derive(Debug)]
pub struct Config;

impl Config {
    pub fn build() -> ResultBox<Config> {
        Ok(Config)
    }
}

pub fn run(config: Config) {
    dbg!(&config);
    let distance_mat = phylip_distance_matrix::read_phylip_distance_matrix(io::stdin().lock())
        .unwrap_or_else(|err| {
            eprintln!("{err}");
            process::exit(1);
        });
    //dbg!(&distance_mat);
    let tree = naive_neighbor_joining::naive_neighbor_joining(distance_mat);
    let graph = &tree.unwrap().tree;
    println!(
        "{:?}",
        Dot::with_config(&graph, &[dot::Config::EdgeNoLabel])
    );

    // let tree = PhyloTree::random(5 + 1);
    // let graph = &tree.tree;
    // println!(
    //     "{:?}",
    //     Dot::with_config(&graph, &[dot::Config::EdgeNoLabel])
    // );
    // let d = phylip_distance_matrix::DistanceMatrix::from(tree);
    // dbg!(&d);
    // let tree = naive_neighbor_joining::naive_neighbor_joining(d);
    // let graph = &tree.unwrap().tree;
    // println!(
    //     "{:?}",
    //     Dot::with_config(&graph, &[dot::Config::EdgeNoLabel])
    // );
}
