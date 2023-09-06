pub mod naive_neighbor_joining;
pub mod phylip_distance_matrix;
pub mod phylogenetic_tree;
use petgraph::{
    algo,
    dot::{self, Dot},
};

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
    // let original_tree = PhyloTree::random(6);
    // println!(
    //     "{:?}",
    //     Dot::with_config(&original_tree.tree, &[dot::Config::EdgeNoLabel])
    // );
    // let d = phylip_distance_matrix::DistanceMatrix::from(original_tree.clone());
    // println!("{:?}", d);

    // let tree = naive_neighbor_joining::naive_neighbor_joining(d).unwrap();
    // println!(
    //     "{:?}",
    //     Dot::with_config(&tree.tree, &[dot::Config::EdgeNoLabel])
    // );
    // dbg!(petgraph::algo::is_isomorphic(
    //     &original_tree.tree,
    //     &tree.tree
    // ));
}
