pub mod naive_neighbor_joining;
pub mod newick;
pub mod phylip_distance_matrix;
pub mod phylogenetic_tree;

use std::{error, io, process};

use crate::newick::to_newick;

type ResultBox<T> = std::result::Result<T, Box<dyn error::Error>>;

#[derive(Debug)]
pub struct Config;

impl Config {
    pub fn build() -> ResultBox<Config> {
        Ok(Config)
    }
}

pub fn run(_config: Config) {
    //dbg!(&config);
    let distance_mat = phylip_distance_matrix::read_phylip_distance_matrix(io::stdin().lock())
        .unwrap_or_else(|err| {
            eprintln!("{err}");
            process::exit(1);
        });
    //dbg!(&distance_mat);
    let tree = naive_neighbor_joining::naive_neighbor_joining(distance_mat);
    let graph = &tree.unwrap().tree;
    println!("{}", to_newick(graph));
}
