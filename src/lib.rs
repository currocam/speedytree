pub mod naive_neighbor_joining;
pub mod newick;
pub mod phylip_distance_matrix;
pub mod phylogenetic_tree;
use petgraph::dot::{self, Dot};

use std::{error, io, process};

use crate::{newick::to_newick, phylogenetic_tree::PhyloTree};

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
    dbg!(&graph);
    println!(
        "{:?}",
        Dot::with_config(&graph, &[dot::Config::EdgeNoLabel])
    );

    let mut tree = PhyloTree::random(5);
    // Set weights to 0.1, 0.2, 0.3, 0.4, 0.5, 0.6
    for edge in tree.tree.edge_indices() {
        tree.tree[edge] = (edge.index() as f64 + 1.0) / 10.0;
    }
    let graph = &tree.tree;
    dbg!(&graph);
    println!("{}", to_newick(graph));
}
