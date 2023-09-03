pub mod distances;
pub mod neighbors;
pub mod newick;

use ndarray::{ArrayBase, Dim, OwnedRepr};
use std::{error, io, process};

use crate::{
    distances::{distance_matrix_from_stdin, DistanceMatrix},
    neighbors::neighbor_joining,
    newick::PhylogeneticTree,
};

type ResultBox<T> = std::result::Result<T, Box<dyn error::Error>>;
type M = ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>;

#[derive(Debug)]
pub struct Config;

impl Config {
    pub fn build() -> ResultBox<Config> {
        Ok(Config)
    }
}

pub fn run(config: Config) {
    dbg!(&config);
    let distance_mat: DistanceMatrix =
        distance_matrix_from_stdin(io::stdin()).unwrap_or_else(|err| {
            eprintln!("{err}");
            process::exit(1);
        });
    dbg!(&distance_mat);
    let tree = neighbor_joining(&distance_mat).unwrap_or_else(|err| {
        eprintln!("{err}");
        process::exit(1);
    });
    dbg!(&tree);
    let newick = tree.into_newick().unwrap_or_else(|| {
        eprintln!("Empty tree");
        process::exit(1);
    });
    println!("{}", newick);
}
