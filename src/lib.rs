#![feature(allocator_api)]
#![feature(btreemap_alloc)]
mod distances;
mod naive_nj;
mod newick;
pub mod property_tests;
pub mod rapid_nj;

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

// Define Enum for Algorithm
#[derive(Debug)]
enum Algorithm {
    Naive,
    Rapid,
}

#[derive(Debug)]
pub struct Config {
    algo: Algorithm,
}

impl Config {
    pub fn build(mut args: impl Iterator<Item = String>) -> ResultBox<Config> {
        // Let match the algorithm, if not specified, use Naive
        let error_msg = "Usage: nj [-r|--rapid] [-n|--naive] < input.phy > output.nwk";
        args.next(); // Skip the first argument
        let algo = match args.next() {
            Some(algo) => match algo.as_str() {
                "-r" | "--rapid" => Algorithm::Rapid,
                "-n" | "--naive" => Algorithm::Naive,
                _ => {
                    return Err(From::from(format!(
                        "Invalid algorithm: {}. \n{}",
                        algo.as_str(),
                        error_msg
                    )))
                }
            },
            None => Algorithm::Naive,
        };
        Ok(Config { algo })
    }
}

pub fn run(config: Config) {
    //dbg!(&config);
    let reader = io::stdin().lock();
    let d = DistanceMatrix::build_from_phylip(reader).unwrap_or_else(|err| {
        eprintln!("{err}");
        process::exit(1);
    });

    // Let match the algorithm
    //dbg!(&d);
    let d = match config.algo {
        Algorithm::Naive => naive_neighbor_joining(d),
        Algorithm::Rapid => rapid_nj(d),
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
