#![feature(allocator_api)]
#![feature(btreemap_alloc)]
pub mod algo;
mod distances;
mod naive_nj;
mod newick;
pub mod property_tests;
pub mod rapid_nj;

use algo::neighbor_joining;

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
    RapidNJ,
    Hybrid,
}

#[derive(Debug)]
pub struct Config {
    algo: Algorithm,
    threads: usize,
}

impl Config {
    pub fn build(mut args: impl Iterator<Item = String>) -> ResultBox<Config> {
        // Let match the algorithm, if not specified, use Naive
        let error_msg = "Usage: nj [-r|--rapid] [-n|--naive] < input.phy > output.nwk";
        args.next(); // Skip the first argument
        let algo = match args.next() {
            Some(algo) => match algo.as_str() {
                "-r" | "--rapid" => Algorithm::RapidNJ,
                "-n" | "--naive" => Algorithm::Naive,
                "-h" | "--hybrid" => Algorithm::Hybrid,
                _ => {
                    return Err(From::from(format!(
                        "Invalid algorithm: {}. \n{}",
                        algo.as_str(),
                        error_msg
                    )))
                }
            },
            None => Algorithm::Hybrid,
        };

        // Parse the number of threads, format -j 2 or -j2
        let threads = match args.next() {
            Some(threads) => match threads.as_str() {
                "-j" => match args.next() {
                    Some(threads) => match threads.parse::<usize>() {
                        Ok(threads) => threads,
                        Err(_) => {
                            return Err(From::from(format!(
                                "Invalid number of threads: {}. \n{}",
                                threads, error_msg
                            )))
                        }
                    },
                    None => {
                        return Err(From::from(format!(
                            "Invalid number of threads: {}. \n{}",
                            threads, error_msg
                        )))
                    }
                },
                _ => {
                    return Err(From::from(format!(
                        "Invalid number of threads: {}. \n{}",
                        threads, error_msg
                    )))
                }
            },
            None => 1,
        };
        Ok(Config { algo, threads })
    }
}

pub fn run(config: Config) {
    //dbg!(&config);
    rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()
        .unwrap();

    let reader = io::stdin().lock();
    let d = DistanceMatrix::build_from_phylip(reader).unwrap_or_else(|err| {
        eprintln!("{err}");
        process::exit(1);
    });

    // Let match the algorithm
    //dbg!(&d);
    let d = match config.algo {
        Algorithm::Naive => naive_neighbor_joining(d),
        Algorithm::RapidNJ => rapid_nj(d),
        Algorithm::Hybrid => {
            let n = d.size();
            neighbor_joining(d, n - n / 5)
        }
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
