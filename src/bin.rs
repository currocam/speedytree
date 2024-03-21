extern crate speedytree;
use clap::Parser;
/// # speedytree
/// `speedytree` is a command line tool for quickly creating a directory tree.
/// It is a Rust implementation of the `tree` command line tool.
/// It is intended to be a drop-in replacement for the `tree` command.
/// It is not intended to be a complete implementation of the `tree` command.
/// It is intended to be a fast implementation of the `tree` command.
use hybrid_nj::neighbor_joining;
use speedytree::hybrid_nj;

use speedytree::distances::DistanceMatrix;
use speedytree::naive_nj::canonical_neighbor_joining;
use speedytree::newick::to_newick;
use speedytree::rapid_nj::rapid_nj;
use std::{
    error,
    io::{self, Write},
    process,
};
type ResultBox<T> = std::result::Result<T, Box<dyn error::Error>>;

/// Define the configuration of the program
/// It contains the algorithm to use and the number of threads to use
///
#[derive(Debug)]
pub struct Config {
    pub(crate) algo: Algorithm,
    pub(crate) threads: usize,
    pub(crate) chunk_size: usize,
    pub(crate) naive_percentage: usize,
}

impl Config {
    /// Build the configuration from the command line arguments
    pub fn build(args: Args) -> ResultBox<Config> {
        // Let match the algorithm, if not specified, use Naive
        let algo = if args.naive {
            Algorithm::Naive
        } else if args.rapidnj {
            Algorithm::RapidNJ
        } else {
            Algorithm::Hybrid
        };
        let cores = args.cores;
        let chunk_size = args.chunk_size;
        // If chunk size is 0, error
        if chunk_size == 0 {
            return Err("Chunk size cannot be 0".into());
        }
        let naive_percentage = args.naive_percentage;
        // If naive percentage is 0, error
        if naive_percentage == 0 {
            return Err("Naive percentage cannot be 0".into());
        }
        // If naive percentage is 100, error
        if naive_percentage == 100 {
            return Err("Naive percentage cannot be 100".into());
        }
        Ok(Config {
            algo,
            threads: cores,
            chunk_size,
            naive_percentage,
        })
    }
}

/// Define the command line arguments
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    /// Use the rapidnj heuristic
    #[arg(long, conflicts_with = "hybrid", conflicts_with = "naive")]
    rapidnj: bool,
    /// Use the naive algorithm
    #[arg(long, conflicts_with = "hybrid", conflicts_with = "rapidnj")]
    naive: bool,
    /// Use the hybrid heuristic
    #[arg(long, conflicts_with = "rapidnj", conflicts_with = "naive")]
    hybrid: bool,
    /// Number of cores to use
    /// Default: 1
    #[arg(short, long, default_value = "1")]
    cores: usize,
    /// Chunk size to be handled by each thread
    /// Default: 30
    #[arg(long, default_value = "30", conflicts_with = "naive")]
    chunk_size: usize,
    /// Percentage of the matrix to be handled by the naive algorithm
    /// Default: 90
    #[arg(
        long,
        default_value = "90",
        conflicts_with = "naive",
        conflicts_with = "rapidnj"
    )]
    naive_percentage: usize,
}

/// Available algorithms in the program
#[derive(Debug, Clone)]
pub enum Algorithm {
    /// Naive neighbor joining
    Naive,
    /// Rapid neighbor joining
    RapidNJ,
    /// Hybrid neighbor joining
    Hybrid,
}
/// Main function of the crate
pub fn run(config: Config) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()
        .unwrap();

    let reader = io::stdin().lock();
    let d = DistanceMatrix::read_from_phylip(reader).unwrap_or_else(|err| {
        eprintln!("{err}");
        process::exit(1);
    });

    let d = match config.algo {
        Algorithm::Naive => canonical_neighbor_joining(d),
        Algorithm::RapidNJ => rapid_nj(d, config.chunk_size),
        Algorithm::Hybrid => {
            let naive_steps = d.size() * config.naive_percentage / 100;
            dbg!(naive_steps);
            neighbor_joining(d, naive_steps, config.chunk_size)
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

fn main() {
    let args = Args::parse();
    //dbg!(&args);
    let config = Config::build(args).unwrap_or_else(|err| {
        eprintln!("Problem parsing arguments: {err}");
        process::exit(1);
    });
    run(config);
}
