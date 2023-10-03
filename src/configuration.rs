use crate::Algorithm;

use super::ResultBox;
use clap::Parser;
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
