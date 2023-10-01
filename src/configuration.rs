use super::Algorithm;
use super::ResultBox;
/// Define the configuration of the program
/// It contains the algorithm to use and the number of threads to use
///
#[derive(Debug)]
pub struct Config {
    pub(crate) algo: Algorithm,
    pub(crate) threads: usize,
}

impl Config {
    /// Build the configuration from the command line arguments
    /// # Examples
    ///
    /// ```
    /// let args = ["speedytree", "-r", "-j", "2"].iter().map(|s| s.to_string());
    /// let config = speedytree::configuration::Config::build(args).unwrap();
    /// println!("{:?}", config);
    /// ```
    pub fn build(mut args: impl Iterator<Item = String>) -> ResultBox<Config> {
        // Let match the algorithm, if not specified, use Naive
        let error_msg =
            "Usage: nj [-r|--rapid, -n|--naive, -h|--hybrid] [-j 1] < input.phy > output.nwk";
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
