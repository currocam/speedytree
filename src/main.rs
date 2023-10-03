use clap::Parser;
/// # speedytree
/// `speedytree` is a command line tool for quickly creating a directory tree.
/// It is a Rust implementation of the `tree` command line tool.
/// It is intended to be a drop-in replacement for the `tree` command.
/// It is not intended to be a complete implementation of the `tree` command.
/// It is intended to be a fast implementation of the `tree` command.
use speedytree::{
    configuration::{Args, Config},
    run,
};
use std::process;
fn main() {
    let args = Args::parse();
    //dbg!(&args);
    let config = Config::build(args).unwrap_or_else(|err| {
        eprintln!("Problem parsing arguments: {err}");
        process::exit(1);
    });
    run(config);
}
