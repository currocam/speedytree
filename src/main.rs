use birc_rapidnj::{run, Config};
use std::process;

fn main() {
    let args = std::env::args();
    //dbg!(&args);
    let config = Config::build(args).unwrap_or_else(|err| {
        eprintln!("Problem parsing arguments: {err}");
        process::exit(1);
    });
    run(config);
}
