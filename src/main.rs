use birc_rapidnj::{run, Config};
use std::process;

fn main() {
    let config = Config::build().unwrap_or_else(|err| {
        eprintln!("Problem parsing arguments: {err}");
        process::exit(1);
    });
    run(config);
}
