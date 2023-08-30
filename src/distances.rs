use ndarray::Array2;
use seq_io::fasta::{Reader, Record};
use std::io::Stdin;

use crate::{ResultBox, M};

#[derive(Debug)]
pub struct DistanceMatrix {
    pub names: Vec<String>,
    matrix: M,
}

pub fn distance_matrix_from_stdin(input: Stdin) -> ResultBox<DistanceMatrix> {
    let mut reader = Reader::new(input);
    let mut n = 0;
    let mut names = Vec::new();
    while let Some(record) = reader.next() {
        let record = record?;
        names.push(record.id()?.to_owned());
        n += 1;
    }
    let matrix = Array2::<f64>::zeros((n, n));
    Ok(DistanceMatrix { names, matrix })
}
