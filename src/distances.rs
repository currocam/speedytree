use ndarray::Array2;
use seq_io::fasta::{OwnedRecord, Reader, Record};
use std::io::Stdin;

use crate::{ResultBox, M};

#[derive(Debug)]
pub struct DistanceMatrix {
    pub names: Vec<String>,
    matrix: M,
}

pub fn distance_matrix_from_stdin(input: Stdin) -> ResultBox<DistanceMatrix> {
    let mut reader = Reader::new(input);
    let records: Result<Vec<_>, _> = reader.records().collect();
    let records = records?;
    let n = records.len();

    let mut names = Vec::with_capacity(n);
    for record in &records {
        let id = record.id()?;
        names.push(id.to_owned());
    }

    let mut matrix = Array2::<f64>::zeros((n, n));
    for i in 0..n {
        for j in i..n {
            let d = distance(records[i].seq(), records[j].seq());
            matrix[[i, j]] = d;
            matrix[[j, i]] = d;
        }
    }
    Ok(DistanceMatrix { names, matrix })
}

fn distance(record_i: &[u8], record_j: &[u8]) -> f64 {
    let mut cost = 0;
    for (char_i, char_j) in record_i.iter().zip(record_j.iter()) {
        if *char_i == b'-' || *char_j == b'-' {
            continue;
        }
        if char_i != char_j {
            cost += 1;
        }
    }
    cost as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_distance_empty() {
        let (record_one, record_two) = (vec![], vec![]);
        assert_eq!(distance(&record_one, &record_two), 0.0);
        let (record_one, record_two) = (vec![], vec![b'A']);
        assert_eq!(distance(&record_one, &record_two), 0.0);
    }
    #[test]
    fn test_happy_path() {
        let (record_one, record_two) = (vec![b'A', b'A', b'-'], vec![b'A', b'A', b'C']);
        assert_eq!(distance(&record_one, &record_two), 0.0);
        let (record_one, record_two) = (vec![b'A', b'C', b'-'], vec![b'A', b'A', b'C']);
        assert_eq!(distance(&record_one, &record_two), 1.0);
    }
    #[test]
    fn test_aminoacid() {
        let (record_one, record_two) = (vec![b'M', b'M', b'L'], vec![b'L', b'L', b'M']);
        assert_eq!(distance(&record_one, &record_two), 3.0);
    }
}
