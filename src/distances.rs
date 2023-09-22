use crate::ResultBox;
use std::io::{self};

#[derive(Debug)]
pub struct DistanceMatrix {
    pub matrix: Vec<Vec<f64>>,
    pub names: Vec<String>,
}

impl DistanceMatrix {
    pub fn build_from_phylip<R>(mut reader: R) -> ResultBox<DistanceMatrix>
    where
        R: io::BufRead,
    {
        let mut line = String::new();
        reader.read_line(&mut line)?;
        let n = line.trim().parse::<usize>()?;
        // Read the next n lines to get the names of the sequences (first word), and parse the vector
        let mut names = Vec::with_capacity(n);
        let mut matrix = Vec::with_capacity(n);
        for _ in 0..n {
            line.clear();
            reader.read_line(&mut line)?;
            let mut words = line.split_whitespace();
            names.push(words.next().unwrap().to_string());
            matrix.push(
                words
                    .map(|s| s.parse::<f64>().unwrap())
                    .collect::<Vec<f64>>(),
            );
        }
        Ok(DistanceMatrix { matrix, names })
    }
    pub fn size(&self) -> usize {
        self.matrix.len()
    }
}

// Test
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_distance_matrix_from_stdin() {
        // create test case
        let input = "4
A 0.0 5.0 9.0 9.0
B 5.0 0.0 10.0 10.0
C 9.0 10.0 0.0 8.0
D 9.0 10.0 8.0 0.0
"
        .as_bytes();
        // run function
        let distance_matrix = DistanceMatrix::build_from_phylip::<&[u8]>(&input[..]).unwrap();
        // check result
        assert_eq!(
            distance_matrix.matrix,
            vec![
                vec![0.0, 5.0, 9.0, 9.0],
                vec![5.0, 0.0, 10.0, 10.0],
                vec![9.0, 10.0, 0.0, 8.0],
                vec![9.0, 10.0, 8.0, 0.0],
            ]
        );
        assert_eq!(
            distance_matrix.names,
            vec![
                "A".to_string(),
                "B".to_string(),
                "C".to_string(),
                "D".to_string()
            ]
        );
    }
}
