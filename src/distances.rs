use rand::seq::SliceRandom;

use crate::ResultBox;
use std::io::{self};
/// Distance matrix struct
#[derive(Debug, Clone)]
pub struct DistanceMatrix {
    /// Distance matrix as a vector of vectors
    pub matrix: Vec<Vec<f64>>,
    /// Names of the sequences
    pub names: Vec<String>,
}

/// Distance matrix from a phylip file
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
    /// Size of the distance matrix
    pub fn size(&self) -> usize {
        self.matrix.len()
    }
    /// Permutate the distance matrix for testing purposes
    pub fn permutate(&mut self) {
        let mut rng = rand::thread_rng();
        let mut perm = (0..self.size()).collect::<Vec<usize>>();
        perm.shuffle(&mut rng);
        let mut new_matrix = vec![vec![0.0; self.size()]; self.size()];
        for i in 0..self.size() {
            for j in 0..self.size() {
                new_matrix[i][j] = self.matrix[perm[i]][perm[j]];
            }
        }
        self.matrix = new_matrix;
        let mut new_names = vec![String::new(); self.size()];
        for i in 0..self.size() {
            new_names[i] = self.names[perm[i]].clone();
        }
        self.names = new_names;
    }
    /// Example from Wikipedia, https://en.wikipedia.org/wiki/Neighbor_joining
    pub fn wikipedia_example() -> DistanceMatrix {
        DistanceMatrix {
            matrix: vec![
                vec![0.0, 5.0, 9.0, 9.0, 8.0],
                vec![5.0, 0.0, 10.0, 10.0, 9.0],
                vec![9.0, 10.0, 0.0, 8.0, 7.0],
                vec![9.0, 10.0, 8.0, 0.0, 3.0],
                vec![8.0, 9.0, 7.0, 3.0, 0.0],
            ],
            names: vec![
                "A".to_string(),
                "B".to_string(),
                "C".to_string(),
                "D".to_string(),
                "E".to_string(),
            ],
        }
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
        let distance_matrix = DistanceMatrix::build_from_phylip::<&[u8]>(input).unwrap();
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
