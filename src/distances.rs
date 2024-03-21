use crate::ResultBox;
use std::io::{self};
/// Distance matrix data structure
#[derive(Debug, Clone)]
pub struct DistanceMatrix {
    /// Distance matrix
    pub matrix: Vec<Vec<f64>>,
    /// Names of the taxa
    pub names: Vec<String>,
}

/// Distance matrix from a [PHYLIP](https://phylipweb.github.io/phylip/) file
impl DistanceMatrix {
    pub fn read_from_phylip<R>(mut reader: R) -> ResultBox<DistanceMatrix>
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
            names.push(words.next().expect("Valid Phylip format").to_string());
            matrix.push(
                words
                    .map(|s| s.parse::<f64>().expect("Valid Phylip format"))
                    .collect::<Vec<f64>>(),
            );
        }
        Ok(DistanceMatrix { matrix, names })
    }
    /// Size of the distance matrix
    pub fn size(&self) -> usize {
        self.matrix.len()
    }
    pub fn build(matrix: Vec<Vec<f64>>, names: Vec<String>) -> ResultBox<DistanceMatrix> {
        if matrix.len() != names.len() {
            return Err("Matrix and names have different lengths".into());
        }
        Ok(DistanceMatrix { matrix, names })
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
        let distance_matrix = DistanceMatrix::read_from_phylip::<&[u8]>(input).unwrap();
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
