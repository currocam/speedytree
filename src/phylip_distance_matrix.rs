use std::io::{self};

use petgraph::algo;

use crate::{phylogenetic_tree::PhyloTree, ResultBox};

#[derive(Debug)]
pub struct DistanceMatrix {
    pub matrix: Vec<Vec<f64>>,
    pub names: Vec<String>,
}

impl From<PhyloTree> for DistanceMatrix {
    fn from(value: PhyloTree) -> Self {
        // Ge keys from hash
        let n_leaves = value.leaves.len();
        let mut matrix = vec![vec![0.0; n_leaves]; n_leaves];
        let mut names = Vec::with_capacity(n_leaves);
        dbg!(&value.leaves);
        for i_leaf in 0..n_leaves {
            let node = value.leaves.get(&i_leaf).unwrap();
            names.push(value.tree[*node].to_owned());
        }
        for i in 0..n_leaves {
            for j in i..n_leaves {
                let node_i = value.leaves.get(&i).unwrap();
                let node_j = value.leaves.get(&j).unwrap();
                let path = algo::astar(
                    &value.tree,
                    *node_i,          // start
                    |n| n == *node_j, // is_goal
                    |e| *e.weight(),  // edge_cost
                    |_| 0.0,          // estimate_cost
                );
                matrix[i][j] = path.unwrap().0;
                matrix[j][i] = matrix[i][j];
            }
        }

        DistanceMatrix { matrix, names }
    }
}

pub fn read_phylip_distance_matrix<R>(mut reader: R) -> ResultBox<DistanceMatrix>
where
    R: io::BufRead,
{
    // Read the first line to get the number of sequences
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
        let distance_matrix = read_phylip_distance_matrix::<&[u8]>(&input[..]).unwrap();
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
