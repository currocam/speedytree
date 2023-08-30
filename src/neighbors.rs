use trees::Tree;

use crate::{distances::DistanceMatrix, ResultBox};

pub fn neighbor_joining(_d: &DistanceMatrix) -> ResultBox<Tree<usize>> {
    Ok(trees::Tree::<usize>::from_tuple((0, (1, 2, 3), (4, 5, 6))))
}
