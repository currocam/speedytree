mod algorithm;
mod data;
pub use algorithm::neighbor_joining;

#[cfg(test)]
mod tests {
    use crate::distances::DistanceMatrix;

    use super::*;
    #[test]
    fn test_example_wikipedia() {
        let d = DistanceMatrix::wikipedia_example();
        let phylo = neighbor_joining(d, 4, 1);
        assert!(phylo.is_ok());
        let tree = phylo.unwrap();
        let mut node_indices = tree.node_indices();
        let (a, b, c, d, e) = (
            node_indices.next().unwrap(),
            node_indices.next().unwrap(),
            node_indices.next().unwrap(),
            node_indices.next().unwrap(),
            node_indices.next().unwrap(),
        );
        // Get internal nodes
        let u = tree.neighbors(a).next().unwrap();
        let w = tree.neighbors(e).next().unwrap();
        let v = tree.neighbors(c).next().unwrap();
        // Create array of node from, node to and distance
        let expected = [
            (a, u, 2.0),
            (b, u, 3.0),
            (e, w, 1.0),
            (d, w, 2.0),
            (c, v, 4.0),
            (v, u, 3.0),
            (v, w, 2.0),
        ];
        for (a, b, dist) in expected.into_iter() {
            assert_eq!(tree.edge_weight(tree.find_edge(a, b).unwrap()), Some(&dist));
        }
    }
}
