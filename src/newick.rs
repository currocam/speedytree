use fixedbitset::FixedBitSet;
use petgraph::stable_graph::NodeIndex;

use crate::Tree;

fn format_edge_float<'a>(
    t: &Tree,
    node: NodeIndex,
    parent: NodeIndex,
    buffer: &'a mut dtoa::Buffer,
) -> &'a str {
    let e: petgraph::stable_graph::EdgeIndex = t.find_edge(node, parent).unwrap();
    // Stop using format! for performance, create a string directly
    buffer.format_finite(*t.edge_weight(e).unwrap())
}

pub fn to_newick(t: &Tree) -> String {
    // Find a node with 3 children
    let mut buffer = dtoa::Buffer::new();
    let root = root(t).unwrap();
    let mut visited = FixedBitSet::with_capacity(t.node_count());
    fn inner(
        t: &Tree,
        visited: &mut FixedBitSet,
        node: NodeIndex,
        parent: NodeIndex,
        buffer: &mut dtoa::Buffer,
    ) -> String {
        let mut newick = String::new();
        visited.insert(node.index());
        // If leaf
        if !t[node].is_empty() {
            // Stop using format! for performance, create a string directly
            let mut output = t[node].to_owned();
            output.push(':');
            output.push_str(format_edge_float(t, node, parent, buffer));
            return output;
        }
        // If internal node
        let mut children: Vec<NodeIndex> = t.neighbors(node).collect();
        // Filter out visited nodes
        children.retain(|child| !visited.contains(child.index()));
        children.reverse();
        let n_children = children.len();
        newick.push('(');
        for (index, child) in children.into_iter().enumerate() {
            if !visited.contains(child.index()) {
                let child_newick = inner(t, visited, child, node, buffer);
                newick.push_str(&child_newick);
                if index < n_children - 1 {
                    newick.push(',');
                }
            }
        }
        newick.push(')');

        if t.find_edge(node, parent).is_some() {
            newick.push(':');
            newick.push_str(format_edge_float(t, node, parent, buffer));
        }
        newick
    }
    let mut output = inner(t, &mut visited, root, root, &mut buffer);
    output.push(';');
    output
}

fn root(t: &Tree) -> Option<NodeIndex> {
    let mut root = None;
    for node in t.node_indices() {
        if t.neighbors(node).count() == 3 {
            root = Some(node);
            break;
        }
    }
    root
}

// If the node has children, create a sub-tree representation

#[cfg(test)]
mod tests {
    use petgraph::{graph, stable_graph::NodeIndex};

    use super::*;
    #[test]
    fn test_non_binary_tree() {
        let graph = &mut petgraph::Graph::<String, f64, petgraph::Undirected>::new_undirected();
        let a = graph.add_node("A".to_string());
        let b = graph.add_node("B".to_string());
        let c = graph.add_node("C".to_string());
        let internal = graph.add_node("".to_string());
        // Extend the graph
        graph.extend_with_edges(&[
            (a, internal, 100.1),
            (b, internal, 20.2),
            (c, internal, 77757.1),
        ]);
        let newick = to_newick(graph);
        assert_eq!(newick, "(A:100.1,B:20.2,C:77757.1);");
    }
    #[test]
    fn test_binary_tree() {
        let graph = &mut petgraph::Graph::<String, f64, petgraph::Undirected>::new_undirected();
        let a = graph.add_node("A".to_string());
        let b = graph.add_node("B".to_string());
        let c = graph.add_node("C".to_string());
        let d = graph.add_node("D".to_string());
        let internal_v: NodeIndex = graph.add_node("".to_string());
        let internal_u: NodeIndex = graph.add_node("".to_string());
        // Extend the graph
        graph.extend_with_edges(&[
            (a, internal_u, 0.1),
            (b, internal_u, 0.2),
            (c, internal_v, 0.3),
            (d, internal_v, 0.4),
            (internal_u, internal_v, 0.5),
        ]);
        let newick = to_newick(graph);
        assert_eq!(newick, "(C:0.3,D:0.4,(A:0.1,B:0.2):0.5);");
    }
}
