use trees::Tree;

pub fn newick_format(t: Tree<usize>, n: Vec<String>) {
    println!("{}", t.to_string())
}
