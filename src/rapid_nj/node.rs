#[derive(Debug, Clone)]
pub struct Node {
    pub index: usize,
    pub value: f64,
}

impl Node {
    pub fn new(index: usize, value: f64) -> Self {
        Self { index, value }
    }
}
// Define comparison for Node, so only the value is compared
impl PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}
impl Eq for Node {}

impl PartialOrd for Node {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Node {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.value.partial_cmp(&other.value).unwrap()
    }
}
