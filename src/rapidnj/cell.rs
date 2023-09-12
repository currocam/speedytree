use std::cmp::Ordering;

#[derive(Debug, Clone, Copy)]
pub struct NJCell {
    pub index: usize,
    pub value: f64,
}

impl NJCell {
    pub fn new(index: usize, value: f64) -> Self {
        Self { index, value }
    }
}

impl PartialEq for NJCell {
    fn eq(&self, other: &Self) -> bool {
        self.index == other.index
    }
}

impl Eq for NJCell {}

impl PartialOrd for NJCell {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match self.value.partial_cmp(&other.value) {
            Some(Ordering::Equal) => self.index.partial_cmp(&other.index),
            result => result,
        }
    }
}

impl Ord for NJCell {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_cell_definition() {
        let cell1 = NJCell::new(0, 1.0);
        let cell2 = NJCell::new(0, 1.0);
        let cell3 = NJCell::new(2, 2.0);
        let cell4 = NJCell::new(1, 3.0);
        assert_eq!(cell1, cell2);
        assert!(cell4 > cell3)
    }
}
