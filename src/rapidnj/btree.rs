use std::collections::BTreeSet;

use super::cell::NJCell;
pub struct NJRow {
    pub val: BTreeSet<NJCell>,
}

impl NJRow {
    pub fn new() -> Self {
        let val = BTreeSet::new();
        Self { val }
    }

    fn first(&self) -> Option<&NJCell> {
        self.val.first()
    }

    pub fn from(value: Vec<f64>, offset: usize) -> Self {
        let iter = value
            .iter()
            .enumerate()
            .map(|(i, v)| NJCell::new(i + offset, *v));
        let val = BTreeSet::from_iter(iter);
        NJRow { val }
    }

    pub fn replace(&mut self, old: NJCell, new: NJCell) {
        self.val.remove(&old);
        self.val.insert(new);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_binary_tree_set_operations() {
        let row = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let mut njrow = NJRow::from(row.clone(), 0);
        let row_rev = row.iter().rev().cloned().collect::<Vec<f64>>();
        let mut njrow_rev = NJRow::from(row_rev.clone(), 0);

        assert_eq!(njrow.first().unwrap().index, 0);
        assert_eq!(njrow.first().unwrap().value, 0.0);
        assert_eq!(njrow_rev.first().unwrap().index, 5);
        assert_eq!(njrow_rev.first().unwrap().value, 0.0);

        // Test replace
        njrow.replace(NJCell::new(0, 0.0), NJCell::new(0, 10.0));
        njrow_rev.replace(NJCell::new(5, 0.0), NJCell::new(5, 10.0));
        assert_eq!(njrow.first().unwrap().index, 1);
        assert_eq!(njrow.first().unwrap().value, 1.0);
        assert_eq!(njrow.val.last().unwrap().index, 0);

        assert_eq!(njrow_rev.first().unwrap().index, 4);
        assert_eq!(njrow_rev.first().unwrap().value, 1.0);
    }
}
