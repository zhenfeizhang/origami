//! Implementation of MinRoot hash
//!

use std::marker::PhantomData;

use ark_ff::PrimeField;

#[derive(Clone, Debug)]
pub struct MinRootHasher<F: PrimeField> {
    pub(crate) alpha: u64,
    #[cfg(feature = "inspection")]
    pub(crate) vec_x: Vec<F>,
    #[cfg(feature = "inspection")]
    pub(crate) vec_y: Vec<F>,
    _phantom: PhantomData<F>,
}

impl<F: PrimeField> MinRootHasher<F> {
    pub fn new() -> Self {
        Self {
            alpha: 5,
            _phantom: PhantomData,
            #[cfg(feature = "inspection")]
            vec_x: vec![],
            #[cfg(feature = "inspection")]
            vec_y: vec![],
        }
    }

    fn iterate_once(&mut self, cur_x: &F, cur_y: &F, indexer: usize) -> (F, F) {
        let next_x = (*cur_x + *cur_y).pow(&[self.alpha]);
        let next_y = *cur_x + F::from(indexer as u64);
        #[cfg(feature = "inspection")]
        {
            self.vec_x.push(next_x);
            self.vec_y.push(next_y)
        }
        (next_x, next_y)
    }

    pub fn hash(&mut self, x0: &F, y0: &F, iteration: usize) -> (F, F) {
        let mut cur_x = *x0;
        let mut cur_y = *y0;
        #[cfg(feature = "inspection")]
        {
            self.vec_x.push(cur_x);
            self.vec_y.push(cur_y)
        }

        for indexer in 1..=iteration {
            (cur_x, cur_y) = self.iterate_once(&cur_x, &cur_y, indexer)
        }
        (cur_x, cur_y)
    }

    #[cfg(feature = "inspection")]
    pub fn display(&self) {
        for (i, (x, y)) in self.vec_x.iter().zip(self.vec_y.iter()).enumerate() {
            println!("{}-th iter: x {} y {}", i, x, y)
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bn254::Fr;
    #[test]
    fn test_min_root_hash() {
        let x = Fr::from(1);
        let y = Fr::from(2);
        let mut hasher = MinRootHasher::<Fr>::new();
        let _res = hasher.hash(&x, &y, 10);
        #[cfg(feature = "inspection")]
        hasher.display()
    }
}
