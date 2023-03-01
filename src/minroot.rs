//! Implementation of MinRoot hash
//!

use std::marker::PhantomData;

use ark_ff::{BigInteger256, PrimeField};

pub trait MinRootParam: PrimeField {
    const ALPHA: u64;
    // 1/alpha mod p-1
    const ALPHA_INV: <Self as PrimeField>::BigInt;
}

#[derive(Clone, Debug)]
pub struct MinRootHasher<F: MinRootParam> {
    pub(crate) vec_x: Vec<F>,
    pub(crate) vec_y: Vec<F>,
    pub(crate) vec_indexer: Vec<F>,
    _phantom: PhantomData<F>,
}

impl<F: MinRootParam> MinRootHasher<F> {
    pub fn new() -> Self {
        Self {
            _phantom: PhantomData,
            vec_x: vec![],
            vec_y: vec![],
            vec_indexer: vec![],
        }
    }

    fn iterate_once(&mut self, cur_x: &F, cur_y: &F, indexer: usize) -> (F, F) {
        let indexer = F::from(indexer as u64);
        let next_x = (*cur_x + *cur_y).pow(F::ALPHA_INV);
        let next_y = *cur_x + indexer;

        self.vec_x.push(next_x);
        self.vec_y.push(next_y);
        self.vec_indexer.push(indexer);

        (next_x, next_y)
    }

    pub fn hash(&mut self, x0: &F, y0: &F, iteration: usize) -> (F, F) {
        let mut cur_x = *x0;
        let mut cur_y = *y0;

        self.vec_x.push(cur_x);
        self.vec_y.push(cur_y);

        for indexer in 1..=iteration {
            (cur_x, cur_y) = self.iterate_once(&cur_x, &cur_y, indexer)
        }
        (cur_x, cur_y)
    }

    pub fn display(&self) {
        for (i, (indexer, (x, y))) in self
            .vec_indexer
            .iter()
            .zip(self.vec_x.iter().zip(self.vec_y.iter()))
            .enumerate()
        {
            println!("{}-th iter: indexer {}; x {}; y {}", i, indexer, x, y)
        }
    }
}

impl MinRootParam for ark_bn254::Fr {
    const ALPHA: u64 = 5;
    // 1/alpha mod p-1
    const ALPHA_INV: <Self as PrimeField>::BigInt = BigInteger256([
        0xcfe7f7a98ccccccd,
        0x535cb9d394945a0d,
        0x93736af8679aad17,
        0x26b6a528b427b354,
    ]);
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bn254::Fr;
    use ark_ff::{Field, UniformRand};
    use ark_std::test_rng;

    #[test]
    fn test_alpha_inv() {
        let mut rng = test_rng();
        for _ in 0..10 {
            let tmp = Fr::rand(&mut rng);
            assert_eq!(tmp.pow(&[Fr::ALPHA]).pow(Fr::ALPHA_INV), tmp);
        }
    }

    #[test]
    fn test_min_root_hash() {
        let x = Fr::from(1);
        let y = Fr::from(2);

        let mut hasher = MinRootHasher::<Fr>::new();
        let _res = hasher.hash(&x, &y, 10);

        hasher.display()
    }
}
