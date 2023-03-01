use ark_ff::{Field, PrimeField};

use crate::minroot::MinRootParam;

fn fold<F: MinRootParam>(witness_1: &[F], error_1: &[F], witness_2: &[F], error_2: &[F]) {
    let num_constraints = witness_1.len();
    assert_eq!(num_constraints, error_1.len());
    assert_eq!(num_constraints, error_2.len());
    assert_eq!(num_constraints, witness_2.len());
    for i in 0..num_constraints - 2 {
        assert_eq!(
            witness_1[i] + witness_1[i + 1] - witness_1[i + 2].pow(&[5]) + error_1[i],
            F::zero()
        );
        assert_eq!(
            witness_2[i] + witness_2[i + 1] - witness_2[i + 2].pow(&[5]) + error_2[i],
            F::zero()
        );
    }

    let randomizer = F::from(2u64);
    let mut witness_res = vec![];
    let mut error_res = vec![];
    for i in 0..num_constraints - 2 {
        witness_res.push(randomizer * witness_1[i] + witness_2[i]);
        let delta = compute_delta(&(randomizer * witness_1[i + 2]), &witness_2[i + 2]);
        error_res.push(error_1[i] + error_2[i] + delta);
    }

    for i in 0..num_constraints - 2 {
        println!(
            "{}-th constraint: {}",
            i,
            witness_res[i] + witness_res[i + 1] - witness_res[i + 2].pow(&[5]) + error_res[i],
        );

        // assert_eq!(
        //     witness_res[i] + witness_res[i + 1] - witness_res[i + 2].pow(&[5]) + error_res[i],
        //     F::zero()
        // );
    }
}

// TODO: improve this code
#[inline]
fn compute_delta<F: PrimeField>(w1: &F, w2: &F) -> F {
    let five = F::from(5u64);
    let ten = F::from(10u64);

    five * w1.pow(&[4]) * w2
        + ten * w1.pow(&[3]) * w2.pow(&[2])
        + ten * w1.pow(&[2]) * w2.pow(&[3])
        + five * w1 * w2.pow(&[4])
}

#[cfg(test)]
mod test {
    use crate::minroot::MinRootHasher;

    use super::*;
    use ark_bn254::Fr;
    use ark_ff::{UniformRand, Zero};
    use ark_poly::{Polynomial, Radix2EvaluationDomain};
    use ark_std::{rand::RngCore, test_rng};
    #[test]
    fn test_fold() {
        let mut rng = test_rng();
        for _ in 0..10 {
            // let iter = rng.next_u32() % 100 + 10;
            let iter = 10;
            let x1 = Fr::rand(&mut rng);
            let y1 = Fr::rand(&mut rng);
            let x2 = Fr::rand(&mut rng);
            let y2 = Fr::rand(&mut rng);

            let mut hasher1 = MinRootHasher::<Fr>::new();
            let mut hasher2 = MinRootHasher::<Fr>::new();
            let _res = hasher1.hash(&x1, &y1, iter as usize);
            let _res = hasher2.hash(&x2, &y2, iter as usize);
            fold(
                &hasher1.vec_x[..iter],
                &hasher1.vec_indexer[..iter],
                &hasher2.vec_x[..iter],
                &hasher2.vec_indexer[..iter],
            )
        }
    }
}
