use ark_ff::{Field, PrimeField};

use crate::minroot::MinRootParam;

fn fold<F: MinRootParam>(
    witness_1: &[F],
    error_1: &[F],
    witness_2: &[F],
    error_2: &[F],
) -> (Vec<F>, Vec<F>) {
    // sanity checks on all the inputs
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

    // for now we use a fixed randomizer for random linear combination
    let randomizer = F::from(2u64);
    // w3 = r * w1 + w2
    let mut witness_res = vec![];
    // e2 = delta + r * e1 + e2 + (r^5-r) w_13^5
    let mut error_res = vec![];

    for i in 0..num_constraints {
        witness_res.push(randomizer * witness_1[i] + witness_2[i]);
    }
    for i in 0..num_constraints - 2 {
        // compute the actual data for the rest rows
        let delta = compute_delta_new(&witness_1[i + 2], &witness_2[i + 2], &randomizer);
        error_res.push(
            randomizer * error_1[i]
                + error_2[i]
                + delta
                + (randomizer.pow(&[5]) - randomizer) * witness_1[i + 2].pow(&[5u64]),
        );
    }
    //  pad the last two rows; the errors here are not used
    error_res.push(F::zero());
    error_res.push(F::zero());

    (witness_res, error_res)
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

// TODO: improve this code
fn compute_delta_new<F: PrimeField>(w1: &F, w2: &F, r: &F) -> F {
    let five = F::from(5u64);
    let ten = F::from(10u64);
    let w1 = *r * *w1;
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
            let iter = rng.next_u32() as usize % 100 + 10;
            // let iter = 10;
            let x1 = Fr::rand(&mut rng);
            let y1 = Fr::rand(&mut rng);
            let x2 = Fr::rand(&mut rng);
            let y2 = Fr::rand(&mut rng);
            let x3 = Fr::rand(&mut rng);
            let y3 = Fr::rand(&mut rng);
            let x4 = Fr::rand(&mut rng);
            let y4 = Fr::rand(&mut rng);

            let mut hasher1 = MinRootHasher::<Fr>::new();
            let mut hasher2 = MinRootHasher::<Fr>::new();
            let mut hasher3 = MinRootHasher::<Fr>::new();
            let mut hasher4 = MinRootHasher::<Fr>::new();
            let _res = hasher1.hash(&x1, &y1, iter as usize);
            let _res = hasher2.hash(&x2, &y2, iter as usize);
            let _res = hasher3.hash(&x3, &y3, iter as usize);
            let _res = hasher4.hash(&x4, &y4, iter as usize);
            let (x12, e12) = fold(
                &hasher1.vec_x[..iter],
                &hasher1.vec_indexer[..iter],
                &hasher2.vec_x[..iter],
                &hasher2.vec_indexer[..iter],
            );
            assert!(MinRootHasher::<Fr>::check(&x12, &e12));
            let (x34, e34) = fold(
                &hasher1.vec_x[..iter],
                &hasher1.vec_indexer[..iter],
                &hasher2.vec_x[..iter],
                &hasher2.vec_indexer[..iter],
            );
            assert!(MinRootHasher::<Fr>::check(&x34, &e34));
            let (x1234, e1234) = fold(&x12, &e12, &x34, &e34);
            assert!(MinRootHasher::<Fr>::check(&x1234, &e1234));
        }
    }
}
