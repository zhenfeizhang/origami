//! IOP component for MinRoot relation.
//! Given a list of vectors $X := x0, \dots, x_N$, generate an IOP that
//! X satisfies the MinRoot relation:
//!  x_i + x_{i+1} - x_{i+2}^\alpha + i + 1 = 0

use ark_ff::PrimeField;
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Radix2EvaluationDomain,
    UVPolynomial,
};

/// Input a MinRoot sequence of length $l$; output a polynomial $h[x]$ where
/// $h[omega^i] = 0$ for $0 <= i < l$
fn compute_polynomial_h<F: PrimeField>(x_i: &[F]) -> DensePolynomial<F> {
    assert!(x_i.len() > 2);
    let domain_size = x_i.len() + 2;
    let domain = Radix2EvaluationDomain::<F>::new(domain_size).unwrap();
    let coset_domain = GeneralEvaluationDomain::<F>::new(domain.size() * 5).unwrap();

    let first = x_i[0];
    let second = x_i[1];

    let mut x_i = x_i.to_vec();
    x_i.insert(0, F::zero());
    x_i.insert(0, F::zero());

    let n = x_i.len();

    // w[X]
    let witness_poly = domain.ifft(&x_i);
    // w[omega X]
    let witness_rotated = {
        x_i.remove(0);
        domain.ifft(&x_i)
    };
    // w[omega^2 X]
    let witness_rotated_twice = {
        x_i.remove(0);
        domain.ifft(&x_i)
    };
    // q(x)
    // q(x) = ifft(1, 2, 3, 4...)
    let mut q_i = (1..=n).map(|i| F::from(i as u64)).collect::<Vec<F>>();
    // front-pad the first two elements of q(x) with
    // x[0]^5 and x[1]^5 - x[0]
    // to ensure the first two rows are also 0
    q_i.insert(0, second.pow(&[5]) - first);
    q_i.insert(0, first.pow(&[5]));
    let selector_poly = domain.ifft(&q_i);

    // h(x) = w(x) + w(omega x) -w(omega^2 x)^\alpha + q(x)
    let h = {
        let mut evals = vec![];
        let witness_poly_coset = coset_domain.coset_fft(&witness_poly);
        let witness_poly_rotated_coset = coset_domain.coset_fft(&witness_rotated);
        let witness_poly_rotated_twice_coset = coset_domain.coset_fft(&witness_rotated_twice);
        let selector_poly_coset = coset_domain.coset_fft(&selector_poly);
        for i in 0..coset_domain.size() {
            evals.push(
                witness_poly_coset[i] + witness_poly_rotated_coset[i]
                    - witness_poly_rotated_twice_coset[i].pow(&[5])
                    + selector_poly_coset[i],
            )
        }
        coset_domain.coset_ifft(&evals)
    };

    DensePolynomial::from_coefficients_vec(h)
}

#[cfg(test)]
mod test {
    use crate::minroot::MinRootHasher;

    use super::*;
    use ark_bn254::Fr;
    use ark_ff::{UniformRand, Zero};
    use ark_poly::Polynomial;
    use ark_std::{rand::RngCore, test_rng};
    #[test]
    fn test_iop() {
        let mut rng = test_rng();
        for _ in 0..10 {
            let iter = rng.next_u32() % 100 + 10;
            let x = Fr::rand(&mut rng);
            let y = Fr::rand(&mut rng);

            let mut hasher = MinRootHasher::<Fr>::new();
            let _res = hasher.hash(&x, &y, iter as usize);

            let domain_size = iter + 2;
            let domain = Radix2EvaluationDomain::<Fr>::new(domain_size as usize).unwrap();

            let h = compute_polynomial_h(&hasher.vec_x);
            for i in 0..iter {
                assert_eq!(h.evaluate(&domain.element(i as usize)), Fr::zero())
            }
        }
    }
}
