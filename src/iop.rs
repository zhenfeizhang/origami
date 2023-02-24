//! IOP component for MinRoot relation.
//! Given a list of vectors $X := x0, \dots, x_N$, generate an IOP that
//! X satisfies the MinRoot relation:
//!  x_i + x_{i+1} - x_{i+2}^\alpha + i + 1 = 0

use ark_ff::PrimeField;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Radix2EvaluationDomain, univariate::DensePolynomial, UVPolynomial, Polynomial};

fn prover_work<F: PrimeField>(x_i: &[F]) {
    let n = x_i.len();
    let domain = Radix2EvaluationDomain::<F>::new(16).unwrap();
    let coset_domain = GeneralEvaluationDomain::<F>::new(128).unwrap();

    let witness_poly = domain.ifft(x_i);

    let witness_rotated = {
        let mut x_i = x_i.to_vec();
        x_i.insert(0, F::zero());

        for (i, e) in x_i.iter().enumerate() {
            println!("x_{}: {}", i, e);
        }
        domain.ifft(&x_i)
    };

    let witness_rotated_twice = {
        let mut x_i = x_i.to_vec();
        x_i.insert(0, F::zero());
        x_i.insert(0, F::zero());

        for (i, e) in x_i.iter().enumerate() {
            println!("x_{}: {}", i, e);
        }
        domain.ifft(&x_i)
    };

    // q(x) = ifft(1,2,3,4...)
    let q_i = (1..=n).map(|i| F::from(i as u64)).collect::<Vec<F>>();
    let selector_poly = domain.ifft(&q_i);

    // h(x) = w(x) + w(omega x) -w(omega^2 x)^\alpha + q(x)
    let h = {
        let mut evals = vec![];
        let witness_poly_coset = coset_domain.coset_fft(&witness_poly);
        let witness_poly_rotated_coset = coset_domain.coset_fft(&witness_rotated);
        let witness_poly_rotated_twice_coset = coset_domain.coset_fft(&witness_rotated_twice);
        let selector_poly_coset = coset_domain.coset_fft(&selector_poly);
        for i in 0..coset_domain.size(){
            println!("{} {} {} {}", i,  witness_poly_coset[i], witness_poly_rotated_coset[i], witness_poly_rotated_twice_coset[i]);
        }
        for i in 0..coset_domain.size() {
            evals.push(
                witness_poly_coset[i] + witness_poly_rotated_coset[i]
                    - witness_poly_rotated_twice_coset[i].pow(&[5])
                    //+ selector_poly_coset[i],
            )
        }
        coset_domain.coset_ifft(&evals)
    };

    for (i, e) in h.iter().enumerate() {
        println!("h_{}: {}", i, e);
    }

    let h_poly = DensePolynomial::from_coefficients_vec(h);

    for i in 0..16{
        println!("h_{}: {}", i, h_poly.evaluate(&domain.element(i)));
    }

    // let res = h_poly.evaluate_over_domain(coset_domain);
    
    // for (i, e) in res.evals.iter().enumerate() {
    //     println!("res_{}: {}", i, e);
    // }

    
    

}
#[cfg(test)]
mod test {
    use crate::minroot::MinRootHasher;

    use super::*;
    use ark_bn254::Fr;
    #[test]
    fn test_iop() {
        let x = Fr::from(1);
        let y = Fr::from(2);
        let mut hasher = MinRootHasher::<Fr>::new();
        let _res = hasher.hash(&x, &y, 10);
        #[cfg(feature = "inspection")]
        {
            hasher.display();
            prover_work( &hasher.vec_x)
        }

    }
}
