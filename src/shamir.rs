use crate::field::OleField;
use rand;
use crate::poly::{poly_from_roots, euclid_division};

pub fn share<F: OleField>(secret: &F, n: u64, rho: u64, omega: &F) -> Vec<F>{
    let mut rng = rand::thread_rng();
    let mut coeffs: Vec<F> = (0..rho).map(|_| F::random(&mut rng)).collect();
    coeffs.resize_with(n as usize, F::zero);
    coeffs[0] = *secret;
    F::fft3(&mut coeffs, omega);
    return coeffs
}

pub fn reconstruct<F: OleField>(indices: &[u64], shares: &[F], n: u64, rho: u64, omega: &F) -> F {
    assert_eq!(indices.len(), shares.len());
    let mut i = 0usize;
    let mut points_with_error: Vec<F> = (0..n).map(|j| {
        if (i < rho as usize) && (j == indices[i]) {
            i += 1;
            shares[i-1]
        } else {
            F::one()
        }
    }).collect();
    let roots: Vec<F> = indices.iter().map(|idx| omega.pow([*idx as u64])).collect();
    F::fft3_inverse(&mut points_with_error, omega);
    let b = poly_from_roots(&roots);
    let (_, mut r) = euclid_division(&points_with_error, &b);
    return r[0]
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::seq::IteratorRandom;
    use crate::field::{OleField, Fp};
    use ff::Field;

    #[test]
    fn test_share_reconstruct() {
        let mut rng = rand::thread_rng();
        let secret = Fp::random(&mut rng);
        let n = 3u64.pow(7);
        let rho = n - 2u64.pow(8);
        let omega = Fp::beta().pow([9]);

        let shares = share(&secret, n, rho, &omega);
        let mut indices: Vec<u64> = (0..n).choose_multiple(&mut rng, rho as usize);
        indices.sort();
        let mut myshares: Vec<Fp> = indices.iter().map(|i| shares[*i as usize]).collect();
        let reconstructed = reconstruct(&indices, &myshares, n, rho, &omega);

        // println!("secret: {:?}\nreconstr: {:?}", secret, reconstructed);

        assert_eq!(reconstructed, secret)
    }
}
