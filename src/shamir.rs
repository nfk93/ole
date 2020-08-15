use crate::field::OleField;
use crate::poly::{euclid_division, poly_from_roots};
use rand;

pub fn share<F: OleField>(secret: &F, n: u64, rho: u64, omega: &F) -> Vec<F> {
    let mut rng = rand::thread_rng();
    let mut coeffs: Vec<F> = (0..rho).map(|_| F::random(&mut rng)).collect();
    coeffs.resize_with(n as usize, F::zero);
    coeffs[0] = *secret;
    F::fft3(&mut coeffs, omega);
    return coeffs;
}

pub fn reconstruct<F: OleField>(indices: &[usize], shares: &[F], n: u64, rho: u64, omega: &F) -> F {
    assert_eq!(indices.len(), shares.len());
    let mut i = 0usize;
    let mut points_with_error: Vec<F> = (0..n)
        .map(|j| {
            if (i < rho as usize) && (j == (indices[i] as u64)) {
                i += 1;
                shares[i - 1]
            } else {
                F::one()
            }
        })
        .collect();
    let roots: Vec<F> = indices.iter().map(|idx| omega.pow([*idx as u64])).collect();
    F::fft3_inverse(&mut points_with_error, omega);
    let b = poly_from_roots(&roots);
    let (_, r) = euclid_division(&points_with_error, &b);
    return r[0];
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::{Fp, OleField};
    use ff::Field;
    use rand::seq::IteratorRandom;

    #[test]
    fn test_share_reconstruct() {
        let mut rng = rand::thread_rng();
        let secret = Fp::random(&mut rng);
        let n = 3u64.pow(7);
        let rho = n - 2u64.pow(8);
        let omega = Fp::beta();

        let shares = share(&secret, n, rho, &omega);
        let mut indices: Vec<usize> = (0..(n as usize)).choose_multiple(&mut rng, rho as usize);
        indices.sort();
        let myshares: Vec<Fp> = indices.iter().map(|i| shares[*i]).collect();
        let reconstructed = reconstruct(&indices, &myshares, n, rho, &omega);

        // println!("secret: {:?}\nreconstr: {:?}", secret, reconstructed);

        assert_eq!(reconstructed, secret)
    }
}
