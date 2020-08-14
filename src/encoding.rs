use crate::field::{Fp, OleField};
use crate::fft::{fft3_in_place, fft2_in_place, fft2_inverse, fft3_inverse};
use crate::poly::{poly_from_roots, euclid_division};
use ff::Field;
use rand::{CryptoRng, Rng, seq::IteratorRandom};

pub fn decode_reed_solomon<F: OleField>(points: &mut [F], pos: &[usize]) -> Vec<F> {
    fft3_inverse(points, &F::beta());
    let roots: Vec<F> = pos.iter().map(|idx| F::beta().pow([*idx as u64])).collect();
    let b = poly_from_roots(&roots);
    let (_, mut r) = euclid_division(points, &b);
    fft2_in_place(&mut r, &F::alpha().pow([2u64]));
    return r
}

pub fn encode_reed_solomon<F: OleField, Crng: CryptoRng + Rng>(x: &[F], rng: &mut Crng) -> (Vec<F>, Vec<usize>) {
    let mut alpha2 = F::alpha();
    alpha2.square();
    let mut xs = x.to_vec();
    // let padding = 0..(F::A/2 - xs.len()).map(|_| F::random(&mut rng)).collect();
    xs.resize_with(F::A/2, || F::random(rng));

    fft2_inverse(&mut xs, &alpha2);
    // let padding = vec!(F::zero(); F::B - xs.len());
    xs.resize_with(F::B, F::zero);

    fft3_in_place(&mut xs, &F::beta());

    let mut pos = (0..(F::B as usize)).choose_multiple(rng, F::A/2);
    pos.sort();
    let pos_len = pos.len();
    let mut j = 0;
    for (i, x) in xs.iter_mut().enumerate() {
        if j < pos_len && i == pos[j] {
            j += 1;
            continue;
        }
        else {
            *x = F::random(rng);
        }
    }

    return (xs, pos)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode() {
        let mut rng = rand::thread_rng();
        let points: Vec<Fp> = (0..500).map(|_| Fp::random(&mut rng)).collect();

        let (mut encoded, pos) = encode_reed_solomon(&points, &mut rng);
        let decoded = decode_reed_solomon(&mut encoded, &pos);

        for (idx, p) in points.iter().enumerate() {
            assert_eq!(*p, decoded[idx]);
        }
    }
}
