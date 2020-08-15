use crate::field::OleField;
use crate::poly::{euclid_division, poly_from_roots};
use rand::{seq::IteratorRandom, CryptoRng, Rng};

pub fn decode_reed_solomon<F: OleField>(points: &mut [F], pos: &[usize]) -> Vec<F> {
    F::fft3_inverse(points, &F::beta());
    let roots: Vec<F> = pos.iter().map(|idx| F::beta().pow([*idx as u64])).collect();
    let b = poly_from_roots(&roots);
    let (_, r) = euclid_division(points, &b);
    return r;
}

pub fn encode_reed_solomon<F: OleField, Crng: CryptoRng + Rng>(
    x: &[F],
    rng: &mut Crng,
) -> (Vec<F>, Vec<F>, Vec<usize>) {
    let pos = pick_indices(F::A, F::B, rng);

    // let mut x_padded = pad_every_other(x, rng);
    // for (i, x_) in x.iter().enumerate() {
    //     assert_eq!(*x_, x_padded[2*i]);
    // }
    let mut x_padded = x.to_vec();
    x_padded.resize_with(F::A / 2, || F::random(rng));
    F::fft2_inverse(&mut x_padded, &F::alpha().pow([2]));
    let x_poly = x_padded.to_vec();
    x_padded.resize_with(F::B, F::zero);
    F::fft3(&mut x_padded, &F::beta());

    let pos_len = pos.len();
    let mut j = 0;
    for (i, x) in x_padded.iter_mut().enumerate() {
        if j < pos_len && i == pos[j] {
            j += 1;
            continue;
        } else {
            *x = F::random(rng);
        }
    }

    return (x_padded, x_poly, pos);
}

pub fn pad_every_other<F: OleField, Crng: CryptoRng + Rng>(input: &[F], rng: &mut Crng) -> Vec<F> {
    let result: Vec<F> = (0..input.len() * 2)
        .map(|i| {
            if i % 2 == 0 {
                input[i / 2]
            } else {
                F::random(rng)
            }
        })
        .collect();
    return result;
}

pub fn pick_indices<Crng: CryptoRng + Rng>(l: usize, n: usize, rng: &mut Crng) -> Vec<usize> {
    let mut pos = (0..n).choose_multiple(rng, l);
    pos.sort();
    return pos;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::Fp;
    use crate::poly;
    use ff::Field;

    #[test]
    fn test_pad_every_other() {
        let mut rng = rand::thread_rng();
        let points: Vec<Fp> = (0..100).map(|_| Fp::random(&mut rng)).collect();

        let padded = pad_every_other(&points, &mut rng);
        for (i, point) in points.iter().enumerate() {
            assert_eq!(*point, padded[i * 2]);
        }
    }

    #[test]
    fn test_encode_decode() {
        let mut rng = rand::thread_rng();
        let points: Vec<Fp> = (0..(Fp::A / 2) - 1).map(|_| Fp::random(&mut rng)).collect();

        let (mut encoded, poly, pos) = encode_reed_solomon(&points, &mut rng);
        for (i, p) in points.iter().enumerate() {
            assert_eq!(*p, poly::horner(&poly, &Fp::alpha().pow([(2 * i) as u64])));
        }
        assert!(encoded.len() == Fp::B);
        assert!(poly.len() == Fp::A / 2);
        assert!(pos.len() == Fp::A);

        let mut decoded_poly = decode_reed_solomon(&mut encoded, &pos);
        println!("len decoded: {}", decoded_poly.len());
        Fp::fft2(&mut decoded_poly, &Fp::alpha().pow([2]));

        for (idx, p) in points.iter().enumerate() {
            assert_eq!(*p, decoded_poly[idx]);
        }
    }

    #[test]
    fn test_encode_decode_linear_operations() {
        let mut rng = rand::thread_rng();
        let points: Vec<Fp> = (0..Fp::A / 2).map(|_| Fp::random(&mut rng)).collect();

        let (mut encoded, _poly, pos) = encode_reed_solomon(&points, &mut rng);

        let mut a: Vec<Fp> = (0..Fp::A / 2).map(|_| Fp::random(&mut rng)).collect();
        let a_copy = a.to_vec();
        Fp::fft2_inverse(&mut a, &Fp::alpha().pow([2]));

        let b_copy: Vec<Fp> = (0..Fp::A / 2).map(|_| Fp::random(&mut rng)).collect();
        let mut b = pad_every_other(&b_copy, &mut rng);
        Fp::fft2_inverse(&mut b, &Fp::alpha());

        for i in 0..Fp::A / 2 {
            assert_eq!(
                a_copy[i],
                poly::horner(&a, &Fp::alpha().pow([(2 * i) as u64]))
            );
            assert_eq!(
                b_copy[i],
                poly::horner(&b, &Fp::alpha().pow([(2 * i) as u64]))
            );
        }

        a.resize_with(Fp::B, Fp::zero);
        Fp::fft3(&mut a, &Fp::beta());
        b.resize_with(Fp::B, Fp::zero);
        Fp::fft3(&mut b, &Fp::beta());

        // Modify the encoding as in the protocol
        let mut j = 0;
        for (idx, x) in encoded.iter_mut().enumerate() {
            x.mul_assign(&a[idx]);
            x.add_assign(&b[idx]);
            if j < pos.len() && idx != pos[j] {
                x.add_assign(&Fp::random(&mut rng));
            } else {
                j += 1;
            }
        }

        let mut decoded_poly = decode_reed_solomon(&mut encoded, &pos);
        // let decoded_copy = decoded_poly.to_vec();
        println!("len decoded: {}", decoded_poly.len());
        Fp::fft2(&mut decoded_poly, &Fp::alpha());

        for (idx, p) in points.iter().enumerate() {
            let mut expected = *p;
            expected.mul_assign(&a_copy[idx]);
            expected.add_assign(&b_copy[idx]);
            assert_eq!(expected, decoded_poly[2 * idx]);
        }
    }
}
