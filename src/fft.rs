use ff::{PrimeField, Field};

// Slow reference implementation, use fft2_in_place
pub fn fft2<F: PrimeField>(a_coeffs: &[F], alpha: &F) -> Vec<F> {
    let l = a_coeffs.len();

    if l == 1 {
        return a_coeffs.to_vec();
    }

    // split A into B and C such that A(x) = B(x^2) + x * C(x^2)
    let b_coeffs: Vec<F> = a_coeffs.iter().step_by(2).map(|x| *x).collect();
    let c_coeffs: Vec<F> = a_coeffs.iter().skip(1).step_by(2).map(|x| *x).collect();

    // apply recursively
    let mut alpha_2 = alpha.clone();
    alpha_2.square();
    let b_values = fft2(&b_coeffs, &alpha_2);
    let c_values = fft2(&c_coeffs, &alpha_2);

    // combine subresults
    let mut a_values = vec![F::zero(); l];
    for i in 0..(l/2) {
        let x = alpha.pow(&[i as u64]);
        a_values[i] = c_values[i];
        a_values[i].mul_assign(&x);
        a_values[i].add_assign(&b_values[i]);

        let j = i + l/2;
        let x = alpha.pow(&[j as u64]);
        a_values[j] = c_values[i];
        a_values[j].mul_assign(&x);
        a_values[j].add_assign(&b_values[i]);
    }

    return a_values
}

// Slow reference implementation, use fft3_in_place
pub fn fft3<F: PrimeField>(coeffs: &[F], beta: &F) -> Vec<F> {
    let l = coeffs.len();
    if l == 1 {
        return coeffs.to_vec();
    }

    let a_coeffs: Vec<F> = coeffs.iter().step_by(3).map(|x| *x).collect();
    let b_coeffs: Vec<F> = coeffs.iter().skip(1).step_by(3).map(|x| *x).collect();
    let c_coeffs: Vec<F> = coeffs.iter().skip(2).step_by(3).map(|x| *x).collect();

    let mut beta3 = beta.clone();
    beta3.square();
    beta3.mul_assign(&beta);

    let a_vals = fft3(&a_coeffs, &beta3);
    let b_vals = fft3(&b_coeffs, &beta3);
    let c_vals = fft3(&c_coeffs, &beta3);

    let mut result = vec![F::zero(); l]; // could use unsafe unitialized arrays
    for i in 0..(l/3) {
        let mut f = |j| {
            let x = beta.pow(&[j as u64]);
            let x2 = x.pow(&[2 as u64]);
            let mut v = c_vals[i];
            v.mul_assign(&x2);
            let mut bx = b_vals[i];
            bx.mul_assign(&x);
            v.add_assign(&bx);
            v.add_assign(&a_vals[i]);
            result[j] = v;
        };
        f(i);
        f(i + l/3);
        f(i + (2*l)/3);
    }

    return result;
}

// only works when coeffs.len() = 2^k for some k, and alpha is a 2^k'th rooth of unity
pub fn fft2_in_place<F: PrimeField>(coeffs: &mut [F], alpha: &F) {
    digit_reverse_swap(coeffs, 2);
    let n = coeffs.len();
    let mut distance = 1usize;
    while distance < n {
        let mut factor = F::one();
        let factor_multiplier = alpha.pow([(n/distance/2) as u64]);
        for k in 0..distance {
            for j in (0..n).step_by(2*distance) {
                let mut x = coeffs[j+k];
                let mut y = coeffs[j+k+distance];
                y.mul_assign(&factor);

                coeffs[j+k].add_assign(&y);
                x.sub_assign(&y);
                coeffs[j+k+distance] = x;
            }
            factor.mul_assign(&factor_multiplier);
        }
        distance <<= 1;
    }
}

pub fn fft3_in_place<F: PrimeField>(coeffs: &mut [F], beta: &F) {
    digit_reverse_swap(coeffs, 3);
    let n = coeffs.len();

    let beta_1 = beta.pow([(n/3) as u64]);
    let mut beta_2 = beta_1;
    beta_2.square();

    let mut distance = 1usize;
    while distance < n {
        let mut factor = F::one();
        let factor_multiplier = beta.pow([(n/distance/3) as u64]);
        for k in 0..distance {
            let mut factor2 = factor;
            factor2.square();
            for j in (0..n).step_by(3*distance) {
                let x = coeffs[j+k];
                let mut y = coeffs[j+k+distance];
                let mut z = coeffs[j+k+2*distance];
                y.mul_assign(&factor);
                z.mul_assign(&factor2);

                coeffs[j+k].add_assign(&y);
                coeffs[j+k].add_assign(&z);

                for i in 1..3 {
                    y.mul_assign(&beta_1);
                    z.mul_assign(&beta_2);
                    coeffs[j+k+i*distance] = x;
                    coeffs[j+k+i*distance].add_assign(&y);
                    coeffs[j+k+i*distance].add_assign(&z);
                }
            }
            factor.mul_assign(&factor_multiplier);
        }
        distance = distance*3;
    }
}

// TODO: dont use from_str hack
pub fn fft3_inverse<F: Field + PrimeField>(points: &mut [F], beta: &F) {
    let len_inv = F::from_str(&points.len().to_string()).unwrap().inverse().unwrap();
    let beta_inv = beta.inverse().unwrap();
    fft3_in_place::<F>(points, &beta_inv);
    points.iter_mut().for_each(|coeff| coeff.mul_assign(&len_inv));
}

// TODO: dont use from_str hack
pub fn fft2_inverse<F: Field + PrimeField>(points: &mut [F], alpha: &F) {
    let len_inv = F::from_str(&points.len().to_string()).unwrap().inverse().unwrap();
    let alpha_inv = alpha.inverse().unwrap();
    fft2_in_place::<F>(points, &alpha_inv);
    points.iter_mut().for_each(|coeff| coeff.mul_assign(&len_inv));
}

// does digit-reversal permutation of the data in the given base
// precondition: data.len() = base^m
pub fn digit_reverse_swap<T: Sized + Copy>(data: &mut[T], base: usize) {
    let n = data.len();
    let (n1, p_odd) = calc_n1(base, n);
    let seed = seed_table(base, n, n1);
    for i in 0..(n1-1) {
        for j in (i+1)..n1 {
            let temp = data[i + seed[j]];
            data[i + seed[j]] = data[seed[i] + j];
            data[seed[i] + j] = temp;
            if p_odd {
                for z in 1..base {
                    let temp = data[i + seed[j] + (z*n1)];
                    data[i + seed[j] + (z*n1)] = data[seed[i] + j + (z*n1)];
                    data[seed[i] + j + (z*n1)] = temp;
                }
            }
        }
    }
}

fn calc_n1(base: usize, n: usize) -> (usize, bool) {
    let mut k = 0usize;
    let mut t = 1usize;
    while t < n {
        k += 1;
        t = t*base;
    }
    if t != n {
        panic!(format!("Cant do digit-reversal reordering. Data size {} is not on the form {}^{}", n, base, k));
    }

    if k.trailing_zeros() == 0 {
        return (base.pow(((k-1)/2) as u32), true)
    } else {
        return (base.pow((k/2) as u32), false)
    }
}

fn seed_table(base: usize, n: usize, n1: usize) -> Vec<usize> {
    let mut r = vec!(0usize; n1);
    r[1] = n/base;
    for i in 2..base {
        r[i] = r[i-1] + r[1];
    }
    for j in 1..(n1/base) {
        r[base*j] = r[j]/base;
        for k in 1..base {
            r[base*j + k] = r[base*j] + r[k];
        }
    }
    return r;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::horner;
    use crate::field::{Fp, OleField};
    use rand;

    #[test]
    fn test_fft2() {
        let mut rng = rand::thread_rng();
        let coeffs: Vec<Fp> = (0..2048).map(|_| Fp::random(&mut rng)).collect();
        let points = fft2(&coeffs, &Fp::alpha());
        for i in 0..1024 {
            assert_eq!(points[i], horner(&coeffs, &Fp::alpha().pow(&[i as u64])));
        }
    }

    #[test]
    fn test_fft2_in_place() {
        let mut rng = rand::thread_rng();
        let coeffs: Vec<Fp> = (0..1024).map(|_| Fp::random(&mut rng)).collect();
        let mut points = coeffs.to_vec();
        fft2_in_place(&mut points, &Fp::alpha());
        for i in 0..1024 {
            let actual = points[i];
            let expected = horner(&coeffs, &Fp::alpha().pow(&[i as u64]));
            assert!(actual == expected,
                    format!("point {} is incorrect\n\
                             \tfound:    {:?}\n\
                             \texpected: {:?}", i, actual, expected));
        }

        let coeffs: Vec<Fp> = (0..512).map(|_| Fp::random(&mut rng)).collect();
        let mut points = coeffs.to_vec();
        fft2_in_place(&mut points, &Fp::alpha().pow([2u64]));
        for i in 0..512 {
            let actual = points[i];
            let expected = horner(&coeffs, &Fp::alpha().pow(&[2*i as u64]));
            assert!(actual == expected,
                    format!("point {} is incorrect\n\
                             \tfound:    {:?}\n\
                             \texpected: {:?}", i, actual, expected));
        }
    }

    #[test]
    fn test_fft2_inverse() {
        let mut rng = rand::thread_rng();
        let points: Vec<Fp> = (0..1024).map(|_| Fp::random(&mut rng)).collect();
        let mut coeffs = points.to_vec();
        fft2_inverse(&mut coeffs, &Fp::alpha());
        points.iter().enumerate().for_each(|(idx, point)| {
            assert_eq!(*point, horner(&coeffs, &Fp::alpha().pow(&[idx as u64])));
        })
    }

    #[test]
    fn test_fft3() {
        let mut rng = rand::thread_rng();
        let coeffs: Vec<Fp> = (0..(Fp::B)).map(|_| Fp::random(&mut rng)).collect();
        let points = fft3(&coeffs, &Fp::beta());
        for i in 0..100 {
            assert_eq!(points[i], horner(&coeffs, &Fp::beta().pow(&[i as u64])));
        }
    }

    #[test]
    fn test_fft3_in_place() {
        let mut rng = rand::thread_rng();
        let coeffs: Vec<Fp> = (0..(Fp::B)).map(|_| Fp::random(&mut rng)).collect();
        let mut points = coeffs.to_vec();
        fft3_in_place(&mut points, &Fp::beta());
        for i in 0..102 {
            let actual = points[i];
            let expected = horner(&coeffs, &Fp::beta().pow(&[i as u64]));
            assert!(actual == expected,
                    format!("point {} is incorrect\n\
                             \tfound:    {:?}\n\
                             \texpected: {:?}", i, actual, expected));
        }

        let coeffs: Vec<Fp> = (0..(Fp::B/9)).map(|_| Fp::random(&mut rng)).collect();
        let mut points = coeffs.to_vec();
        fft3_in_place(&mut points, &Fp::beta().pow([9]));
        for i in 0..102 {
            let actual = points[i];
            let expected = horner(&coeffs, &Fp::beta().pow(&[9*i as u64]));
            assert!(actual == expected,
                    format!("point {} is incorrect\n\
                             \tfound:    {:?}\n\
                             \texpected: {:?}", i, actual, expected));
        }
    }

    #[test]
    fn test_fft3_inverse() {
        let mut rng = rand::thread_rng();
        let points: Vec<Fp> = (0..(Fp::B)).map(|_| Fp::random(&mut rng)).collect();
        let mut coeffs = points.to_vec();
        fft3_inverse(&mut coeffs, &Fp::beta());
        points.iter().take(100).enumerate().for_each(|(idx, point)| {
            assert_eq!(*point, horner(&coeffs, &Fp::beta().pow(&[idx as u64])));
        })
    }

    #[test]
    fn test_digit_reverse_swap() {
        let mut a: Vec<usize> = (0..(3usize.pow(4))).collect();
        digit_reverse_swap(&mut a, 3usize);

        assert_eq!(a[7], 45);
        assert_eq!(a[2], 54);
        assert_eq!(a[57], 11);
        assert_eq!(a[30], 10);

        let mut a: Vec<usize> = (0..(3usize.pow(5))).collect();
        digit_reverse_swap(&mut a, 3usize);

        assert_eq!(a[140], 196);
        assert_eq!(a[9], 9);
        assert_eq!(a[225], 17);
    }

    #[test]
    fn test_seed_table() {
        let n = 3usize.pow(4);
        let n1 = 3usize.pow(2);
        let base = 3usize;

        let table = seed_table(base, n, n1);
        assert_eq!(table[1], 3usize.pow(3));
        assert_eq!(table[7], 45);
        assert_eq!(table[2], 54);
        assert_eq!(table[3], 9);
        assert_eq!(table[6], 18);
        assert_eq!(table[8], 72);
    }

    #[test]
    fn test_fft_fft_inverse_identity() {
        let mut rng = rand::thread_rng();
        let coeffs: Vec<Fp> = (0..(Fp::B)).map(|_| Fp::random(&mut rng)).collect();
        let mut coeffs_clone = coeffs.to_vec();
        fft3_in_place(&mut coeffs_clone, &Fp::beta());
        fft3_inverse(&mut coeffs_clone, &Fp::beta());
        for (c_actual, c_expected) in coeffs_clone.iter().zip(&coeffs) {
            assert_eq!(*c_actual, *c_expected);
        }

        let coeffs: Vec<Fp> = (0..(Fp::A)).map(|_| Fp::random(&mut rng)).collect();
        let mut coeffs_clone = coeffs.to_vec();
        fft2_in_place(&mut coeffs_clone, &Fp::alpha());
        fft2_inverse(&mut coeffs_clone, &Fp::alpha());
        for (c_actual, c_expected) in coeffs_clone.iter().zip(&coeffs) {
            assert_eq!(*c_actual, *c_expected);
        }
    }
}
