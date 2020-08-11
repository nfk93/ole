use ff::{Field, PrimeField, PrimeFieldDecodingError, PrimeFieldRepr};
use crate::ole::OleField;

// Prime q = 152137607412117916810699707336809121793, with bit size 127
// q = (138 * 54697345034152330060240659727 * 20155392) + 1
// factors of q-1: 2^11 * 3^10 * 23 * 54697345034152330060240659727
// generator of multiplicative group: 5
#[derive(PrimeField)]
#[PrimeFieldModulus = "152137607412117916810699707336809121793"]
#[PrimeFieldGenerator = "5"]
pub struct Fp(FpRepr);

impl OleField for Fp {
    const A: usize = 1024;
    const B: usize = 19683;

    const alpha_gen: Self = Fp(FpRepr([0xc06ef38a81bc942a, 0x0a90ee143aa2da39]));
    const beta_gen: Self = Fp(FpRepr([0x7eaba55f322cb079, 0x5bdbc004d5e45ace]));

    fn interpolate_alpha(ys: &[Self]) -> &[Self] {
        unimplemented!();
    }
    fn interpolate_beta(ys: &[Self]) -> &[Self] {
        unimplemented!();
    }
    fn multipoint_evaluate_alpha(a: &[Self]) -> &[Self] {
        unimplemented!();
    }
    fn multipoint_evaluate_beta(a: &[Self]) -> &[Self] {
        unimplemented!();
    }
}

// Slow reference implementation, use fft3_in_place
pub fn fft3(coeffs: &[Fp], beta: &Fp) -> Vec<Fp> {
    let l = coeffs.len();
    if l == 1 {
        return coeffs.to_vec();
    }

    let a_coeffs: Vec<Fp> = coeffs.iter().step_by(3).map(|x| *x).collect();
    let b_coeffs: Vec<Fp> = coeffs.iter().skip(1).step_by(3).map(|x| *x).collect();
    let c_coeffs: Vec<Fp> = coeffs.iter().skip(2).step_by(3).map(|x| *x).collect();

    let mut beta3 = beta.clone();
    beta3.square();
    beta3.mul_assign(&beta);

    let a_vals = fft3(&a_coeffs, &beta3);
    let b_vals = fft3(&b_coeffs, &beta3);
    let c_vals = fft3(&c_coeffs, &beta3);

    let mut result = vec![Fp::zero(); l]; // could use unsafe unitialized arrays
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

pub fn fft3_inverse(points: &mut [Fp], beta: &Fp) {
    let len_inv = Fp::from_str("152129878020896450676634234608701779969").unwrap();
    let beta_inv = beta.inverse().unwrap();
    fft3_in_place(points, &beta_inv);
    points.iter_mut().for_each(|coeff| coeff.mul_assign(&len_inv));
}

// Slow reference implementation, use fft2_in_place
pub fn fft2(a_coeffs: &[Fp], alpha: &Fp) -> Vec<Fp> {
    let l = a_coeffs.len();

    if l == 1 {
        return a_coeffs.to_vec();
    }

    // split A into B and C such that A(x) = B(x^2) + x * C(x^2)
    let b_coeffs: Vec<Fp> = a_coeffs.iter().step_by(2).map(|x| *x).collect();
    let c_coeffs: Vec<Fp> = a_coeffs.iter().skip(1).step_by(2).map(|x| *x).collect();

    // apply recursively
    let mut alpha_2 = alpha.clone();
    alpha_2.square();
    let b_values = fft2(&b_coeffs, &alpha_2);
    let c_values = fft2(&c_coeffs, &alpha_2);

    // combine subresults
    let mut a_values = vec![Fp::zero(); l];
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

pub fn fft2_inverse(points: &mut [Fp], alpha: &Fp) {
    let len_inv = Fp::from_str("151989035529879520407564258403863019135").unwrap();
    let alpha_inv = alpha.inverse().unwrap();
    fft2_in_place(points, &alpha_inv);
    points.iter_mut().for_each(|coeff| coeff.mul_assign(&len_inv));
}

// only works when coeffs.len() = 2^k for some k, and alpha is a 2^k'th rooth of unity
pub fn fft2_in_place(coeffs: &mut [Fp], alpha: &Fp) {
    digit_reverse_swap(coeffs, 2);
    let n = coeffs.len();
    let mut distance = 1usize;
    let mut iter = 0usize;
    while distance < n {
        let mut factor = Fp::one();
        let mut factor_multiplier = alpha.pow([(n/distance/2) as u64]);
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
        iter += 1;
    }
}

pub fn fft3_in_place(coeffs: &mut [Fp], beta: &Fp) {
    digit_reverse_swap(coeffs, 3);
    let n = coeffs.len();

    let mut beta_1 = beta.pow([(n/3) as u64]);
    let mut beta_2 = beta_1;
    beta_2.square();

    let mut distance = 1usize;
    let mut iter = 0usize;
    while distance < n {
        let mut factor = Fp::one();
        let mut factor_multiplier = beta.pow([(n/distance/3) as u64]);
        for k in 0..distance {
            let mut factor2 = factor;
            factor2.square();
            for j in (0..n).step_by(3*distance) {
                let mut x = coeffs[j+k];
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
        iter += 1;
    }
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

fn horner(coefficients: &[Fp], variable: &Fp) -> Fp {
    coefficients.iter().rev().fold(Fp::zero(), |acc, coeff| {
        let mut a = acc.clone();
        a.mul_assign(variable);
        a.add_assign(coeff);
        return a
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand;

    // const P_MINUS_ONE_FACTORS: [u64; 4] = [2, 3, 23, 54697345034152330060240659727];

    #[test]
    fn test_horner() {
        let coeffs = [
            Fp::from_repr(FpRepr::from(17)).unwrap(),
            Fp::from_repr(FpRepr::from(12)).unwrap(),
            Fp::from_repr(FpRepr::from(19)).unwrap()];
        assert_eq!(Fp::from_str("288944").unwrap(), horner(&coeffs, &Fp::from_str("123").unwrap()))
    }


    #[test]
    fn test_alpha() {
        assert_eq!(Fp::from_str("47236001212243800275399077178102119783").unwrap(), Fp::alpha_gen);
        assert_eq!(Fp::one(), Fp::alpha_gen.pow(&[2u64.pow(10)]));
        assert_eq!( Fp::from_str("152137607412117916810699707336809121792").unwrap(),
                    Fp::alpha_gen.pow(&[2u64.pow(9)]));
    }

    #[test]
    fn test_beta() {
        assert_eq!(Fp::from_str("68906203517960419541774993486852188831").unwrap(), Fp::beta_gen);
        assert_eq!(Fp::one(), Fp::beta_gen.pow(&[3u64.pow(9)]));
        assert_eq!( Fp::from_str("15255265600908596712800720452291744384").unwrap(),
                    Fp::beta_gen.pow(&[3u64.pow(8)]));
    }

    #[test]
    fn test_fft2() {
        let mut rng = rand::thread_rng();
        let coeffs: Vec<Fp> = (0..2048).map(|_| Fp::random(&mut rng)).collect();
        let points = fft2(&coeffs, &Fp::alpha_gen);
        for i in 0..1024 {
            assert_eq!(points[i], horner(&coeffs, &Fp::alpha_gen.pow(&[i as u64])));
        }
    }

    #[test]
    fn test_fft2_in_place() {
        let mut rng = rand::thread_rng();
        let coeffs: Vec<Fp> = (0..1024).map(|_| Fp::random(&mut rng)).collect();
        let mut points = coeffs.to_vec();
        fft2_in_place(&mut points, &Fp::alpha_gen);
        for i in 0..1024 {
            let actual = points[i];
            let expected = horner(&coeffs, &Fp::alpha_gen.pow(&[i as u64]));
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
        fft2_inverse(&mut coeffs, &Fp::alpha_gen);
        points.iter().enumerate().for_each(|(idx, point)| {
            assert_eq!(*point, horner(&coeffs, &Fp::alpha_gen.pow(&[idx as u64])));
        })
    }

    #[test]
    fn test_fft3() {
        let mut rng = rand::thread_rng();
        let coeffs: Vec<Fp> = (0..(Fp::B)).map(|_| Fp::random(&mut rng)).collect();
        let points = fft3(&coeffs, &Fp::beta_gen);
        for i in 0..100 {
            assert_eq!(points[i], horner(&coeffs, &Fp::beta_gen.pow(&[i as u64])));
        }
    }

    #[test]
    fn test_fft3_in_place() {
        let mut rng = rand::thread_rng();
        let coeffs: Vec<Fp> = (0..(Fp::B)).map(|_| Fp::random(&mut rng)).collect();
        let mut points = coeffs.to_vec();
        fft3_in_place(&mut points, &Fp::beta_gen);
        for i in 0..1024 {
            let actual = points[i];
            let expected = horner(&coeffs, &Fp::beta_gen.pow(&[i as u64]));
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
        fft3_inverse(&mut coeffs, &Fp::beta_gen);
        points.iter().take(100).enumerate().for_each(|(idx, point)| {
            assert_eq!(*point, horner(&coeffs, &Fp::beta_gen.pow(&[idx as u64])));
        })
    }

    // #[test]
    // fn test_bit_reverse_copy() {
    //     let mut rng = rand::thread_rng();
    //     let a: Vec<Fp> = (0..256).map(|_| Fp::random(&mut rng)).collect();
    //     let bit_reversed = bit_reverse_copy_u16(&a, 8);
    //     assert_eq!(a[0], bit_reversed[0]);
    //     assert_eq!(a[146], bit_reversed[73]);
    //     assert_eq!(a[160], bit_reversed[5]);
    // }

    // #[test]
    // fn test_reverse_bits() {
    //     assert_eq!(809u16, reverse_bits_u16(595u16, 10));
    //     assert_eq!(8183u16, reverse_bits_u16(7679u16, 13));
    //     assert_eq!(7679u16, reverse_bits_u16(8183u16, 13));
    //     assert_eq!(11u16, reverse_bits_u16(52u16, 6));
    // }

    #[test]
    fn test_digit_reverse_swap() {
        let mut rng = rand::thread_rng();
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
}
