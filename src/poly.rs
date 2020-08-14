use ff::Field;

pub fn horner<F: Field>(coefficients: &[F], variable: &F) -> F {
    coefficients.iter().rev().fold(F::zero(), |acc, coeff| {
        let mut a = acc.clone();
        a.mul_assign(variable);
        a.add_assign(coeff);
        return a;
    })
}

// input: polynomials a and b in coefficient representation
// outputs: (q, r) such that q*b + r = a, with deg(r) < deg(b)
// precondition: leading coefficient of a and b (i.e. a[a.len() -1]) is nonzero
pub fn euclid_division<F: Field>(a: &[F], b: &[F]) -> (Vec<F>, Vec<F>) {
    let mut deg_r = a.len() - 1;
    let mut q = vec![F::zero(); deg_r + 1];
    let mut r = a.to_vec();
    let d = b.len() - 1; // degree of b
    let c = b[d].inverse().unwrap(); // c = lc(b)^-1
    while deg_r >= d {
        let mut s = r[deg_r]; // s = lc(r)
        s.mul_assign(&c); // s = lc(r)/lc(b)
        q[deg_r - d].add_assign(&s);
        for i in 0..(d + 1) {
            let mut sb_i = s;
            sb_i.mul_assign(&b[i]);
            r[deg_r - d + i].sub_assign(&sb_i);
        }
        'update_deg_r: loop {
            if r[deg_r].is_zero() {
                deg_r -= 1;
            } else {
                r.truncate(deg_r + 1);
                break 'update_deg_r;
            }
        }
    }
    return (q, r);
}

#[cfg(test)]
fn naive_poly_mult<F: Field>(a: &[F], b: &[F]) -> Vec<F> {
    let mut prod = vec![F::zero(); a.len() + b.len() - 1];

    for i in 0..a.len() {
        for j in 0..b.len() {
            let mut ab = a[i];
            ab.mul_assign(&b[j]);
            prod[i + j].add_assign(&ab);
        }
    }

    return prod;
}

// adds polynomial b to polynomial a.
// precondition deg(a) > deg(b)
fn poly_add<F: Field>(a: &mut [F], b: &[F]) {
    for (idx, coeff) in b.iter().enumerate() {
        a[idx].add_assign(coeff);
    }
}

pub fn poly_from_roots<F: Field>(a: &[F]) -> Vec<F> {
    let mut result = vec![F::zero(); a.len()];
    for i in 0..result.len() {
        let mut a_neg = a[i];
        a_neg.negate();
        result[i] = F::one();
        for j in (1..i + 1).rev() {
            result[j].mul_assign(&a_neg);
            let rj_1 = result[j - 1];
            result[j].add_assign(&rj_1);
        }
        result[0].mul_assign(&a_neg)
    }
    result.push(F::one());
    return result;
}

// Performs lagrangian interpolation at zero
// Precondition: xs.len() == ys.len()
pub fn lagrangian_interpolation<F: Field>(xs: &[F], ys: &[F]) -> F {
    let mut result = F::zero();
    for i in 0..xs.len() {
        // let mut xi_neg = xs[i];
        // xi_neg.negate();
        let mut numer = F::one();
        let mut denom = F::one();
        for (j, xj) in xs.iter().enumerate() {
            if j != i {
                numer.mul_assign(xj);
                let mut xj_minus_xi = *xj;
                // xj_minus_xi.add_assign(&xi_neg);
                xj_minus_xi.sub_assign(&xs[i]);
                denom.mul_assign(&xj_minus_xi);
            }
        }
        let mut yi_lagrange = ys[i];
        yi_lagrange.mul_assign(&numer);
        yi_lagrange.mul_assign(&denom.inverse().unwrap());
        result.add_assign(&yi_lagrange);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fft::fft2;
    use crate::field::{Fp, FpRepr, OleField};
    use ff::PrimeField;
    use rand;
    use rand::seq::IteratorRandom;

    #[test]
    fn test_lagrangian_interpolation() {
        let mut rng = rand::thread_rng();
        let mut coeffs: Vec<Fp> = (0..Fp::A / 2).map(|_| Fp::random(&mut rng)).collect();
        coeffs.resize_with(Fp::A, Fp::zero);

        let mut ys = coeffs.to_vec();
        Fp::fft2(&mut ys, &Fp::alpha());

        let xs = (0..Fp::A).choose_multiple(&mut rng, Fp::A / 2);
        let ys_: Vec<Fp> = xs.iter().map(|i| ys[*i]).collect();
        let xs_: Vec<Fp> = xs.iter().map(|i| Fp::alpha().pow([*i as u64])).collect();

        let result = lagrangian_interpolation(&xs_, &ys_);
        assert_eq!(result, coeffs[0]);
    }

    #[test]
    fn test_poly_add() {
        let mut rng = rand::thread_rng();
        let a: Vec<Fp> = (0..101).map(|_| Fp::random(&mut rng)).collect();
        let b: Vec<Fp> = (0..22).map(|_| Fp::random(&mut rng)).collect();
        let mut c = a.clone();
        poly_add(&mut c, &b);
        for _ in 0..1000 {
            let point = Fp::random(&mut rng);
            let mut a_ = horner(&a, &point);
            let b_ = horner(&b, &point);
            let c_ = horner(&c, &point);

            a_.add_assign(&b_);
            assert_eq!(a_, c_);
        }
    }

    #[test]
    fn test_naive_poly_mult() {
        let mut rng = rand::thread_rng();
        let a: Vec<Fp> = (0..101).map(|_| Fp::random(&mut rng)).collect();
        let b: Vec<Fp> = (0..22).map(|_| Fp::random(&mut rng)).collect();
        let prod = naive_poly_mult(&a, &b);
        for _ in 0..1000 {
            let point = Fp::random(&mut rng);
            let mut a_ = horner(&a, &point);
            let b_ = horner(&b, &point);
            let c_ = horner(&prod, &point);

            a_.mul_assign(&b_);
            assert_eq!(a_, c_);
        }
    }

    #[test]
    fn test_euclid_division() {
        let mut rng = rand::thread_rng();
        let a: Vec<Fp> = (0..83).map(|_| Fp::random(&mut rng)).collect();
        let b: Vec<Fp> = (0..22).map(|_| Fp::random(&mut rng)).collect();
        let (q, r) = euclid_division(&a, &b);
        let mut actual = naive_poly_mult(&q, &b);
        println!("{:?}", q);
        println!("{:?}", r);
        poly_add(&mut actual, &r);

        for (idx, coeff) in a.iter().enumerate() {
            assert_eq!(actual[idx], *coeff);
        }
    }

    #[test]
    fn test_horner() {
        let coeffs = [
            Fp::from_repr(FpRepr::from(17)).unwrap(),
            Fp::from_repr(FpRepr::from(12)).unwrap(),
            Fp::from_repr(FpRepr::from(19)).unwrap(),
        ];
        assert_eq!(
            Fp::from_str("288944").unwrap(),
            horner(&coeffs, &Fp::from_str("123").unwrap())
        )
    }

    #[test]
    fn test_poly_from_roots() {
        let mut rng = rand::thread_rng();
        let roots: Vec<Fp> = (0..100).map(|_| Fp::random(&mut rng)).collect();
        let coeffs = poly_from_roots(&roots);
        for root in roots {
            assert_eq!(horner(&coeffs, &root), Fp::zero());
        }

        let roots = [
            Fp::from_repr(FpRepr::from(1231)).unwrap(),
            Fp::from_repr(FpRepr::from(2)).unwrap(),
            Fp::from_repr(FpRepr::from(17)).unwrap(),
        ];
        let coeffs = poly_from_roots(&roots);
        assert_eq!(
            horner(&coeffs, &Fp::from_repr(FpRepr::from(1111)).unwrap()),
            Fp::from_str("152137607412117916810699707336663532273").unwrap()
        );
        assert_eq!(
            horner(&coeffs, &Fp::from_str("213131").unwrap()),
            Fp::from_str("9624661948301400").unwrap()
        );
    }
}
