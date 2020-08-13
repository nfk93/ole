use ff::{Field, PrimeField, PrimeFieldDecodingError, PrimeFieldRepr};
use crate::fft;

// Prime q = 152137607412117916810699707336809121793, with bit size 127
// q = (138 * 54697345034152330060240659727 * 20155392) + 1
// factors of q-1: 2^11 * 3^10 * 23 * 54697345034152330060240659727
// generator of multiplicative group: 5
#[derive(PrimeField)]
#[PrimeFieldModulus = "152137607412117916810699707336809121793"]
#[PrimeFieldGenerator = "5"]
pub struct Fp(FpRepr);

// A field with special properties that enable interpolation and evaluation at a set of
// predetermined alphas and betas in the field using FFT.
pub trait OleField: PrimeField {
    const ALPHA: Self; // generator of order A multiplicative subgroup
    const A: usize; // order of ALPHA

    const BETA: Self; // generator of order B multiplicative subgroup
    const B: usize; // order of BETA

    // In-place radix2 DIT FFT
    // After inputting coeffs with coeffs.len() = n = 2^k we have coeffs[i] = f(alpha^i)
    // where f(x) = coeffs[0] + coeffs[1]*x + ... + coeffs[n-1]*x^(n-1)
    // Precondition: alpha^n = 1
    fn fft2(coeffs: &mut [Self], alpha: &Self);

    // In-place radix3 DIT FFT
    // After inputting coeffs with coeffs.len() = n = 3^k we have coeffs[i] = f(beta^i)
    // where f(x) = coeffs[0] + coeffs[1]*x + ... + coeffs[n-1]*x^(n-1)
    // Precondition: beta^n = 1
    fn fft3(coeffs: &mut [Self], beta: &Self);

    // Inverse radix2 DIT FFT
    // Inverse of fft2
    fn fft2_inverse(ys: &mut [Self], alpha: &Self);

    // Inverse radix3 DIT FFT
    // Inverse of fft3
    fn fft3_inverse(ys: &mut [Self], beta: &Self);
}

impl OleField for Fp {
    const A: usize = 1024;
    const B: usize = 19683;

    const ALPHA: Self = Fp(FpRepr([0xc06ef38a81bc942a, 0x0a90ee143aa2da39]));
    const BETA: Self = Fp(FpRepr([0x7eaba55f322cb079, 0x5bdbc004d5e45ace]));

    fn fft2(coeffs: &mut [Self], alpha: &Self) {
        fft::fft2_in_place(coeffs, alpha);
    }

    fn fft3(coeffs: &mut [Self], beta: &Self) {
        fft::fft3_in_place(coeffs, beta);
    }

    fn fft2_inverse(ys: &mut [Self], alpha: &Self) {
        fft::fft2_inverse(ys, alpha);
    }

    fn fft3_inverse(ys: &mut [Self], beta: &Self) {
        fft::fft3_inverse(ys, beta);
    }
}
