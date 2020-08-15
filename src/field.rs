use crate::fft;
use ff::{Field, PrimeField, PrimeFieldDecodingError, PrimeFieldRepr};
use scuttlebutt::Block;

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
pub trait OleField: PrimeField + Field + From<Block> {
    fn alpha() -> Self; // generator of order A multiplicative subgroup
    const A: usize; // order of alpha()

    fn beta() -> Self; // generator of order B multiplicative subgroup
    const B: usize; // order of beta()

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

    fn to_block(self) -> Block;
}

impl OleField for Fp {
    const A: usize = 256usize; //2usize.pow(8);
    const B: usize = 2187; //3usize.pow(7);

    fn alpha() -> Self {
        Fp(FpRepr([0xc06ef38a81bc942a, 0x0a90ee143aa2da39])).pow([4])
    }

    fn beta() -> Self {
        Fp(FpRepr([0x7eaba55f322cb079, 0x5bdbc004d5e45ace])).pow([9])
    }

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

    fn to_block(self) -> Block {
        unsafe { std::mem::transmute((self.0).0) }
    }
}

impl From<Block> for Fp {
    #[inline]
    fn from(m: Block) -> Fp {
        unsafe { Fp(FpRepr(*(&m as *const _ as *const [u64; 2]))) }
    }
}

// impl From<Fp> for Block {
//     #[inline]
//     fn from(Fp(FpRepr(data)): Fp) -> Block {
//         unsafe { std::mem::transmute(data) }
//     }
// }

#[cfg(test)]
mod tests {
    use super::*;
    use rand::thread_rng;

    #[test]
    fn test_fp_block_conversion() {
        let mut rng = thread_rng();
        for _ in 0..1000 {
            let v = Fp::random(&mut rng);
            let v_ = v.clone();
            let block: Block = v.to_block();
            let v_from_block: Fp = block.into();

            assert_eq!(v_, v_from_block);
        }
    }
}
