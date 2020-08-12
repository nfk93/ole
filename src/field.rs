use ff::{Field, PrimeField, PrimeFieldDecodingError, PrimeFieldRepr};

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
    const A: usize;
    const B: usize;
    // pub const G: usize;

    const ALPHA: Self;
    const BETA: Self;
    // const gammas: [Self; Self::G];

    // Interpolates the degree A-1 polynomial through points (xs[i], ys[i]) for i=0,...,A-1
    // Returns the coefficient of the polynomial a_0 + a_1*x + ... + a_{A-1}*x^{A-1} in the form
    // interpolate(ys)[i] = a_i
    fn interpolate_alpha(ys: &[Self]) -> &[Self];

    // Interpolates the degree B-1 polynomial through points (xs[i], ys[i]) for i=0,...,B-1
    // Returns the coefficient of the polynomial a_0 + a_1*x + ... + a_{B-1}*x^{B-1} in the form
    // interpolate(ys)[i] = a_i
    fn interpolate_beta(ys: &[Self]) -> &[Self];

    // // Interpolates the degree G-1 polynomial through points (xs[i], ys[i]) for i=0,...,G-1
    // // Returns the coefficient of the polynomial a_0 + a_1*x + ... + a_{G-1}*x^{G-1} in the form
    // // interpolate(ys)[i] = a_i
    // fn interpolate_gamma(ys: &[Self]) -> [Self; Self::G];

    // Evaluates the polynomial a[0] + a[1]*x + ... + a[K-1]*x^{K-1} at the points alphas
    // Evaluate(a)[i] is the polynomial evaluated in alphas[i]
    fn multipoint_evaluate_alpha(a: &[Self]) -> &[Self];

    // Evaluates the polynomial a[0] + a[1]*x + ... + a[K-1]*x^{K-1} at the points betas
    // Evaluate(a)[i] is the polynomial evaluated in betas[i]
    fn multipoint_evaluate_beta(a: &[Self]) -> &[Self];

    // // Evaluates the polynomial a[0] + a[1]*x + ... + a[K-1]*x^{K-1} at the points gammas
    // // Evaluate(a)[i] is the polynomial evaluated in gammas[i]
    // fn multipoint_evaluate_gamma(a: &[Self]) -> [Self; Self::G];
}

impl OleField for Fp {
    const A: usize = 1024;
    const B: usize = 19683;

    const ALPHA: Self = Fp(FpRepr([0xc06ef38a81bc942a, 0x0a90ee143aa2da39]));
    const BETA: Self = Fp(FpRepr([0x7eaba55f322cb079, 0x5bdbc004d5e45ace]));

    fn interpolate_alpha(_ys: &[Self]) -> &[Self] {
        unimplemented!();
    }
    fn interpolate_beta(_ys: &[Self]) -> &[Self] {
        unimplemented!();
    }
    fn multipoint_evaluate_alpha(_a: &[Self]) -> &[Self] {
        unimplemented!();
    }
    fn multipoint_evaluate_beta(_a: &[Self]) -> &[Self] {
        unimplemented!();
    }
}
