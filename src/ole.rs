use crate::error::{OleError};
use ff::Field;
// use ocelot::ot::Receiver as OTReceiver;
// use ocelot::ot::Sender as OTSender;
use rand::{CryptoRng, Rng};
use scuttlebutt::channel::AbstractChannel;

pub trait Sender<F: OleField>
where
    Self: Sized,
{
    // builds a pool of OTs via OT-extension, to be used for OLE
    fn init<C: AbstractChannel, Crng: CryptoRng + Rng>(
        alpha: &[F],
        beta: &[F],
        channel: &mut C,
        rng: &mut Crng,
    ) -> Result<Self, OleError>;

    fn input<Crng: CryptoRng + Rng>(self, a: &[F], b: &[F], rng: &mut Crng) -> Result<(), OleError> {

        // check that input has the correct size
        let l = a.len();
        if l > F::A {
            return Err(OleError::TODOError)
        }

        let mut padded_a = a.to_vec();
        for i in l..F::A {
            padded_a.push(F::random(rng));
        }

        let mut padded_b = b.to_vec();
        for i in l..F::A {
            padded_b.push(F::random(rng));
        }

        let a_poly = F::interpolate_alpha(&padded_a);
        let b_poly = F::interpolate_alpha(&padded_b);

        Ok(())
    }
}

pub trait Receiver<F: OleField>
where
    Self: Sized,
{
    // builds a pool of OTs via OT-extension, to be used for OLE
    fn init<C: AbstractChannel, Crng: CryptoRng + Rng>(
        alpha: &[F],
        beta: &[F],
        channel: &mut C,
        rng: &mut Crng,
    ) -> Result<Self, OleError>;

    fn input(self, x: &[F]) -> Result<Vec<F>, OleError>;
}

// A field with special properties that enable interpolation and evaluation at a set of
// predetermined alphas and betas in the field using FFT.
pub trait OleField: Field {
    const A: usize;
    const B: usize;
    // pub const G: usize;

    const alpha_gen: Self;
    const beta_gen: Self;
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
