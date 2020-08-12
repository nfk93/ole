use crate::error::{OleError};
use crate::field::*;

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

    fn input<Crng: CryptoRng + Rng>(self, _a: &[F], _b: &[F], _rng: &mut Crng) -> Result<(), OleError> {

        // // check that input has the correct size
        // let l = a.len();
        // if l > F::A {
        //     return Err(OleError::TODOError)
        // }
        //
        // let mut padded_a = a.to_vec();
        // for i in l..F::A {
        //     padded_a.push(F::random(rng));
        // }
        //
        // let mut padded_b = b.to_vec();
        // for i in l..F::A {
        //     padded_b.push(F::random(rng));
        // }
        //
        // let a_poly = F::interpolate_alpha(&padded_a);
        // let b_poly = F::interpolate_alpha(&padded_b);

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
