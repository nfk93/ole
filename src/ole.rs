use crate::error::{OleError};
use crate::field::*;
use crate::shamir;
use sha2::{Sha256, Digest};
use ocelot::ot::{Receiver as OTReceiver, Sender as OTSender};
use rand::{CryptoRng, Rng};
use scuttlebutt::channel::AbstractChannel;
use ff::{Field, PrimeField, PrimeFieldRepr};

pub trait Sender<F: OleField>
where
    Self: Sized,
{
    // builds a pool of OTs via OT-extension, to be used for OLE
    fn init<C: AbstractChannel, Crng: CryptoRng + Rng>(
        channel: &mut C,
        rng: &mut Crng,
    ) -> Result<Self, OleError>;

    fn input<Crng: CryptoRng + Rng>(self, a: &[F], b: &[F], rng: &mut Crng) -> Result<(), OleError>;
    // {


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
//
    //     Ok(())
    // }
}

pub struct OleSenderImpl<OT: OTReceiver> {
    ot: OT
}

impl<F: OleField, OT: OTReceiver> Sender<F> for OleSenderImpl<OT> {
    fn init<C: AbstractChannel, Crng: CryptoRng + Rng>(
        channel: &mut C,
        rng: &mut Crng,
    ) -> Result<Self, OleError> {
        let ot = OT::init(channel, rng)?;
        Ok(Self{
            ot: ot,
        })
    }

    fn input<Crng: CryptoRng + Rng>(self, a: &[F], b: &[F], rng: &mut Crng) -> Result<(), OleError> {
        assert_eq!(a.len(), b.len());
        let mut rng = rand::thread_rng();

        let mask = (0..F::B).map(|_| F::random(&mut rng));
        let secret = F::random(&mut rng);
        let shares = shamir::share(&secret, F::B as u64, (F::B - F::A) as u64, &F::beta());

        let mut hasher = Sha256::new();
        secret.into_repr().write_be(&mut hasher);
        // let mut secret_bytes = secret.into_repr().as_ref().to_le_bytes();
        // hasher.update(&secret_bytes);
        let com = hasher.finalize();
    }
}

pub trait Receiver<F: OleField>
where
    Self: Sized,
{
    // builds a pool of OTs via OT-extension, to be used for OLE
    fn init<C: AbstractChannel, Crng: CryptoRng + Rng>(
        channel: &mut C,
        rng: &mut Crng,
    ) -> Result<Self, OleError>;

    fn input(self, x: &[F]) -> Result<Vec<F>, OleError>;
}

#[cfg(tests)]
mod tests {
    fn test_smoke() {
        let mut rng = rand::thread_rng();
        let channel
        let a = Fp::random(&mut rng);

        let bytes =
    }
}
