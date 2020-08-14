use crate::error::OleError;
use crate::field::*;
use crate::shamir;
use ff::{Field, PrimeField, PrimeFieldRepr};
use ocelot::ot::{KosReceiver, KosSender, Receiver as OTReceiver, Sender as OTSender};
use rand::{CryptoRng, Rng};
use scuttlebutt::channel::{AbstractChannel, Channel};
use sha2::{Digest, Sha256};

pub trait Sender
where
    Self: Sized,
{
    // builds a pool of OTs via OT-extension, to be used for OLE
    fn init<C: AbstractChannel, Crng: CryptoRng + Rng>(
        channel: &mut C,
        rng: &mut Crng,
    ) -> Result<Self, OleError>;

    fn input<F: OleField, C: AbstractChannel, Crng: CryptoRng + Rng>(
        self,
        a: &[F],
        b: &[F],
        c: &mut C,
        rng: &mut Crng,
    ) -> Result<(), OleError>;
    // {
    //
    //
    //     // check that input has the correct size
    //     let l = a.len();
    //     if l > F::A {
    //         return Err(OleError::TODOError)
    //     }
    //
    //     let mut padded_a = a.to_vec();
    //     for i in l..F::A {
    //         padded_a.push(F::random(rng));
    //     }
    //
    //     let mut padded_b = b.to_vec();
    //     for i in l..F::A {
    //         padded_b.push(F::random(rng));
    //     }
    //
    //     let a_poly = F::interpolate_alpha(&padded_a);
    //     // let b_poly = F::interpolate_alpha(&padded_b);
    //
    //     Ok(())
    // }
}

pub struct OleSender<OT: OTSender> {
    ot: OT,
}

impl<OT: OTSender> Sender for OleSender<OT> {
    fn init<C: AbstractChannel, Crng: CryptoRng + Rng>(
        channel: &mut C,
        rng: &mut Crng,
    ) -> Result<Self, OleError> {
        let ot = OT::init(channel, rng)?;
        Ok(Self { ot: ot })
    }

    fn input<F: OleField, C: AbstractChannel, Crng: CryptoRng + Rng>(
        self,
        a: &[F],
        b: &[F],
        channel: &mut C,
        rng: &mut Crng,
    ) -> Result<(), OleError> {
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
        println!("secret: {:?}\ncommit: {:?}", &secret, &com);

        return Err(OleError::TODOError);
    }
}

pub trait Receiver
where
    Self: Sized,
{
    // builds a pool of OTs via OT-extension, to be used for OLE
    fn init<C: AbstractChannel, Crng: CryptoRng + Rng>(
        channel: &mut C,
        rng: &mut Crng,
    ) -> Result<Self, OleError>;

    fn input<F: OleField, C: AbstractChannel, Crng: CryptoRng + Rng>(
        self,
        x: &[F],
        channel: &mut C,
        rng: &mut Crng,
    ) -> Result<Vec<F>, OleError>;
}

pub struct OleReceiver<OT: OTReceiver> {
    ot: OT,
}

impl<OT: OTReceiver> Receiver for OleReceiver<OT> {
    fn init<C: AbstractChannel, Crng: CryptoRng + Rng>(
        channel: &mut C,
        rng: &mut Crng,
    ) -> Result<Self, OleError> {
        let ot = OT::init(channel, rng)?;
        Ok(Self { ot: ot })
    }

    fn input<F: OleField, C: AbstractChannel, Crng: CryptoRng + Rng>(
        self,
        x: &[F],
        channel: &mut C,
        rng: &mut Crng,
    ) -> Result<Vec<F>, OleError> {
        let mut rng = rand::thread_rng();

        let mask = (0..F::B).map(|_| F::random(&mut rng));
        let secret = F::random(&mut rng);
        let shares = shamir::share(&secret, F::B as u64, (F::B - F::A) as u64, &F::beta());

        let mut hasher = Sha256::new();
        secret.into_repr().write_be(&mut hasher);
        // let mut secret_bytes = secret.into_repr().as_ref().to_le_bytes();
        // hasher.update(&secret_bytes);
        let com = hasher.finalize();
        println!("secret: {:?}\ncommit: {:?}", &secret, &com);

        return Err(OleError::TODOError);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // use std::io::{BufReader, BufWriter};
    use rand;
    use rand::rngs::ThreadRng;
    use scuttlebutt::{AesRng, Block, Channel};
    use std::{
        fmt::Display,
        io::prelude::*,
        io::{BufReader, BufWriter},
        os::unix::net::UnixStream,
        sync::{Arc, Mutex},
    };

    #[test]
    fn test_ole() {}

    #[test]
    fn test_ot() {
        let ninputs = 128;
        let m0s: Vec<Block> = (0..ninputs).map(|_| rand::random::<Block>()).collect();
        let m1s: Vec<Block> = (0..ninputs).map(|_| rand::random::<Block>()).collect();
        let bs: Vec<bool> = (0..ninputs).map(|_| rand::random::<bool>()).collect();
        let m0s_ = m0s.clone();
        let m1s_ = m1s.clone();
        let (sender, receiver) = UnixStream::pair().unwrap();

        let handle = std::thread::spawn(move || {
            let mut rng = rand::thread_rng();
            let reader = BufReader::new(sender.try_clone().unwrap());
            let writer = BufWriter::new(sender);
            let mut channel = Channel::new(reader, writer);
            let mut ot = KosSender::init(&mut channel, &mut rng).unwrap();
            let ms = m0s
                .into_iter()
                .zip(m1s.into_iter())
                .collect::<Vec<(Block, Block)>>();
            ot.send(&mut channel, &ms, &mut rng).unwrap();
        });

        let mut rng = rand::thread_rng();
        let reader = BufReader::new(receiver.try_clone().unwrap());
        let writer = BufWriter::new(receiver);
        let mut channel = Channel::new(reader, writer);
        let mut ot = KosReceiver::init(&mut channel, &mut rng).unwrap();
        let result = ot.receive(&mut channel, &bs, &mut rng).unwrap();
        handle.join().unwrap();
        for j in 0..ninputs {
            assert_eq!(result[j], if bs[j] { m1s_[j] } else { m0s_[j] });
        }
    }
}
