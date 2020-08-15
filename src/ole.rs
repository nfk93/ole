use crate::error::OleError;
use crate::field::OleField;
use crate::shamir;
use crate::poly;
use crate::encoding;
use ff::{Field, PrimeField, PrimeFieldRepr};
use ocelot::ot::{KosReceiver, KosSender, Receiver as OTReceiver, Sender as OTSender};
use rand::{CryptoRng, Rng, seq::IteratorRandom};
use scuttlebutt::{Block, Channel, channel::AbstractChannel};
use sha2::{Digest, Sha256};
// use itertools::interleave;

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
        &mut self,
        a: &[F],
        b: &[F],
        c: &mut C,
        rng: &mut Crng,
    ) -> Result<(), OleError>;
}

pub struct OleSender {
    ot: KosSender,
}

impl Sender for OleSender {
    fn init<C: AbstractChannel, Crng: CryptoRng + Rng>(
        channel: &mut C,
        rng: &mut Crng,
    ) -> Result<Self, OleError> {
        let ot = KosSender::init(channel, rng)?;
        Ok(Self { ot: ot })
    }

    fn input<F: OleField, C: AbstractChannel, Crng: CryptoRng + Rng>(
        &mut self,
        a: &[F],
        b: &[F],
        channel: &mut C,
        rng: &mut Crng,
    ) -> Result<(), OleError> {
        assert_eq!(a.len(), b.len());

        let mask: Vec<F> = (0..F::B).map(|_| F::random(rng)).collect();
        let secret = F::random(rng);
        let shares: Vec<F> = shamir::share(&secret, F::B as u64, (F::B - F::A) as u64, &F::beta());

        let mut hasher = Sha256::new();
        secret.into_repr().write_be(&mut hasher);

        let com = hasher.finalize();
        channel.write_bytes(com.as_ref())?;
        channel.flush()?;

        let shares_blocks: Vec<Block> = shares.iter().map(|share| share.to_block()).collect();
        let mask_blocks: Vec<Block> = mask.iter().map(|mask_elem| mask_elem.to_block()).collect();
        let ot_input = shares_blocks
            .into_iter()
            .zip(mask_blocks.into_iter())
            .collect::<Vec<(Block, Block)>>();
        // let ot_input: Vec<(Block, Block)> = shares.into_iter().zip(mask.into_iter()).map(|(share, mask): (&F, &F)| {
        //     (share.to_block(), mask.to_block())
        // }).collect();
        // let ot_input: Vec<(Block, Block)> = shares.iter().zip(&mask).collect();
        self.ot.send(channel, ot_input.as_slice(), rng)?;

        let v_blocks = channel.read_blocks(F::B)?;

        let mut alpha2 = F::alpha();
        alpha2.square();

        let mut a_poly = a.to_vec();
        a_poly.resize_with(F::A/2, || F::random(rng));
        F::fft2_inverse(&mut a_poly, &alpha2);
        // a_poly.resize_with(F::B, F::zero);

        let mut b_poly = encoding::pad_every_other(&b, rng);
        b_poly.resize_with(F::A, || F::random(rng));
        for (idx, b_) in b.iter().enumerate() {
            assert_eq!(b_poly[2*idx], *b_);
        }
        F::fft2_inverse(&mut b_poly, &F::alpha());
        // b_poly.resize_with(F::B, F::zero);

        let mut a_vals = a_poly.to_vec();
        a_vals.resize_with(F::B, F::zero);
        F::fft3(&mut a_vals, &F::beta());

        let mut b_vals = b_poly.to_vec();
        b_vals.resize_with(F::B, F::zero);
        F::fft3(&mut b_vals, &F::beta());

        for (i, v) in v_blocks.iter().enumerate() {
            a_vals[i].mul_assign(&(*v).into());
            a_vals[i].add_assign(&b_vals[i]);
            a_vals[i].add_assign(&mask[i]);
        }

        for w in a_vals.iter() {
            channel.write_block(&w.to_block())?;
        }
        channel.flush()?;

        let zr = (channel.read_block()?).into();
        let a_zr = poly::horner(&a_poly, &zr);
        let b_zr = poly::horner(&b_poly, &zr);
        let zs = F::random(rng);

        channel.write_block(&a_zr.to_block())?;
        channel.write_block(&b_zr.to_block())?;
        channel.write_block(&zs.to_block())?;
        channel.flush()?;



        return Ok(())
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
        &mut self,
        x: &[F],
        channel: &mut C,
        rng: &mut Crng,
    ) -> Result<Vec<F>, OleError>;
}

pub struct OleReceiver {
    ot: KosReceiver,
}

impl Receiver for OleReceiver {
    fn init<C: AbstractChannel, Crng: CryptoRng + Rng>(
        channel: &mut C,
        rng: &mut Crng,
    ) -> Result<Self, OleError> {
        let ot = KosReceiver::init(channel, rng)?;
        Ok(Self { ot: ot })
    }

    fn input<F: OleField, C: AbstractChannel, Crng: CryptoRng + Rng>(
        &mut self,
        x: &[F],
        channel: &mut C,
        rng: &mut Crng,
    ) -> Result<Vec<F>, OleError> {
        let mut com = [0u8; 32];
        channel.read_bytes(&mut com);

        let t = F::A/2;
        assert!(x.len() <= t);
        let (mut encoded, x_poly, indices) = encoding::encode_reed_solomon(&x, rng);

        let mut share_indices = vec!();
        let mut j = 0;
        for i in 0..F::B {
            if (j < t as usize) && (i == indices[j]) {
                j += 1;
            } else {
                share_indices.push(i);
            }
        }
        assert_eq!(share_indices.len(),  F::B - t);

        let mut j = 0;
        let choices: Vec<bool> = (0..F::B).map(|i| {
            if (j < t as usize) && (i == indices[j]) {
                j+=1;
                true
            } else {
                false
            }
        }).collect();
        let vals = self.ot.receive(channel, &choices, rng)?;

        let mask: Vec<F> = indices.iter().map(|i| vals[*i as usize].into()).collect();
        let mut shares: Vec<F> = share_indices.iter().map(|i| vals[*i as usize].into()).collect();
        let secret = shamir::reconstruct(&share_indices, &shares, F::B as u64, (F::B - t) as u64, &F::beta());

        // check that commitment is correct
        let mut hasher = Sha256::new();
        secret.into_repr().write_be(&mut hasher);
        let com_check = hasher.finalize();
        assert_eq!(com_check.as_slice(), com);

        // // Make RS encoding of xs
        // let mut alpha2 = F::alpha();
        // alpha2.square();
        // let mut xs = x.to_vec();
        // xs.resize_with(F::A / 2, || F::random(rng));
        // // perform convolution to evaluate in F::B points
        // F::fft2_inverse(&mut xs, &alpha2);
        // let x_poly = xs.to_vec(); // save the polynomial X
        // xs.resize_with(F::B, F::zero);
        // F::fft3(&mut xs, &F::beta());
        // let mut j = 0;
        // for (i, x) in xs.iter_mut().enumerate() {
        //     if j < t && i == indices[j] {
        //         j += 1;
        //         continue;
        //     } else {
        //         *x = F::random(rng);
        //     }
        // }
        // Send to Sender

        for v in encoded.iter() {
            channel.write_block(&v.to_block())?;
        }
        channel.flush()?;

        let mut w_blocks = channel.read_blocks(F::B)?;
        let mut ws: Vec<F> = w_blocks.iter().map(|w| (*w).into()).collect();
        for (i, ti) in indices.iter().zip(mask) {
            ws[*i].sub_assign(&ti);
        }
        let y_poly = encoding::decode_reed_solomon(&mut ws, &indices);

        let zr = F::random(rng);
        channel.write_block(&zr.to_block())?;
        channel.flush()?;

        let a_zr: F = (channel.read_block()?).into();
        let b_zr: F = (channel.read_block()?).into();
        let zs: F = (channel.read_block()?).into();

        let mut x_zr = poly::horner(&x_poly, &zr);
        let y_zr = poly::horner(&y_poly, &zr);
        x_zr.mul_assign(&a_zr);
        x_zr.add_assign(&b_zr);
        assert_eq!(x_zr, y_zr);

        return Ok(vec!(F::one()));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // use std::io::{BufReader, BufWriter};
    use rand;
    use std::{
        fmt::Display,
        io::prelude::*,
        io::{BufReader, BufWriter},
        os::unix::net::UnixStream,
        sync::{Arc, Mutex},
    };
    use crate::field::Fp;

    #[test]
    fn test_ole() {
        let mut rng = rand::thread_rng();
        let (sender, receiver) = UnixStream::pair().unwrap();

        let handle = std::thread::spawn(move || {
            let mut rng = rand::thread_rng();
            let reader = BufReader::new(receiver.try_clone().unwrap());
            let writer = BufWriter::new(receiver);
            let mut channel = Channel::new(reader, writer);
            let mut olereceiver = OleReceiver::init(&mut channel, &mut rng).unwrap();
            let v = olereceiver.input(&[Fp::one()], &mut channel, &mut rng).unwrap();
            println!("result: {:?}", v)
        });

        let reader = BufReader::new(sender.try_clone().unwrap());
        let writer = BufWriter::new(sender);
        let mut channel = Channel::new(reader, writer);

        let mut olesender = OleSender::init(&mut channel, &mut rng).unwrap();
        olesender.input(&[Fp::one()], &[Fp::one()], &mut channel, &mut rng).unwrap();
        handle.join();
    }

    #[test]
    fn test_channel() {
        let n = 1000u64;
        let (sender, receiver) = UnixStream::pair().unwrap();
        let xs: Vec<Fp> = (0..n).map(|i| Fp::alpha().pow([i])).collect();

        let handle = std::thread::spawn(move || {
            let xs: Vec<Fp> = (0..n).map(|i| Fp::alpha().pow([i])).collect();
            let reader = BufReader::new(receiver.try_clone().unwrap());
            let writer = BufWriter::new(receiver);
            let mut channel = Channel::new(reader, writer);
            for x in xs {
                let x_received: Fp = (channel.read_block().unwrap()).into();
                assert_eq!(x, x_received);
            }
        });

        let reader = BufReader::new(sender.try_clone().unwrap());
        let writer = BufWriter::new(sender);
        let mut channel = Channel::new(reader, writer);
        for x in xs.iter() {
            channel.write_block(&x.to_block()).unwrap();
        }
        channel.flush().unwrap();
        handle.join();
    }

    #[test]
    // test copied from ocelot's test suite
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
