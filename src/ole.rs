use crate::encoding;
use crate::error::OleError;
use crate::field::OleField;
use crate::poly;
use crate::shamir;
use ff::PrimeFieldRepr;
use ocelot::ot::{KosReceiver, KosSender, Receiver as OTReceiver, Sender as OTSender};
use rand::{CryptoRng, Rng};
use scuttlebutt::{channel::AbstractChannel, Block};
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
        secret.into_repr().write_be(&mut hasher)?;

        let com = hasher.finalize();
        channel.write_bytes(com.as_ref())?;
        channel.flush()?;

        let shares_blocks = shares.iter().map(|share| share.to_block());
        let mask_blocks = mask.iter().map(|mask_elem| mask_elem.to_block());
        let ot_input: Vec<(Block, Block)> = shares_blocks.zip(mask_blocks).collect();
        self.ot.send(channel, ot_input.as_slice(), rng)?;

        let secret_from_receiver = (channel.read_block()?).into();
        assert_eq!(secret, secret_from_receiver);

        let v_blocks = channel.read_blocks(F::B)?;

        let mut alpha2 = F::alpha();
        alpha2.square();

        let mut a_poly = a.to_vec();
        a_poly.resize_with(F::A / 2, || F::random(rng));
        F::fft2_inverse(&mut a_poly, &alpha2);
        let mut a_vals = a_poly.to_vec();
        a_vals.resize_with(F::B, F::zero);
        F::fft3(&mut a_vals, &F::beta());

        let mut b_poly = encoding::pad_every_other(&b, rng);
        b_poly.resize_with(F::A, || F::random(rng));
        F::fft2_inverse(&mut b_poly, &F::alpha());
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

        let mut ax_b = poly::horner(&a_poly, &zs);
        let b_zs = poly::horner(&b_poly, &zs);
        let x_zs: F = (channel.read_block()?).into();
        let y_zs: F = (channel.read_block()?).into();
        ax_b.mul_assign(&x_zs);
        ax_b.add_assign(&b_zs);
        assert_eq!(ax_b, y_zs);

        return Ok(());
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
        channel.read_bytes(&mut com)?;

        let (encoded, x_poly, indices) = encoding::encode_reed_solomon(&x, rng);

        let mut share_indices = vec![];
        let mut j = 0;
        for i in 0..F::B {
            if (j < F::A as usize) && (i == indices[j]) {
                j += 1;
            } else {
                share_indices.push(i);
            }
        }

        let mut j = 0;
        let choices: Vec<bool> = (0..F::B)
            .map(|i| {
                if (j < F::A as usize) && (i == indices[j]) {
                    j += 1;
                    true
                } else {
                    false
                }
            })
            .collect();
        let vals = self.ot.receive(channel, &choices, rng)?;

        let mask: Vec<F> = indices.iter().map(|i| vals[*i as usize].into()).collect();
        let shares: Vec<F> = share_indices
            .iter()
            .map(|i| vals[*i as usize].into())
            .collect();

        // check that commitment is correct
        let secret = shamir::reconstruct(
            &share_indices,
            &shares,
            F::B as u64,
            (F::B - F::A) as u64,
            &F::beta(),
        );
        let mut hasher = Sha256::new();
        secret.into_repr().write_be(&mut hasher)?;
        let com_check = hasher.finalize();
        assert_eq!(com_check.as_slice(), com);
        channel.write_block(&secret.to_block())?;
        channel.flush()?;

        for v in encoded.iter() {
            channel.write_block(&v.to_block())?;
        }
        channel.flush()?;

        let w_blocks = channel.read_blocks(F::B)?;
        let mut ws: Vec<F> = w_blocks.iter().map(|w| (*w).into()).collect();
        for (i, ti) in indices.iter().zip(&mask) {
            ws[*i].sub_assign(&ti);
        }

        let mut y_poly = encoding::decode_reed_solomon(&mut ws, &indices);
        assert!(y_poly.len() == F::A);

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

        let x_zs = poly::horner(&x_poly, &zs);
        let y_zs = poly::horner(&y_poly, &zs);
        channel.write_block(&x_zs.to_block())?;
        channel.write_block(&y_zs.to_block())?;
        channel.flush()?;

        F::fft2(&mut y_poly, &F::alpha());
        let result = y_poly.into_iter().step_by(2).collect::<Vec<F>>();
        // let result: Vec<F> = (0..F::A/2).map(|i| poly::horner(&y_poly, &F::alpha().pow([(2*i) as u64]))).collect();
        return Ok(result);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::Fp;
    use ff::Field;
    use rand;
    use scuttlebutt::channel::Channel;
    use std::{
        io::{BufReader, BufWriter},
        os::unix::net::UnixStream,
    };

    #[test]
    fn test_ole() {
        let mut rng = rand::thread_rng();
        let (sender, receiver) = UnixStream::pair().unwrap();
        let b: Vec<Fp> = (0..Fp::A / 2).map(|_| Fp::random(&mut rng)).collect();
        let a: Vec<Fp> = (0..Fp::A / 2).map(|_| Fp::random(&mut rng)).collect();
        // let a: Vec<Fp> = (0..Fp::A/2).map(|i| Fp::one()).collect();
        // let b: Vec<Fp> = (0..Fp::A/2).map(|i| Fp::one()).collect();

        let a_copy = a.to_vec();
        let b_copy = b.to_vec();
        let handle = std::thread::spawn(move || {
            let mut rng = rand::thread_rng();
            let reader = BufReader::new(sender.try_clone().unwrap());
            let writer = BufWriter::new(sender);
            let mut channel = Channel::new(reader, writer);

            let mut olesender = OleSender::init(&mut channel, &mut rng).unwrap();
            olesender
                .input(&a_copy, &b_copy, &mut channel, &mut rng)
                .unwrap();
        });

        let reader = BufReader::new(receiver.try_clone().unwrap());
        let writer = BufWriter::new(receiver);
        let mut channel = Channel::new(reader, writer);
        let mut olereceiver = OleReceiver::init(&mut channel, &mut rng).unwrap();
        let x: Vec<Fp> = (0..Fp::A / 2).map(|_| Fp::random(&mut rng)).collect();
        let result = olereceiver.input(&x, &mut channel, &mut rng).unwrap();
        handle.join().unwrap();

        for i in 0..Fp::A / 2 {
            let mut expected = x[i];
            expected.mul_assign(&a[i]);
            expected.add_assign(&b[i]);
            println!(
                "expected:\n{:?} * {:?} + {:?} = {:?}",
                a[i], x[i], b[i], expected
            );
            println!("received: {:?}", result[i]);
            assert_eq!(result[i], expected);
        }
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
        handle.join().unwrap();
    }

    #[test]
    // test copied from ocelot's test suite
    fn test_ot() {
        let mut rng = rand::thread_rng();

        let ninputs = 128;
        let m0s: Vec<Fp> = (0..ninputs).map(|_| Fp::random(&mut rng)).collect();
        let m1s: Vec<Fp> = (0..ninputs).map(|_| Fp::random(&mut rng)).collect();
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
                .map(|(m0, m1)| (m0.to_block(), m1.to_block()))
                .collect::<Vec<(Block, Block)>>();
            ot.send(&mut channel, &ms, &mut rng).unwrap();
        });

        let reader = BufReader::new(receiver.try_clone().unwrap());
        let writer = BufWriter::new(receiver);
        let mut channel = Channel::new(reader, writer);
        let mut ot = KosReceiver::init(&mut channel, &mut rng).unwrap();
        let result_ = ot.receive(&mut channel, &bs, &mut rng).unwrap();
        let result: Vec<Fp> = result_.iter().map(|v| (*v).into()).collect();
        handle.join().unwrap();
        for j in 0..ninputs {
            assert_eq!(result[j], if bs[j] { m1s_[j] } else { m0s_[j] });
        }
    }
}
