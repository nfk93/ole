#[macro_use]
extern crate criterion;

use criterion::Criterion;
use ff::Field;
use ole::fft::{
    digit_reverse_swap, fft2, fft2_in_place, fft2_inverse, fft3, fft3_in_place, fft3_inverse,
};
use ole::field::{Fp, OleField};
use ole::ole::{OleReceiver, OleSender, Receiver, Sender};
use ole::poly::{euclid_division, lagrangian_interpolation, poly_from_roots};
use ole::shamir::{reconstruct, share};
use rand;
use rand::seq::IteratorRandom;
use scuttlebutt::{channel::AbstractChannel, Channel};
use std::{
    io::{BufReader, BufWriter},
    os::unix::net::UnixStream,
};

pub fn bench_fft2_in_place(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let mut coeffs: Vec<Fp> = (0..1024).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function("fft2_in_place, 1024 points", move |b| {
        b.iter(|| fft2_in_place(&mut coeffs, &Fp::alpha()))
    });
}

pub fn bench_fft2_out_of_place(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let coeffs: Vec<Fp> = (0..1024).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function("fft2, 1024 points", move |b| {
        b.iter(|| fft2(&coeffs, &Fp::alpha()))
    });
}

pub fn bench_fft3_in_place(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let mut coeffs: Vec<Fp> = (0..Fp::B).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function(&format!("fft3_in_place, {} points", Fp::B), move |b| {
        b.iter(|| fft3_in_place(&mut coeffs, &Fp::beta()))
    });
}

pub fn bench_fft3_out_of_place(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let coeffs: Vec<Fp> = (0..Fp::B).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function(&format!("fft3, {} points", Fp::B), move |b| {
        b.iter(|| fft3(&coeffs, &Fp::beta()))
    });
}

pub fn bench_fft2_inverse(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let mut points: Vec<Fp> = (0..Fp::A).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function(&format!("fft2_inverse, {} points", Fp::A), move |b| {
        b.iter(|| fft2_inverse(&mut points, &Fp::beta()))
    });
}

pub fn bench_fft3_inverse(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let mut points: Vec<Fp> = (0..Fp::B).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function(&format!("fft3_inverse, {} points", Fp::B), move |b| {
        b.iter(|| fft3_inverse(&mut points, &Fp::beta()))
    });
}

pub fn bench_digit_reverse_swap(c: &mut Criterion) {
    let mut rng = rand::thread_rng();

    let n = 2usize.pow(11 as u32);
    let mut data: Vec<Fp> = (0..n).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function(
        &format!("digit_reverse_swap, size {}, base 2", n),
        move |b| b.iter(|| digit_reverse_swap(&mut data, 2)),
    );

    let n = 2usize.pow(22 as u32);
    let mut data: Vec<Fp> = (0..n).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function(
        &format!("digit_reverse_swap, size {}, base 2", n),
        move |b| b.iter(|| digit_reverse_swap(&mut data, 2)),
    );

    let n = 3usize.pow(10 as u32);
    let mut data: Vec<Fp> = (0..n).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function(
        &format!("digit_reverse_swap, size {}, base 3", n),
        move |b| b.iter(|| digit_reverse_swap(&mut data, 3)),
    );

    let n = 3usize.pow(15 as u32);
    let mut data: Vec<Fp> = (0..n).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function(
        &format!("digit_reverse_swap, size {}, base 3", n),
        move |b| b.iter(|| digit_reverse_swap(&mut data, 3)),
    );
}

pub fn bench_poly_from_roots(c: &mut Criterion) {
    let mut rng = rand::thread_rng();

    let n = 10000;
    let mut data: Vec<Fp> = (0..n).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function(&format!("poly_from_roots1, size {}", n), move |b| {
        b.iter(|| poly_from_roots(&mut data))
    });

    let n = 3usize.pow(7) - 2usize.pow(8);
    let mut data: Vec<Fp> = (0..n).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function(&format!("poly_from_roots2, size {}", n), move |b| {
        b.iter(|| poly_from_roots(&mut data))
    });
}

pub fn bench_euclid_division(c: &mut Criterion) {
    let mut rng = rand::thread_rng();

    let na = 3usize.pow(9 as u32);
    let nb = 2usize.pow(11 as u32);
    let mut a: Vec<Fp> = (0..na).map(|_| Fp::random(&mut rng)).collect();
    let mut b: Vec<Fp> = (0..nb).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function(
        &format!("euclid_division, size a = {}, size b = {}", na, nb),
        move |b_| b_.iter(|| euclid_division(&mut a, &b)),
    );

    let na = Fp::B;
    let n = Fp::A;
    let mut a: Vec<Fp> = (0..na).map(|_| Fp::random(&mut rng)).collect();
    let mut b: Vec<Fp> = (0..nb).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function(
        &format!("euclid_division, size a = {}, size b = {}", na, nb),
        move |b_| b_.iter(|| euclid_division(&mut a, &b)),
    );
}

pub fn bench_lagrange(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let n = Fp::B - Fp::A;
    let ys: Vec<Fp> = (0..n).map(|_| Fp::random(&mut rng)).collect();
    let xs: Vec<Fp> = (0..n).map(|_| Fp::random(&mut rng)).collect();

    c.bench_function(
        &format!("lagrangian_interpolation in 0, {} points", n),
        move |b_| b_.iter(|| lagrangian_interpolation(&xs, &ys)),
    );
}

pub fn bench_share(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let secret = Fp::random(&mut rng);
    let n = Fp::B;
    let rho = Fp::B - Fp::A;
    let omega = Fp::beta();

    c.bench_function(&format!("share, n = {}, rho = {}", n, rho), move |b_| {
        b_.iter(|| share(&secret, n as u64, rho as u64, &omega))
    });
}

pub fn bench_reconstruct(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let secret = Fp::random(&mut rng);
    let n = Fp::B;
    let rho = Fp::B - Fp::A;
    let omega = Fp::beta();

    let shares = share(&secret, n as u64, rho as u64, &omega);
    let mut indices: Vec<usize> = (0..n).choose_multiple(&mut rng, rho);
    indices.sort();
    let mut myshares: Vec<Fp> = indices.iter().map(|i| shares[*i]).collect();

    c.bench_function(
        &format!("reconstruct, n = {}, rho = {}", n, rho),
        move |b_| b_.iter(|| reconstruct(&indices, &myshares, n as u64, rho as u64, &omega)),
    );
}

fn run_ole_bench<F: OleField>(n: usize, a: Vec<F>, b: Vec<F>, x: Vec<F>) {
    let (sender, receiver) = UnixStream::pair().unwrap();
    let handle = std::thread::spawn(move || {
        let mut rng = rand::thread_rng();
        let reader = BufReader::new(sender.try_clone().unwrap());
        let writer = BufWriter::new(sender);
        let mut sender_channel = Channel::new(reader, writer);
        let mut olesender = OleSender::init(&mut sender_channel, &mut rng).unwrap();
        for _ in 0..n {
            olesender
                .input(&a.clone(), &b.clone(), &mut sender_channel, &mut rng)
                .unwrap();
        }
    });
    let mut rng = rand::thread_rng();
    let reader = BufReader::new(receiver.try_clone().unwrap());
    let writer = BufWriter::new(receiver);
    let mut receiver_channel = Channel::new(reader, writer);
    let mut olereceiver = OleReceiver::init(&mut receiver_channel, &mut rng).unwrap();
    for _ in 0..n {
        let result = olereceiver
            .input(&x, &mut receiver_channel, &mut rng)
            .unwrap();
    }
    handle.join();
}

pub fn bench_ole_send_receive(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let a: Vec<Fp> = (0..Fp::A / 2).map(|i| Fp::random(&mut rng)).collect();
    let b: Vec<Fp> = (0..Fp::A / 2).map(|i| Fp::random(&mut rng)).collect();
    let x: Vec<Fp> = (0..Fp::A / 2).map(|i| Fp::random(&mut rng)).collect();
    let n = 1;
    c.bench_function(&format!("ole128 * {}", n), move |b_| {
        b_.iter(|| {
            run_ole_bench(n, a.clone(), b.clone(), x.clone());
        })
    });
    let a: Vec<Fp> = (0..Fp::A / 2).map(|i| Fp::random(&mut rng)).collect();
    let b: Vec<Fp> = (0..Fp::A / 2).map(|i| Fp::random(&mut rng)).collect();
    let x: Vec<Fp> = (0..Fp::A / 2).map(|i| Fp::random(&mut rng)).collect();
    let n = 2;
    c.bench_function(&format!("ole128 * {}", n), move |b_| {
        b_.iter(|| {
            run_ole_bench(n, a.clone(), b.clone(), x.clone());
        })
    });
    let a: Vec<Fp> = (0..Fp::A / 2).map(|i| Fp::random(&mut rng)).collect();
    let b: Vec<Fp> = (0..Fp::A / 2).map(|i| Fp::random(&mut rng)).collect();
    let x: Vec<Fp> = (0..Fp::A / 2).map(|i| Fp::random(&mut rng)).collect();
    let n = 3;
    c.bench_function(&format!("ole128 * {}", n), move |b_| {
        b_.iter(|| {
            run_ole_bench(n, a.clone(), b.clone(), x.clone());
        })
    });
}

pub fn bench_ole_init(c: &mut Criterion) {
    c.bench_function("ole init", move |b_| {
        b_.iter(|| {
            let (sender, receiver) = UnixStream::pair().unwrap();
            let handle = std::thread::spawn(move || {
                let mut rng = rand::thread_rng();
                let reader = BufReader::new(sender.try_clone().unwrap());
                let writer = BufWriter::new(sender);
                let mut sender_channel = Channel::new(reader, writer);
                let olesender = OleSender::init(&mut sender_channel, &mut rng).unwrap();
            });
            let mut rng = rand::thread_rng();
            let reader = BufReader::new(receiver.try_clone().unwrap());
            let writer = BufWriter::new(receiver);
            let mut receiver_channel = Channel::new(reader, writer);
            let mut olereceiver = OleReceiver::init(&mut receiver_channel, &mut rng).unwrap();
            handle.join();
        })
    });
}

criterion_group!(bench_ole, bench_ole_send_receive, bench_ole_init);
criterion_group!(bench_ss, bench_share, bench_reconstruct);
criterion_group!(
    bench_poly,
    bench_poly_from_roots,
    bench_euclid_division,
    bench_lagrange
);
criterion_group!(
    bench_fft2,
    bench_fft2_in_place,
    bench_fft2_inverse,
    bench_fft2_out_of_place
);
criterion_group!(
    bench_fft3,
    bench_fft3_in_place,
    bench_fft3_inverse,
    bench_fft3_out_of_place
);
criterion_group!(bench_digit_reverse, bench_digit_reverse_swap);
criterion_main!(
    bench_fft2,
    bench_fft3,
    bench_digit_reverse,
    bench_poly,
    bench_ss,
    bench_ole
);
