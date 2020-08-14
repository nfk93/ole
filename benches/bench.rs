#[macro_use]
extern crate criterion;

use ole::fft::{fft2_in_place, fft2_inverse, fft2, fft3_in_place, fft3_inverse, fft3, digit_reverse_swap};
use ole::field::{Fp, OleField};
use ole::poly::{euclid_division, poly_from_roots, lagrangian_interpolation};
use ff::Field;
use rand;
use rand::seq::IteratorRandom;
use criterion::Criterion;
use ole::shamir::{share, reconstruct};

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
    c.bench_function(&format!("digit_reverse_swap, size {}, base 2", n), move |b| {
        b.iter(|| digit_reverse_swap(&mut data, 2))
    });

    let n = 2usize.pow(22 as u32);
    let mut data: Vec<Fp> = (0..n).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function(&format!("digit_reverse_swap, size {}, base 2", n), move |b| {
        b.iter(|| digit_reverse_swap(&mut data, 2))
    });

    let n = 3usize.pow(10 as u32);
    let mut data: Vec<Fp> = (0..n).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function(&format!("digit_reverse_swap, size {}, base 3", n), move |b| {
        b.iter(|| digit_reverse_swap(&mut data, 3))
    });

    let n = 3usize.pow(15 as u32);
    let mut data: Vec<Fp> = (0..n).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function(&format!("digit_reverse_swap, size {}, base 3", n), move |b| {
        b.iter(|| digit_reverse_swap(&mut data, 3))
    });
}

pub fn bench_poly_from_roots(c: &mut Criterion) {
    let mut rng = rand::thread_rng();

    let n = 2usize.pow(11 as u32);
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
    c.bench_function(&format!("euclid_division1, size a = {}, size b = {}", na, nb), move |b_| {
        b_.iter(|| euclid_division(&mut a, &b))
    });

    let na = 3usize.pow(9 as u32);
    let n = 3usize.pow(9) - 2usize.pow(11);
    let mut a: Vec<Fp> = (0..na).map(|_| Fp::random(&mut rng)).collect();
    let mut b: Vec<Fp> = (0..nb).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function(&format!("euclid_division2, size a = {}, size b = {}", na, nb), move |b_| {
        b_.iter(|| euclid_division(&mut a, &b))
    });
}

pub fn bench_lagrange(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let n = 3u64.pow(7) - 2u64.pow(8);
    let ys: Vec<Fp> = (0..n).map(|_| Fp::random(&mut rng)).collect();
    let xs: Vec<Fp> = (0..n).map(|_| Fp::random(&mut rng)).collect();

    c.bench_function(&format!("lagrangian_interpolation in 0, {} points", n), move |b_| {
        b_.iter(|| lagrangian_interpolation(&xs, &ys))
    });
}

pub fn bench_share(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let secret = Fp::random(&mut rng);
    let n = 3u64.pow(7);
    let rho = n - 2u64.pow(8);
    let omega = Fp::beta().pow([9]);

    c.bench_function(&format!("share, n = {}, rho = {}", n, rho), move |b_| {
        b_.iter(|| share(&secret, n, rho, &omega))
    });
}

pub fn bench_reconstruct(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let secret = Fp::random(&mut rng);
    let n = 3u64.pow(7);
    let rho = n - 2u64.pow(8);
    let omega = Fp::beta().pow([9]);

    let shares = share(&secret, n, rho, &omega);
    let mut indices: Vec<u64> = (0..n).choose_multiple(&mut rng, rho as usize);
    indices.sort();
    let mut myshares: Vec<Fp> = indices.iter().map(|i| shares[*i as usize]).collect();

    c.bench_function(&format!("reconstruct, n = {}, rho = {}", n, rho), move |b_| {
        b_.iter(|| reconstruct(&indices, &myshares, n, rho, &omega))
    });
}

criterion_group!(bench_ss, bench_share, bench_reconstruct);
criterion_group!(bench_poly, bench_poly_from_roots, bench_euclid_division, bench_lagrange);
criterion_group!(bench_fft2, bench_fft2_in_place, bench_fft2_inverse, bench_fft2_out_of_place);
criterion_group!(bench_fft3, bench_fft3_in_place, bench_fft3_inverse, bench_fft3_out_of_place);
criterion_group!(bench_digit_reverse, bench_digit_reverse_swap);
criterion_main!(bench_fft2, bench_fft3, bench_digit_reverse, bench_poly, bench_ss);
