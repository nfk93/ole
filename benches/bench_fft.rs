#[macro_use]
extern crate criterion;

use ole::field::{Fp, fft2_in_place, fft2_inverse, fft2, fft3_in_place, fft3_inverse, fft3, digit_reverse_swap};
use ole::ole::OleField;
use ff::Field;
use rand;
use criterion::Criterion;

pub fn bench_fft2_in_place(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let mut coeffs: Vec<Fp> = (0..1024).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function("fft2_in_place, 1024 points", move |b| {
        b.iter(|| fft2_in_place(&mut coeffs, &Fp::alpha_gen))
    });
}

pub fn bench_fft2_out_of_place(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let coeffs: Vec<Fp> = (0..1024).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function("fft2, 1024 points", move |b| {
        b.iter(|| fft2(&coeffs, &Fp::alpha_gen))
    });
}

pub fn bench_fft3_in_place(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let mut coeffs: Vec<Fp> = (0..Fp::B).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function(&format!("fft3_in_place, {} points", Fp::B), move |b| {
        b.iter(|| fft3_in_place(&mut coeffs, &Fp::beta_gen))
    });
}

pub fn bench_fft3_out_of_place(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let coeffs: Vec<Fp> = (0..Fp::B).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function(&format!("fft3, {} points", Fp::B), move |b| {
        b.iter(|| fft3(&coeffs, &Fp::beta_gen))
    });
}

pub fn bench_fft2_inverse(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let mut points: Vec<Fp> = (0..Fp::A).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function(&format!("fft2_inverse, {} points", Fp::A), move |b| {
        b.iter(|| fft2_inverse(&mut points, &Fp::beta_gen))
    });
}

pub fn bench_fft3_inverse(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let mut points: Vec<Fp> = (0..Fp::B).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function(&format!("fft3_inverse, {} points", Fp::B), move |b| {
        b.iter(|| fft3_inverse(&mut points, &Fp::beta_gen))
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

criterion_group!(bench_fft2, bench_fft2_in_place, bench_fft2_inverse, bench_fft2_out_of_place);
criterion_group!(bench_fft3, bench_fft3_in_place, bench_fft3_inverse, bench_fft3_out_of_place);
criterion_group!(bench_digit_reverse, bench_digit_reverse_swap);
criterion_main!(bench_fft2, bench_fft3, bench_digit_reverse);
