#[macro_use]
extern crate criterion;

use ole::field::{Fp, fft2_in_place_u16, fft2_forward_alpha};
use ole::ole::OleField;
use ff::Field;
use rand;
use criterion::Criterion;

pub fn bench_fft2_in_place(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let coeffs: Vec<Fp> = (0..1024).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function("fft2_in_place, 1024 points", move |b| {
        b.iter(|| fft2_in_place(&coeffs, &Fp::alpha_gen))
    });
}

pub fn bench_fft2_out_of_place(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let coeffs: Vec<Fp> = (0..1024).map(|_| Fp::random(&mut rng)).collect();
    c.bench_function("fft2, 1024 points", move |b| {
        b.iter(|| fft2(&coeffs, &Fp::alpha_gen))
    });
}

criterion_group!(bench_fft2, bench_fft2_in_place, bench_fft2_out_of_place);
criterion_main!(bench_fft2);
