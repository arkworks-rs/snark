#[macro_use]
extern crate criterion;

#[macro_use]
extern crate bench_utils;

use algebra::fft::{
    get_best_evaluation_domain, BasicRadix2Domain, DensePolynomial, EvaluationDomain,
};
use algebra::{fields::tweedle::Fr, PrimeField, UniformRand};

use std::{
    fs::File,
    path::Path,
    time::{SystemTime, UNIX_EPOCH},
};

use criterion::{BatchSize, BenchmarkId, Criterion};
const DATA_PATH: &str = "./coeffs_tweedle";

fn save_data<F: PrimeField>(num_coeffs: usize) {
    let mut fs = File::create(DATA_PATH).unwrap();
    let rng = &mut rand::thread_rng();

    for _ in 0..num_coeffs {
        let elem: F = UniformRand::rand(rng);
        match elem.write(&mut fs) {
            Ok(_) => {}
            Err(msg) => {
                panic!("Cannot save coeffs to file: {}", msg)
            }
        }
    }
}

fn load_data<F: PrimeField>(samples: usize) -> Vec<F> {
    if !Path::new(DATA_PATH).exists() {
        save_data::<F>(1 << 23);
    }

    let mut fs = File::open(DATA_PATH).unwrap();
    let mut a: Vec<F> = Vec::with_capacity(samples);

    for _i in 0..samples {
        let elem1 = F::read(&mut fs).unwrap();
        a.push(elem1);
    }
    a
}

fn bench_ffts<F: PrimeField, D: EvaluationDomain<F>>(
    c: &mut Criterion,
    num_coeffs: usize,
    name: &'static str,
) {
    let mut group = c.benchmark_group(name);

    // We expect the num_coeffs input to be a compatible size for the domain.
    let domain = get_best_evaluation_domain::<F>(num_coeffs).unwrap();
    let domain_size = domain.size();
    assert_eq!(num_coeffs, domain_size);

    group.bench_with_input(
        BenchmarkId::from_parameter(num_coeffs),
        &num_coeffs,
        |b, _samples| {
            b.iter_batched(
                || {
                    let a: Vec<F> = load_data(num_coeffs);
                    a
                },
                |a| {
                    add_to_trace!(
                        || format!("****************{}*******************", domain_size),
                        || format!(
                            "--->START TIMESTAMP: {:?}",
                            SystemTime::now()
                                .duration_since(UNIX_EPOCH)
                                .unwrap()
                                .as_secs()
                        )
                    );

                    domain.fft(&mut a.as_slice());

                    add_to_trace!(
                        || format!("****************{}*******************", domain_size),
                        || format!(
                            "--->END TIMESTAMP: {:?}",
                            SystemTime::now()
                                .duration_since(UNIX_EPOCH)
                                .unwrap()
                                .as_secs()
                        )
                    );
                },
                BatchSize::PerIteration,
            );
        },
    );
}

fn bench_fft_tweedle(c: &mut Criterion) {
    bench_ffts::<Fr, BasicRadix2Domain<Fr>>(c, 1 << 14, "radix-2 FFT - 2^14 - tweedle");
    bench_ffts::<Fr, BasicRadix2Domain<Fr>>(c, 1 << 15, "radix-2 FFT - 2^15 - tweedle");
    bench_ffts::<Fr, BasicRadix2Domain<Fr>>(c, 1 << 16, "radix-2 FFT - 2^16 - tweedle");
    bench_ffts::<Fr, BasicRadix2Domain<Fr>>(c, 1 << 17, "radix-2 FFT - 2^17 - tweedle");
    bench_ffts::<Fr, BasicRadix2Domain<Fr>>(c, 1 << 18, "radix-2 FFT - 2^18 - tweedle");
    bench_ffts::<Fr, BasicRadix2Domain<Fr>>(c, 1 << 19, "radix-2 FFT - 2^19 - tweedle");
    bench_ffts::<Fr, BasicRadix2Domain<Fr>>(c, 1 << 20, "radix-2 FFT - 2^20 - tweedle");
    bench_ffts::<Fr, BasicRadix2Domain<Fr>>(c, 1 << 21, "radix-2 FFT - 2^21 - tweedle");
    bench_ffts::<Fr, BasicRadix2Domain<Fr>>(c, 1 << 22, "radix-2 FFT - 2^22 - tweedle");
    bench_ffts::<Fr, BasicRadix2Domain<Fr>>(c, 1 << 23, "radix-2 FFT - 2^23 - tweedle");
}

fn bench_iffts<F: PrimeField, D: EvaluationDomain<F>>(
    c: &mut Criterion,
    num_coeffs: usize,
    name: &'static str,
) {
    let mut group = c.benchmark_group(name);

    // We expect the num_coeffs input to be a compatible size for the domain.
    let domain = get_best_evaluation_domain::<F>(num_coeffs).unwrap();
    let domain_size = domain.size();
    assert_eq!(num_coeffs, domain_size);

    group.bench_with_input(
        BenchmarkId::from_parameter(num_coeffs),
        &num_coeffs,
        |b, _samples| {
            b.iter_batched(
                || {
                    let a: Vec<F> = load_data(num_coeffs);
                    a
                },
                |mut a| {
                    add_to_trace!(
                        || format!("****************{}*******************", domain_size),
                        || format!(
                            "--->START TIMESTAMP: {:?}",
                            SystemTime::now()
                                .duration_since(UNIX_EPOCH)
                                .unwrap()
                                .as_secs()
                        )
                    );

                    domain.ifft(&mut a);

                    add_to_trace!(
                        || format!("****************{}*******************", domain_size),
                        || format!(
                            "--->END TIMESTAMP: {:?}",
                            SystemTime::now()
                                .duration_since(UNIX_EPOCH)
                                .unwrap()
                                .as_secs()
                        )
                    );
                },
                BatchSize::PerIteration,
            );
        },
    );
}

fn bench_ifft_tweedle(c: &mut Criterion) {
    bench_iffts::<Fr, BasicRadix2Domain<Fr>>(c, 1 << 14, "radix-2 iFFT - 2^14 - tweedle");
    bench_iffts::<Fr, BasicRadix2Domain<Fr>>(c, 1 << 15, "radix-2 iFFT - 2^15 - tweedle");
    bench_iffts::<Fr, BasicRadix2Domain<Fr>>(c, 1 << 16, "radix-2 iFFT - 2^16 - tweedle");
    bench_iffts::<Fr, BasicRadix2Domain<Fr>>(c, 1 << 17, "radix-2 iFFT - 2^17 - tweedle");
    bench_iffts::<Fr, BasicRadix2Domain<Fr>>(c, 1 << 18, "radix-2 iFFT - 2^18 - tweedle");
    bench_iffts::<Fr, BasicRadix2Domain<Fr>>(c, 1 << 19, "radix-2 iFFT - 2^19 - tweedle");
    bench_iffts::<Fr, BasicRadix2Domain<Fr>>(c, 1 << 20, "radix-2 iFFT - 2^20 - tweedle");
    bench_iffts::<Fr, BasicRadix2Domain<Fr>>(c, 1 << 21, "radix-2 iFFT - 2^21 - tweedle");
    bench_iffts::<Fr, BasicRadix2Domain<Fr>>(c, 1 << 22, "radix-2 iFFT - 2^22 - tweedle");
    bench_iffts::<Fr, BasicRadix2Domain<Fr>>(c, 1 << 23, "radix-2 iFFT - 2^23 - tweedle");
}

fn bench_dense_poly_muls<F: PrimeField, D: EvaluationDomain<F>>(
    c: &mut Criterion,
    num_degree: usize,
    name: &'static str,
) {
    // Per benchmark setup
    let rng = &mut rand::thread_rng();

    c.bench_function(name, move |bencher| {
        let p1 = DensePolynomial::<F>::rand(num_degree, rng);
        let p2 = DensePolynomial::<F>::rand(num_degree, rng);

        bencher.iter(|| {
            add_to_trace!(
                || format!("****************{}*******************", num_degree),
                || format!(
                    "--->START TIMESTAMP: {:?}",
                    SystemTime::now()
                        .duration_since(UNIX_EPOCH)
                        .unwrap()
                        .as_secs()
                )
            );

            let _ab = (&p1) * (&p2);

            add_to_trace!(
                || format!("****************{}*******************", num_degree),
                || format!(
                    "--->END TIMESTAMP: {:?}",
                    SystemTime::now()
                        .duration_since(UNIX_EPOCH)
                        .unwrap()
                        .as_secs()
                )
            );
        })
    });
}

fn bench_dense_poly_mul_tweedle(c: &mut Criterion) {
    bench_dense_poly_muls::<Fr, BasicRadix2Domain<Fr>>(
        c,
        1 << 14,
        "radix-2 DensePolynomial::mul - 2^14 - tweedle",
    );
    bench_dense_poly_muls::<Fr, BasicRadix2Domain<Fr>>(
        c,
        1 << 15,
        "radix-2 DensePolynomial::mul - 2^15 - tweedle",
    );
    bench_dense_poly_muls::<Fr, BasicRadix2Domain<Fr>>(
        c,
        1 << 16,
        "radix-2 DensePolynomial::mul - 2^16 - tweedle",
    );
    bench_dense_poly_muls::<Fr, BasicRadix2Domain<Fr>>(
        c,
        1 << 17,
        "radix-2 DensePolynomial::mul - 2^17 - tweedle",
    );
    bench_dense_poly_muls::<Fr, BasicRadix2Domain<Fr>>(
        c,
        1 << 18,
        "radix-2 DensePolynomial::mul - 2^18 - tweedle",
    );
    bench_dense_poly_muls::<Fr, BasicRadix2Domain<Fr>>(
        c,
        1 << 19,
        "radix-2 DensePolynomial::mul - 2^19 - tweedle",
    );
    bench_dense_poly_muls::<Fr, BasicRadix2Domain<Fr>>(
        c,
        1 << 20,
        "radix-2 DensePolynomial::mul - 2^20 - tweedle",
    );
    bench_dense_poly_muls::<Fr, BasicRadix2Domain<Fr>>(
        c,
        1 << 21,
        "radix-2 DensePolynomial::mul - 2^21 - tweedle",
    );
    bench_dense_poly_muls::<Fr, BasicRadix2Domain<Fr>>(
        c,
        1 << 22,
        "radix-2 DensePolynomial::mul - 2^22 - tweedle",
    );
    bench_dense_poly_muls::<Fr, BasicRadix2Domain<Fr>>(
        c,
        1 << 23,
        "radix-2 DensePolynomial::mul - 2^23 - tweedle",
    );
}

fn bench_dense_poly_div_by_vanishing_poly<F: PrimeField, D: EvaluationDomain<F>>(
    c: &mut Criterion,
    num_coeffs: usize,
    name: &'static str,
) {
    // Per benchmark setup
    let rng = &mut rand::thread_rng();

    // We expect the num_coeffs input to be a compatible size for the domain.
    let domain = get_best_evaluation_domain::<F>(num_coeffs).unwrap();
    let domain_size = domain.size();
    assert_eq!(num_coeffs, domain_size);

    c.bench_function(name, move |bencher| {
        let p = DensePolynomial::<F>::rand(4 * num_coeffs, rng);

        bencher.iter(|| {
            add_to_trace!(
                || format!("****************{}*******************", num_coeffs),
                || format!(
                    "--->START TIMESTAMP: {:?}",
                    SystemTime::now()
                        .duration_since(UNIX_EPOCH)
                        .unwrap()
                        .as_secs()
                )
            );

            let _ans1 = p.divide_by_vanishing_poly(&domain.clone());

            add_to_trace!(
                || format!("****************{}*******************", num_coeffs),
                || format!(
                    "--->END TIMESTAMP: {:?}",
                    SystemTime::now()
                        .duration_since(UNIX_EPOCH)
                        .unwrap()
                        .as_secs()
                )
            );
        })
    });
}

fn bench_dense_poly_divide_by_vanishing_poly_tweedle(c: &mut Criterion) {
    bench_dense_poly_div_by_vanishing_poly::<Fr, BasicRadix2Domain<Fr>>(
        c,
        1 << 14,
        "radix-2 DensePolynomial::div by vanishing poly - 2^14 - tweedle",
    );
    bench_dense_poly_div_by_vanishing_poly::<Fr, BasicRadix2Domain<Fr>>(
        c,
        1 << 15,
        "radix-2 DensePolynomial::div by vanishing poly - 2^15 - tweedle",
    );
    bench_dense_poly_div_by_vanishing_poly::<Fr, BasicRadix2Domain<Fr>>(
        c,
        1 << 16,
        "radix-2 DensePolynomial::div by vanishing poly - 2^16 - tweedle",
    );
    bench_dense_poly_div_by_vanishing_poly::<Fr, BasicRadix2Domain<Fr>>(
        c,
        1 << 17,
        "radix-2 DensePolynomial::div by vanishing poly - 2^17 - tweedle",
    );
    bench_dense_poly_div_by_vanishing_poly::<Fr, BasicRadix2Domain<Fr>>(
        c,
        1 << 18,
        "radix-2 DensePolynomial::div by vanishing poly - 2^18 - tweedle",
    );
    bench_dense_poly_div_by_vanishing_poly::<Fr, BasicRadix2Domain<Fr>>(
        c,
        1 << 19,
        "radix-2 DensePolynomial::div by vanishing poly - 2^19 - tweedle",
    );
    bench_dense_poly_div_by_vanishing_poly::<Fr, BasicRadix2Domain<Fr>>(
        c,
        1 << 20,
        "radix-2 DensePolynomial::div by vanishing poly - 2^20 - tweedle",
    );
    bench_dense_poly_div_by_vanishing_poly::<Fr, BasicRadix2Domain<Fr>>(
        c,
        1 << 21,
        "radix-2 DensePolynomial::div by vanishing poly - 2^21 - tweedle",
    );
    bench_dense_poly_div_by_vanishing_poly::<Fr, BasicRadix2Domain<Fr>>(
        c,
        1 << 22,
        "radix-2 DensePolynomial::div by vanishing poly - 2^22 - tweedle",
    );
    bench_dense_poly_div_by_vanishing_poly::<Fr, BasicRadix2Domain<Fr>>(
        c,
        1 << 23,
        "radix-2 DensePolynomial::div by vanishing poly - 2^23 - tweedle",
    );
}

criterion_group! {
    name = radix_2_fft;
    config = Criterion::default().sample_size(10);
    targets = bench_fft_tweedle, bench_ifft_tweedle, bench_dense_poly_mul_tweedle, bench_dense_poly_divide_by_vanishing_poly_tweedle
}

criterion_main!(radix_2_fft);
