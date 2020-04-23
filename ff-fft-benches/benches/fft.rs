use rand;

#[macro_use]
extern crate criterion;

use algebra::{mnt4_753::Fr as MNT4Fr, mnt6_753::Fr as MNT6Fr, FftField, UniformRand};
use criterion::Criterion;
use ff_fft::{EvaluationDomain, MixedRadixEvaluationDomain, Radix2EvaluationDomain};

fn bench_groth16_ffts<F: FftField, D: EvaluationDomain<F>>(
    c: &mut Criterion,
    num_coeffs: usize,
    name: &'static str,
) {
    // Per benchmark setup
    let rng = &mut rand::thread_rng();

    // We expect the num_coeffs input to be a compatible size for the domain.
    let domain = D::new(num_coeffs).unwrap();
    let domain_size = domain.size();
    assert_eq!(num_coeffs, domain_size);

    c.bench_function(name, move |bencher| {
        // Per sample setup
        let mut a = Vec::<F>::with_capacity(num_coeffs);
        let mut b = Vec::<F>::with_capacity(num_coeffs);
        let mut c = Vec::<F>::with_capacity(num_coeffs);
        for _ in 0..num_coeffs {
            a.push(UniformRand::rand(rng));
            b.push(UniformRand::rand(rng));
            c.push(UniformRand::rand(rng));
        }

        bencher.iter(|| {
            // Emulate the FFT operations Groth16 performs in a call to `witness_map`.
            domain.ifft_in_place(&mut a);
            domain.ifft_in_place(&mut b);

            domain.coset_fft_in_place(&mut a);
            domain.coset_fft_in_place(&mut b);

            let mut ab = domain.mul_polynomials_in_evaluation_domain(&a, &b);

            domain.ifft_in_place(&mut c);
            domain.coset_fft_in_place(&mut c);

            domain.divide_by_vanishing_poly_on_coset_in_place(&mut ab);

            domain.coset_ifft_in_place(&mut ab);
        })
    });
}

fn bench_groth16_ffts_radix2(c: &mut Criterion) {
    // Choose 2^16 = 65,536 coefficients for the radix-2 FFT.
    // This comes closest to the 51,200 coefficients chosen in the mixed-radix FFT.
    bench_groth16_ffts::<MNT4Fr, Radix2EvaluationDomain<MNT4Fr>>(c, 1 << 16, "radix-2 FFT");
}

fn bench_groth16_ffts_mixed_radix(c: &mut Criterion) {
    // Choose 2^11*5^2 = 51,200 coefficients for the mixed-radix FFT.
    // This is above the maximum of 32,768 coefficients for the radix-2 FFT.
    bench_groth16_ffts::<MNT6Fr, MixedRadixEvaluationDomain<MNT6Fr>>(c, 51200, "mixed-radix FFT");
}

criterion_group! {
    name = radix_2;
    config = Criterion::default().sample_size(10);
    targets = bench_groth16_ffts_radix2
}

criterion_group! {
    name = mixed_radix;
    config = Criterion::default().sample_size(10);
    targets = bench_groth16_ffts_mixed_radix
}

criterion_main!(radix_2, mixed_radix);
