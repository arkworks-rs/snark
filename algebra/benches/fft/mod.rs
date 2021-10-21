use algebra::{
    fields::mnt6753::Fr,
    fft::get_best_evaluation_domain,
    UniformRand,
};
use rand::SeedableRng;
use rand_xorshift::XorShiftRng;

#[bench]
fn bench_basic_domain_fft(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let num_coeffs = 1 << 14;
    let domain = get_best_evaluation_domain::<Fr>(num_coeffs).unwrap();
    let domain_size = domain.size();
    assert_eq!(num_coeffs, domain_size);

    let v_a: Vec<Vec<Fr>> = (0..SAMPLES)
        .map(|_| {
            let mut v = vec![];
            for _ in 0..domain_size {
                v.push(Fr::rand(&mut rng));
            }
            v
        }).collect();

    let v_b: Vec<Vec<Fr>> = (0..SAMPLES)
        .map(|_| {
            let mut v = vec![];
            for _ in 0..domain_size {
                v.push(Fr::rand(&mut rng));
            }
            v
        }).collect();

    let v_c: Vec<Vec<Fr>> = (0..SAMPLES)
        .map(|_| {
            let mut v = vec![];
            for _ in 0..domain_size {
                v.push(Fr::rand(&mut rng));
            }
            v
        }).collect();

    let mut count = 0;
    b.iter(|| {
        count = (count + 1) % SAMPLES;
        let mut a = v_a[count].clone();
        let mut b = v_b[count].clone();
        domain.ifft_in_place(&mut a);
        domain.ifft_in_place(&mut b);
        domain.coset_fft_in_place(&mut a);
        domain.coset_fft_in_place(&mut b);
        let mut ab = domain.mul_polynomials_in_evaluation_domain(&a, &b).unwrap();
        drop(a);
        drop(b);
        let mut c = v_c[count].clone();
        domain.ifft_in_place(&mut c);
        domain.coset_fft_in_place(&mut c);
        domain.divide_by_vanishing_poly_on_coset_in_place(&mut ab);
        domain.coset_ifft_in_place(&mut ab);
    });
}

#[bench]
fn bench_mixed_domain_fft(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let num_coeffs = 10240;
    let domain = get_best_evaluation_domain::<Fr>(num_coeffs).unwrap();
    let domain_size = domain.size();
    assert_eq!(num_coeffs, domain_size);

    let v_a: Vec<Vec<Fr>> = (0..SAMPLES)
        .map(|_| {
            let mut v = vec![];
            for _ in 0..domain_size {
                v.push(Fr::rand(&mut rng));
            }
            v
        }).collect();

    let v_b: Vec<Vec<Fr>> = (0..SAMPLES)
        .map(|_| {
            let mut v = vec![];
            for _ in 0..domain_size {
                v.push(Fr::rand(&mut rng));
            }
            v
        }).collect();

    let v_c: Vec<Vec<Fr>> = (0..SAMPLES)
        .map(|_| {
            let mut v = vec![];
            for _ in 0..domain_size {
                v.push(Fr::rand(&mut rng));
            }
            v
        }).collect();

    let mut count = 0;
    b.iter(|| {
        count = (count + 1) % SAMPLES;
        let mut a = v_a[count].clone();
        let mut b = v_b[count].clone();
        domain.ifft_in_place(&mut a);
        domain.ifft_in_place(&mut b);
        domain.coset_fft_in_place(&mut a);
        domain.coset_fft_in_place(&mut b);
        let mut ab = domain.mul_polynomials_in_evaluation_domain(&a, &b).unwrap();
        drop(a);
        drop(b);
        let mut c = v_c[count].clone();
        domain.ifft_in_place(&mut c);
        domain.coset_fft_in_place(&mut c);
        domain.divide_by_vanishing_poly_on_coset_in_place(&mut ab);
        domain.coset_ifft_in_place(&mut ab);
    });
}