use algebra::UniformRand;
use rand::SeedableRng;
use rand_xorshift::XorShiftRng;

use algebra::{
    biginteger::{BigInteger, BigInteger768 as FrRepr},
    fields::{mnt6753::fr::Fr, Field, PrimeField, SquareRootField},
};
use std::ops::{AddAssign, MulAssign, SubAssign};

#[bench]
fn bench_fr_repr_add_nocarry(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let v: Vec<(FrRepr, FrRepr)> = (0..SAMPLES)
        .map(|_| {
            let mut tmp1 = FrRepr::rand(&mut rng);
            let mut tmp2 = FrRepr::rand(&mut rng);
            // Shave a few bits off to avoid overflow.
            for _ in 0..3 {
                tmp1.div2();
                tmp2.div2();
            }
            (tmp1, tmp2)
        })
        .collect();

    let mut count = 0;
    b.iter(|| {
        let mut tmp = v[count].0;
        tmp.add_nocarry(&v[count].1);
        count = (count + 1) % SAMPLES;
        tmp
    });
}

#[bench]
fn bench_fr_repr_sub_noborrow(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let v: Vec<(FrRepr, FrRepr)> = (0..SAMPLES)
        .map(|_| {
            let tmp1 = FrRepr::rand(&mut rng);
            let mut tmp2 = tmp1;
            // Ensure tmp2 is smaller than tmp1.
            for _ in 0..10 {
                tmp2.div2();
            }
            (tmp1, tmp2)
        })
        .collect();

    let mut count = 0;
    b.iter(|| {
        let mut tmp = v[count].0;
        tmp.sub_noborrow(&v[count].1);
        count = (count + 1) % SAMPLES;
        tmp
    });
}

#[bench]
fn bench_fr_repr_num_bits(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let v: Vec<FrRepr> = (0..SAMPLES).map(|_| FrRepr::rand(&mut rng)).collect();

    let mut count = 0;
    b.iter(|| {
        let tmp = v[count].num_bits();
        count = (count + 1) % SAMPLES;
        tmp
    });
}

#[bench]
fn bench_fr_repr_mul2(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let v: Vec<FrRepr> = (0..SAMPLES).map(|_| FrRepr::rand(&mut rng)).collect();

    let mut count = 0;
    b.iter(|| {
        let mut tmp = v[count];
        tmp.mul2();
        count = (count + 1) % SAMPLES;
        tmp
    });
}

#[bench]
fn bench_fr_repr_div2(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let v: Vec<FrRepr> = (0..SAMPLES).map(|_| FrRepr::rand(&mut rng)).collect();

    let mut count = 0;
    b.iter(|| {
        let mut tmp = v[count];
        tmp.div2();
        count = (count + 1) % SAMPLES;
        tmp
    });
}

#[bench]
fn bench_fr_add_assign(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let v: Vec<(Fr, Fr)> = (0..SAMPLES)
        .map(|_| (Fr::rand(&mut rng), Fr::rand(&mut rng)))
        .collect();

    let mut count = 0;
    b.iter(|| {
        let mut tmp = v[count].0;
        tmp.add_assign(&v[count].1);
        count = (count + 1) % SAMPLES;
        tmp
    });
}

#[bench]
fn bench_fr_sub_assign(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let v: Vec<(Fr, Fr)> = (0..SAMPLES)
        .map(|_| (Fr::rand(&mut rng), Fr::rand(&mut rng)))
        .collect();

    let mut count = 0;
    b.iter(|| {
        let mut tmp = v[count].0;
        tmp.sub_assign(&v[count].1);
        count = (count + 1) % SAMPLES;
        tmp
    });
}

#[bench]
fn bench_fr_mul_assign(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let v: Vec<(Fr, Fr)> = (0..SAMPLES)
        .map(|_| (Fr::rand(&mut rng), Fr::rand(&mut rng)))
        .collect();

    let mut count = 0;
    b.iter(|| {
        let mut tmp = v[count].0;
        tmp.mul_assign(&v[count].1);
        count = (count + 1) % SAMPLES;
        tmp
    });
}

#[bench]
fn bench_fr_double(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let v: Vec<Fr> = (0..SAMPLES).map(|_| Fr::rand(&mut rng)).collect();

    let mut count = 0;
    b.iter(|| {
        let mut tmp = v[count];
        tmp.double_in_place();
        count = (count + 1) % SAMPLES;
        tmp
    });
}

#[bench]
fn bench_fr_square(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let v: Vec<Fr> = (0..SAMPLES).map(|_| Fr::rand(&mut rng)).collect();

    let mut count = 0;
    b.iter(|| {
        let mut tmp = v[count];
        tmp.square_in_place();
        count = (count + 1) % SAMPLES;
        tmp
    });
}

#[bench]
fn bench_fr_inverse(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let v: Vec<Fr> = (0..SAMPLES).map(|_| Fr::rand(&mut rng)).collect();

    let mut count = 0;
    b.iter(|| {
        count = (count + 1) % SAMPLES;
        v[count].inverse()
    });
}

#[bench]
fn bench_fr_negate(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let v: Vec<Fr> = (0..SAMPLES).map(|_| Fr::rand(&mut rng)).collect();

    let mut count = 0;
    b.iter(|| {
        let mut tmp = v[count];
        tmp = -tmp;
        count = (count + 1) % SAMPLES;
        tmp
    });
}

#[bench]
fn bench_fr_sqrt(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let v: Vec<Fr> = (0..SAMPLES)
        .map(|_| {
            let mut tmp = Fr::rand(&mut rng);
            tmp.square_in_place();
            tmp
        })
        .collect();

    let mut count = 0;
    b.iter(|| {
        count = (count + 1) % SAMPLES;
        v[count].sqrt()
    });
}

#[bench]
fn bench_fr_into_repr(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let v: Vec<Fr> = (0..SAMPLES).map(|_| Fr::rand(&mut rng)).collect();

    let mut count = 0;
    b.iter(|| {
        count = (count + 1) % SAMPLES;
        v[count].into_repr()
    });
}

#[bench]
fn bench_fr_from_repr(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let v: Vec<FrRepr> = (0..SAMPLES)
        .map(|_| Fr::rand(&mut rng).into_repr())
        .collect();

    let mut count = 0;
    b.iter(|| {
        count = (count + 1) % SAMPLES;
        Fr::from_repr(v[count])
    });
}

// The two benchmarks below prove that, for a given k, if the evaluation size is greater that 2^k,
// a MixedRadix2Domain is more convenient than a BasicRadix2Domain because the last one increases
// with power of two and will work on a domain of size 2^(k+1), while the latter will work on a
// domain of intermediate size because it increases as 2^n * p^q: the reduced domain size overcomes
// the loss in performance you get from using a MixedRadix2Domain algorithm for FFT.

use algebra::fft::{EvaluationDomain, EvaluationDomainImpl};

#[bench]
fn bench_basic_domain_fft_over_fr(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let num_coeffs = 1 << 14;
    let domain = EvaluationDomain::<Fr>::new(num_coeffs).unwrap();
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
        let mut ab = domain.mul_polynomials_in_evaluation_domain(&a, &b);
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
fn bench_mixed_domain_fft_over_fr(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let num_coeffs = 10240;
    let domain = EvaluationDomain::<Fr>::new(num_coeffs).unwrap();
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
        let mut ab = domain.mul_polynomials_in_evaluation_domain(&a, &b);
        drop(a);
        drop(b);
        let mut c = v_c[count].clone();
        domain.ifft_in_place(&mut c);
        domain.coset_fft_in_place(&mut c);
        domain.divide_by_vanishing_poly_on_coset_in_place(&mut ab);
        domain.coset_ifft_in_place(&mut ab);
    });
}