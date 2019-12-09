use algebra::UniformRand;
use rand::SeedableRng;
use rand_xorshift::XorShiftRng;
use std::ops::{AddAssign, MulAssign, SubAssign};

use algebra::{fields::mnt4753::fq4::Fq4, Field};

#[bench]
fn bench_fq4_add_assign(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let v: Vec<(Fq4, Fq4)> = (0..SAMPLES)
        .map(|_| (Fq4::rand(&mut rng), Fq4::rand(&mut rng)))
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
fn bench_fq4_sub_assign(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let v: Vec<(Fq4, Fq4)> = (0..SAMPLES)
        .map(|_| (Fq4::rand(&mut rng), Fq4::rand(&mut rng)))
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
fn bench_fq4_mul_assign(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let v: Vec<(Fq4, Fq4)> = (0..SAMPLES)
        .map(|_| (Fq4::rand(&mut rng), Fq4::rand(&mut rng)))
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
fn bench_fq4_double(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let v: Vec<Fq4> = (0..SAMPLES).map(|_| Fq4::rand(&mut rng)).collect();

    let mut count = 0;
    b.iter(|| {
        let mut tmp = v[count];
        tmp.double_in_place();
        count = (count + 1) % SAMPLES;
        tmp
    });
}

#[bench]
fn bench_fq4_square(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let v: Vec<Fq4> = (0..SAMPLES).map(|_| Fq4::rand(&mut rng)).collect();

    let mut count = 0;
    b.iter(|| {
        let mut tmp = v[count];
        tmp.square_in_place();
        count = (count + 1) % SAMPLES;
        tmp
    });
}

#[bench]
fn bench_fq4_inverse(b: &mut ::test::Bencher) {
    const SAMPLES: usize = 1000;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let v: Vec<Fq4> = (0..SAMPLES).map(|_| Fq4::rand(&mut rng)).collect();

    let mut count = 0;
    b.iter(|| {
        let tmp = v[count].inverse();
        count = (count + 1) % SAMPLES;
        tmp
    });
}