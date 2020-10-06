#![allow(unused)]
use algebra_core::{
    batch_bucketed_add, batch_verify_in_subgroup,
    biginteger::BigInteger64,
    curves::{AffineCurve, BatchGroupArithmeticSlice, ProjectiveCurve},
    io::Cursor,
    BucketPosition, CanonicalDeserialize, CanonicalSerialize, Field, MontgomeryModelParameters,
    One, PrimeField, SWFlags, SWModelParameters, SerializationError, TEModelParameters,
    UniformRand, Vec, VerificationError, Zero,
};
use rand::{
    distributions::{Distribution, Uniform},
    SeedableRng,
};
use rand_xorshift::XorShiftRng;

use std::ops::Neg;

use crate::tests::helpers::create_pseudo_uniform_random_elems;

use crate::cfg_chunks_mut;
#[cfg(any(feature = "parallel"))]
use rayon::prelude::*;

pub const AFFINE_BATCH_SIZE: usize = 4096;
pub const ITERATIONS: usize = 10;

fn random_addition_test<G: ProjectiveCurve>() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..ITERATIONS {
        let a = G::rand(&mut rng);
        let b = G::rand(&mut rng);
        let c = G::rand(&mut rng);
        let a_affine = a.into_affine();
        let b_affine = b.into_affine();
        let c_affine = c.into_affine();

        // a + a should equal the doubling
        {
            let mut aplusa = a;
            aplusa.add_assign(&a);

            let mut aplusamixed = a;
            aplusamixed.add_assign_mixed(&a.into_affine());

            let mut adouble = a;
            adouble.double_in_place();

            assert_eq!(aplusa, adouble);
            assert_eq!(aplusa, aplusamixed);
        }

        let mut tmp = vec![G::zero(); 6];

        // (a + b) + c
        tmp[0] = (a + &b) + &c;

        // a + (b + c)
        tmp[1] = a + &(b + &c);

        // (a + c) + b
        tmp[2] = (a + &c) + &b;

        // Mixed addition

        // (a + b) + c
        tmp[3] = a_affine.into_projective();
        tmp[3].add_assign_mixed(&b_affine);
        tmp[3].add_assign_mixed(&c_affine);

        // a + (b + c)
        tmp[4] = b_affine.into_projective();
        tmp[4].add_assign_mixed(&c_affine);
        tmp[4].add_assign_mixed(&a_affine);

        // (a + c) + b[G]: BatchArithmetic
        tmp[5] = a_affine.into_projective();
        tmp[5].add_assign_mixed(&c_affine);
        tmp[5].add_assign_mixed(&b_affine);

        // Comparisons
        for i in 0..6 {
            for j in 0..6 {
                if tmp[i] != tmp[j] {
                    println!("{} \n{}", tmp[i], tmp[j]);
                }
                assert_eq!(tmp[i], tmp[j], "Associativity failed {} {}", i, j);
                assert_eq!(
                    tmp[i].into_affine(),
                    tmp[j].into_affine(),
                    "Associativity failed"
                );
            }

            assert!(tmp[i] != a);
            assert!(tmp[i] != b);
            assert!(tmp[i] != c);

            assert!(a != tmp[i]);
            assert!(b != tmp[i]);
            assert!(c != tmp[i]);
        }
    }
}

fn random_multiplication_test<G: ProjectiveCurve>() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..ITERATIONS {
        let mut a = G::rand(&mut rng);
        let mut b = G::rand(&mut rng);
        let a_affine = a.into_affine();
        let b_affine = b.into_affine();

        let s = G::ScalarField::rand(&mut rng);

        // s ( a + b )
        let mut tmp1 = a;
        tmp1.add_assign(&b);
        tmp1.mul_assign(s);

        // sa + sb
        a.mul_assign(s);
        b.mul_assign(s);

        let mut tmp2 = a;
        tmp2.add_assign(&b);

        // Affine multiplication
        let mut tmp3 = a_affine.mul(s.into_repr());
        tmp3.add_assign(&b_affine.mul(s.into_repr()));

        assert_eq!(tmp1, tmp2);
        assert_eq!(tmp1, tmp3);
    }
}

fn random_doubling_test<G: ProjectiveCurve>() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..ITERATIONS {
        let mut a = G::rand(&mut rng);
        let mut b = G::rand(&mut rng);

        // 2(a + b)
        let mut tmp1 = a;
        tmp1.add_assign(&b);
        tmp1.double_in_place();

        // 2a + 2b
        a.double_in_place();
        b.double_in_place();

        let mut tmp2 = a;
        tmp2.add_assign(&b);

        let mut tmp3 = a;
        tmp3.add_assign_mixed(&b.into_affine());

        assert_eq!(tmp1, tmp2);
        assert_eq!(tmp1, tmp3);
    }
}

fn random_negation_test<G: ProjectiveCurve>() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..ITERATIONS {
        let r = G::rand(&mut rng);

        let s = G::ScalarField::rand(&mut rng);
        let sneg = -s;
        assert!((s + &sneg).is_zero());

        let mut t1 = r;
        t1.mul_assign(s);

        let mut t2 = r;
        t2.mul_assign(sneg);

        let mut t3 = t1;
        t3.add_assign(&t2);
        assert!(t3.is_zero());

        let mut t4 = t1;
        t4.add_assign_mixed(&t2.into_affine());
        assert!(t4.is_zero());

        t1 = -t1;
        assert_eq!(t1, t2);
    }
}

fn random_transformation_test<G: ProjectiveCurve>() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..ITERATIONS {
        let g = G::rand(&mut rng);
        let g_affine = g.into_affine();
        let g_projective = g_affine.into_projective();
        assert_eq!(g, g_projective);
    }

    // Batch normalization
    for _ in 0..10 {
        let mut v = (0..ITERATIONS)
            .map(|_| G::rand(&mut rng))
            .collect::<Vec<_>>();

        for i in &v {
            assert!(!i.is_normalized());
        }

        use rand::distributions::{Distribution, Uniform};
        let between = Uniform::from(0..ITERATIONS);
        // Sprinkle in some normalized points
        for _ in 0..5 {
            v[between.sample(&mut rng)] = G::zero();
        }
        for _ in 0..5 {
            let s = between.sample(&mut rng);
            v[s] = v[s].into_affine().into_projective();
        }

        let expected_v = v
            .iter()
            .map(|v| v.into_affine().into_projective())
            .collect::<Vec<_>>();
        G::batch_normalization(&mut v);

        for i in &v {
            assert!(i.is_normalized());
        }

        assert_eq!(v, expected_v);
    }
}

pub fn random_batch_doubling_test<G: ProjectiveCurve>() {
    use algebra_core::curves::models::short_weierstrass_jacobian::{GroupAffine, GroupProjective};
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for j in 0..ITERATIONS {
        let size = std::cmp::min(1 << 7, 1 << (j + 5));
        let mut a = Vec::with_capacity(size);
        let mut b = Vec::with_capacity(size);

        for i in 0..size {
            a.push(G::rand(&mut rng));
            b.push(G::rand(&mut rng));
        }

        let mut c = a.clone();

        let mut a: Vec<G::Affine> = a.iter().map(|p| p.into_affine()).collect();

        a[..].batch_double_in_place(&(0..size).map(|x| x as u32).collect::<Vec<_>>()[..]);

        for p_c in c.iter_mut() {
            *p_c.double_in_place();
        }

        let c: Vec<G::Affine> = c.iter().map(|p| p.into_affine()).collect();

        assert_eq!(a, c);
    }
}

pub fn random_batch_addition_test<G: ProjectiveCurve>() {
    use algebra_core::curves::models::short_weierstrass_jacobian::{GroupAffine, GroupProjective};
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for j in 0..ITERATIONS {
        let size = std::cmp::min(1 << 7, 1 << (j + 5));
        let mut a = Vec::with_capacity(size);
        let mut b = Vec::with_capacity(size);

        for i in 0..size {
            a.push(G::rand(&mut rng));
            b.push(G::rand(&mut rng));
        }

        let mut c = a.clone();
        let mut d = b.clone();

        let mut a: Vec<G::Affine> = a.iter().map(|p| p.into_affine()).collect();
        let mut b: Vec<G::Affine> = b.iter().map(|p| p.into_affine()).collect();

        a[..].batch_add_in_place(
            &mut b[..],
            &(0..size).map(|x| (x as u32, x as u32)).collect::<Vec<_>>()[..],
        );

        for (p_c, p_d) in c.iter_mut().zip(d.iter()) {
            *p_c += *p_d;
        }

        let c: Vec<G::Affine> = c.iter().map(|p| p.into_affine()).collect();

        assert_eq!(a, c);
    }
}

pub fn random_batch_add_doubling_test<G: ProjectiveCurve>() {
    use algebra_core::curves::models::short_weierstrass_jacobian::{GroupAffine, GroupProjective};
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for j in 0..ITERATIONS {
        let size = std::cmp::min(1 << 7, 1 << (j + 5));
        let mut a = Vec::<G>::with_capacity(size);
        let mut b = Vec::<G>::with_capacity(size);

        for i in 0..size {
            a.push(G::rand(&mut rng));
        }

        let mut b = a.clone();
        let mut c = a.clone();
        let mut d = b.clone();

        let mut a: Vec<G::Affine> = a.iter().map(|p| p.into_affine()).collect();
        let mut b: Vec<G::Affine> = b.iter().map(|p| p.into_affine()).collect();

        a[..].batch_add_in_place(
            &mut b[..],
            &(0..size).map(|x| (x as u32, x as u32)).collect::<Vec<_>>()[..],
        );

        for (p_c, p_d) in c.iter_mut().zip(d.iter()) {
            *p_c += *p_d;
        }

        let c: Vec<G::Affine> = c.iter().map(|p| p.into_affine()).collect();

        assert_eq!(a, c);
    }
}

pub fn random_batch_scalar_mul_test<G: ProjectiveCurve>() {
    use algebra_core::curves::models::short_weierstrass_jacobian::{GroupAffine, GroupProjective};
    use std::ops::MulAssign;
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for j in 0..ITERATIONS {
        let size = std::cmp::min(1 << 7, 1 << (j + 4));
        let mut a = Vec::with_capacity(size);
        let mut s = Vec::with_capacity(size);

        for i in 0..size {
            a.push(G::rand(&mut rng));
            s.push(G::ScalarField::rand(&mut rng));
        }

        let mut c = a.clone();
        let mut t = s.clone();

        let mut a: Vec<G::Affine> = a.iter().map(|p| p.into_affine()).collect();

        let mut s: Vec<<G::ScalarField as PrimeField>::BigInt> =
            s.iter().map(|p| p.into_repr()).collect();

        let now = std::time::Instant::now();
        a[..].batch_scalar_mul_in_place::<<G::ScalarField as PrimeField>::BigInt>(&mut s[..], 4);
        println!(
            "Batch affine mul for {} elems: {}us",
            size,
            now.elapsed().as_micros()
        );

        let now = std::time::Instant::now();
        for (p_c, s_t) in c.iter_mut().zip(t.iter()) {
            p_c.mul_assign(*s_t);
        }
        println!(
            "Proj mul for {} elems: {}us",
            size,
            now.elapsed().as_micros()
        );

        let c: Vec<G::Affine> = c.iter().map(|p| p.into_affine()).collect();

        for (p1, p2) in a.iter().zip(c) {
            assert_eq!(*p1, p2);
        }
    }
}

fn batch_bucketed_add_test<C: AffineCurve>() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    #[cfg(not(feature = "big_n"))]
    const MAX_LOGN: usize = 12;
    #[cfg(feature = "big_n")]
    const MAX_LOGN: usize = 22;

    let random_elems = create_pseudo_uniform_random_elems(&mut rng, MAX_LOGN);

    for i in (MAX_LOGN - 4)..(ITERATIONS / 2 + MAX_LOGN - 4) {
        let n_elems = 1 << i;
        let n_buckets = 1 << (i - 3);

        let mut bucket_assign = Vec::<_>::with_capacity(n_elems);
        let step = Uniform::new(0, n_buckets);

        for i in 0..n_elems {
            bucket_assign.push(BucketPosition {
                bucket: step.sample(&mut rng) as u32,
                position: i as u32,
            });
        }

        let mut res1 = vec![];
        let mut elems_mut = random_elems[0..n_elems].to_vec();
        let now = std::time::Instant::now();
        res1 = batch_bucketed_add::<C>(
            n_buckets,
            &mut elems_mut[..],
            &mut bucket_assign.to_vec()[..],
        );
        println!(
            "batch bucketed add for {} elems: {:?}",
            n_elems,
            now.elapsed().as_micros()
        );

        let mut res2 = vec![C::Projective::zero(); n_buckets];
        let mut elems = random_elems[0..n_elems].to_vec();

        let now = std::time::Instant::now();
        for (&bucket_idx, elem) in bucket_assign.iter().zip(elems) {
            res2[bucket_idx.bucket as usize].add_assign_mixed(&elem);
        }
        println!(
            "bucketed add for {} elems: {:?}",
            n_elems,
            now.elapsed().as_micros()
        );

        let res1: Vec<C::Projective> = res1.iter().map(|&p| p.into()).collect();

        for (i, (p1, p2)) in res1.iter().zip(res2).enumerate() {
            assert_eq!(*p1, p2);
        }
    }
}

macro_rules! batch_verify_test {
    ($P: ident, $GroupAffine: ident, $GroupProjective: ident) => {
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        #[cfg(not(feature = "big_n"))]
        const MAX_LOGN: usize = 14;
        #[cfg(feature = "big_n")]
        const MAX_LOGN: usize = 22;

        const SECURITY_PARAM: usize = 128;
        // Generate pseudorandom group elements
        let random_elems: Vec<$GroupAffine<P>> = create_pseudo_uniform_random_elems(&mut rng, MAX_LOGN);

        let now = std::time::Instant::now();
        let mut non_subgroup_points = Vec::with_capacity(1 << 10);
        while non_subgroup_points.len() < 1 << 10 {
            if let Some(elem) = $GroupAffine::<P>::get_point_from_x($P::BaseField::rand(&mut rng), false)
            {
                // If the cofactor is small, with non-negligible probability the sampled point
                // is in the group, so we should check it isn't. Else we don't waste compute.
                if $P::COFACTOR[1..].iter().all(|&x| x == 0u64) {
                    if !elem.is_in_correct_subgroup_assuming_on_curve() {
                        non_subgroup_points.push(elem);
                    }
                } else {
                    non_subgroup_points.push(elem);
                }
            }
        }
        println!(
            "Generate non-subgroup points: {:?}",
            now.elapsed().as_micros()
        );

        println!("Security Param: {}", SECURITY_PARAM);
        let mut estimated_timing = 0;
        for i in (MAX_LOGN - 4)..(ITERATIONS / 2 + MAX_LOGN - 4) {
            let n_elems = 1 << i;
            println!("n: {}", n_elems);

            if i == MAX_LOGN - 4 {
                let mut tmp_elems_for_naive = random_elems[0..n_elems].to_vec();
                let now = std::time::Instant::now();
                cfg_chunks_mut!(tmp_elems_for_naive, AFFINE_BATCH_SIZE).map(|e| {
                    // Probably could optimise this further: single scalar
                    // We also need to make GLV work with the characteristic
                    let size = e.len();
                    e[..].batch_scalar_mul_in_place::<<<$GroupAffine<P> as AffineCurve>::ScalarField as PrimeField>::BigInt>(
                        &mut vec![<<$GroupAffine<P> as AffineCurve>::ScalarField as PrimeField>::modulus().into(); size][..],
                        4,
                    );
                    e.iter().all(|p| p.is_zero())
                })
                .all(|b| b);

                estimated_timing = now.elapsed().as_micros();
                println!(
                    "Success: In Subgroup. n: {}, time: {} (naive)",
                    n_elems,
                    estimated_timing
                );
            } else {
                estimated_timing *= 2;
                println!(
                    "Estimated timing for n: {}, time: {} (naive)",
                    n_elems,
                    estimated_timing
                );
            }

            let random_location = Uniform::new(0, n_elems);
            let mut tmp_elems = random_elems[0..n_elems].to_vec();

            let now = std::time::Instant::now();
            batch_verify_in_subgroup::<$GroupAffine<P>, XorShiftRng>(&tmp_elems[..], SECURITY_PARAM, &mut rng)
                .expect("Should have verified as correct");
            println!(
                "Success: In Subgroup. n: {}, time: {}",
                n_elems,
                now.elapsed().as_micros()
            );

            for j in 0..10 {
                // Randomly insert random non-subgroup elems
                for k in 0..(1 << j) {
                    tmp_elems[random_location.sample(&mut rng)] = non_subgroup_points[k];
                }
                let now = std::time::Instant::now();
                match batch_verify_in_subgroup::<$GroupAffine<P>, XorShiftRng>(&tmp_elems[..], SECURITY_PARAM, &mut rng) {
                    Ok(_) => assert!(false, "did not detect non-subgroup elems"),
                    _ => assert!(true),
                };
                println!(
                    "Success: Not in subgroup. n: {}, non-subgroup elems: {}, time: {}",
                    n_elems,
                    (1 << (j + 1)) - 1,
                    now.elapsed().as_micros()
                );
            }
        }

        // // We can induce a collision and thus failure to identify non-subgroup elements with the following
        // // for small security parameters. This is a non-deterministic "anti-test" that should fail and cause
        // // panic. It is meant for sanity checking.
        // for j in 0..10000 {
        //     // Randomly insert random non-subgroup elems
        //     if j == 0 {
        //         for _ in 0..(1 << j) {
        //             loop {
        //                 if let Some(non_subgroup_elem) =
        //                     GroupAffine::<P>::get_point_from_x(P::BaseField::rand(&mut rng), false)
        //                 {
        //                     tmp_elems[random_location.sample(&mut rng)] = non_subgroup_elem;
        //                     tmp_elems[random_location.sample(&mut rng) + 1] = non_subgroup_elem.neg();
        //                     break;
        //                 }
        //             }
        //         }
        //     }
        //     let now = std::time::Instant::now();
        //     match batch_verify_in_subgroup::<GroupAffine<P>>(&tmp_elems[..], SECURITY_PARAM) {
        //         Ok(_) => assert!(false, "did not detect non-subgroup elems"),
        //         _ => assert!(true),
        //     };
        //     println!(
        //         "Success: Not in subgroup. n: {}, non-subgroup elems: {}, time: {}",
        //         n_elems,
        //         (1 << (j + 1)) - 1,
        //         now.elapsed().as_micros()
        //     );
        // }
    }
}

fn sw_batch_verify_test<P: SWModelParameters>() {
    use algebra_core::curves::models::short_weierstrass_jacobian::{GroupAffine, GroupProjective};
    batch_verify_test!(P, GroupAffine, GroupProjective);
}

fn te_batch_verify_test<P: TEModelParameters>() {
    use algebra_core::curves::models::twisted_edwards_extended::{GroupAffine, GroupProjective};
    batch_verify_test!(P, GroupAffine, GroupProjective);
}

pub fn curve_tests<G: ProjectiveCurve>() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    // Negation edge case with zero.
    {
        let z = -G::zero();
        assert!(z.is_zero());
    }

    // Doubling edge case with zero.
    {
        let mut z = -G::zero();
        z.double_in_place();
        assert!(z.is_zero());
    }

    // Addition edge cases with zero
    {
        let mut r = G::rand(&mut rng);
        let rcopy = r;
        r.add_assign(&G::zero());
        assert_eq!(r, rcopy);
        r.add_assign_mixed(&G::Affine::zero());
        assert_eq!(r, rcopy);

        let mut z = G::zero();
        z.add_assign(&G::zero());
        assert!(z.is_zero());
        z.add_assign_mixed(&G::Affine::zero());
        assert!(z.is_zero());

        let mut z2 = z;
        z2.add_assign(&r);

        z.add_assign_mixed(&r.into_affine());

        assert_eq!(z, z2);
        assert_eq!(z, r);
    }

    // Transformations
    {
        let a = G::rand(&mut rng);
        let b = a.into_affine().into_projective();
        let c = a
            .into_affine()
            .into_projective()
            .into_affine()
            .into_projective();
        assert_eq!(a, b);
        assert_eq!(b, c);
    }

    // Test COFACTOR and COFACTOR_INV
    {
        let a = G::rand(&mut rng);
        let b = a.into_affine();
        let c = b.mul_by_cofactor_inv().mul_by_cofactor();
        assert_eq!(b, c);
    }

    random_addition_test::<G>();
    random_multiplication_test::<G>();
    random_doubling_test::<G>();
    random_negation_test::<G>();
    random_transformation_test::<G>();
}

pub fn batch_affine_test<G: ProjectiveCurve>() {
    random_batch_doubling_test::<G>();
    random_batch_add_doubling_test::<G>();
    random_batch_addition_test::<G>();
    random_batch_scalar_mul_test::<G>();
    batch_bucketed_add_test::<G::Affine>();
}

pub fn sw_tests<P: SWModelParameters>() {
    #[cfg(feature = "serialisation")]
    sw_curve_serialization_test::<P>();
    #[cfg(feature = "random_bytes")]
    sw_from_random_bytes::<P>();
    // Only check batch verification for non-unit cofactor
    #[cfg(feature = "verify")]
    if !(P::COFACTOR[0] == 1u64 && P::COFACTOR[1..].iter().all(|&x| x == 0u64)) {
        sw_batch_verify_test::<P>();
    }
}

pub fn sw_from_random_bytes<P: SWModelParameters>() {
    use algebra_core::curves::models::short_weierstrass_jacobian::{GroupAffine, GroupProjective};

    let buf_size = GroupAffine::<P>::zero().serialized_size();

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..ITERATIONS {
        let a = GroupProjective::<P>::rand(&mut rng);
        let mut a = a.into_affine();
        {
            let mut serialized = vec![0; buf_size];
            let mut cursor = Cursor::new(&mut serialized[..]);
            a.serialize(&mut cursor).unwrap();

            let mut cursor = Cursor::new(&serialized[..]);
            let p1 = GroupAffine::<P>::deserialize(&mut cursor).unwrap();
            let p2 = GroupAffine::<P>::from_random_bytes(&serialized).unwrap();
            assert_eq!(p1, p2);
        }
    }
}

pub fn sw_curve_serialization_test<P: SWModelParameters>() {
    use algebra_core::curves::models::short_weierstrass_jacobian::{GroupAffine, GroupProjective};

    let buf_size = GroupAffine::<P>::zero().serialized_size();

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..ITERATIONS {
        let a = GroupProjective::<P>::rand(&mut rng);
        let mut a = a.into_affine();
        {
            let mut serialized = vec![0; buf_size];
            let mut cursor = Cursor::new(&mut serialized[..]);
            a.serialize(&mut cursor).unwrap();

            let mut cursor = Cursor::new(&serialized[..]);
            let b = GroupAffine::<P>::deserialize(&mut cursor).unwrap();
            assert_eq!(a, b);
        }

        {
            a.y = -a.y;
            let mut serialized = vec![0; buf_size];
            let mut cursor = Cursor::new(&mut serialized[..]);
            a.serialize(&mut cursor).unwrap();
            let mut cursor = Cursor::new(&serialized[..]);
            let b = GroupAffine::<P>::deserialize(&mut cursor).unwrap();
            assert_eq!(a, b);
        }

        {
            let a = GroupAffine::<P>::zero();
            let mut serialized = vec![0; buf_size];
            let mut cursor = Cursor::new(&mut serialized[..]);
            a.serialize(&mut cursor).unwrap();
            let mut cursor = Cursor::new(&serialized[..]);
            let b = GroupAffine::<P>::deserialize(&mut cursor).unwrap();
            assert_eq!(a, b);
        }

        {
            let a = GroupAffine::<P>::zero();
            let mut serialized = vec![0; buf_size - 1];
            let mut cursor = Cursor::new(&mut serialized[..]);
            a.serialize(&mut cursor).unwrap_err();
        }

        {
            let serialized = vec![0; buf_size - 1];
            let mut cursor = Cursor::new(&serialized[..]);
            GroupAffine::<P>::deserialize(&mut cursor).unwrap_err();
        }

        {
            let mut serialized = vec![0; a.uncompressed_size()];
            let mut cursor = Cursor::new(&mut serialized[..]);
            a.serialize_uncompressed(&mut cursor).unwrap();

            let mut cursor = Cursor::new(&serialized[..]);
            let b = GroupAffine::<P>::deserialize_uncompressed(&mut cursor).unwrap();
            assert_eq!(a, b);
        }

        {
            a.y = -a.y;
            let mut serialized = vec![0; a.uncompressed_size()];
            let mut cursor = Cursor::new(&mut serialized[..]);
            a.serialize_uncompressed(&mut cursor).unwrap();
            let mut cursor = Cursor::new(&serialized[..]);
            let b = GroupAffine::<P>::deserialize_uncompressed(&mut cursor).unwrap();
            assert_eq!(a, b);
        }

        {
            let a = GroupAffine::<P>::zero();
            let mut serialized = vec![0; a.uncompressed_size()];
            let mut cursor = Cursor::new(&mut serialized[..]);
            a.serialize_uncompressed(&mut cursor).unwrap();
            let mut cursor = Cursor::new(&serialized[..]);
            let b = GroupAffine::<P>::deserialize_uncompressed(&mut cursor).unwrap();
            assert_eq!(a, b);
        }
    }
}

pub(crate) fn montgomery_conversion_test<P>()
where
    P: TEModelParameters,
{
    // A = 2 * (a + d) / (a - d)
    let a = P::BaseField::one().double()
        * &(P::COEFF_A + &P::COEFF_D)
        * &(P::COEFF_A - &P::COEFF_D).inverse().unwrap();
    // B = 4 / (a - d)
    let b = P::BaseField::one().double().double() * &(P::COEFF_A - &P::COEFF_D).inverse().unwrap();

    assert_eq!(a, P::MontgomeryModelParameters::COEFF_A);
    assert_eq!(b, P::MontgomeryModelParameters::COEFF_B);
}

pub fn edwards_tests<P: TEModelParameters>()
where
    P::BaseField: PrimeField,
{
    #[cfg(feature = "serialisation")]
    edwards_curve_serialization_test::<P>();
    #[cfg(feature = "random_bytes")]
    edwards_from_random_bytes::<P>();
    // Only check batch verification for non-unit cofactor
    #[cfg(feature = "verify")]
    if !(P::COFACTOR[0] == 1u64 && P::COFACTOR[1..].iter().all(|&x| x == 0u64)) {
        te_batch_verify_test::<P>();
    }
}

pub fn edwards_from_random_bytes<P: TEModelParameters>()
where
    P::BaseField: PrimeField,
{
    use algebra_core::curves::models::twisted_edwards_extended::{GroupAffine, GroupProjective};
    use algebra_core::{to_bytes, ToBytes};

    let buf_size = GroupAffine::<P>::zero().serialized_size();

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..ITERATIONS {
        let a = GroupProjective::<P>::rand(&mut rng);
        let mut a = a.into_affine();
        {
            let mut serialized = vec![0; buf_size];
            let mut cursor = Cursor::new(&mut serialized[..]);
            a.serialize(&mut cursor).unwrap();

            let mut cursor = Cursor::new(&serialized[..]);
            let p1 = GroupAffine::<P>::deserialize(&mut cursor).unwrap();
            let p2 = GroupAffine::<P>::from_random_bytes(&serialized).unwrap();
            assert_eq!(p1, p2);
        }
    }

    for _ in 0..ITERATIONS {
        let mut biginteger =
            <<GroupAffine<P> as AffineCurve>::BaseField as PrimeField>::BigInt::rand(&mut rng);
        let mut bytes = to_bytes![biginteger].unwrap();
        let mut g = GroupAffine::<P>::from_random_bytes(&bytes);
        while g.is_none() {
            bytes.iter_mut().for_each(|i| *i = i.wrapping_sub(1));
            g = GroupAffine::<P>::from_random_bytes(&bytes);
        }
        let _g = g.unwrap();
    }
}

pub fn edwards_curve_serialization_test<P: TEModelParameters>() {
    use algebra_core::curves::models::twisted_edwards_extended::{GroupAffine, GroupProjective};

    let buf_size = GroupAffine::<P>::zero().serialized_size();

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..ITERATIONS {
        let a = GroupProjective::<P>::rand(&mut rng);
        let a = a.into_affine();
        {
            let mut serialized = vec![0; buf_size];
            let mut cursor = Cursor::new(&mut serialized[..]);
            a.serialize(&mut cursor).unwrap();

            let mut cursor = Cursor::new(&serialized[..]);
            let b = GroupAffine::<P>::deserialize(&mut cursor).unwrap();
            assert_eq!(a, b);
        }

        {
            let a = GroupAffine::<P>::zero();
            let mut serialized = vec![0; buf_size];
            let mut cursor = Cursor::new(&mut serialized[..]);
            a.serialize(&mut cursor).unwrap();
            let mut cursor = Cursor::new(&serialized[..]);
            let b = GroupAffine::<P>::deserialize(&mut cursor).unwrap();
            assert_eq!(a, b);
        }

        {
            let a = GroupAffine::<P>::zero();
            let mut serialized = vec![0; buf_size - 1];
            let mut cursor = Cursor::new(&mut serialized[..]);
            a.serialize(&mut cursor).unwrap_err();
        }

        {
            let serialized = vec![0; buf_size - 1];
            let mut cursor = Cursor::new(&serialized[..]);
            GroupAffine::<P>::deserialize(&mut cursor).unwrap_err();
        }

        {
            let mut serialized = vec![0; a.uncompressed_size()];
            let mut cursor = Cursor::new(&mut serialized[..]);
            a.serialize_uncompressed(&mut cursor).unwrap();

            let mut cursor = Cursor::new(&serialized[..]);
            let b = GroupAffine::<P>::deserialize_uncompressed(&mut cursor).unwrap();
            assert_eq!(a, b);
        }

        {
            let a = GroupAffine::<P>::zero();
            let mut serialized = vec![0; a.uncompressed_size()];
            let mut cursor = Cursor::new(&mut serialized[..]);
            a.serialize_uncompressed(&mut cursor).unwrap();
            let mut cursor = Cursor::new(&serialized[..]);
            let b = GroupAffine::<P>::deserialize_uncompressed(&mut cursor).unwrap();
            assert_eq!(a, b);
        }
    }
}
