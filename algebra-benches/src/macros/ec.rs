macro_rules! ec_bench {
    () => {
        #[bench]
        fn bench_g1_rand(b: &mut ::test::Bencher) {
            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
            b.iter(|| G1::rand(&mut rng));
        }

        #[bench]
        fn bench_g1_mul_assign(b: &mut ::test::Bencher) {
            const SAMPLES: usize = 1000;

            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

            let v: Vec<(G1, Fr)> = (0..SAMPLES)
                .map(|_| (G1::rand(&mut rng), Fr::rand(&mut rng)))
                .collect();

            let mut count = 0;
            b.iter(|| {
                let mut tmp = v[count].0;
                tmp *= v[count].1;
                count = (count + 1) % SAMPLES;
                tmp
            });
        }

        #[bench]
        fn bench_g1_add_assign(b: &mut ::test::Bencher) {
            const SAMPLES: usize = 1000;

            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

            let v: Vec<(G1, G1)> = (0..SAMPLES)
                .map(|_| (G1::rand(&mut rng), G1::rand(&mut rng)))
                .collect();

            let mut count = 0;
            b.iter(|| {
                let mut tmp = v[count].0;
                n_fold!(tmp, v, add_assign, count);
                count = (count + 1) % SAMPLES;
                tmp
            });
        }

        #[bench]
        fn bench_g1_add_assign_mixed(b: &mut ::test::Bencher) {
            const SAMPLES: usize = 1000;

            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

            let v: Vec<(G1, G1Affine)> = (0..SAMPLES)
                .map(|_| (G1::rand(&mut rng), G1::rand(&mut rng).into()))
                .collect();

            let mut count = 0;
            b.iter(|| {
                let mut tmp = v[count].0;
                n_fold!(tmp, v, add_assign_mixed, count);
                count = (count + 1) % SAMPLES;
                tmp
            });
        }

        #[bench]
        fn bench_g1_double(b: &mut ::test::Bencher) {
            const SAMPLES: usize = 1000;

            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

            let v: Vec<(G1, G1)> = (0..SAMPLES)
                .map(|_| (G1::rand(&mut rng), G1::rand(&mut rng)))
                .collect();

            let mut count = 0;
            b.iter(|| {
                let mut tmp = v[count].0;
                n_fold!(tmp, double_in_place);
                count = (count + 1) % SAMPLES;
                tmp
            });
        }

        #[bench]
        fn bench_g1_deser(b: &mut ::test::Bencher) {
            use algebra::{CanonicalDeserialize, CanonicalSerialize, ProjectiveCurve};
            const SAMPLES: usize = 1000;

            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

            let mut num_bytes = 0;
            let tmp = G1::rand(&mut rng).into_affine();
            let v: Vec<_> = (0..SAMPLES)
                .flat_map(|_| {
                    let mut bytes = Vec::with_capacity(1000);
                    tmp.serialize(&mut bytes).unwrap();
                    num_bytes = bytes.len();
                    bytes
                })
                .collect();

            let mut count = 0;
            b.iter(|| {
                count = (count + 1) % SAMPLES;
                let index = count * num_bytes;
                G1Affine::deserialize(&v[index..(index + num_bytes)]).unwrap()
            });
        }

        #[bench]
        fn bench_g1_ser(b: &mut ::test::Bencher) {
            use algebra::{CanonicalSerialize, ProjectiveCurve};
            const SAMPLES: usize = 1000;

            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

            let mut v: Vec<_> = (0..SAMPLES).map(|_| G1::rand(&mut rng)).collect();
            let v = G1::batch_normalization_into_affine(v.as_mut_slice());
            let mut bytes = Vec::with_capacity(1000);

            let mut count = 0;
            b.iter(|| {
                let tmp = v[count];
                count = (count + 1) % SAMPLES;
                bytes.clear();
                tmp.serialize(&mut bytes)
            });
        }

        #[bench]
        fn bench_g1_deser_unchecked(b: &mut ::test::Bencher) {
            use algebra::{CanonicalDeserialize, CanonicalSerialize, ProjectiveCurve};
            const SAMPLES: usize = 1000;

            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

            let mut num_bytes = 0;
            let tmp = G1::rand(&mut rng).into_affine();
            let v: Vec<_> = (0..SAMPLES)
                .flat_map(|_| {
                    let mut bytes = Vec::with_capacity(1000);
                    tmp.serialize_unchecked(&mut bytes).unwrap();
                    num_bytes = bytes.len();
                    bytes
                })
                .collect();

            let mut count = 0;
            b.iter(|| {
                count = (count + 1) % SAMPLES;
                let index = count * num_bytes;
                G1Affine::deserialize_unchecked(&v[index..(index + num_bytes)]).unwrap()
            });
        }

        #[bench]
        fn bench_g1_ser_unchecked(b: &mut ::test::Bencher) {
            use algebra::{CanonicalSerialize, ProjectiveCurve};
            const SAMPLES: usize = 1000;

            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

            let mut v: Vec<_> = (0..SAMPLES).map(|_| G1::rand(&mut rng)).collect();
            let v = G1::batch_normalization_into_affine(v.as_mut_slice());
            let mut bytes = Vec::with_capacity(1000);

            let mut count = 0;
            b.iter(|| {
                let tmp = v[count];
                count = (count + 1) % SAMPLES;
                bytes.clear();
                tmp.serialize_unchecked(&mut bytes)
            });
        }

        #[bench]
        fn bench_g2_rand(b: &mut ::test::Bencher) {
            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
            b.iter(|| G2::rand(&mut rng));
        }

        #[bench]
        fn bench_g2_mul_assign(b: &mut ::test::Bencher) {
            const SAMPLES: usize = 1000;

            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

            let v: Vec<(G2, Fr)> = (0..SAMPLES)
                .map(|_| (G2::rand(&mut rng), Fr::rand(&mut rng)))
                .collect();

            let mut count = 0;
            b.iter(|| {
                let mut tmp = v[count].0;
                tmp *= v[count].1;
                count = (count + 1) % SAMPLES;
                tmp
            });
        }

        #[bench]
        fn bench_g2_add_assign(b: &mut ::test::Bencher) {
            const SAMPLES: usize = 1000;

            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

            let v: Vec<(G2, G2)> = (0..SAMPLES)
                .map(|_| (G2::rand(&mut rng), G2::rand(&mut rng)))
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
        fn bench_g2_add_assign_mixed(b: &mut ::test::Bencher) {
            const SAMPLES: usize = 1000;

            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

            let v: Vec<(G2, G2Affine)> = (0..SAMPLES)
                .map(|_| (G2::rand(&mut rng), G2::rand(&mut rng).into()))
                .collect();

            let mut count = 0;
            b.iter(|| {
                let mut tmp = v[count].0;
                tmp.add_assign_mixed(&v[count].1);
                count = (count + 1) % SAMPLES;
                tmp
            });
        }

        #[bench]
        fn bench_g2_double(b: &mut ::test::Bencher) {
            const SAMPLES: usize = 1000;

            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

            let v: Vec<(G2, G2)> = (0..SAMPLES)
                .map(|_| (G2::rand(&mut rng), G2::rand(&mut rng)))
                .collect();

            let mut count = 0;
            b.iter(|| {
                let mut tmp = v[count].0;
                tmp.double_in_place();
                count = (count + 1) % SAMPLES;
                tmp
            });
        }

        #[bench]
        fn bench_g2_deser(b: &mut ::test::Bencher) {
            use algebra::{CanonicalDeserialize, CanonicalSerialize, ProjectiveCurve};
            const SAMPLES: usize = 1000;

            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

            let mut num_bytes = 0;
            let tmp = G2::rand(&mut rng).into_affine();
            let v: Vec<_> = (0..SAMPLES)
                .flat_map(|_| {
                    let mut bytes = Vec::with_capacity(1000);
                    tmp.serialize(&mut bytes).unwrap();
                    num_bytes = bytes.len();
                    bytes
                })
                .collect();

            let mut count = 0;
            b.iter(|| {
                count = (count + 1) % SAMPLES;
                let index = count * num_bytes;
                G2Affine::deserialize(&v[index..(index + num_bytes)]).unwrap()
            });
        }

        #[bench]
        fn bench_g2_ser(b: &mut ::test::Bencher) {
            use algebra::{CanonicalSerialize, ProjectiveCurve};
            const SAMPLES: usize = 1000;

            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

            let mut v: Vec<_> = (0..SAMPLES).map(|_| G2::rand(&mut rng)).collect();
            let v = G2::batch_normalization_into_affine(v.as_mut_slice());
            let mut bytes = Vec::with_capacity(1000);

            let mut count = 0;
            b.iter(|| {
                let tmp = v[count];
                count = (count + 1) % SAMPLES;
                bytes.clear();
                tmp.serialize(&mut bytes)
            });
        }

        #[bench]
        fn bench_g2_deser_unchecked(b: &mut ::test::Bencher) {
            use algebra::{CanonicalDeserialize, CanonicalSerialize, ProjectiveCurve};
            const SAMPLES: usize = 1000;

            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

            let mut num_bytes = 0;
            let tmp = G2::rand(&mut rng).into_affine();
            let v: Vec<_> = (0..SAMPLES)
                .flat_map(|_| {
                    let mut bytes = Vec::with_capacity(1000);
                    tmp.serialize_unchecked(&mut bytes).unwrap();
                    num_bytes = bytes.len();
                    bytes
                })
                .collect();

            let mut count = 0;
            b.iter(|| {
                count = (count + 1) % SAMPLES;
                let index = count * num_bytes;
                G2Affine::deserialize_unchecked(&v[index..(index + num_bytes)]).unwrap()
            });
        }

        #[bench]
        fn bench_g2_ser_unchecked(b: &mut ::test::Bencher) {
            use algebra::{CanonicalSerialize, ProjectiveCurve};
            const SAMPLES: usize = 1000;

            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

            let mut v: Vec<_> = (0..SAMPLES).map(|_| G2::rand(&mut rng)).collect();
            let v = G2::batch_normalization_into_affine(v.as_mut_slice());
            let mut bytes = Vec::with_capacity(1000);

            let mut count = 0;
            b.iter(|| {
                let tmp = v[count];
                count = (count + 1) % SAMPLES;
                bytes.clear();
                tmp.serialize_unchecked(&mut bytes)
            });
        }
    };
}
