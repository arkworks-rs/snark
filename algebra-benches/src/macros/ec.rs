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
                n_fold!(tmp, v, add_assign, count);
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
                n_fold!(tmp, v, add_assign_mixed, count);
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
                n_fold!(tmp, double_in_place);
                count = (count + 1) % SAMPLES;
                tmp
            });
        }
    }
}
