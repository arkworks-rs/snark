macro_rules! pairing_bench {
    ($curve:ident, $pairing_field:ident, $pairing_type:ident) => {
        #[bench]
        fn bench_pairing_miller_loop(b: &mut ::test::Bencher) {
            const SAMPLES: usize = 1000;

            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

            $pairing_type!(v, rng);

            let mut count = 0;
            b.iter(|| {
                let tmp = $curve::miller_loop(&[(v[count].0.clone(), v[count].1.clone())]);
                count = (count + 1) % SAMPLES;
                tmp
            });
        }

        #[bench]
        fn bench_pairing_final_exponentiation(b: &mut ::test::Bencher) {
            const SAMPLES: usize = 1000;

            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

            let v: Vec<$pairing_field> = (0..SAMPLES)
                .map(|_| {
                    (
                        G1Affine::from(G1::rand(&mut rng)).into(),
                        G2Affine::from(G2::rand(&mut rng)).into(),
                    )
                })
                .map(|(p, q)| $curve::miller_loop(&[(p, q)]))
                .collect();

            let mut count = 0;
            b.iter(|| {
                let tmp = $curve::final_exponentiation(&v[count]);
                count = (count + 1) % SAMPLES;
                tmp
            });
        }

        #[bench]
        fn bench_pairing_full(b: &mut ::test::Bencher) {
            const SAMPLES: usize = 1000;

            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

            let v: Vec<(G1, G2)> = (0..SAMPLES)
                .map(|_| (G1::rand(&mut rng), G2::rand(&mut rng)))
                .collect();

            let mut count = 0;
            b.iter(|| {
                let tmp = $curve::pairing(v[count].0, v[count].1);
                count = (count + 1) % SAMPLES;
                tmp
            });
        }
    };
}
