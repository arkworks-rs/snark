mod pairing {
    use rand::{Rand, SeedableRng, XorShiftRng};

    use algebra::{
        curves::{
            bls12::{G1Prepared, G2Prepared},
            bls12_377::{
                Bls12_377, Bls12_377Parameters, G1Affine, G1Projective as G1, G2Affine,
                G2Projective as G2,
            },
            PairingCurve, PairingEngine,
        },
        fields::bls12_377::Fq12,
    };

    #[bench]
    fn bench_pairing_miller_loop(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        let v: Vec<(
            G1Prepared<Bls12_377Parameters>,
            G2Prepared<Bls12_377Parameters>,
        )> = (0..SAMPLES)
            .map(|_| {
                (
                    G1Affine::from(G1::rand(&mut rng)).prepare(),
                    G2Affine::from(G2::rand(&mut rng)).prepare(),
                )
            })
            .collect();

        let mut count = 0;
        b.iter(|| {
            let tmp = Bls12_377::miller_loop(&[(&v[count].0, &v[count].1)]);
            count = (count + 1) % SAMPLES;
            tmp
        });
    }

    #[bench]
    fn bench_pairing_final_exponentiation(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        let v: Vec<Fq12> = (0..SAMPLES)
            .map(|_| {
                (
                    G1Affine::from(G1::rand(&mut rng)).prepare(),
                    G2Affine::from(G2::rand(&mut rng)).prepare(),
                )
            })
            .map(|(ref p, ref q)| Bls12_377::miller_loop(&[(p, q)]))
            .collect();

        let mut count = 0;
        b.iter(|| {
            let tmp = Bls12_377::final_exponentiation(&v[count]);
            count = (count + 1) % SAMPLES;
            tmp
        });
    }

    #[bench]
    fn bench_pairing_full(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        let v: Vec<(G1, G2)> = (0..SAMPLES)
            .map(|_| (G1::rand(&mut rng), G2::rand(&mut rng)))
            .collect();

        let mut count = 0;
        b.iter(|| {
            let tmp = Bls12_377::pairing(v[count].0, v[count].1);
            count = (count + 1) % SAMPLES;
            tmp
        });
    }
}
