use crate::domain::*;
use algebra::{
    bls12_381::{Bls12_381, Fr, G1Projective},
    mnt6_753::{Fr as MNT6Fr, G1Projective as MNT6G1Projective, MNT6_753},
};
use algebra_core::{test_rng, PairingEngine, PrimeField, UniformRand};

// Test multiplying various (low degree) polynomials together and
// comparing with naive evaluations.
#[test]
fn fft_composition() {
    fn test_fft_composition<
        F: PrimeField,
        T: DomainCoeff<F> + UniformRand + core::fmt::Debug + Eq,
        R: rand::Rng,
    >(
        rng: &mut R,
        max_coeffs: usize,
    ) {
        for coeffs in 0..max_coeffs {
            let coeffs = 1 << coeffs;

            let mut v = vec![];
            for _ in 0..coeffs {
                v.push(T::rand(rng));
            }
            let mut v2 = v.clone();

            let domain = EvaluationDomain::<F>::new(coeffs).unwrap();
            domain.ifft_in_place(&mut v2);
            domain.fft_in_place(&mut v2);
            assert_eq!(v, v2, "ifft(fft(.)) != iden");

            domain.fft_in_place(&mut v2);
            domain.ifft_in_place(&mut v2);
            assert_eq!(v, v2, "fft(ifft(.)) != iden");

            domain.coset_ifft_in_place(&mut v2);
            domain.coset_fft_in_place(&mut v2);
            assert_eq!(v, v2, "coset_fft(coset_ifft(.)) != iden");

            domain.coset_fft_in_place(&mut v2);
            domain.coset_ifft_in_place(&mut v2);
            assert_eq!(v, v2, "coset_ifft(coset_fft(.)) != iden");
        }
    }

    let rng = &mut test_rng();

    test_fft_composition::<Fr, Fr, _>(rng, 10);
    test_fft_composition::<Fr, G1Projective, _>(rng, 10);
    test_fft_composition::<MNT6Fr, MNT6Fr, _>(rng, 16);
    test_fft_composition::<MNT6Fr, MNT6G1Projective, _>(rng, 5);
}

#[test]
#[cfg(feature = "parallel")]
fn parallel_fft_consistency() {
    use crate::Vec;
    use core::cmp::min;

    fn test_consistency<E: PairingEngine, R: rand::Rng>(rng: &mut R, max_coeffs: u32) {
        for _ in 0..5 {
            for log_d in 0..max_coeffs {
                let d = 1 << log_d;

                let mut v1 = (0..d).map(|_| E::Fr::rand(rng)).collect::<Vec<_>>();
                let mut v2 = v1.clone();

                let domain = EvaluationDomain::new(v1.len()).unwrap();

                for log_cpus in log_d..min(log_d + 1, 3) {
                    parallel_fft::<E::Fr, E::Fr>(&mut v1, domain.group_gen, log_d, log_cpus);
                    serial_fft::<E::Fr, E::Fr>(&mut v2, domain.group_gen, log_d);

                    assert_eq!(v1, v2);
                }
            }
        }
    }

    let rng = &mut test_rng();

    test_consistency::<Bls12_381, _>(rng, 10);
    test_consistency::<MNT6_753, _>(rng, 16);
}
