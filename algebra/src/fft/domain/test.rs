use crate::curves::{bls12_381::Bls12_381, PairingEngine};
use crate::{domain::*, multicore::*};
use crate::{Field, FpParameters, PrimeField, UniformRand};
use rand;
use rand::Rng;
use std::cmp::min;

// Test multiplying various (low degree) polynomials together and
// comparing with naive evaluations.
#[test]
fn fft_composition() {
    fn test_fft_composition<E: PairingEngine, R: Rng>(rng: &mut R) {
        for coeffs in 0..18 {
            let coeffs = 1 << coeffs;

            let mut v = vec![];
            for _ in 0..coeffs {
                v.push(E::Fr::rand(rng));
            }
            let mut v2 = v.clone();

            let domain = get_best_evaluation_domain::<E::Fr>(coeffs).unwrap();
            v.resize(domain.size(), E::Fr::zero());

            domain.ifft_in_place(&mut v2);
            domain.fft_in_place(&mut v2);
            assert_eq!(v, v2, "ifft(fft(.)) != iden");

            domain.fft_in_place(&mut v2);
            domain.ifft_in_place(&mut v2);
            assert_eq!(v, v2, "fft(ifft(.)) != iden");

            domain.coset_ifft_in_place(&mut v2);
            domain.coset_fft_in_place(&mut v2);
            assert_eq!(v, v2, "coset_ifft(coset_fft(.)) != iden");

            domain.coset_fft_in_place(&mut v2);
            domain.coset_ifft_in_place(&mut v2);
            assert_eq!(v, v2, "coset_ifft(coset_fft(.)) != iden");
        }
    }

    let rng = &mut rand::thread_rng();

    test_fft_composition::<Bls12_381, _>(rng);
}

#[test]
fn fft_consistency() {
    fn test_consistency<E: PairingEngine, R: Rng>(rng: &mut R) {
        let worker = Worker::new();

        for _ in 0..5 {
            for log_d in 0..18 {
                let d = 1 << log_d;

                let mut v1 = (0..d).map(|_| E::Fr::rand(rng)).collect::<Vec<_>>();
                let mut v2 = v1.clone();

                let domain = get_best_evaluation_domain::<E::Fr>(v1.len()).unwrap();

                for log_cpus in log_d..min(log_d + 1, 3) {
                    if log_d < <E::Fr as PrimeField>::Params::TWO_ADICITY {
                        BasicRadix2Domain::parallel_fft(
                            &mut v1,
                            &worker,
                            domain.group_gen(),
                            log_d,
                            log_cpus,
                        );
                        BasicRadix2Domain::serial_fft(&mut v2, domain.group_gen(), log_d);
                    } else {
                        MixedRadix2Domain::mixed_parallel_fft(
                            &mut v1,
                            &worker,
                            domain.group_gen(),
                            log_d,
                            log_cpus,
                        );
                        MixedRadix2Domain::mixed_serial_fft(&mut v2, domain.group_gen(), log_d);
                    }
                    assert_eq!(v1, v2);
                }
            }
        }
    }

    let rng = &mut rand::thread_rng();

    test_consistency::<Bls12_381, _>(rng);
}
