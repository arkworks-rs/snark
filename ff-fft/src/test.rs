use crate::domain::*;
use algebra::{
    bls12_381::{Fr, G1Projective},
    mnt6_753::{Fr as MNT6Fr, G1Projective as MNT6G1Projective},
};
use algebra_core::{test_rng, PrimeField, UniformRand};

// Test multiplying various (low degree) polynomials together and
// comparing with naive evaluations.
#[test]
fn fft_composition() {
    fn test_fft_composition<
        F: PrimeField,
        T: DomainCoeff<F> + UniformRand + core::fmt::Debug + Eq,
        R: rand::Rng,
        D: EvaluationDomain<F>,
    >(
        rng: &mut R,
        max_coeffs: usize,
    ) {
        for coeffs in 0..max_coeffs {
            let coeffs = 1 << coeffs;

            let domain = D::new(coeffs).unwrap();

            let mut v = vec![];
            for _ in 0..coeffs {
                v.push(T::rand(rng));
            }
            // Fill up with zeros.
            v.resize(domain.size(), T::zero());
            let mut v2 = v.clone();

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

    test_fft_composition::<Fr, Fr, _, GeneralEvaluationDomain<Fr>>(rng, 10);
    test_fft_composition::<Fr, G1Projective, _, GeneralEvaluationDomain<Fr>>(rng, 10);
    // This will result in a mixed-radix domain being used.
    test_fft_composition::<MNT6Fr, MNT6Fr, _, MixedRadixEvaluationDomain<MNT6Fr>>(rng, 17);
    test_fft_composition::<MNT6Fr, MNT6G1Projective, _, MixedRadixEvaluationDomain<MNT6Fr>>(rng, 5);
}
