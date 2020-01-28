use crate::{
    curves::models::twisted_edwards_extended::{GroupAffine, GroupProjective},
    fields::Field,
    CanonicalDeserialize, CanonicalSerialize, MontgomeryModelParameters, ProjectiveCurve,
    SerializationError, TEModelParameters, UniformRand,
};
use num_traits::{One, Zero};
use rand::SeedableRng;
use rand_xorshift::XorShiftRng;

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

pub const ITERATIONS: usize = 10;

pub fn edwards_curve_serialization_test<P: TEModelParameters>(buf_size: usize) {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..ITERATIONS {
        let a = GroupProjective::<P>::rand(&mut rng);
        let a = a.into_affine();
        {
            let mut serialized = vec![0; buf_size];
            a.serialize(&[], &mut serialized).unwrap();

            let mut extra_info_buf = [false; 0];
            let b = GroupAffine::<P>::deserialize(&serialized, &mut extra_info_buf).unwrap();
            assert_eq!(a, b);
        }

        {
            let a = GroupAffine::<P>::zero();
            let mut serialized = vec![0; buf_size];
            a.serialize(&[], &mut serialized).unwrap();
            let mut extra_info_buf = [false; 0];
            let b = GroupAffine::<P>::deserialize(&serialized, &mut extra_info_buf).unwrap();
            assert_eq!(a, b);
        }

        {
            let a = GroupAffine::<P>::zero();
            let mut serialized = vec![0; buf_size];
            assert!(if let SerializationError::ExtraInfoWrongSize =
                a.serialize(&[false; 1], &mut serialized).unwrap_err()
            {
                true
            } else {
                false
            });
        }

        {
            let a = GroupAffine::<P>::zero();
            let mut serialized = vec![0; buf_size];
            a.serialize(&[], &mut serialized).unwrap();
            let mut extra_info_buf = [false; 1];
            assert!(if let SerializationError::ExtraInfoWrongSize =
                GroupAffine::<P>::deserialize(&serialized, &mut extra_info_buf).unwrap_err()
            {
                true
            } else {
                false
            });
        }

        {
            let a = GroupAffine::<P>::zero();
            let mut serialized = vec![0; buf_size + 1];
            assert!(if let SerializationError::BufferWrongSize =
                a.serialize(&[false; 0], &mut serialized).unwrap_err()
            {
                true
            } else {
                false
            });
        }

        {
            let serialized = vec![0; buf_size + 1];
            let mut extra_info_buf = [false; 0];
            assert!(if let SerializationError::BufferWrongSize =
                GroupAffine::<P>::deserialize(&serialized, &mut extra_info_buf).unwrap_err()
            {
                true
            } else {
                false
            });
        }
    }
}
