use crate::prelude::*;
use algebra::{Field, PairingEngine};
use core::fmt::Debug;
use r1cs_core::SynthesisError;

/// This module implements pairings for BLS12 bilinear groups.
pub mod bls12;
/// This module implements pairings for MNT4 bilinear groups.
pub mod mnt4;
/// This module implements pairings for MNT6 bilinear groups.
pub mod mnt6;

/// Specifies the constraints for computing a pairing in the yybilinear group `E`.
pub trait PairingVar<E: PairingEngine, ConstraintF: Field = <E as PairingEngine>::Fq> {
    /// An variable representing an element of `G1`.
    /// This is the R1CS equivalent of `E::G1Projective`.
    type G1Var: CurveVar<E::G1Projective, ConstraintF>
        + AllocVar<E::G1Projective, ConstraintF>
        + AllocVar<E::G1Affine, ConstraintF>;

    /// An variable representing an element of `G2`.
    /// This is the R1CS equivalent of `E::G2Projective`.
    type G2Var: CurveVar<E::G2Projective, ConstraintF>
        + AllocVar<E::G2Projective, ConstraintF>
        + AllocVar<E::G2Affine, ConstraintF>;

    /// An variable representing an element of `GT`.
    /// This is the R1CS equivalent of `E::GT`.
    type GTVar: FieldVar<E::Fqk, ConstraintF>;

    /// An variable representing cached precomputation  that can speed up pairings computations.
    /// This is the R1CS equivalent of `E::G1Prepared`.
    type G1PreparedVar: ToBytesGadget<ConstraintF>
        + AllocVar<E::G1Prepared, ConstraintF>
        + Clone
        + Debug;
    /// An variable representing cached precomputation  that can speed up pairings computations.
    /// This is the R1CS equivalent of `E::G2Prepared`.
    type G2PreparedVar: ToBytesGadget<ConstraintF>
        + AllocVar<E::G2Prepared, ConstraintF>
        + Clone
        + Debug;

    /// Computes a multi-miller loop between elements
    /// of `p` and `q`.
    fn miller_loop(
        p: &[Self::G1PreparedVar],
        q: &[Self::G2PreparedVar],
    ) -> Result<Self::GTVar, SynthesisError>;

    /// Computes a final exponentiation over `p`.
    fn final_exponentiation(p: &Self::GTVar) -> Result<Self::GTVar, SynthesisError>;

    /// Computes a pairing over `p` and `q`.
    #[tracing::instrument(target = "r1cs")]
    fn pairing(
        p: Self::G1PreparedVar,
        q: Self::G2PreparedVar,
    ) -> Result<Self::GTVar, SynthesisError> {
        let tmp = Self::miller_loop(&[p], &[q])?;
        Self::final_exponentiation(&tmp)
    }

    /// Computes a product of pairings over the elements in `p` and `q`.
    #[must_use]
    #[tracing::instrument(target = "r1cs")]
    fn product_of_pairings(
        p: &[Self::G1PreparedVar],
        q: &[Self::G2PreparedVar],
    ) -> Result<Self::GTVar, SynthesisError> {
        let miller_result = Self::miller_loop(p, q)?;
        Self::final_exponentiation(&miller_result)
    }

    /// Performs the precomputation to generate `Self::G1PreparedVar`.
    fn prepare_g1(q: &Self::G1Var) -> Result<Self::G1PreparedVar, SynthesisError>;

    /// Performs the precomputation to generate `Self::G2PreparedVar`.
    fn prepare_g2(q: &Self::G2Var) -> Result<Self::G2PreparedVar, SynthesisError>;
}

#[cfg(test)]
pub(crate) mod tests {
    use crate::{prelude::*, Vec};
    use algebra::{
        test_rng, BitIteratorLE, Field, PairingEngine, PrimeField, ProjectiveCurve, UniformRand,
    };
    use r1cs_core::{ConstraintSystem, SynthesisError};

    #[allow(dead_code)]
    pub(crate) fn bilinearity_test<E: PairingEngine, P: PairingVar<E>>(
    ) -> Result<(), SynthesisError>
    where
        for<'a> &'a P::G1Var: GroupOpsBounds<'a, E::G1Projective, P::G1Var>,
        for<'a> &'a P::G2Var: GroupOpsBounds<'a, E::G2Projective, P::G2Var>,
        for<'a> &'a P::GTVar: FieldOpsBounds<'a, E::Fqk, P::GTVar>,
    {
        let cs = ConstraintSystem::<E::Fq>::new_ref();

        let mut rng = test_rng();
        let a = E::G1Projective::rand(&mut rng);
        let b = E::G2Projective::rand(&mut rng);
        let s = E::Fr::rand(&mut rng);

        let mut sa = a;
        sa *= s;
        let mut sb = b;
        sb *= s;

        let a_g = P::G1Var::new_witness(cs.clone(), || Ok(a.into_affine()))?;
        let b_g = P::G2Var::new_witness(cs.clone(), || Ok(b.into_affine()))?;
        let sa_g = P::G1Var::new_witness(cs.clone(), || Ok(sa.into_affine()))?;
        let sb_g = P::G2Var::new_witness(cs.clone(), || Ok(sb.into_affine()))?;

        let mut preparation_num_constraints = cs.num_constraints();
        let a_prep_g = P::prepare_g1(&a_g)?;
        let b_prep_g = P::prepare_g2(&b_g)?;
        preparation_num_constraints = cs.num_constraints() - preparation_num_constraints;
        println!(
            "Preparation num constraints: {}",
            preparation_num_constraints
        );

        let sa_prep_g = P::prepare_g1(&sa_g)?;
        let sb_prep_g = P::prepare_g2(&sb_g)?;

        let (ans1_g, ans1_n) = {
            let ml_constraints = cs.num_constraints();
            let ml_g = P::miller_loop(&[sa_prep_g], &[b_prep_g.clone()])?;
            println!(
                "ML num constraints: {}",
                cs.num_constraints() - ml_constraints
            );
            let fe_constraints = cs.num_constraints();
            let ans_g = P::final_exponentiation(&ml_g)?;
            println!(
                "FE num constraints: {}",
                cs.num_constraints() - fe_constraints
            );
            let ans_n = E::pairing(sa, b);
            (ans_g, ans_n)
        };

        let (ans2_g, ans2_n) = {
            let ans_g = P::pairing(a_prep_g.clone(), sb_prep_g)?;
            let ans_n = E::pairing(a, sb);
            (ans_g, ans_n)
        };

        let (ans3_g, ans3_n) = {
            let s_iter = BitIteratorLE::without_trailing_zeros(s.into_repr())
                .map(Boolean::constant)
                .collect::<Vec<_>>();

            let mut ans_g = P::pairing(a_prep_g, b_prep_g)?;
            let mut ans_n = E::pairing(a, b);
            ans_n = ans_n.pow(s.into_repr());
            ans_g = ans_g.pow_le(&s_iter)?;

            (ans_g, ans_n)
        };

        ans1_g.enforce_equal(&ans2_g)?;
        ans2_g.enforce_equal(&ans3_g)?;

        assert_eq!(ans1_g.value()?, ans1_n, "Failed native test 1");
        assert_eq!(ans2_g.value()?, ans2_n, "Failed native test 2");
        assert_eq!(ans3_g.value()?, ans3_n, "Failed native test 3");

        assert_eq!(ans1_n, ans2_n, "Failed ans1_native == ans2_native");
        assert_eq!(ans2_n, ans3_n, "Failed ans2_native == ans3_native");
        assert_eq!(ans1_g.value()?, ans3_g.value()?, "Failed ans1 == ans3");
        assert_eq!(ans1_g.value()?, ans2_g.value()?, "Failed ans1 == ans2");
        assert_eq!(ans2_g.value()?, ans3_g.value()?, "Failed ans2 == ans3");

        if !cs.is_satisfied().unwrap() {
            println!("Unsatisfied: {:?}", cs.which_is_unsatisfied());
        }

        assert!(cs.is_satisfied().unwrap(), "cs is not satisfied");
        Ok(())
    }
}
