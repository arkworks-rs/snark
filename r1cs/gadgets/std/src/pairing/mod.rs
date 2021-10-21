use crate::prelude::*;
use algebra::{Field, PairingEngine};
use r1cs_core::{ConstraintSystem, SynthesisError};
use std::fmt::Debug;

pub mod bls12;
pub mod bn;
pub mod mnt4;
pub mod mnt6;

pub trait PairingGadget<PairingE: PairingEngine, ConstraintF: Field> {
    type G1Gadget: GroupGadget<PairingE::G1Projective, ConstraintF>;
    type G2Gadget: GroupGadget<PairingE::G2Projective, ConstraintF>;

    type G1PreparedGadget: ToBytesGadget<ConstraintF> + Clone + Debug;

    type G2PreparedGadget: ToBytesGadget<ConstraintF> + Clone + Debug;

    type GTGadget: FieldGadget<PairingE::Fqk, ConstraintF> + Clone;

    fn miller_loop<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        p: &[Self::G1PreparedGadget],
        q: &[Self::G2PreparedGadget],
    ) -> Result<Self::GTGadget, SynthesisError>;

    fn final_exponentiation<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        p: &Self::GTGadget,
    ) -> Result<Self::GTGadget, SynthesisError>;

    fn pairing<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        p: Self::G1PreparedGadget,
        q: Self::G2PreparedGadget,
    ) -> Result<Self::GTGadget, SynthesisError> {
        let tmp = Self::miller_loop(cs.ns(|| "miller loop"), &[p], &[q])?;
        Self::final_exponentiation(cs.ns(|| "final_exp"), &tmp)
    }

    /// Computes a product of pairings.
    #[must_use]
    fn product_of_pairings<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        p: &[Self::G1PreparedGadget],
        q: &[Self::G2PreparedGadget],
    ) -> Result<Self::GTGadget, SynthesisError> {
        let miller_result = Self::miller_loop(&mut cs.ns(|| "Miller loop"), p, q)?;
        Self::final_exponentiation(&mut cs.ns(|| "Final Exp"), &miller_result)
    }

    fn prepare_g1<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        q: &Self::G1Gadget,
    ) -> Result<Self::G1PreparedGadget, SynthesisError>;

    fn prepare_g2<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        q: &Self::G2Gadget,
    ) -> Result<Self::G2PreparedGadget, SynthesisError>;
}

#[cfg(test)]
pub(crate) mod tests {
    use crate::{bits::boolean::Boolean, prelude::*, test_constraint_system::TestConstraintSystem};
    use algebra::{BitIterator, Field, Group, PairingEngine, PrimeField, UniformRand};
    use r1cs_core::ConstraintSystem;
    use rand;
    use rand::thread_rng;

    #[allow(dead_code)]
    pub(crate) fn bilinearity_test<
        E: PairingEngine,
        ConstraintF: Field,
        P: PairingGadget<E, ConstraintF>,
    >() {
        let mut cs = TestConstraintSystem::<ConstraintF>::new();

        let mut rng = &mut thread_rng();
        let a = E::G1Projective::rand(&mut rng);
        let b = E::G2Projective::rand(&mut rng);
        let s = E::Fr::rand(&mut rng);

        let sa = a.mul(&s);
        let sb = b.mul(&s);

        let a_g = P::G1Gadget::alloc(&mut cs.ns(|| "a"), || Ok(a)).unwrap();
        let b_g = P::G2Gadget::alloc(&mut cs.ns(|| "b"), || Ok(b)).unwrap();
        let sa_g = P::G1Gadget::alloc(&mut cs.ns(|| "sa"), || Ok(sa)).unwrap();
        let sb_g = P::G2Gadget::alloc(&mut cs.ns(|| "sb"), || Ok(sb)).unwrap();

        let a_prep_g = P::prepare_g1(&mut cs.ns(|| "a_prep"), &a_g).unwrap();
        let b_prep_g = P::prepare_g2(&mut cs.ns(|| "b_prep"), &b_g).unwrap();

        let sa_prep_g = P::prepare_g1(&mut cs.ns(|| "sa_prep"), &sa_g).unwrap();
        let sb_prep_g = P::prepare_g2(&mut cs.ns(|| "sb_prep"), &sb_g).unwrap();

        let (ans1_g, ans1_n) = {
            let ans_g = P::pairing(cs.ns(|| "pair(sa, b)"), sa_prep_g, b_prep_g.clone()).unwrap();
            let ans_n = E::pairing(sa, b).unwrap();
            (ans_g, ans_n)
        };

        let (ans2_g, ans2_n) = {
            let ans_g = P::pairing(cs.ns(|| "pair(a, sb)"), a_prep_g.clone(), sb_prep_g).unwrap();
            let ans_n = E::pairing(a, sb).unwrap();
            (ans_g, ans_n)
        };

        let (ans3_g, ans3_n) = {
            let s_iter = BitIterator::new(s.into_repr())
                .map(Boolean::constant)
                .collect::<Vec<_>>();

            let mut ans_g = P::pairing(cs.ns(|| "pair(a, b)"), a_prep_g, b_prep_g).unwrap();
            let mut ans_n = E::pairing(a, b).unwrap();
            ans_n = ans_n.pow(s.into_repr());
            ans_g = ans_g.pow(cs.ns(|| "pow"), &s_iter).unwrap();

            (ans_g, ans_n)
        };

        assert_eq!(ans1_n, ans2_n, "Failed ans1_native == ans2_native");
        assert_eq!(ans2_n, ans3_n, "Failed ans2_native == ans3_native");
        assert_eq!(
            ans1_g.get_value(),
            ans3_g.get_value(),
            "Failed ans1 == ans3"
        );
        assert_eq!(
            ans1_g.get_value(),
            ans2_g.get_value(),
            "Failed ans1 == ans2"
        );
        assert_eq!(
            ans2_g.get_value(),
            ans3_g.get_value(),
            "Failed ans2 == ans3"
        );

        ans1_g
            .enforce_equal(&mut cs.ns(|| "ans1 == ans2?"), &ans2_g)
            .unwrap();
        ans2_g
            .enforce_equal(&mut cs.ns(|| "ans2 == ans3?"), &ans3_g)
            .unwrap();

        assert_eq!(ans1_g.get_value().unwrap(), ans1_n, "Failed native test 1");
        assert_eq!(ans2_g.get_value().unwrap(), ans2_n, "Failed native test 2");
        assert_eq!(ans3_g.get_value().unwrap(), ans3_n, "Failed native test 3");

        if !cs.is_satisfied() {
            println!("Unsatisfied: {:?}", cs.which_is_unsatisfied());
        }
        assert!(cs.is_satisfied(), "cs is not satisfied");
    }
}
