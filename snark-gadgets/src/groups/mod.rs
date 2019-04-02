use crate::{
    bits::boolean::Boolean,
    utils::{AllocGadget, CondSelectGadget, EqGadget, NEqGadget, ToBitsGadget, ToBytesGadget},
};
use algebra::{Group, PairingEngine};
use snark::{ConstraintSystem, SynthesisError};

use std::{borrow::Borrow, fmt::Debug};

pub mod curves;

pub use self::curves::{
    short_weierstrass::bls12,
    twisted_edwards::{edwards_sw6, jubjub},
};

pub trait GroupGadget<G: Group, E: PairingEngine>:
    Sized
    + ToBytesGadget<E>
    + NEqGadget<E>
    + EqGadget<E>
    + ToBitsGadget<E>
    + CondSelectGadget<E>
    + AllocGadget<G, E>
    + Clone
    + Debug
{
    type Value: Debug;
    type Variable;

    fn get_value(&self) -> Option<Self::Value>;

    fn get_variable(&self) -> Self::Variable;

    fn zero<CS: ConstraintSystem<E>>(cs: CS) -> Result<Self, SynthesisError>;

    fn add<CS: ConstraintSystem<E>>(&self, cs: CS, other: &Self) -> Result<Self, SynthesisError>;

    fn sub<CS: ConstraintSystem<E>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        let neg_other = other.negate(cs.ns(|| "Negate other"))?;
        self.add(cs.ns(|| "Self - other"), &neg_other)
    }

    fn add_constant<CS: ConstraintSystem<E>>(
        &self,
        cs: CS,
        other: &G,
    ) -> Result<Self, SynthesisError>;

    fn sub_constant<CS: ConstraintSystem<E>>(
        &self,
        mut cs: CS,
        other: &G,
    ) -> Result<Self, SynthesisError> {
        let neg_other = -(*other);
        self.add_constant(cs.ns(|| "Self - other"), &neg_other)
    }

    fn double_in_place<CS: ConstraintSystem<E>>(&mut self, cs: CS) -> Result<(), SynthesisError>;

    fn negate<CS: ConstraintSystem<E>>(&self, cs: CS) -> Result<Self, SynthesisError>;

    /// Inputs must be specified in *little-endian* form.
    /// If the addition law is incomplete for the identity element,
    /// `result` must not be the identity element.
    fn mul_bits<'a, CS: ConstraintSystem<E>>(
        &self,
        mut cs: CS,
        result: &Self,
        bits: impl Iterator<Item = &'a Boolean>,
    ) -> Result<Self, SynthesisError> {
        let mut power = self.clone();
        let mut result = result.clone();
        for (i, bit) in bits.enumerate() {
            let new_encoded = result.add(&mut cs.ns(|| format!("Add {}-th power", i)), &power)?;
            result = Self::conditionally_select(
                &mut cs.ns(|| format!("Select {}", i)),
                bit.borrow(),
                &new_encoded,
                &result,
            )?;
            power.double_in_place(&mut cs.ns(|| format!("{}-th Doubling", i)))?;
        }
        Ok(result)
    }

    fn precomputed_base_scalar_mul<'a, CS, I, B>(
        &mut self,
        mut cs: CS,
        scalar_bits_with_base_powers: I,
    ) -> Result<(), SynthesisError>
    where
        CS: ConstraintSystem<E>,
        I: Iterator<Item = (B, &'a G)>,
        B: Borrow<Boolean>,
        G: 'a,
    {
        for (i, (bit, base_power)) in scalar_bits_with_base_powers.enumerate() {
            let new_encoded = self.add_constant(
                &mut cs.ns(|| format!("Add {}-th base power", i)),
                &base_power,
            )?;
            *self = Self::conditionally_select(
                &mut cs.ns(|| format!("Conditional Select {}", i)),
                bit.borrow(),
                &new_encoded,
                &self,
            )?;
        }
        Ok(())
    }

    fn precomputed_base_multiscalar_mul<'a, CS, T, I, B>(
        mut cs: CS,
        bases: &[B],
        scalars: I,
    ) -> Result<Self, SynthesisError>
    where
        CS: ConstraintSystem<E>,
        T: 'a + ToBitsGadget<E> + ?Sized,
        I: Iterator<Item = &'a T>,
        B: Borrow<[G]>,
    {
        let mut result = Self::zero(&mut cs.ns(|| "Declare Result"))?;
        // Compute ∏(h_i^{m_i}) for all i.
        for (i, (bits, base_powers)) in scalars.zip(bases).enumerate() {
            let base_powers = base_powers.borrow();
            let bits = bits.to_bits(&mut cs.ns(|| format!("Convert Scalar {} to bits", i)))?;
            result.precomputed_base_scalar_mul(
                cs.ns(|| format!("Chunk {}", i)),
                bits.iter().zip(base_powers),
            )?;
        }
        Ok(result)
    }

    fn cost_of_add() -> usize;

    fn cost_of_double() -> usize;
}

#[cfg(test)]
mod test {
    use algebra::PairingEngine;
    use snark::ConstraintSystem;

    use crate::{
        groups::GroupGadget, test_constraint_system::TestConstraintSystem, utils::AllocGadget,
    };
    use algebra::groups::Group;
    use rand;

    pub(crate) fn group_test<
        E: PairingEngine,
        G: Group,
        GG: GroupGadget<G, E>,
        CS: ConstraintSystem<E>,
    >(
        cs: &mut CS,
        a: GG,
        b: GG,
    ) {
        let zero = GG::zero(cs.ns(|| "Zero")).unwrap();
        assert_eq!(zero, zero);

        // a == a
        assert_eq!(a, a);
        // a + 0 = a
        assert_eq!(a.add(cs.ns(|| "a_plus_zero"), &zero).unwrap(), a);
        // a - 0 = a
        assert_eq!(a.sub(cs.ns(|| "a_minus_zero"), &zero).unwrap(), a);
        // a - a = 0
        assert_eq!(a.sub(cs.ns(|| "a_minus_a"), &a).unwrap(), zero);
        // a + b = b + a
        let a_b = a.add(cs.ns(|| "a_plus_b"), &b).unwrap();
        let b_a = b.add(cs.ns(|| "b_plus_a"), &a).unwrap();
        assert_eq!(a_b, b_a);
        // (a + b) + a = a + (b + a)
        let ab_a = a_b.add(&mut cs.ns(|| "a_b_plus_a"), &a).unwrap();
        let a_ba = a.add(&mut cs.ns(|| "a_plus_b_a"), &b_a).unwrap();
        assert_eq!(ab_a, a_ba);
        // a.double() = a + a
        let a_a = a.add(cs.ns(|| "a + a"), &a).unwrap();
        let mut a2 = a.clone();
        a2.double_in_place(cs.ns(|| "2a")).unwrap();
        assert_eq!(a2, a_a);
        // b.double() = b + b
        let mut b2 = b.clone();
        b2.double_in_place(cs.ns(|| "2b")).unwrap();
        let b_b = b.add(cs.ns(|| "b + b"), &b).unwrap();
        assert_eq!(b2, b_b);

        let _ = a.to_bytes(&mut cs.ns(|| "ToBytes")).unwrap();
        let _ = a.to_bytes_strict(&mut cs.ns(|| "ToBytes Strict")).unwrap();

        let _ = b.to_bytes(&mut cs.ns(|| "b ToBytes")).unwrap();
        let _ = b
            .to_bytes_strict(&mut cs.ns(|| "b ToBytes Strict"))
            .unwrap();
    }

    #[test]
    fn jubjub_group_gadgets_test() {
        use crate::groups::jubjub::JubJubGadget;
        use algebra::curves::{bls12_381::Bls12_381, jubjub::JubJubProjective};

        let mut cs = TestConstraintSystem::<Bls12_381>::new();

        let a: JubJubProjective = rand::random();
        let b: JubJubProjective = rand::random();

        let a = JubJubGadget::alloc(&mut cs.ns(|| "generate_a"), || Ok(a)).unwrap();
        let b = JubJubGadget::alloc(&mut cs.ns(|| "generate_b"), || Ok(b)).unwrap();
        group_test::<_, JubJubProjective, _, _>(&mut cs.ns(|| "GroupTest(a, b)"), a, b);
    }
}
