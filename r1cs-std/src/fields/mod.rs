// use std::ops::{Mul, MulAssign};
use algebra::Field;
use r1cs_core::{ConstraintSystem, SynthesisError};
use std::fmt::Debug;

use crate::prelude::*;

pub mod fp;
pub mod fp12;
pub mod fp2;
pub mod fp6_3over2;

pub mod bls12_377;
pub mod edwards_bls12;
pub mod edwards_sw6;
pub mod jubjub;

pub trait FieldGadget<F: Field, ConstraintF: Field>:
    Sized
    + Clone
    + EqGadget<ConstraintF>
    + NEqGadget<ConstraintF>
    + ConditionalEqGadget<ConstraintF>
    + ToBitsGadget<ConstraintF>
    + AllocGadget<F, ConstraintF>
    + ToBytesGadget<ConstraintF>
    + CondSelectGadget<ConstraintF>
    + TwoBitLookupGadget<ConstraintF, TableConstant = F>
    + ThreeBitCondNegLookupGadget<ConstraintF, TableConstant = F>
    + Debug
{
    type Variable: Clone + Debug;

    fn get_value(&self) -> Option<F>;

    fn get_variable(&self) -> Self::Variable;

    fn zero<CS: ConstraintSystem<ConstraintF>>(_: CS) -> Result<Self, SynthesisError>;

    fn one<CS: ConstraintSystem<ConstraintF>>(_: CS) -> Result<Self, SynthesisError>;

    fn conditionally_add_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _: CS,
        _: &Boolean,
        _: F,
    ) -> Result<Self, SynthesisError>;

    fn add<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _: CS,
        _: &Self,
    ) -> Result<Self, SynthesisError>;

    fn add_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        cs: CS,
        other: &Self,
    ) -> Result<&mut Self, SynthesisError> {
        *self = self.add(cs, other)?;
        Ok(self)
    }

    fn double<CS: ConstraintSystem<ConstraintF>>(&self, cs: CS) -> Result<Self, SynthesisError> {
        self.add(cs, &self)
    }

    fn double_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        cs: CS,
    ) -> Result<&mut Self, SynthesisError> {
        *self = self.double(cs)?;
        Ok(self)
    }

    fn sub<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _: CS,
        _: &Self,
    ) -> Result<Self, SynthesisError>;

    fn sub_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        cs: CS,
        other: &Self,
    ) -> Result<&mut Self, SynthesisError> {
        *self = self.sub(cs, other)?;
        Ok(self)
    }

    fn negate<CS: ConstraintSystem<ConstraintF>>(&self, _: CS) -> Result<Self, SynthesisError>;

    #[inline]
    fn negate_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        cs: CS,
    ) -> Result<&mut Self, SynthesisError> {
        *self = self.negate(cs)?;
        Ok(self)
    }

    fn mul<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _: CS,
        _: &Self,
    ) -> Result<Self, SynthesisError>;

    fn mul_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        cs: CS,
        other: &Self,
    ) -> Result<&mut Self, SynthesisError> {
        *self = self.mul(cs, other)?;
        Ok(self)
    }

    fn square<CS: ConstraintSystem<ConstraintF>>(&self, cs: CS) -> Result<Self, SynthesisError> {
        self.mul(cs, &self)
    }

    fn square_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        cs: CS,
    ) -> Result<&mut Self, SynthesisError> {
        *self = self.square(cs)?;
        Ok(self)
    }

    fn mul_equals<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        result: &Self,
    ) -> Result<(), SynthesisError> {
        let actual_result = self.mul(cs.ns(|| "calc_actual_result"), other)?;
        result.enforce_equal(&mut cs.ns(|| "test_equals"), &actual_result)
    }

    fn square_equals<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        result: &Self,
    ) -> Result<(), SynthesisError> {
        let actual_result = self.square(cs.ns(|| "calc_actual_result"))?;
        result.enforce_equal(&mut cs.ns(|| "test_equals"), &actual_result)
    }

    fn add_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _: CS,
        _: &F,
    ) -> Result<Self, SynthesisError>;

    fn add_constant_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        cs: CS,
        other: &F,
    ) -> Result<&mut Self, SynthesisError> {
        *self = self.add_constant(cs, other)?;
        Ok(self)
    }

    fn sub_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        fe: &F,
    ) -> Result<Self, SynthesisError> {
        self.add_constant(cs, &(-(*fe)))
    }

    fn sub_constant_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        cs: CS,
        other: &F,
    ) -> Result<&mut Self, SynthesisError> {
        self.add_constant_in_place(cs, &(-(*other)))
    }

    fn mul_by_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _: CS,
        _: &F,
    ) -> Result<Self, SynthesisError>;

    fn mul_by_constant_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        cs: CS,
        other: &F,
    ) -> Result<&mut Self, SynthesisError> {
        *self = self.mul_by_constant(cs, other)?;
        Ok(self)
    }

    fn inverse<CS: ConstraintSystem<ConstraintF>>(&self, _: CS) -> Result<Self, SynthesisError>;

    fn frobenius_map<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _: CS,
        power: usize,
    ) -> Result<Self, SynthesisError>;

    fn frobenius_map_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        cs: CS,
        power: usize,
    ) -> Result<&mut Self, SynthesisError> {
        *self = self.frobenius_map(cs, power)?;
        Ok(self)
    }

    /// Accepts as input a list of bits which, when interpreted in big-endian
    /// form, are a scalar.
    #[inline]
    fn pow<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        bits: &[Boolean],
    ) -> Result<Self, SynthesisError> {
        let mut res = Self::one(cs.ns(|| "Alloc result"))?;
        for (i, bit) in bits.iter().enumerate() {
            res = res.square(cs.ns(|| format!("Double {}", i)))?;
            let tmp = res.mul(cs.ns(|| format!("Add {}-th base power", i)), self)?;
            res = Self::conditionally_select(
                cs.ns(|| format!("Conditional Select {}", i)),
                bit,
                &tmp,
                &res,
            )?;
        }
        Ok(res)
    }

    fn cost_of_mul() -> usize;

    fn cost_of_inv() -> usize;
}

#[cfg(test)]
mod test {
    use rand::{self, thread_rng, SeedableRng};
    use rand_xorshift::XorShiftRng;

    use crate::{prelude::*, test_constraint_system::TestConstraintSystem};
    use algebra::{BitIterator, Field, UniformRand};
    use r1cs_core::ConstraintSystem;

    fn field_test<
        FE: Field,
        ConstraintF: Field,
        F: FieldGadget<FE, ConstraintF>,
        CS: ConstraintSystem<ConstraintF>,
    >(
        mut cs: CS,
        a: F,
        b: F,
    ) {
        let a_native = a.get_value().unwrap();
        let b_native = b.get_value().unwrap();

        let zero = F::zero(cs.ns(|| "zero")).unwrap();
        let zero_native = zero.get_value().unwrap();
        zero.enforce_equal(&mut cs.ns(|| "zero_equals?"), &zero)
            .unwrap();
        assert_eq!(zero, zero);

        let one = F::one(cs.ns(|| "one")).unwrap();
        let one_native = one.get_value().unwrap();
        assert_eq!(one, one);
        one.enforce_equal(&mut cs.ns(|| "one_equals?"), &one)
            .unwrap();
        assert_ne!(one, zero);

        let one_dup = zero.add(cs.ns(|| "zero_plus_one"), &one).unwrap();
        one_dup
            .enforce_equal(&mut cs.ns(|| "one_plus_zero_equals"), &one)
            .unwrap();
        assert_eq!(one_dup, one);

        let two = one.add(cs.ns(|| "one_plus_one"), &one).unwrap();
        two.enforce_equal(&mut cs.ns(|| "two_equals?"), &two)
            .unwrap();
        assert_eq!(two, two);
        assert_ne!(zero, two);
        assert_ne!(one, two);

        // a == a
        assert_eq!(a, a);

        // a + 0 = a
        let a_plus_zero = a.add(cs.ns(|| "a_plus_zero"), &zero).unwrap();
        assert_eq!(a_plus_zero, a);
        assert_eq!(a_plus_zero.get_value().unwrap(), a_native);
        a_plus_zero
            .enforce_equal(&mut cs.ns(|| "a_plus_zero_equals?"), &a)
            .unwrap();

        // a - 0 = a
        let a_minus_zero = a.sub(cs.ns(|| "a_minus_zero"), &zero).unwrap();
        assert_eq!(a_minus_zero, a);
        assert_eq!(a_minus_zero.get_value().unwrap(), a_native);
        a_minus_zero
            .enforce_equal(&mut cs.ns(|| "a_minus_zero_equals?"), &a)
            .unwrap();

        // a - a = 0
        let a_minus_a = a.sub(cs.ns(|| "a_minus_a"), &a).unwrap();
        assert_eq!(a_minus_a, zero);
        assert_eq!(a_minus_a.get_value().unwrap(), zero_native);
        a_minus_a
            .enforce_equal(&mut cs.ns(|| "a_minus_a_equals?"), &zero)
            .unwrap();

        // a + b = b + a
        let a_b = a.add(cs.ns(|| "a_plus_b"), &b).unwrap();
        let b_a = b.add(cs.ns(|| "b_plus_a"), &a).unwrap();
        assert_eq!(a_b, b_a);
        assert_eq!(a_b.get_value().unwrap(), a_native + &b_native);
        a_b.enforce_equal(&mut cs.ns(|| "a+b == b+a"), &b_a)
            .unwrap();

        // (a + b) + a = a + (b + a)
        let ab_a = a_b.add(cs.ns(|| "a_b_plus_a"), &a).unwrap();
        let a_ba = a.add(cs.ns(|| "a_plus_b_a"), &b_a).unwrap();
        assert_eq!(ab_a, a_ba);
        assert_eq!(ab_a.get_value().unwrap(), a_native + &b_native + &a_native);
        ab_a.enforce_equal(&mut cs.ns(|| "a+b + a == a+ b+a"), &a_ba)
            .unwrap();

        let b_times_a_plus_b = a_b.mul(cs.ns(|| "b * (a + b)"), &b).unwrap();
        let b_times_b_plus_a = b_a.mul(cs.ns(|| "b * (b + a)"), &b).unwrap();
        assert_eq!(b_times_b_plus_a, b_times_a_plus_b);
        assert_eq!(
            b_times_a_plus_b.get_value().unwrap(),
            b_native * &(b_native + &a_native)
        );
        assert_eq!(
            b_times_a_plus_b.get_value().unwrap(),
            (b_native + &a_native) * &b_native
        );
        assert_eq!(
            b_times_a_plus_b.get_value().unwrap(),
            (a_native + &b_native) * &b_native
        );
        b_times_b_plus_a
            .enforce_equal(&mut cs.ns(|| "b*(a+b) == b * (b+a)"), &b_times_a_plus_b)
            .unwrap();

        // a * 0 = 0
        assert_eq!(a.mul(cs.ns(|| "a_times_zero"), &zero).unwrap(), zero);

        // a * 1 = a
        assert_eq!(a.mul(cs.ns(|| "a_times_one"), &one).unwrap(), a);
        assert_eq!(
            a.mul(cs.ns(|| "a_times_one2"), &one)
                .unwrap()
                .get_value()
                .unwrap(),
            a_native * &one_native
        );

        // a * b = b * a
        let ab = a.mul(cs.ns(|| "a_times_b"), &b).unwrap();
        let ba = b.mul(cs.ns(|| "b_times_a"), &a).unwrap();
        assert_eq!(ab, ba);
        assert_eq!(ab.get_value().unwrap(), a_native * &b_native);

        // (a * b) * a = a * (b * a)
        let ab_a = ab.mul(cs.ns(|| "ab_times_a"), &a).unwrap();
        let a_ba = a.mul(cs.ns(|| "a_times_ba"), &ba).unwrap();
        assert_eq!(ab_a, a_ba);
        assert_eq!(ab_a.get_value().unwrap(), a_native * &b_native * &a_native);

        let aa = a.mul(cs.ns(|| "a * a"), &a).unwrap();
        let a_squared = a.square(cs.ns(|| "a^2")).unwrap();
        a_squared
            .enforce_equal(&mut cs.ns(|| "a^2 == a*a"), &aa)
            .unwrap();
        assert_eq!(aa, a_squared);
        assert_eq!(aa.get_value().unwrap(), a_native.square());

        let aa = a
            .mul_by_constant(cs.ns(|| "a * a via mul_by_const"), &a.get_value().unwrap())
            .unwrap();
        a_squared
            .enforce_equal(&mut cs.ns(|| "a^2 == a*a via mul_by_const"), &aa)
            .unwrap();
        assert_eq!(aa, a_squared);
        assert_eq!(aa.get_value().unwrap(), a_native.square());

        let a_b2 = a
            .add_constant(cs.ns(|| "a + b via add_const"), &b.get_value().unwrap())
            .unwrap();
        a_b.enforce_equal(&mut cs.ns(|| "a + b == a + b via add_const"), &a_b2)
            .unwrap();
        assert_eq!(a_b, a_b2);

        let a_inv = a.inverse(cs.ns(|| "a_inv")).unwrap();
        a_inv
            .mul_equals(cs.ns(|| "check_equals"), &a, &one)
            .unwrap();
        assert_eq!(
            a_inv.get_value().unwrap(),
            a.get_value().unwrap().inverse().unwrap()
        );
        assert_eq!(a_inv.get_value().unwrap(), a_native.inverse().unwrap());
        // a * a * a = a^3
        let bits = BitIterator::new([0x3])
            .map(Boolean::constant)
            .collect::<Vec<_>>();
        assert_eq!(
            a_native * &(a_native * &a_native),
            a.pow(cs.ns(|| "test_pow"), &bits)
                .unwrap()
                .get_value()
                .unwrap()
        );

        // a * a * a = a^3
        let mut constants = [FE::zero(); 4];
        for c in &mut constants {
            *c = UniformRand::rand(&mut thread_rng());
            println!("Current c[i]: {:?}", c);
        }
        let bits = [Boolean::constant(false), Boolean::constant(true)];
        let lookup_result =
            F::two_bit_lookup(cs.ns(|| "Lookup"), &bits, constants.as_ref()).unwrap();
        assert_eq!(lookup_result.get_value().unwrap(), constants[2]);

        let negone: FE = UniformRand::rand(&mut thread_rng());

        let n = F::alloc(&mut cs.ns(|| "alloc new var"), || Ok(negone)).unwrap();
        let _ = n.to_bytes(&mut cs.ns(|| "ToBytes")).unwrap();
        let _ = n.to_bytes_strict(&mut cs.ns(|| "ToBytes Strict")).unwrap();

        let ab_false = a
            .conditionally_add_constant(
                cs.ns(|| "Add bool with coeff false"),
                &Boolean::constant(false),
                b_native,
            )
            .unwrap();
        assert_eq!(ab_false.get_value().unwrap(), a_native);
        let ab_true = a
            .conditionally_add_constant(
                cs.ns(|| "Add bool with coeff true"),
                &Boolean::constant(true),
                b_native,
            )
            .unwrap();
        assert_eq!(ab_true.get_value().unwrap(), a_native + &b_native);
    }

    fn random_frobenius_tests<
        FE: Field,
        ConstraintF: Field,
        F: FieldGadget<FE, ConstraintF>,
        CS: ConstraintSystem<ConstraintF>,
    >(
        mut cs: CS,
        maxpower: usize,
    ) {
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
        for i in 0..=maxpower {
            let mut a = FE::rand(&mut rng);
            let mut a_gadget = F::alloc(cs.ns(|| format!("a_gadget_{:?}", i)), || Ok(a)).unwrap();
            a_gadget = a_gadget
                .frobenius_map(cs.ns(|| format!("frob_map_{}", i)), i)
                .unwrap();
            a.frobenius_map(i);

            assert_eq!(a_gadget.get_value().unwrap(), a);
        }
    }

    #[test]
    fn bls12_377_field_gadgets_test() {
        use crate::fields::bls12_377::{Fq12Gadget, Fq2Gadget, Fq6Gadget, FqGadget};
        use algebra::fields::bls12_377::{Fq, Fq12, Fq2, Fq6};

        let mut cs = TestConstraintSystem::<Fq>::new();

        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        let a = FqGadget::alloc(&mut cs.ns(|| "generate_a"), || Ok(Fq::rand(&mut rng))).unwrap();
        let b = FqGadget::alloc(&mut cs.ns(|| "generate_b"), || Ok(Fq::rand(&mut rng))).unwrap();
        field_test(cs.ns(|| "test_fq"), a, b);
        if !cs.is_satisfied() {
            println!("{:?}", cs.which_is_unsatisfied().unwrap());
        }

        let c = Fq2Gadget::alloc(&mut cs.ns(|| "generate_c"), || Ok(Fq2::rand(&mut rng))).unwrap();
        let d = Fq2Gadget::alloc(&mut cs.ns(|| "generate_d"), || Ok(Fq2::rand(&mut rng))).unwrap();
        field_test(cs.ns(|| "test_fq2"), c, d);
        random_frobenius_tests::<Fq2, _, Fq2Gadget, _>(cs.ns(|| "test_frob_fq2"), 13);
        if !cs.is_satisfied() {
            println!("{:?}", cs.which_is_unsatisfied().unwrap());
        }

        let a = Fq6Gadget::alloc(&mut cs.ns(|| "generate_e"), || Ok(Fq6::rand(&mut rng))).unwrap();
        let b = Fq6Gadget::alloc(&mut cs.ns(|| "generate_f"), || Ok(Fq6::rand(&mut rng))).unwrap();
        field_test(cs.ns(|| "test_fq6"), a, b);
        random_frobenius_tests::<Fq6, _, Fq6Gadget, _>(cs.ns(|| "test_frob_fq6"), 13);
        if !cs.is_satisfied() {
            println!("{:?}", cs.which_is_unsatisfied().unwrap());
        }

        let c =
            Fq12Gadget::alloc(&mut cs.ns(|| "generate_g"), || Ok(Fq12::rand(&mut rng))).unwrap();
        let d =
            Fq12Gadget::alloc(&mut cs.ns(|| "generate_h"), || Ok(Fq12::rand(&mut rng))).unwrap();
        field_test(cs.ns(|| "test_fq12"), c, d);
        random_frobenius_tests::<Fq12, _, Fq12Gadget, _>(cs.ns(|| "test_frob_fq12"), 13);
        if !cs.is_satisfied() {
            println!("Here!");
            println!("{:?}", cs.which_is_unsatisfied().unwrap());
        }

        assert!(cs.is_satisfied());
    }

    #[test]
    fn jubjub_field_gadgets_test() {
        use crate::fields::jubjub::FqGadget;
        use algebra::fields::jubjub::fq::Fq;

        let mut cs = TestConstraintSystem::<Fq>::new();

        let mut rng = thread_rng();

        let a = FqGadget::alloc(&mut cs.ns(|| "generate_a"), || Ok(Fq::rand(&mut rng))).unwrap();
        let b = FqGadget::alloc(&mut cs.ns(|| "generate_b"), || Ok(Fq::rand(&mut rng))).unwrap();
        field_test(cs.ns(|| "test_fq"), a, b);
        if !cs.is_satisfied() {
            println!("{:?}", cs.which_is_unsatisfied().unwrap());
        }
        assert!(cs.is_satisfied());
    }

    #[test]
    fn edwards_field_gadgets_test() {
        use crate::fields::edwards_bls12::FqGadget;
        use algebra::fields::edwards_bls12::fq::Fq;

        let mut cs = TestConstraintSystem::<Fq>::new();

        let mut rng = thread_rng();

        let a = FqGadget::alloc(&mut cs.ns(|| "generate_a"), || Ok(Fq::rand(&mut rng))).unwrap();
        let b = FqGadget::alloc(&mut cs.ns(|| "generate_b"), || Ok(Fq::rand(&mut rng))).unwrap();
        field_test(cs.ns(|| "test_fq"), a, b);
        if !cs.is_satisfied() {
            println!("{:?}", cs.which_is_unsatisfied().unwrap());
        }
        assert!(cs.is_satisfied());
    }
}
