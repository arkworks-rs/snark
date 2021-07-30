// use std::ops::{Mul, MulAssign};
use algebra::{Field, BitIterator};
use r1cs_core::{ConstraintSystem, SynthesisError};
use std::fmt::Debug;

use crate::{prelude::*, Assignment};

pub mod fp;
pub mod fp12;
pub mod fp2;
pub mod fp3;
pub mod fp4;
pub mod fp6_3over2;
pub mod fp6_2over3;
pub mod quadratic_extension;
pub mod cubic_extension;

#[cfg(feature = "nonnative")]
pub mod nonnative;

pub trait FieldGadget<F: Field, ConstraintF: Field>:
    Sized
    + Clone
    + EqGadget<ConstraintF>
    + ToBitsGadget<ConstraintF>
    + AllocGadget<F, ConstraintF>
    + ConstantGadget<F, ConstraintF>
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

    fn inverse<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Self, SynthesisError> {
        let one = Self::one(&mut cs.ns(|| "one"))?;
        let inverse = Self::alloc(&mut cs.ns(|| "alloc inverse"), || {
            self.get_value().and_then(|val| val.inverse()).get()
        })?;
        self.mul_equals(cs.ns(|| "check inv"), &inverse, &one)?;
        Ok(inverse)
    }

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
        for (i, bit) in bits.into_iter().enumerate() {
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

    /// Computes `self^S`, where S is interpreted as an little-endian
    /// u64-decomposition of an integer.
    fn pow_by_constant<S: AsRef<[u64]>, CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        exp: S
    ) -> Result<Self, SynthesisError> {
        let mut res = Self::one(cs.ns(|| "alloc result"))?;
        let mut found_one = false;

        for i in BitIterator::new(exp) {

            // Skip leading zeros
            if !found_one {
                if i {
                    found_one = true;
                } else {
                    continue;
                }
            }
            res.square_in_place(cs.ns(|| format!("square_{}", i)))?;

            if i {
                res.mul_in_place(cs.ns(|| format!("multiply_{}", i)), self)?;
            }
        }
        Ok(res)
    }

    fn cost_of_mul() -> usize;

    fn cost_of_mul_equals() -> usize;

    fn cost_of_inv() -> usize;
}

#[cfg(test)]
pub(crate) mod tests {
    use rand::{self, thread_rng, SeedableRng, Rng};
    use rand_xorshift::XorShiftRng;

    use crate::{prelude::*, test_constraint_system::TestConstraintSystem, fields::fp::FpGadget};
    use algebra::{BitIterator, Field, UniformRand, PrimeField, leading_zeros};
    use r1cs_core::ConstraintSystem;

    #[allow(dead_code)]
    pub(crate) fn field_test<
        FE: Field,
        ConstraintF: Field,
        F: FieldGadget<FE, ConstraintF>,
    >()
    {
        let mut cs = TestConstraintSystem::<ConstraintF>::new();

        let mut rng = &mut thread_rng();

        let a_native = FE::rand(&mut rng);
        let b_native = FE::rand(&mut rng);
        let a = F::alloc(&mut cs.ns(|| "generate_a"), || Ok(a_native)).unwrap();
        let b = F::alloc(&mut cs.ns(|| "generate_b"), || Ok(b_native)).unwrap();

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
            .map(|bit| Boolean::constant(bit))
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

    #[allow(dead_code)]
    pub(crate) fn frobenius_tests<
        FE: Field,
        ConstraintF: Field,
        F: FieldGadget<FE, ConstraintF>,
    >(maxpower: usize)
    {
        let mut cs = TestConstraintSystem::<ConstraintF>::new();
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
        for i in 0..(maxpower + 1) {
            let mut a = FE::rand(&mut rng);
            let mut a_gadget = F::alloc(cs.ns(|| format!("a_gadget_{:?}", i)), || Ok(a)).unwrap();
            a_gadget = a_gadget
                .frobenius_map(cs.ns(|| format!("frob_map_{}", i)), i)
                .unwrap();
            a.frobenius_map(i);

            assert_eq!(a_gadget.get_value().unwrap(), a);
        }
    }

    #[allow(dead_code)]
    pub(crate) fn from_bits_fp_gadget_test<ConstraintF: PrimeField>()
    {
        let mut rng = thread_rng();
        let mut cs = TestConstraintSystem::<ConstraintF>::new();

        // Sample a random field element with bit length MODULUS_BITS - 1
        // (Because `from_bits` pack only up until MODULUS_BITS - 1 bits)
        let (f, leading_zeros) = loop {
            let val = ConstraintF::rand(&mut rng);
            let zeros = leading_zeros(val.write_bits().as_slice());
            if zeros > 1 {
                break (val, zeros as usize)
            }
        };

        //Positive case
        let f_g_bits = Vec::<Boolean>::alloc(
            cs.ns(|| "alloc f bits"),
            || Ok(f.write_bits()[leading_zeros..].to_vec())
        ).unwrap();
        let f_g = FpGadget::<ConstraintF>::from_bits(
            cs.ns(|| "pack f_g_bits"),
            f_g_bits.as_slice()
        ).unwrap();
        assert_eq!(f, f_g.get_value().unwrap());
        assert!(cs.is_satisfied());

        //Let's alter one random bit and check that the cs is not satisfied anymore
        let random_bit: usize = rng.gen_range(leading_zeros..f_g_bits.len());
        let prev_value = f_g_bits[random_bit].get_value().unwrap();
        let new_value = if prev_value {
            ConstraintF::zero()
        } else {
            ConstraintF::one()
        };
        cs.set(format!("alloc f bits/value_{}/boolean", random_bit).as_ref(), new_value);
        assert!(!cs.is_satisfied());
        assert_eq!("pack f_g_bits/packing constraint", cs.which_is_unsatisfied().unwrap());

        //Let's change the value of the packed variable and check that the cs is not satisfied anymore

        //Bringing back the modified bit's value to its original one
        let prev_value = if prev_value {ConstraintF::one()} else {ConstraintF::zero()};
        cs.set(format!("alloc f bits/value_{}/boolean", random_bit).as_ref(), prev_value);
        assert!(cs.is_satisfied()); //Situation should be back to positive case

        //Modify packed value
        cs.set(format!("pack f_g_bits/variable/alloc").as_ref(), ConstraintF::rand(&mut rng));
        assert!(!cs.is_satisfied());
        assert_eq!("pack f_g_bits/packing constraint", cs.which_is_unsatisfied().unwrap());
    }

    #[allow(dead_code)]
    pub(crate) fn bit_fp_gadgets_test<ConstraintF: PrimeField>(){
        use crate::algebra::FpParameters;

        let mut rng = thread_rng();
        let mut cs = TestConstraintSystem::<ConstraintF>::new();

        //Native to_bits test
        let a = ConstraintF::rand(&mut rng);
        let a_g = FpGadget::<ConstraintF>::alloc(
            cs.ns(|| "alloc a"),
            || Ok(a),
        ).unwrap();

        let a_bits = a.write_bits();
        let a_g_bits = a_g.to_bits(cs.ns(|| "a_to_bits")).unwrap();
        assert_eq!(
            a_bits,
            a_g_bits.iter().map(|b| b.get_value().unwrap()).collect::<Vec<_>>(),
        );

        //Native from_bits test
        //Let's cut off one bit from both in order to be able to use from_bits gadget
        let a_bits = a_bits[1..].to_vec();
        let a_g_bits = a_g_bits[1..].as_ref();

        let a_read = ConstraintF::read_bits(a_bits).unwrap();
        let a_g_read = FpGadget::<ConstraintF>::from_bits(
            cs.ns(|| "read a_g"),
            a_g_bits,
        ).unwrap();

        assert_eq!(a_read, a_g_read.get_value().unwrap());

        //to_bits_with_length_restriction test
        let (b, leading_zeros) = loop {
            let val = ConstraintF::rand(&mut rng);
            let zeros = leading_zeros(val.write_bits().as_slice());
            if zeros >= 3 {
                break (val, zeros)
            }
        };

        let b_g = FpGadget::<ConstraintF>::alloc(
            cs.ns(|| "alloc b"),
            || Ok(b),
        ).unwrap();

        //Positive case
        let b_g_restricted_bits = b_g.to_bits_with_length_restriction(
            cs.ns(|| "serialize with length restriction"),
            leading_zeros as usize,
        ).unwrap();

        assert_eq!(b_g_restricted_bits.len() as u32, ConstraintF::Params::MODULUS_BITS - leading_zeros);

        //Of course we should be able to reconstruct the original field element
        let b_g_read = FpGadget::<ConstraintF>::from_bits(
            cs.ns(|| "read b_g_restricted"),
            b_g_restricted_bits.as_slice(),
        ).unwrap();
        b_g.enforce_equal(cs.ns(|| "should pass"), &b_g_read).unwrap();
        assert!(cs.is_satisfied());

        //If we cut off more bits we will reconstruct a different field element
        let bad_b_g_bits = b_g.to_bits_with_length_restriction(
            cs.ns(|| "serialize bad with length restriction"),
            leading_zeros as usize + 1,
        ).unwrap();

        let bad_b_g_read = FpGadget::<ConstraintF>::from_bits(
            cs.ns(|| "read bad_b_g_restricted"),
            bad_b_g_bits.as_slice(),
        ).unwrap();
        b_g.enforce_equal(cs.ns(|| "should not pass"), &bad_b_g_read).unwrap();
        assert!(!cs.is_satisfied());
    }

    #[allow(dead_code)]
    pub(crate) fn equ_verdict_fp_gadget_test<ConstraintF: PrimeField>()
    {
        let mut rng = thread_rng();
        let a = ConstraintF::rand(&mut rng);

        //Case a == b
        {
            let mut cs = TestConstraintSystem::<ConstraintF>::new();

            let a_gadget = FpGadget::<ConstraintF>::alloc(
                cs.ns(|| "alloc a"),
                || Ok(a)
            ).unwrap();

            //If a == b then v = True
            let b_gadget = FpGadget::<ConstraintF>::alloc(
                cs.ns(|| "alloc b"),
                || Ok(a.clone())
            ).unwrap();

            let v = a_gadget.is_eq(cs.ns(|| "a == b"), &b_gadget).unwrap();
            v.enforce_equal(cs.ns(|| " v == True"), &Boolean::constant(true)).unwrap();
            assert!(cs.is_satisfied());

            //If a == b but the prover maliciously witness v as False, cs will not be satisfied
            cs.set("a == b/alloc verdict/boolean", ConstraintF::zero());
            assert!(!cs.is_satisfied());
            assert_eq!("a == b/1 - v = c * (x - y)", cs.which_is_unsatisfied().unwrap());

            //If a == b the prover can freely choose c without invalidating any constraint
            cs.set("a == b/alloc verdict/boolean", ConstraintF::one()); //Let's bring back v to True
            assert!(cs.is_satisfied()); //Situation should be back to positive case
            cs.set("a == b/alloc c/alloc", ConstraintF::rand(&mut rng)); //Let's choose a random v
            assert!(cs.is_satisfied());
        }

        //Case a != b
        {
            let mut cs = TestConstraintSystem::<ConstraintF>::new();

            let a_gadget = FpGadget::<ConstraintF>::alloc(
                cs.ns(|| "alloc a"),
                || Ok(a)
            ).unwrap();

            //If a != b then v = False
            let b_gadget = FpGadget::<ConstraintF>::alloc(
                cs.ns(|| "alloc b"),
                || Ok(ConstraintF::rand(&mut rng))
            ).unwrap();

            let v = a_gadget.is_eq(cs.ns(|| "a != b"), &b_gadget).unwrap();
            v.enforce_equal(cs.ns(|| " v == False"), &Boolean::constant(false)).unwrap();
            assert!(cs.is_satisfied());

            //If a != b but the prover maliciously witness v as True, cs will not be satisfied
            cs.set("a != b/alloc verdict/boolean", ConstraintF::one());
            assert!(!cs.is_satisfied());
            assert_eq!("a != b/0 = v * (x - y)/conditional_equals", cs.which_is_unsatisfied().unwrap());

            //If a != b the prover is forced to choose c as 1/(a-b)
            cs.set("a != b/alloc verdict/boolean", ConstraintF::zero()); //Let's bring back v to False
            assert!(cs.is_satisfied()); //Situation should be back to normal
            cs.set("a != b/alloc c/alloc", ConstraintF::rand(&mut rng)); //Let's choose a random c
            assert!(!cs.is_satisfied());
            assert_eq!("a != b/1 - v = c * (x - y)", cs.which_is_unsatisfied().unwrap());
        }
    }
}