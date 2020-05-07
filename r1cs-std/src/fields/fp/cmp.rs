use crate::{
    boolean::Boolean,
    fields::{fp::FpGadget, FieldGadget},
    ToBitsGadget,
};
use algebra::PrimeField;
use core::cmp::Ordering;
use r1cs_core::{ConstraintSystem, SynthesisError};

impl<F: PrimeField> FpGadget<F> {
    /// This function enforces the ordering between `self` and `b`. The
    /// constraint system will not be satisfied otherwise. If `self` should
    /// also be checked for equality, e.g. `a <= b` instead of `a < b`, set
    /// `should_also_check_quality` to `true`. This variant verifies `a` and `b`
    /// are `<= (p-1)/2`.
    pub fn enforce_cmp<CS: ConstraintSystem<F>>(
        &self,
        mut cs: CS,
        b: &FpGadget<F>,
        ordering: Ordering,
        should_also_check_equality: bool,
    ) -> Result<(), SynthesisError> {
        let (left, right) = Self::process_cmp_inputs(
            cs.ns(|| "process cmp inputs"),
            &self,
            b,
            ordering,
            should_also_check_equality,
        )?;
        Self::enforce_smaller_than(cs.ns(|| "enforce smaller than"), &left, &right)
    }

    /// This function enforces the ordering between `self` and `b`. The
    /// constraint system will not be satisfied otherwise. If `self` should
    /// also be checked for equality, e.g. `a <= b` instead of `a < b`, set
    /// `should_also_check_quality` to `true`. This variant assumes `a` and `b`
    /// are `<= (p-1)/2` and does not generate constraints to verify that.
    pub fn enforce_cmp_unchecked<CS: ConstraintSystem<F>>(
        &self,
        mut cs: CS,
        b: &FpGadget<F>,
        ordering: Ordering,
        should_also_check_equality: bool,
    ) -> Result<(), SynthesisError> {
        let (left, right) = Self::process_cmp_inputs(
            cs.ns(|| "process cmp inputs"),
            &self,
            b,
            ordering,
            should_also_check_equality,
        )?;
        Self::enforce_smaller_than_unchecked(cs.ns(|| "enforce smaller than"), &left, &right)
    }

    /// This function checks the ordering between `self` and `b`. It outputs a
    /// `Boolean` that contains the result - `1` if true, `0` otherwise. The
    /// constraint system will be satisfied in any case. If `self` should
    /// also be checked for equality, e.g. `a <= b` instead of `a < b`, set
    /// `should_also_check_quality` to `true`. This variant verifies `a` and `b`
    /// are `<= (p-1)/2`.
    pub fn is_cmp<CS: ConstraintSystem<F>>(
        &self,
        mut cs: CS,
        b: &FpGadget<F>,
        ordering: Ordering,
        should_also_check_equality: bool,
    ) -> Result<Boolean, SynthesisError> {
        let (left, right) = Self::process_cmp_inputs(
            cs.ns(|| "process cmp inputs"),
            &self,
            b,
            ordering,
            should_also_check_equality,
        )?;
        Self::is_smaller_than(cs.ns(|| "enforce smaller than"), &left, &right)
    }

    /// This function checks the ordering between `self` and `b`. It outputs a
    /// `Boolean` that contains the result - `1` if true, `0` otherwise. The
    /// constraint system will be satisfied in any case. If `self` should
    /// also be checked for equality, e.g. `a <= b` instead of `a < b`, set
    /// `should_also_check_quality` to `true`. This variant assumes `a` and `b`
    /// are `<= (p-1)/2` and does not generate constraints to verify that.
    pub fn is_cmp_unchecked<CS: ConstraintSystem<F>>(
        &self,
        mut cs: CS,
        b: &FpGadget<F>,
        ordering: Ordering,
        should_also_check_equality: bool,
    ) -> Result<Boolean, SynthesisError> {
        let (left, right) = Self::process_cmp_inputs(
            cs.ns(|| "process cmp inputs"),
            &self,
            b,
            ordering,
            should_also_check_equality,
        )?;
        Self::is_smaller_than_unchecked(cs.ns(|| "enforce smaller than"), &left, &right)
    }

    fn process_cmp_inputs<CS: ConstraintSystem<F>>(
        mut cs: CS,
        a: &FpGadget<F>,
        b: &FpGadget<F>,
        ordering: Ordering,
        should_also_check_equality: bool,
    ) -> Result<(FpGadget<F>, FpGadget<F>), SynthesisError> {
        let left;
        let right;
        match ordering {
            Ordering::Less => {
                left = a;
                right = b;
            }
            Ordering::Greater => {
                left = b;
                right = a;
            }
            Ordering::Equal => {
                return Err(SynthesisError::Unsatisfiable);
            }
        };
        let right_for_check = if should_also_check_equality {
            right.add_constant(cs.ns(|| "plus one"), &F::one())?
        } else {
            right.clone()
        };

        Ok((left.clone(), right_for_check))
    }

    // Helper function to enforce `a <= (p-1)/2`.
    pub fn enforce_smaller_or_equal_than_mod_minus_one_div_two<CS: ConstraintSystem<F>>(
        mut cs: CS,
        a: &FpGadget<F>,
    ) -> Result<(), SynthesisError> {
        let a_bits = a.to_bits(cs.ns(|| "a to bits"))?;
        Boolean::enforce_smaller_or_equal_than::<_, _, F, _>(
            cs.ns(|| "enforce smaller than modulus minus one div two"),
            &a_bits,
            F::modulus_minus_one_div_two(),
        )?;

        Ok(())
    }

    /// Helper function to check `a < b` and output a result bit. This function
    /// verifies `a` and `b` are `<= (p-1)/2`.
    fn is_smaller_than<CS: ConstraintSystem<F>>(
        mut cs: CS,
        a: &FpGadget<F>,
        b: &FpGadget<F>,
    ) -> Result<Boolean, SynthesisError> {
        Self::enforce_smaller_or_equal_than_mod_minus_one_div_two(cs.ns(|| "check a in range"), a)?;
        Self::enforce_smaller_or_equal_than_mod_minus_one_div_two(cs.ns(|| "check b in range"), b)?;
        Self::is_smaller_than_unchecked(cs.ns(|| "enforce smaller than"), a, b)
    }

    /// Helper function to check `a < b` and output a result bit. This function
    /// assumes `a` and `b` are `<= (p-1)/2` and does not generate constraints
    /// to verify that.
    fn is_smaller_than_unchecked<CS: ConstraintSystem<F>>(
        mut cs: CS,
        a: &FpGadget<F>,
        b: &FpGadget<F>,
    ) -> Result<Boolean, SynthesisError> {
        let two = F::one() + F::one();
        let d0 = a.sub(cs.ns(|| "a - b"), b)?;
        let d = d0.mul_by_constant(cs.ns(|| "mul 2"), &two)?;
        let d_bits = d.to_bits(cs.ns(|| "d to bits"))?;
        let d_bits_len = d_bits.len();
        Ok(d_bits[d_bits_len - 1])
    }

    /// Helper function to enforce `a < b`. This function verifies `a` and `b`
    /// are `<= (p-1)/2`.
    fn enforce_smaller_than<CS: ConstraintSystem<F>>(
        mut cs: CS,
        a: &FpGadget<F>,
        b: &FpGadget<F>,
    ) -> Result<(), SynthesisError> {
        Self::enforce_smaller_or_equal_than_mod_minus_one_div_two(cs.ns(|| "check a in range"), a)?;
        Self::enforce_smaller_or_equal_than_mod_minus_one_div_two(cs.ns(|| "check b in range"), b)?;
        Self::enforce_smaller_than_unchecked(cs.ns(|| "enforce smaller than"), a, b)
    }

    /// Helper function to enforce `a < b`. This function assumes `a` and `b`
    /// are `<= (p-1)/2` and does not generate constraints to verify that.
    fn enforce_smaller_than_unchecked<CS: ConstraintSystem<F>>(
        mut cs: CS,
        a: &FpGadget<F>,
        b: &FpGadget<F>,
    ) -> Result<(), SynthesisError> {
        let is_smaller_than = Self::is_smaller_than_unchecked(cs.ns(|| "is smaller than"), a, b)?;
        cs.enforce(
            || "enforce smaller than",
            |_| is_smaller_than.lc(CS::one(), F::one()),
            |lc| lc + (F::one(), CS::one()),
            |lc| lc + (F::one(), CS::one()),
        );

        Ok(())
    }
}

#[cfg(test)]
mod test {
    use rand::{Rng, SeedableRng};
    use rand_xorshift::XorShiftRng;
    use std::cmp::Ordering;

    use crate::{
        alloc::AllocGadget, fields::fp::FpGadget, test_constraint_system::TestConstraintSystem,
    };
    use algebra::{bls12_381::Fr, PrimeField, UniformRand};
    use r1cs_core::ConstraintSystem;

    #[test]
    fn test_cmp() {
        let mut rng = &mut XorShiftRng::from_seed([
            0x5d, 0xbe, 0x62, 0x59, 0x8d, 0x31, 0x3d, 0x76, 0x32, 0x37, 0xdb, 0x17, 0xe5, 0xbc,
            0x06, 0x54,
        ]);
        fn rand_in_range<R: Rng>(rng: &mut R) -> Fr {
            let pminusonedivtwo = Fr::from_repr(Fr::modulus_minus_one_div_two());
            let mut r;
            loop {
                r = Fr::rand(rng);
                if r <= pminusonedivtwo {
                    break;
                }
            }
            r
        }
        for i in 0..10 {
            let mut cs = TestConstraintSystem::<Fr>::new();
            let a = rand_in_range(&mut rng);
            let a_var = FpGadget::<Fr>::alloc(cs.ns(|| "a"), || Ok(a)).unwrap();
            let b = rand_in_range(&mut rng);
            let b_var = FpGadget::<Fr>::alloc(cs.ns(|| "b"), || Ok(b)).unwrap();

            match a.cmp(&b) {
                Ordering::Less => {
                    a_var
                        .enforce_cmp(cs.ns(|| "smaller than test"), &b_var, Ordering::Less, false)
                        .unwrap();
                    a_var
                        .enforce_cmp(
                            cs.ns(|| "smaller than test 2"),
                            &b_var,
                            Ordering::Less,
                            true,
                        )
                        .unwrap();
                }
                Ordering::Greater => {
                    a_var
                        .enforce_cmp(
                            cs.ns(|| "smaller than test"),
                            &b_var,
                            Ordering::Greater,
                            false,
                        )
                        .unwrap();
                    a_var
                        .enforce_cmp(
                            cs.ns(|| "smaller than test 2"),
                            &b_var,
                            Ordering::Greater,
                            true,
                        )
                        .unwrap();
                }
                _ => {}
            }

            if i == 0 {
                println!("number of constraints: {}", cs.num_constraints());
            }
            assert!(cs.is_satisfied());
        }

        for _i in 0..10 {
            let mut cs = TestConstraintSystem::<Fr>::new();
            let a = rand_in_range(&mut rng);
            let a_var = FpGadget::<Fr>::alloc(cs.ns(|| "a"), || Ok(a)).unwrap();
            let b = rand_in_range(&mut rng);
            let b_var = FpGadget::<Fr>::alloc(cs.ns(|| "b"), || Ok(b)).unwrap();

            match b.cmp(&a) {
                Ordering::Less => {
                    a_var
                        .enforce_cmp(cs.ns(|| "smaller than test"), &b_var, Ordering::Less, false)
                        .unwrap();
                    a_var
                        .enforce_cmp(
                            cs.ns(|| "smaller than test 2"),
                            &b_var,
                            Ordering::Less,
                            true,
                        )
                        .unwrap();
                }
                Ordering::Greater => {
                    a_var
                        .enforce_cmp(
                            cs.ns(|| "smaller than test"),
                            &b_var,
                            Ordering::Greater,
                            false,
                        )
                        .unwrap();
                    a_var
                        .enforce_cmp(
                            cs.ns(|| "smaller than test 2"),
                            &b_var,
                            Ordering::Greater,
                            true,
                        )
                        .unwrap();
                }
                _ => {}
            }

            assert!(!cs.is_satisfied());
        }

        for _i in 0..10 {
            let mut cs = TestConstraintSystem::<Fr>::new();
            let a = rand_in_range(&mut rng);
            let a_var = FpGadget::<Fr>::alloc(cs.ns(|| "a"), || Ok(a)).unwrap();
            a_var
                .enforce_cmp(cs.ns(|| "smaller than test"), &a_var, Ordering::Less, false)
                .unwrap();

            assert!(!cs.is_satisfied());
        }

        for _i in 0..10 {
            let mut cs = TestConstraintSystem::<Fr>::new();
            let a = rand_in_range(&mut rng);
            let a_var = FpGadget::<Fr>::alloc(cs.ns(|| "a"), || Ok(a)).unwrap();
            a_var
                .enforce_cmp(
                    cs.ns(|| "smaller than or equal to test"),
                    &a_var,
                    Ordering::Less,
                    true,
                )
                .unwrap();

            assert!(cs.is_satisfied());
        }
    }
}
