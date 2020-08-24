use crate::{
    boolean::Boolean,
    fields::{fp::FpVar, FieldVar},
    prelude::*,
    ToBitsGadget,
};
use algebra::PrimeField;
use core::cmp::Ordering;
use r1cs_core::{lc, SynthesisError, Variable};

impl<F: PrimeField> FpVar<F> {
    /// This function enforces the ordering between `self` and `other`. The
    /// constraint system will not be satisfied otherwise. If `self` should
    /// also be checked for equality, e.g. `self <= other` instead of `self < other`, set
    /// `should_also_check_quality` to `true`. This variant verifies `self` and `other`
    /// are `<= (p-1)/2`.
    pub fn enforce_cmp(
        &self,
        other: &FpVar<F>,
        ordering: Ordering,
        should_also_check_equality: bool,
    ) -> Result<(), SynthesisError> {
        let (left, right) = self.process_cmp_inputs(other, ordering, should_also_check_equality)?;
        left.enforce_smaller_than(&right)
    }

    /// This function enforces the ordering between `self` and `other`. The
    /// constraint system will not be satisfied otherwise. If `self` should
    /// also be checked for equality, e.g. `self <= other` instead of `self < other`, set
    /// `should_also_check_quality` to `true`. This variant assumes `self` and `other`
    /// are `<= (p-1)/2` and does not generate constraints to verify that.
    pub fn enforce_cmp_unchecked(
        &self,
        other: &FpVar<F>,
        ordering: Ordering,
        should_also_check_equality: bool,
    ) -> Result<(), SynthesisError> {
        let (left, right) = self.process_cmp_inputs(other, ordering, should_also_check_equality)?;
        left.enforce_smaller_than_unchecked(&right)
    }

    /// This function checks the ordering between `self` and `other`. It outputs self
    /// `Boolean` that contains the result - `1` if true, `0` otherwise. The
    /// constraint system will be satisfied in any case. If `self` should
    /// also be checked for equality, e.g. `self <= other` instead of `self < other`, set
    /// `should_also_check_quality` to `true`. This variant verifies `self` and `other`
    /// are `<= (p-1)/2`.
    pub fn is_cmp(
        &self,
        other: &FpVar<F>,
        ordering: Ordering,
        should_also_check_equality: bool,
    ) -> Result<Boolean<F>, SynthesisError> {
        let (left, right) = self.process_cmp_inputs(other, ordering, should_also_check_equality)?;
        left.is_smaller_than(&right)
    }

    /// This function checks the ordering between `self` and `other`. It outputs a
    /// `Boolean` that contains the result - `1` if true, `0` otherwise. The
    /// constraint system will be satisfied in any case. If `self` should
    /// also be checked for equality, e.g. `self <= other` instead of `self < other`, set
    /// `should_also_check_quality` to `true`. This variant assumes `self` and `other`
    /// are `<= (p-1)/2` and does not generate constraints to verify that.
    pub fn is_cmp_unchecked(
        &self,
        other: &FpVar<F>,
        ordering: Ordering,
        should_also_check_equality: bool,
    ) -> Result<Boolean<F>, SynthesisError> {
        let (left, right) = self.process_cmp_inputs(other, ordering, should_also_check_equality)?;
        left.is_smaller_than_unchecked(&right)
    }

    fn process_cmp_inputs(
        &self,
        other: &Self,
        ordering: Ordering,
        should_also_check_equality: bool,
    ) -> Result<(Self, Self), SynthesisError> {
        let (left, right) = match ordering {
            Ordering::Less => (self, other),
            Ordering::Greater => (other, self),
            Ordering::Equal => Err(SynthesisError::Unsatisfiable)?,
        };
        let right_for_check = if should_also_check_equality {
            right + F::one()
        } else {
            right.clone()
        };

        Ok((left.clone(), right_for_check))
    }

    // Helper function to enforce `self <= (p-1)/2`.
    pub fn enforce_smaller_or_equal_than_mod_minus_one_div_two(
        &self,
    ) -> Result<(), SynthesisError> {
        let _ = Boolean::enforce_smaller_or_equal_than(
            &self.to_bits()?,
            F::modulus_minus_one_div_two(),
        )?;
        Ok(())
    }

    /// Helper function to check `self < other` and output a result bit. This function
    /// verifies `self` and `other` are `<= (p-1)/2`.
    fn is_smaller_than(&self, other: &FpVar<F>) -> Result<Boolean<F>, SynthesisError> {
        self.enforce_smaller_or_equal_than_mod_minus_one_div_two()?;
        other.enforce_smaller_or_equal_than_mod_minus_one_div_two()?;
        self.is_smaller_than_unchecked(other)
    }

    /// Helper function to check `self < other` and output a result bit. This function
    /// assumes `self` and `other` are `<= (p-1)/2` and does not generate constraints
    /// to verify that.
    fn is_smaller_than_unchecked(&self, other: &FpVar<F>) -> Result<Boolean<F>, SynthesisError> {
        Ok((self - other).double()?.to_bits()?.last().unwrap().clone())
    }

    /// Helper function to enforce `self < other`. This function verifies `self` and `other`
    /// are `<= (p-1)/2`.
    fn enforce_smaller_than(&self, other: &FpVar<F>) -> Result<(), SynthesisError> {
        self.enforce_smaller_or_equal_than_mod_minus_one_div_two()?;
        other.enforce_smaller_or_equal_than_mod_minus_one_div_two()?;
        self.enforce_smaller_than_unchecked(other)
    }

    /// Helper function to enforce `self < other`. This function assumes `self` and `other`
    /// are `<= (p-1)/2` and does not generate constraints to verify that.
    fn enforce_smaller_than_unchecked(&self, other: &FpVar<F>) -> Result<(), SynthesisError> {
        let cs = [self, other].cs().unwrap();
        let is_smaller_than = self.is_smaller_than_unchecked(other)?;
        let lc_one = lc!() + Variable::One;
        cs.enforce_constraint(is_smaller_than.lc(), lc_one.clone(), lc_one)
    }
}

#[cfg(test)]
mod test {
    use rand::{Rng, SeedableRng};
    use rand_xorshift::XorShiftRng;
    use std::cmp::Ordering;

    use crate::{alloc::AllocVar, fields::fp::FpVar};
    use algebra::{bls12_381::Fr, PrimeField, UniformRand};
    use r1cs_core::ConstraintSystem;

    #[test]
    fn test_cmp() {
        let mut rng = &mut XorShiftRng::from_seed([
            0x5d, 0xbe, 0x62, 0x59, 0x8d, 0x31, 0x3d, 0x76, 0x32, 0x37, 0xdb, 0x17, 0xe5, 0xbc,
            0x06, 0x54,
        ]);
        fn rand_in_range<R: Rng>(rng: &mut R) -> Fr {
            let pminusonedivtwo: Fr = Fr::modulus_minus_one_div_two().into();
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
            let cs = ConstraintSystem::<Fr>::new_ref();
            let a = rand_in_range(&mut rng);
            let a_var = FpVar::<Fr>::new_witness(cs.ns("a"), || Ok(a)).unwrap();
            let b = rand_in_range(&mut rng);
            let b_var = FpVar::<Fr>::new_witness(cs.ns("b"), || Ok(b)).unwrap();

            match a.cmp(&b) {
                Ordering::Less => {
                    a_var.enforce_cmp(&b_var, Ordering::Less, false).unwrap();
                    a_var.enforce_cmp(&b_var, Ordering::Less, true).unwrap();
                }
                Ordering::Greater => {
                    a_var.enforce_cmp(&b_var, Ordering::Greater, false).unwrap();
                    a_var.enforce_cmp(&b_var, Ordering::Greater, true).unwrap();
                }
                _ => {}
            }

            if i == 0 {
                println!("number of constraints: {}", cs.num_constraints());
            }
            assert!(cs.is_satisfied().unwrap());
        }

        for _i in 0..10 {
            let cs = ConstraintSystem::<Fr>::new_ref();
            let a = rand_in_range(&mut rng);
            let a_var = FpVar::<Fr>::new_witness(cs.ns("a"), || Ok(a)).unwrap();
            let b = rand_in_range(&mut rng);
            let b_var = FpVar::<Fr>::new_witness(cs.ns("b"), || Ok(b)).unwrap();

            match b.cmp(&a) {
                Ordering::Less => {
                    a_var.enforce_cmp(&b_var, Ordering::Less, false).unwrap();
                    a_var.enforce_cmp(&b_var, Ordering::Less, true).unwrap();
                }
                Ordering::Greater => {
                    a_var.enforce_cmp(&b_var, Ordering::Greater, false).unwrap();
                    a_var.enforce_cmp(&b_var, Ordering::Greater, true).unwrap();
                }
                _ => {}
            }

            assert!(!cs.is_satisfied().unwrap());
        }

        for _i in 0..10 {
            let cs = ConstraintSystem::<Fr>::new_ref();
            let a = rand_in_range(&mut rng);
            let a_var = FpVar::<Fr>::new_witness(cs.ns("a"), || Ok(a)).unwrap();
            a_var.enforce_cmp(&a_var, Ordering::Less, false).unwrap();

            assert!(!cs.is_satisfied().unwrap());
        }

        for _i in 0..10 {
            let cs = ConstraintSystem::<Fr>::new_ref();
            let a = rand_in_range(&mut rng);
            let a_var = FpVar::<Fr>::new_witness(cs.ns("a"), || Ok(a)).unwrap();
            a_var.enforce_cmp(&a_var, Ordering::Less, true).unwrap();
            assert!(cs.is_satisfied().unwrap());
        }
    }
}
