use crate::{
    boolean::Boolean,
    fields::{fp::FpGadget, FieldGadget},
    ToBitsGadget,
};
use algebra::PrimeField;
use core::{cmp::Ordering, marker::PhantomData};
use r1cs_core::{ConstraintSystem, SynthesisError};

pub struct CmpGadget<ConstraintF: PrimeField> {
    constraint_field_type: PhantomData<ConstraintF>,
}

impl<ConstraintF: PrimeField> CmpGadget<ConstraintF> {
    fn process_cmp_inputs<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        a: &FpGadget<ConstraintF>,
        b: &FpGadget<ConstraintF>,
        ordering: Ordering,
        should_also_check_equality: bool,
    ) -> Result<(FpGadget<ConstraintF>, FpGadget<ConstraintF>), SynthesisError> {
        let left;
        let right;
        match ordering {
            Ordering::Less => {
                left = a;
                right = b;
            },
            Ordering::Greater => {
                left = b;
                right = a;
            },
            Ordering::Equal => {
                return Err(SynthesisError::Unsatisfiable);
            },
        };
        let right_for_check = if should_also_check_equality {
            right.add_constant(cs.ns(|| "plus one"), &ConstraintF::one())?
        } else {
            right.clone()
        };

        Ok((left.clone(), right_for_check))
    }

    fn check_smaller_than_mod_minus_one_div_two<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        a: &FpGadget<ConstraintF>,
    ) -> Result<(), SynthesisError> {
        let a_bits = a.to_bits(cs.ns(|| "a to bits"))?;
        Boolean::enforce_smaller_or_equal_than::<_, _, ConstraintF, _>(
            cs.ns(|| "enforce smaller than modulus minus one div two"),
            &a_bits,
            ConstraintF::modulus_minus_one_div_two(),
        )?;

        Ok(())
    }

    /// this function verifies a and b are <= (p-1)/2
    pub fn enforce_cmp<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        a: &FpGadget<ConstraintF>,
        b: &FpGadget<ConstraintF>,
        ordering: Ordering,
        should_also_check_equality: bool,
    ) -> Result<(), SynthesisError> {
        let (left, right) = Self::process_cmp_inputs(
            cs.ns(|| "process cmp inputs"),
            a,
            b,
            ordering,
            should_also_check_equality,
        )?;
        Self::enforce_smaller_than_unchecked(cs.ns(|| "enforce smaller than"), &left, &right)
    }

    /// this function assumes a and b are known to be <= (p-1)/2
    pub fn enforce_cmp_unchecked<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        a: &FpGadget<ConstraintF>,
        b: &FpGadget<ConstraintF>,
        ordering: Ordering,
        should_also_check_equality: bool,
    ) -> Result<(), SynthesisError> {
        let (left, right) = Self::process_cmp_inputs(
            cs.ns(|| "process cmp inputs"),
            a,
            b,
            ordering,
            should_also_check_equality,
        )?;
        Self::enforce_smaller_than(cs.ns(|| "enforce smaller than"), &left, &right)
    }

    /// this function verifies a and b are <= (p-1)/2
    pub fn is_cmp<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        a: &FpGadget<ConstraintF>,
        b: &FpGadget<ConstraintF>,
        ordering: Ordering,
        should_also_check_equality: bool,
    ) -> Result<Boolean, SynthesisError> {
        let (left, right) = Self::process_cmp_inputs(
            cs.ns(|| "process cmp inputs"),
            a,
            b,
            ordering,
            should_also_check_equality,
        )?;
        Self::is_smaller_than(cs.ns(|| "enforce smaller than"), &left, &right)
    }

    /// this function assumes a and b are known to be <= (p-1)/2
    pub fn is_cmp_unchecked<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        a: &FpGadget<ConstraintF>,
        b: &FpGadget<ConstraintF>,
        ordering: Ordering,
        should_also_check_equality: bool,
    ) -> Result<Boolean, SynthesisError> {
        let (left, right) = Self::process_cmp_inputs(
            cs.ns(|| "process cmp inputs"),
            a,
            b,
            ordering,
            should_also_check_equality,
        )?;
        Self::is_smaller_than_unchecked(cs.ns(|| "enforce smaller than"), &left, &right)
    }

    /// this function verifies a and b are <= (p-1)/2
    fn is_smaller_than<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        a: &FpGadget<ConstraintF>,
        b: &FpGadget<ConstraintF>,
    ) -> Result<Boolean, SynthesisError> {
        Self::check_smaller_than_mod_minus_one_div_two(cs.ns(|| "check a in range"), a)?;
        Self::check_smaller_than_mod_minus_one_div_two(cs.ns(|| "check b in range"), b)?;
        Self::is_smaller_than_unchecked(cs.ns(|| "enforce smaller than"), a, b)
    }

    /// this function assumes a and b are known to be <= (p-1)/2
    fn is_smaller_than_unchecked<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        a: &FpGadget<ConstraintF>,
        b: &FpGadget<ConstraintF>,
    ) -> Result<Boolean, SynthesisError> {
        let two = ConstraintF::one() + ConstraintF::one();
        let d0 = a.sub(cs.ns(|| "a - b"), b)?;
        let d = d0.mul_by_constant(cs.ns(|| "mul 2"), &two)?;
        let d_bits = d.to_bits_strict(cs.ns(|| "d to bits"))?;
        let d_bits_len = d_bits.len();
        Ok(d_bits[d_bits_len - 1])
    }

    /// this function verifies a and b are <= (p-1)/2
    fn enforce_smaller_than<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        a: &FpGadget<ConstraintF>,
        b: &FpGadget<ConstraintF>,
    ) -> Result<(), SynthesisError> {
        Self::check_smaller_than_mod_minus_one_div_two(cs.ns(|| "check a in range"), a)?;
        Self::check_smaller_than_mod_minus_one_div_two(cs.ns(|| "check b in range"), b)?;
        Self::enforce_smaller_than_unchecked(cs.ns(|| "enforce smaller than"), a, b)
    }

    /// this function assumes a and b are known to be <= (p-1)/2
    fn enforce_smaller_than_unchecked<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        a: &FpGadget<ConstraintF>,
        b: &FpGadget<ConstraintF>,
    ) -> Result<(), SynthesisError> {
        let is_smaller_than = Self::is_smaller_than_unchecked(cs.ns(|| "is smaller than"), a, b)?;
        cs.enforce(
            || "enforce smaller than",
            |_| is_smaller_than.lc(CS::one(), ConstraintF::one()),
            |lc| lc + (ConstraintF::one(), CS::one()),
            |lc| lc + (ConstraintF::one(), CS::one()),
        );

        Ok(())
    }
}

#[cfg(test)]
mod test {
    use rand::{Rng, SeedableRng};
    use rand_xorshift::XorShiftRng;
    use std::cmp::Ordering;

    use super::CmpGadget;
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
                    CmpGadget::<Fr>::enforce_cmp(
                        cs.ns(|| "smaller than test"),
                        &a_var,
                        &b_var,
                        Ordering::Less,
                        false,
                    )
                    .unwrap();
                    CmpGadget::<Fr>::enforce_cmp(
                        cs.ns(|| "smaller than test 2"),
                        &a_var,
                        &b_var,
                        Ordering::Less,
                        true,
                    )
                    .unwrap();
                },
                Ordering::Greater => {
                    CmpGadget::<Fr>::enforce_cmp(
                        cs.ns(|| "smaller than test"),
                        &a_var,
                        &b_var,
                        Ordering::Greater,
                        false,
                    )
                    .unwrap();
                    CmpGadget::<Fr>::enforce_cmp(
                        cs.ns(|| "smaller than test 2"),
                        &a_var,
                        &b_var,
                        Ordering::Greater,
                        true,
                    )
                    .unwrap();
                },
                _ => {},
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
                    CmpGadget::<Fr>::enforce_cmp(
                        cs.ns(|| "smaller than test"),
                        &a_var,
                        &b_var,
                        Ordering::Less,
                        false,
                    )
                    .unwrap();
                    CmpGadget::<Fr>::enforce_cmp(
                        cs.ns(|| "smaller than test 2"),
                        &a_var,
                        &b_var,
                        Ordering::Less,
                        true,
                    )
                    .unwrap();
                },
                Ordering::Greater => {
                    CmpGadget::<Fr>::enforce_cmp(
                        cs.ns(|| "smaller than test"),
                        &a_var,
                        &b_var,
                        Ordering::Greater,
                        false,
                    )
                    .unwrap();
                    CmpGadget::<Fr>::enforce_cmp(
                        cs.ns(|| "smaller than test 2"),
                        &a_var,
                        &b_var,
                        Ordering::Greater,
                        true,
                    )
                    .unwrap();
                },
                _ => {},
            }

            assert!(!cs.is_satisfied());
        }

        for _i in 0..10 {
            let mut cs = TestConstraintSystem::<Fr>::new();
            let a = rand_in_range(&mut rng);
            let a_var = FpGadget::<Fr>::alloc(cs.ns(|| "a"), || Ok(a)).unwrap();
            CmpGadget::<Fr>::enforce_cmp(
                cs.ns(|| "smaller than test"),
                &a_var,
                &a_var,
                Ordering::Less,
                false,
            )
            .unwrap();

            assert!(!cs.is_satisfied());
        }

        for _i in 0..10 {
            let mut cs = TestConstraintSystem::<Fr>::new();
            let a = rand_in_range(&mut rng);
            let a_var = FpGadget::<Fr>::alloc(cs.ns(|| "a"), || Ok(a)).unwrap();
            CmpGadget::<Fr>::enforce_cmp(
                cs.ns(|| "smaller than or equal to test"),
                &a_var,
                &a_var,
                Ordering::Less,
                true,
            )
            .unwrap();

            assert!(cs.is_satisfied());
        }
    }
}
