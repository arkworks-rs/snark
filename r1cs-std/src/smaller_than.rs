use crate::{
    boolean::Boolean,
    fields::{fp::FpGadget, FieldGadget},
    ToBitsGadget,
};
use algebra::PrimeField;
use r1cs_core::{ConstraintSystem, SynthesisError};
use std::marker::PhantomData;

pub struct SmallerThanGadget<ConstraintF: PrimeField> {
    constraint_field_type: PhantomData<ConstraintF>,
}

impl<ConstraintF: PrimeField> SmallerThanGadget<ConstraintF> {
    // the function assumes a and b are known to be <= (p-1)/2
    pub fn is_smaller_than<CS: ConstraintSystem<ConstraintF>>(
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

    // the function assumes a and b are known to be <= (p-1)/2
    pub fn is_smaller_than_or_equal_to<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        a: &FpGadget<ConstraintF>,
        b: &FpGadget<ConstraintF>,
    ) -> Result<Boolean, SynthesisError> {
        let b_plus_one = b.add_constant(cs.ns(|| "plus one"), &ConstraintF::one())?;
        Self::is_smaller_than(cs.ns(|| "is smaller than"), a, &b_plus_one)
    }

    // the function assumes a and b are known to be <= (p-1)/2
    pub fn enforce_smaller_than<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        a: &FpGadget<ConstraintF>,
        b: &FpGadget<ConstraintF>,
    ) -> Result<(), SynthesisError> {
        let is_smaller_than = Self::is_smaller_than(cs.ns(|| "is smaller than"), a, b)?;
        cs.enforce(
            || "enforce smaller than",
            |_| is_smaller_than.lc(CS::one(), ConstraintF::one()),
            |lc| lc + (ConstraintF::one(), CS::one()),
            |lc| lc + (ConstraintF::one(), CS::one()),
        );

        Ok(())
    }

    // the function assumes a and b are known to be <= (p-1)/2
    pub fn enforce_smaller_than_or_equal_to<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        a: &FpGadget<ConstraintF>,
        b: &FpGadget<ConstraintF>,
    ) -> Result<(), SynthesisError> {
        let b_plus_one = b.add_constant(cs.ns(|| "plus one"), &ConstraintF::one())?;
        Self::enforce_smaller_than(cs.ns(|| "enforce smaller than"), a, &b_plus_one)
    }

    pub fn enforce_smaller_than_strict<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        a: &FpGadget<ConstraintF>,
        b: &FpGadget<ConstraintF>,
    ) -> Result<(), SynthesisError> {
        let a_bits = a.to_bits(cs.ns(|| "a to bits"))?;
        Boolean::enforce_smaller_or_equal_than::<_, _, ConstraintF, _>(
            cs.ns(|| "enforce a smaller than modulus minus one div two"),
            &a_bits,
            ConstraintF::modulus_minus_one_div_two(),
        )?;
        let b_bits = b.to_bits(cs.ns(|| "b to bits"))?;
        Boolean::enforce_smaller_or_equal_than::<_, _, ConstraintF, _>(
            cs.ns(|| "enforce b smaller than modulus minus one div two"),
            &b_bits,
            ConstraintF::modulus_minus_one_div_two(),
        )?;

        let is_smaller_than = Self::is_smaller_than(cs.ns(|| "is smaller than"), a, b)?;
        cs.enforce(
            || "enforce smaller than",
            |_| is_smaller_than.lc(CS::one(), ConstraintF::one()),
            |lc| lc + (ConstraintF::one(), CS::one()),
            |lc| lc + (ConstraintF::one(), CS::one()),
        );

        Ok(())
    }

    pub fn enforce_smaller_than_or_equal_to_strict<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        a: &FpGadget<ConstraintF>,
        b: &FpGadget<ConstraintF>,
    ) -> Result<(), SynthesisError> {
        let b_plus_one = b.add_constant(cs.ns(|| "plus one"), &ConstraintF::one())?;
        Self::enforce_smaller_than_strict(cs.ns(|| "enforce smaller than strict"), a, &b_plus_one)
    }
}

#[cfg(test)]
mod test {
    use rand::{Rng, SeedableRng};
    use rand_xorshift::XorShiftRng;
    use std::cmp::Ordering;

    use super::SmallerThanGadget;
    use crate::{
        alloc::AllocGadget, fields::fp::FpGadget, test_constraint_system::TestConstraintSystem,
    };
    use algebra::{bls12_381::Fr, PrimeField, UniformRand};
    use r1cs_core::ConstraintSystem;

    #[test]
    fn test_smaller_than() {
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
                    SmallerThanGadget::<Fr>::enforce_smaller_than_strict(
                        cs.ns(|| "smaller than test"),
                        &a_var,
                        &b_var,
                    )
                    .unwrap();
                    SmallerThanGadget::<Fr>::enforce_smaller_than_or_equal_to_strict(
                        cs.ns(|| "smaller than test 2"),
                        &a_var,
                        &b_var,
                    )
                    .unwrap();
                },
                Ordering::Greater => {
                    SmallerThanGadget::<Fr>::enforce_smaller_than_strict(
                        cs.ns(|| "smaller than test"),
                        &b_var,
                        &a_var,
                    )
                    .unwrap();
                    SmallerThanGadget::<Fr>::enforce_smaller_than_or_equal_to_strict(
                        cs.ns(|| "smaller than test 2"),
                        &b_var,
                        &a_var,
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
                    SmallerThanGadget::<Fr>::enforce_smaller_than_strict(
                        cs.ns(|| "smaller than test"),
                        &a_var,
                        &b_var,
                    )
                    .unwrap();
                    SmallerThanGadget::<Fr>::enforce_smaller_than_or_equal_to_strict(
                        cs.ns(|| "smaller than test 2"),
                        &a_var,
                        &b_var,
                    )
                    .unwrap();
                },
                Ordering::Greater => {
                    SmallerThanGadget::<Fr>::enforce_smaller_than_strict(
                        cs.ns(|| "smaller than test"),
                        &b_var,
                        &a_var,
                    )
                    .unwrap();
                    SmallerThanGadget::<Fr>::enforce_smaller_than_or_equal_to(
                        cs.ns(|| "smaller than test 2"),
                        &b_var,
                        &a_var,
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
            SmallerThanGadget::<Fr>::enforce_smaller_than_strict(
                cs.ns(|| "smaller than test"),
                &a_var,
                &a_var,
            )
            .unwrap();

            assert!(!cs.is_satisfied());
        }

        for _i in 0..10 {
            let mut cs = TestConstraintSystem::<Fr>::new();
            let a = rand_in_range(&mut rng);
            let a_var = FpGadget::<Fr>::alloc(cs.ns(|| "a"), || Ok(a)).unwrap();
            SmallerThanGadget::<Fr>::enforce_smaller_than_or_equal_to(
                cs.ns(|| "smaller than or equal to test"),
                &a_var,
                &a_var,
            )
            .unwrap();

            assert!(cs.is_satisfied());
        }
    }
}
