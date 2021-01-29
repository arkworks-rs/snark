use crate::{prelude::*, fields::{
    fp::FpGadget,
    nonnative::{
        params::get_params,
        reduce::{bigint_to_constraint_field, limbs_to_bigint, Reducer},
        nonnative_field_gadget::NonNativeFieldGadget,
    }
}, overhead, FromGadget};
use algebra::fields::{FpParameters, PrimeField};
use r1cs_core::{ConstraintSystem, SynthesisError};
use std::{
    marker::PhantomData,
    vec::Vec,
};
use num_bigint::BigUint;

#[derive(Debug)]
#[must_use]
pub struct NonNativeFieldMulResultGadget<SimulationF: PrimeField, ConstraintF: PrimeField> {
    /// Limbs of the intermediate representations
    pub limbs: Vec<FpGadget<ConstraintF>>,
    /// The cumulative num of additions
    pub prod_of_num_of_additions: ConstraintF,
    #[doc(hidden)]
    pub simulation_phantom: PhantomData<SimulationF>,
}

impl<SimulationF: PrimeField, ConstraintF: PrimeField> FromGadget<&NonNativeFieldGadget<SimulationF, ConstraintF>, ConstraintF>
for NonNativeFieldMulResultGadget<SimulationF, ConstraintF>
{
    fn from<CS: ConstraintSystem<ConstraintF>>(
        other: &NonNativeFieldGadget<SimulationF, ConstraintF>,
        cs: CS
    ) -> Result<Self, SynthesisError> {
        let params = get_params(SimulationF::size_in_bits(), ConstraintF::size_in_bits());

        let mut limbs = other.limbs.clone();
        limbs.reverse();
        limbs.resize(2 * params.num_limbs - 1, FpGadget::<ConstraintF>::zero(cs)?);
        limbs.reverse();

        let prod_of_num_of_additions = other.num_of_additions_over_normal_form + &ConstraintF::one();

        Ok(Self {
            limbs,
            prod_of_num_of_additions,
            simulation_phantom: PhantomData,
        })
    }
}

impl<SimulationF: PrimeField, ConstraintF: PrimeField> NonNativeFieldMulResultGadget<SimulationF, ConstraintF>
{
    /// Get the value of the multiplication result
    pub fn value(&self) -> Result<SimulationF, SynthesisError> {
        let params = get_params(SimulationF::size_in_bits(), ConstraintF::size_in_bits());

        let p_representations =
            NonNativeFieldGadget::<SimulationF, ConstraintF>::get_limbs_representations_from_big_integer(
                &<SimulationF as PrimeField>::Params::MODULUS
            )?;
        let p_bigint = limbs_to_bigint(params.bits_per_limb, &p_representations);

        let mut limbs_values = Vec::<ConstraintF>::new();
        for limb in self.limbs.iter() {
            limbs_values.push(limb.get_value().unwrap_or_default());
        }
        let value_bigint = limbs_to_bigint(params.bits_per_limb, &limbs_values);

        let res = bigint_to_constraint_field::<SimulationF>(&(value_bigint % p_bigint));
        Ok(res)
    }

    /// Constraints for reducing the result of a multiplication mod p, to get an original representation.
    pub fn reduce<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<NonNativeFieldGadget<SimulationF, ConstraintF>, SynthesisError>
    {
        let params = get_params(SimulationF::size_in_bits(), ConstraintF::size_in_bits());

        // Step 1: get p
        let p_representations =
            NonNativeFieldGadget::<SimulationF, ConstraintF>::get_limbs_representations_from_big_integer(
                &<SimulationF as PrimeField>::Params::MODULUS,
            )?;
        let p_bigint = limbs_to_bigint(params.bits_per_limb, &p_representations);

        let mut p_gadget_limbs = Vec::new();
        for (i, limb) in p_representations.iter().enumerate() {
            p_gadget_limbs.push(FpGadget::<ConstraintF>::from_value(
                cs.ns(|| format!("hardcode limb {}", i)),
                limb
            ));
        }
        let p_gadget = NonNativeFieldGadget::<SimulationF, ConstraintF> {
            limbs: p_gadget_limbs,
            num_of_additions_over_normal_form: ConstraintF::one(),
            is_in_the_normal_form: false,
            simulation_phantom: PhantomData,
        };

        // Step 2: compute surfeit
        let surfeit = overhead!(self.prod_of_num_of_additions + ConstraintF::one()) + 1 + 1;

        // Step 3: allocate k
        let k_bits = {
            let mut res = Vec::new();

            let mut limbs_values = Vec::<ConstraintF>::new();
            for limb in self.limbs.iter() {
                limbs_values.push(limb.get_value().unwrap_or_default());
            }

            let value_bigint = limbs_to_bigint(params.bits_per_limb, &limbs_values);
            let mut k_cur = value_bigint / p_bigint; // drops the remainder

            let total_len = SimulationF::size_in_bits() + surfeit;

            for i in 0..total_len {
                res.push(Boolean::alloc(
                    cs.ns(|| format!("alloc k bit {}", i)),
                    || {
                    Ok(&k_cur % 2u64 == BigUint::from(1u64))
                })?);
                k_cur /= 2u64; // drops the remainder
            }
            res
        };

        let k_limbs = {
            let zero = FpGadget::zero(cs.ns(|| "hardcode zero for k_limbs"))?;
            let mut limbs = Vec::new();

            let mut k_bits_cur = k_bits.clone();

            for i in 0..params.num_limbs {
                let this_limb_size = if i != params.num_limbs - 1 {
                    params.bits_per_limb
                } else {
                    k_bits.len() - (params.num_limbs - 1) * params.bits_per_limb
                };

                let this_limb_bits = k_bits_cur[0..this_limb_size].to_vec();
                k_bits_cur = k_bits_cur[this_limb_size..].to_vec();

                let mut limb = zero.clone();
                let mut cur = ConstraintF::one();

                for (j, bit) in this_limb_bits.iter().enumerate() {
                    limb = limb.conditionally_add_constant(
                        cs.ns(|| format!("add bit {} for limb {}", j, i)),
                        bit,
                        cur
                    )?;
                    cur.double_in_place();
                }
                limbs.push(limb);
            }

            limbs.reverse();
            limbs
        };

        let k_gadget = NonNativeFieldGadget::<SimulationF, ConstraintF> {
            limbs: k_limbs,
            num_of_additions_over_normal_form: self.prod_of_num_of_additions,
            is_in_the_normal_form: false,
            simulation_phantom: PhantomData,
        };

        let r_gadget = NonNativeFieldGadget::<SimulationF, ConstraintF>::alloc(
            cs.ns(|| "alloc r"),
            || Ok(self.value()?),
        )?;

        let params = get_params(SimulationF::size_in_bits(), ConstraintF::size_in_bits());

        // Step 1: reduce `self` and `other` if necessary
        let mut prod_limbs = Vec::new();
        let zero = FpGadget::<ConstraintF>::zero(cs.ns(|| "hardcode zero for step 1"))?;

        for _ in 0..2 * params.num_limbs - 1 {
            prod_limbs.push(zero.clone());
        }

        for i in 0..params.num_limbs {
            for j in 0..params.num_limbs {
                let mul_result = p_gadget.limbs[i].mul(
                    cs.ns(|| format!("mul_result = p_gadget.limbs[{}] * k_gadget.limbs[{}]", i, j)),
                    &k_gadget.limbs[j]
                )?;
                prod_limbs[i + j] = prod_limbs[i + j].add(
                    cs.ns(|| format!("prod_limbs[{},{}] = prod_limbs[{},{}] + mul_result", i, j, i, j)),
                    &mul_result
                )?;
            }
        }

        let mut kp_plus_r_gadget = Self {
            limbs: prod_limbs,
            prod_of_num_of_additions: (p_gadget.num_of_additions_over_normal_form
                + ConstraintF::one())
                * (k_gadget.num_of_additions_over_normal_form + ConstraintF::one()),
            simulation_phantom: PhantomData,
        };

        let kp_plus_r_limbs_len = kp_plus_r_gadget.limbs.len();
        for (i, limb) in r_gadget.limbs.iter().rev().enumerate() {
            kp_plus_r_gadget.limbs[kp_plus_r_limbs_len - 1 - i].add_in_place(
                cs.ns(|| format!("kp_plus_r_gadget.limbs[{}] + r_gadget.limbs_rev[{}]", kp_plus_r_limbs_len - 1 - i, i)),
                limb
            )?;
        }

        Reducer::<SimulationF, ConstraintF>::group_and_check_equality(
            cs.ns(|| "group and check equality"),
            surfeit,
            2 * params.bits_per_limb,
            params.bits_per_limb,
            &self.limbs,
            &kp_plus_r_gadget.limbs,
        )?;

        Ok(r_gadget)
    }

    /// Add unreduced elements.
    pub fn add<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self
    ) -> Result<Self, SynthesisError> {
        let mut new_limbs = Vec::new();

        for (i, (l1, l2)) in self.limbs.iter().zip(other.limbs.iter()).enumerate() {
            let new_limb = l1.add(cs.ns(|| format!("l1_{} + l2_{}", i, i)), l2)?;
            new_limbs.push(new_limb);
        }

        Ok(Self {
            limbs: new_limbs,
            prod_of_num_of_additions: self.prod_of_num_of_additions
                + other.prod_of_num_of_additions,
            simulation_phantom: PhantomData,
        })
    }

    /// Add native constant elem
    pub fn add_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &SimulationF
    ) -> Result<Self, SynthesisError> {
        let mut other_limbs =
            NonNativeFieldGadget::<SimulationF, ConstraintF>::get_limbs_representations(other)?;
        other_limbs.reverse();

        let mut new_limbs = Vec::new();

        for (i, limb) in self.limbs.iter().rev().enumerate() {
            if i < other_limbs.len() {
                let new_limb = limb.add_constant(
                    cs.ns(|| format!("limb_{} + other_limb_{}", i, i)),
                    &other_limbs[i]
                )?;
                new_limbs.push(new_limb);
            } else {
                new_limbs.push((*limb).clone());
            }
        }

        new_limbs.reverse();

        Ok(Self {
            limbs: new_limbs,
            prod_of_num_of_additions: self.prod_of_num_of_additions + ConstraintF::one(),
            simulation_phantom: PhantomData,
        })
    }
}
