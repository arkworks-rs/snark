use algebra::Field;
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::prelude::*;
use r1cs_std::boolean::AllocatedBit;

use crate::merkle_tree::*;
use crate::crh::{FixedLengthCRH, FixedLengthCRHGadget};

use std::borrow::Borrow;

pub struct MerkleTreePathGadget<P, HGadget, ConstraintF>
where
    P: MerkleTreeConfig,
    HGadget: FixedLengthCRHGadget<P::H, ConstraintF>,
    ConstraintF: Field,
{
    path:    Vec<(HGadget::OutputGadget, HGadget::OutputGadget)>,
}

impl<P, CRHGadget, ConstraintF> MerkleTreePathGadget<P, CRHGadget, ConstraintF>
where
    P: MerkleTreeConfig,
    ConstraintF: Field,
    CRHGadget: FixedLengthCRHGadget<P::H, ConstraintF>,
{
    pub fn check_membership<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        parameters: &CRHGadget::ParametersGadget,
        root: &CRHGadget::OutputGadget,
        leaf: impl ToBytesGadget<ConstraintF>,
    ) -> Result<(), SynthesisError> {
    self.conditionally_check_membership(
            cs,
            parameters,
            root,
            leaf,
            &Boolean::Constant(true),
        )
    }

    pub fn conditionally_check_membership<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        parameters: &CRHGadget::ParametersGadget,
        root: &CRHGadget::OutputGadget,
        leaf: impl ToBytesGadget<ConstraintF>,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        assert_eq!(self.path.len(), P::HEIGHT - 1);
        // Check that the hash of the given leaf matches the leaf hash in the membership
        // proof.
        let leaf_bits = leaf.to_bytes(&mut cs.ns(|| "leaf_to_bytes"))?;
        let leaf_hash = CRHGadget::check_evaluation_gadget(
            cs.ns(|| "check_evaluation_gadget"),
            parameters,
            &leaf_bits,
        )?;

        // Check if leaf is one of the bottom-most siblings.
        let leaf_is_left = AllocatedBit::alloc(&mut cs.ns(|| "leaf_is_left"), || {
            Ok(leaf_hash == self.path[0].0)
        })?
        .into();
        CRHGadget::OutputGadget::conditional_enforce_equal_or(
            &mut cs.ns(|| "check_leaf_is_left"),
            &leaf_is_left,
            &leaf_hash,
            &self.path[0].0,
            &self.path[0].1,
            should_enforce,
        )?;

        // Check levels between leaf level and root.
        let mut previous_hash = leaf_hash;
        for (i, &(ref left_hash, ref right_hash)) in self.path.iter().enumerate() {
            // Check if the previous_hash matches the correct current hash.
            let previous_is_left =
                AllocatedBit::alloc(&mut cs.ns(|| format!("previous_is_left_{}", i)), || {
                    Ok(&previous_hash == left_hash)
                })?
                .into();

            CRHGadget::OutputGadget::conditional_enforce_equal_or(
                &mut cs.ns(|| format!("check_equals_which_{}", i)),
                &previous_is_left,
                &previous_hash,
                left_hash,
                right_hash,
                should_enforce,
            )?;

            previous_hash = hash_inner_node_gadget::<P::H, CRHGadget, ConstraintF, _>(
                &mut cs.ns(|| format!("hash_inner_node_{}", i)),
                parameters,
                left_hash,
                right_hash,
            )?;
        }

        root.conditional_enforce_equal(
            &mut cs.ns(|| "root_is_last"),
            &previous_hash,
            should_enforce,
        )
    }
}

pub(crate) fn hash_inner_node_gadget<H, HG, ConstraintF, CS>(
    mut cs: CS,
    parameters: &HG::ParametersGadget,
    left_child: &HG::OutputGadget,
    right_child: &HG::OutputGadget,
) -> Result<HG::OutputGadget, SynthesisError>
where
    ConstraintF: Field,
    CS: ConstraintSystem<ConstraintF>,
    H: FixedLengthCRH,
    HG: FixedLengthCRHGadget<H, ConstraintF>,
{
    let left_bytes = left_child.to_bytes(&mut cs.ns(|| "left_to_bytes"))?;
    let right_bytes = right_child.to_bytes(&mut cs.ns(|| "right_to_bytes"))?;
    let mut bytes = left_bytes;
    bytes.extend_from_slice(&right_bytes);

    HG::check_evaluation_gadget(cs, parameters, &bytes)
}

impl<P, HGadget, ConstraintF> AllocGadget<MerkleTreePath<P>, ConstraintF>
    for MerkleTreePathGadget<P, HGadget, ConstraintF>
where
    P: MerkleTreeConfig,
    HGadget: FixedLengthCRHGadget<P::H, ConstraintF>,
    ConstraintF: Field,
{
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<MerkleTreePath<P>>,
    {
        let mut path = Vec::new();
        for (i, &(ref l, ref r)) in value_gen()?.borrow().path.iter().enumerate() {
            let l_hash =
                HGadget::OutputGadget::alloc(&mut cs.ns(|| format!("l_child_{}", i)), || {
                    Ok(l.clone())
                })?;
            let r_hash =
                HGadget::OutputGadget::alloc(&mut cs.ns(|| format!("r_child_{}", i)), || {
                    Ok(r.clone())
                })?;
            path.push((l_hash, r_hash));
        }
        Ok(MerkleTreePathGadget {
            path,
        })
    }

    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<MerkleTreePath<P>>,
    {
        let mut path = Vec::new();
        for (i, &(ref l, ref r)) in value_gen()?.borrow().path.iter().enumerate() {
            let l_hash = HGadget::OutputGadget::alloc_input(
                &mut cs.ns(|| format!("l_child_{}", i)),
                || Ok(l.clone()),
            )?;
            let r_hash = HGadget::OutputGadget::alloc_input(
                &mut cs.ns(|| format!("r_child_{}", i)),
                || Ok(r.clone()),
            )?;
            path.push((l_hash, r_hash));
        }

        Ok(MerkleTreePathGadget {
            path,
        })
    }
}

#[cfg(test)]
mod test {
    use std::rc::Rc;

    use crate::{
        crh::{
            pedersen::{PedersenCRH, PedersenWindow},
            pedersen::constraints::PedersenCRHGadget, 
            FixedLengthCRH,
            FixedLengthCRHGadget,
        },
        merkle_tree::*,
    };
    use algebra::{
        curves::jubjub::JubJubAffine as JubJub,
        fields::jubjub::fq::Fq,
    };
    use rand::SeedableRng;
    use rand_xorshift::XorShiftRng;
    use r1cs_core::ConstraintSystem;

    use super::*;
    use r1cs_std::{
        groups::curves::twisted_edwards::jubjub::JubJubGadget,
        test_constraint_system::TestConstraintSystem,
    };

    #[derive(Clone)]
    pub(super) struct Window4x256;
    impl PedersenWindow for Window4x256 {
        const WINDOW_SIZE: usize = 4;
        const NUM_WINDOWS: usize = 256;
    }

    type H = PedersenCRH<JubJub, Window4x256>;
    type HG = PedersenCRHGadget<JubJub, Fq, JubJubGadget>;

    struct JubJubMerkleTreeParams;
    
    impl MerkleTreeConfig for JubJubMerkleTreeParams {
        const HEIGHT: usize = 32;
        type H = H;
    }

    type JubJubMerkleTree = MerkleHashTree<JubJubMerkleTreeParams>;

    fn generate_merkle_tree(leaves: &[[u8; 30]], use_bad_root: bool) -> () {
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);

        let crh_parameters = Rc::new(H::setup(&mut rng).unwrap());
        let tree = JubJubMerkleTree::new(crh_parameters.clone(), leaves).unwrap();
        let root = tree.root();
        let mut satisfied = true;
        for (i, leaf) in leaves.iter().enumerate() {
            let mut cs = TestConstraintSystem::<Fq>::new();
            let proof = tree.generate_proof(i, &leaf).unwrap();
            assert!(proof.verify(&crh_parameters, &root, &leaf).unwrap());

            // Allocate Merkle Tree Root
            let root = <HG as FixedLengthCRHGadget<H, _>>::OutputGadget::alloc(&mut cs.ns(|| format!("new_digest_{}", i)), || {
                if use_bad_root {
                    Ok(<H as FixedLengthCRH>::Output::default())
                } else {
                    Ok(root)
                }
            })
            .unwrap();

            let constraints_from_digest = cs.num_constraints();
            println!("constraints from digest: {}", constraints_from_digest);

            // Allocate Parameters for CRH
            let crh_parameters =
                <HG as FixedLengthCRHGadget<H, Fq>>::ParametersGadget::alloc(
                    &mut cs.ns(|| format!("new_parameters_{}", i)),
                    || Ok(crh_parameters.clone()),
                )
                .unwrap();

            let constraints_from_parameters = cs.num_constraints() - constraints_from_digest;
            println!(
                "constraints from parameters: {}",
                constraints_from_parameters
            );

            // Allocate Leaf
            let leaf_g = UInt8::constant_vec(leaf);
                
            let constraints_from_leaf =
                cs.num_constraints() - constraints_from_parameters - constraints_from_digest;
            println!("constraints from leaf: {}", constraints_from_leaf);

            // Allocate Merkle Tree Path
            let cw = MerkleTreePathGadget::<_, HG, _>::alloc(&mut cs.ns(|| format!("new_witness_{}", i)), || {
                Ok(proof)
            })
            .unwrap();

            let constraints_from_path = cs.num_constraints()
                - constraints_from_parameters
                - constraints_from_digest
                - constraints_from_leaf;
            println!("constraints from path: {}", constraints_from_path);
            let leaf_g: &[UInt8] = leaf_g.as_slice();
            cw.check_membership(
                &mut cs.ns(|| format!("new_witness_check_{}", i)),
                &crh_parameters,
                &root,
                &leaf_g,
            )
            .unwrap();
            if !cs.is_satisfied() {
                satisfied = false;
                println!(
                    "Unsatisfied constraint: {}",
                    cs.which_is_unsatisfied().unwrap()
                );
            }
            let setup_constraints = constraints_from_leaf
                + constraints_from_digest
                + constraints_from_parameters
                + constraints_from_path;
            println!(
                "number of constraints: {}",
                cs.num_constraints() - setup_constraints
            );
        }

        assert!(satisfied);
    }

    #[test]
    fn good_root_test() {
        let mut leaves = Vec::new();
        for i in 0..4u8 {
            let input = [i ; 30];
            leaves.push(input);
        }
        generate_merkle_tree(&leaves, false);
    }

    #[should_panic]
    #[test]
    fn bad_root_test() {
        let mut leaves = Vec::new();
        for i in 0..4u8 {
            let input = [i ; 30];
            leaves.push(input);
        }
        generate_merkle_tree(&leaves, true);
    }
}
