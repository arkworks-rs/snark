use algebra::{Field, PrimeField};
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::prelude::*;

use primitives::{crh::FixedLengthCRH, merkle_tree::*, FieldBasedHash};
use crate::{FixedLengthCRHGadget, FieldBasedHashGadget};

use std::borrow::Borrow;
use r1cs_std::fields::fp::FpGadget;

pub mod field_based_mht;

pub trait FieldBasedMerkleTreePathGadget<
    P: FieldBasedMerkleTreePath<H = H>,
    H:  FieldBasedHash<Data = ConstraintF>,
    HGadget: FieldBasedHashGadget<H, ConstraintF>,
    ConstraintF: PrimeField,
>: AllocGadget<P, ConstraintF> + ConstantGadget<P, ConstraintF> + EqGadget<ConstraintF> + Clone
where
{
    /// Return the length of the `self` path.
    fn length(&self) -> usize;

    /// Enforce that the root reconstructed from `self` and `leaf` is equal to
    /// `expected_root`.
    fn check_membership<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        expected_root: &HGadget::DataGadget,
        leaf: &HGadget::DataGadget,
    ) -> Result<(), SynthesisError> {
        self.conditionally_check_membership(cs, expected_root, leaf, &Boolean::Constant(true))
    }

    /// Enforce that the root reconstructed from `self` and `leaf` is equal to
    /// `expected_root` if `should_enforce` is True, otherwise enforce nothing.
    fn conditionally_check_membership<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        expected_root: &HGadget::DataGadget,
        leaf: &HGadget::DataGadget,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError>
    {
        let root = self.enforce_root_from_leaf(
            cs.ns(|| "reconstruct root"),
            leaf
        )?;

        root.conditional_enforce_equal(
            &mut cs.ns(|| "root_is_last"),
            expected_root,
            should_enforce,
        )
    }

    /// Enforce correct reconstruction of the root of the Merkle Tree
    /// from `self` path and `leaf`.
    fn enforce_root_from_leaf<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        leaf: &HGadget::DataGadget,
    ) -> Result<HGadget::DataGadget, SynthesisError>;

    /// Given a field element `leaf_index` representing the position of a leaf in a
    /// Merkle Tree, enforce that the leaf index corresponding to `self` path is the
    /// same of `leaf_index`.
    fn enforce_leaf_index<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        leaf_index: &FpGadget<ConstraintF>,
    ) -> Result<(), SynthesisError>
    {
        self.conditionally_enforce_leaf_index(cs, leaf_index, &Boolean::Constant(true))
    }


    /// Given a field element `leaf_index` representing the position of a leaf in a
    /// Merkle Tree, enforce that the leaf index corresponding to `self` path is the
    /// same of `leaf_index` if `should_enforce` is True, otherwise enforce nothing.
    fn conditionally_enforce_leaf_index<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        leaf_index: &FpGadget<ConstraintF>,
        should_enforce: &Boolean
    ) -> Result<(), SynthesisError>;
}

pub struct MerkleTreePathGadget<P, HGadget, ConstraintF>
where
    P: MerkleTreeConfig,
    HGadget: FixedLengthCRHGadget<P::H, ConstraintF>,
    ConstraintF: Field,
{
    path: Vec<(HGadget::OutputGadget, Boolean)>,
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
        self.conditionally_check_membership(cs, parameters, root, leaf, &Boolean::Constant(true))
    }

    pub fn conditionally_check_membership<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        parameters: &CRHGadget::ParametersGadget,
        root: &CRHGadget::OutputGadget,
        leaf: impl ToBytesGadget<ConstraintF>,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        if self.path.len() != P::HEIGHT {
            return Err(SynthesisError::Other(format!(
                "Path length must be equal to height. Path len: {}, Height: {}",
                self.path.len(),
                P::HEIGHT
            ).to_owned()));
        }

        // Check that the hash of the given leaf matches the leaf hash in the membership
        // proof.
        let leaf_bits = leaf.to_bytes(&mut cs.ns(|| "leaf_to_bytes"))?;
        let leaf_hash = CRHGadget::check_evaluation_gadget(
            cs.ns(|| "check_evaluation_gadget"),
            parameters,
            &leaf_bits,
        )?;

        // Check levels between leaf level and root.
        let mut previous_hash = leaf_hash;
        for (i, &(ref sibling_hash, ref direction)) in self.path.iter().enumerate() {

            //Select left hash based on direction
            let lhs = CRHGadget::OutputGadget::conditionally_select(cs.ns(|| format!("Choose left hash {}", i)),
                                                                direction,
                                                                &sibling_hash,
                                                                &previous_hash
            )?;

            //Select right hash based on direction
            let rhs = CRHGadget::OutputGadget::conditionally_select(cs.ns(|| format!("Choose right hash {}", i)),
                                                                direction,
                                                                &previous_hash,
                                                                &sibling_hash
            )?;

            previous_hash = hash_inner_node_gadget::<P::H, CRHGadget, ConstraintF, _>(
                &mut cs.ns(|| format!("hash_inner_node_{}", i)),
                parameters,
                &lhs,
                &rhs,
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
        for (i, &(ref sibling, ref d)) in value_gen()?.borrow().path.iter().enumerate() {
            let sibling_hash =
                HGadget::OutputGadget::alloc(&mut cs.ns(|| format!("sibling_hash_{}", i)), || {
                    Ok(sibling)
                })?;
            let direction =
                Boolean::alloc(&mut cs.ns(|| format!("direction_bit_{}", i)), || {
                    Ok(d)
                })?;
            path.push((sibling_hash, direction));
        }
        Ok(MerkleTreePathGadget { path })
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
        for (i, &(ref sibling, ref d)) in value_gen()?.borrow().path.iter().enumerate() {
            let sibling_hash =
                HGadget::OutputGadget::alloc_input(&mut cs.ns(|| format!("sibling_hash_{}", i)), || {
                    Ok(sibling)
                })?;
            let direction =
                Boolean::alloc_input(&mut cs.ns(|| format!("direction_bit_{}", i)), || {
                    Ok(d)
                })?;
            path.push((sibling_hash, direction));
        }
        Ok(MerkleTreePathGadget { path })
    }
}

#[cfg(test)]
mod test {
    use std::rc::Rc;

    use primitives::{
        crh::{
            FixedLengthCRH,
            pedersen::PedersenWindow,
            injective_map::{PedersenCRHCompressor, TECompressor},
        },
        merkle_tree::*,
    };
    use crate::crh::{
        FixedLengthCRHGadget,
        injective_map::{PedersenCRHCompressorGadget, TECompressorGadget},
    };
    use algebra::{curves::jubjub::JubJubAffine as JubJub, fields::jubjub::fq::Fq};
    use r1cs_core::ConstraintSystem;
    use rand::SeedableRng;
    use rand_xorshift::XorShiftRng;
    use super::*;
    use r1cs_std::{
        instantiated::jubjub::JubJubGadget,
        test_constraint_system::TestConstraintSystem,
    };

    #[derive(Clone)]
    pub(super) struct Window4x128;
    impl PedersenWindow for Window4x128 {
        const WINDOW_SIZE: usize = 4;
        const NUM_WINDOWS: usize = 128;
    }

    type H = PedersenCRHCompressor<JubJub, TECompressor, Window4x128>;
    type HG = PedersenCRHCompressorGadget<JubJub, TECompressor, Fq, JubJubGadget, TECompressorGadget>;

    struct JubJubMerkleTreeParams;

    impl MerkleTreeConfig for JubJubMerkleTreeParams {
        const HEIGHT: usize = 3;
        type H = H;
    }

    type JubJubMerkleTree = MerkleHashTree<JubJubMerkleTreeParams>;

    fn generate_merkle_tree(leaves: &[[u8; 8]], use_bad_root: bool) -> bool {
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);

        let crh_parameters = Rc::new(H::setup(&mut rng).unwrap());
        let tree = JubJubMerkleTree::new(crh_parameters.clone(), leaves).unwrap();
        let root = tree.root().unwrap();
        let mut satisfied = true;
        for (i, leaf) in leaves.iter().enumerate() {
            let mut cs = TestConstraintSystem::<Fq>::new();
            let proof = tree.generate_proof(i, &leaf).unwrap();
            assert!(proof.verify(&crh_parameters, &root, &leaf).unwrap());

            // Allocate Merkle Tree Root
            let root = <HG as FixedLengthCRHGadget<H, _>>::OutputGadget::alloc(
                &mut cs.ns(|| format!("new_digest_{}", i)),
                || {
                    if use_bad_root {
                        Ok(<H as FixedLengthCRH>::Output::default())
                    } else {
                        Ok(root)
                    }
                },
            )
                .unwrap();

            let constraints_from_digest = cs.num_constraints();
            println!("constraints from digest: {}", constraints_from_digest);

            // Allocate Parameters for CRH
            let crh_parameters = <HG as FixedLengthCRHGadget<H, Fq>>::ParametersGadget::alloc(
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
            let cw = MerkleTreePathGadget::<_, HG, _>::alloc(
                &mut cs.ns(|| format!("new_witness_{}", i)),
                || Ok(proof),
            )
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

        satisfied
    }

    #[test]
    fn good_root_test() {

        //Test #leaves << 2^HEIGHT
        let mut leaves = Vec::new();
        for i in 0..2u8 {
            let input = [i; 8];
            leaves.push(input);
        }
        assert!(generate_merkle_tree(&leaves, false));

        //Test #leaves = 2^HEIGHT - 1
        let mut leaves = Vec::new();
        for i in 0..4u8 {
            let input = [i; 8];
            leaves.push(input);
        }
        assert!(generate_merkle_tree(&leaves, false));

        //Test #leaves = 2^HEIGHT
        let mut leaves = Vec::new();
        for i in 0..8u8 {
            let input = [i; 8];
            leaves.push(input);
        }
        assert!(generate_merkle_tree(&leaves, false));
    }

    #[test]
    fn bad_root_test() {

        //Test #leaves << 2^HEIGHT
        let mut leaves = Vec::new();
        for i in 0..2u8 {
            let input = [i; 8];
            leaves.push(input);
        }
        assert!(!generate_merkle_tree(&leaves, true));

        //Test #leaves = 2^HEIGHT - 1
        let mut leaves = Vec::new();
        for i in 0..4u8 {
            let input = [i; 8];
            leaves.push(input);
        }
        assert!(!generate_merkle_tree(&leaves, true));

        //Test #leaves = 2^HEIGHT
        let mut leaves = Vec::new();
        for i in 0..8u8 {
            let input = [i; 8];
            leaves.push(input);
        }
        assert!(!generate_merkle_tree(&leaves, true));
    }
}