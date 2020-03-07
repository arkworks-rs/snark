use algebra_core::Field;
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::{boolean::AllocatedBit, prelude::*, uint64::UInt64};

use crate::{
    crh::{FixedLengthCRH, FixedLengthCRHGadget},
    merkle_sparse_tree::*,
};

use core::borrow::Borrow;

pub struct MerkleTreePathGadget<P, HGadget, ConstraintF>
where
    P: MerkleTreeConfig,
    HGadget: FixedLengthCRHGadget<P::H, ConstraintF>,
    ConstraintF: Field,
{
    path: Vec<(HGadget::OutputGadget, HGadget::OutputGadget)>,
}

pub struct MerkleTreeTwoPathsGadget<P, HGadget, ConstraintF>
where
    P: MerkleTreeConfig,
    HGadget: FixedLengthCRHGadget<P::H, ConstraintF>,
    ConstraintF: Field,
{
    old_path: Vec<(HGadget::OutputGadget, HGadget::OutputGadget)>,
    new_path: Vec<(HGadget::OutputGadget, HGadget::OutputGadget)>,
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

    pub fn check_membership_with_index<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        parameters: &CRHGadget::ParametersGadget,
        root: &CRHGadget::OutputGadget,
        leaf: impl ToBytesGadget<ConstraintF>,
        index: &UInt64,
    ) -> Result<(), SynthesisError> {
        self.conditionally_check_membership_with_index(
            cs,
            parameters,
            root,
            leaf,
            index,
            &Boolean::Constant(true),
        )
    }

    pub fn conditionally_check_membership_with_index<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        parameters: &CRHGadget::ParametersGadget,
        root: &CRHGadget::OutputGadget,
        leaf: impl ToBytesGadget<ConstraintF>,
        index: &UInt64,
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

        // Check levels between leaf level and root.
        let mut previous_hash = leaf_hash;
        let index_bits = index.to_bits_le();
        for (i, &(ref left_hash, ref right_hash)) in self.path.iter().enumerate() {
            // Check if the previous_hash matches the correct current hash.
            let previous_is_left = index_bits[i].not();

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

impl<P, CRHGadget, ConstraintF> MerkleTreeTwoPathsGadget<P, CRHGadget, ConstraintF>
where
    P: MerkleTreeConfig,
    ConstraintF: Field,
    CRHGadget: FixedLengthCRHGadget<P::H, ConstraintF>,
{
    pub fn check_update<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        parameters: &CRHGadget::ParametersGadget,
        old_root: &CRHGadget::OutputGadget,
        new_root: &CRHGadget::OutputGadget,
        new_leaf: impl ToBytesGadget<ConstraintF>,
        index: &UInt64,
    ) -> Result<(), SynthesisError> {
        self.conditionally_check_update(
            cs,
            parameters,
            old_root,
            new_root,
            new_leaf,
            index,
            &Boolean::Constant(true),
        )
    }

    pub fn conditionally_check_update<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        parameters: &CRHGadget::ParametersGadget,
        old_root: &CRHGadget::OutputGadget,
        new_root: &CRHGadget::OutputGadget,
        new_leaf: impl ToBytesGadget<ConstraintF>,
        index: &UInt64,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        assert_eq!(self.old_path.len(), P::HEIGHT - 1);
        assert_eq!(self.new_path.len(), P::HEIGHT - 1);
        // Check that the hash of the given leaf matches the leaf hash in the membership
        // proof.
        let new_leaf_bits = new_leaf.to_bytes(&mut cs.ns(|| "leaf_to_bytes"))?;
        let new_leaf_hash = CRHGadget::check_evaluation_gadget(
            cs.ns(|| "check_evaluation_gadget"),
            parameters,
            &new_leaf_bits,
        )?;

        // Check levels between leaf level and root of the new tree.
        let mut previous_hash = new_leaf_hash;
        let index_bits = index.to_bits_le();
        for (i, &(ref left_hash, ref right_hash)) in self.new_path.iter().enumerate() {
            // Check if the previous_hash matches the correct current hash.
            let previous_is_left = index_bits[i].not();

            CRHGadget::OutputGadget::conditional_enforce_equal_or(
                &mut cs.ns(|| format!("new_path_check_equals_which_{}", i)),
                &previous_is_left,
                &previous_hash,
                left_hash,
                right_hash,
                should_enforce,
            )?;

            previous_hash = hash_inner_node_gadget::<P::H, CRHGadget, ConstraintF, _>(
                &mut cs.ns(|| format!("new_path_hash_inner_node_{}", i)),
                parameters,
                left_hash,
                right_hash,
            )?;
        }

        new_root.conditional_enforce_equal(
            &mut cs.ns(|| "new_root_is_last"),
            &previous_hash,
            should_enforce,
        )?;

        let mut old_path_iter = self.old_path.iter();
        let old_path_first_entry = old_path_iter.next().unwrap();

        previous_hash = hash_inner_node_gadget::<P::H, CRHGadget, ConstraintF, _>(
            &mut cs.ns(|| format!("hash_leaf_above_level_inner_node")),
            parameters,
            &old_path_first_entry.0,
            &old_path_first_entry.1,
        )?;

        let mut current_loc = 1;
        loop {
            let pair = old_path_iter.next();

            match pair {
                Some((left_hash, right_hash)) => {
                    // Check if the previous_hash matches the correct current hash.
                    let previous_is_left = index_bits[current_loc].not();

                    CRHGadget::OutputGadget::conditional_enforce_equal_or(
                        &mut cs.ns(|| format!("old_path_check_equals_which_{}", current_loc)),
                        &previous_is_left,
                        &previous_hash,
                        left_hash,
                        right_hash,
                        should_enforce,
                    )?;

                    previous_hash = hash_inner_node_gadget::<P::H, CRHGadget, ConstraintF, _>(
                        &mut cs.ns(|| format!("old_path_hash_inner_node_{}", current_loc)),
                        parameters,
                        left_hash,
                        right_hash,
                    )?;
                },
                None => break,
            }
            current_loc += 1;
        }

        old_path_iter = self.old_path.iter();
        for (i, &(ref left_hash, ref right_hash)) in self.new_path.iter().enumerate() {
            // Check if the previous_hash matches the correct current hash.
            let previous_is_left = index_bits[i].not();
            let previous_is_right = previous_is_left.not();

            let old_path_corresponding_entry = old_path_iter.next().unwrap();

            right_hash.conditional_enforce_equal(
                &mut cs.ns(|| format!("check_copath_matching_case_left_{}", i)),
                &old_path_corresponding_entry.1,
                &previous_is_left,
            )?;

            left_hash.conditional_enforce_equal(
                &mut cs.ns(|| format!("check_copath_matching_case_right_{}", i)),
                &old_path_corresponding_entry.0,
                &previous_is_right,
            )?;
        }

        old_root.conditional_enforce_equal(
            &mut cs.ns(|| "old_root_is_last"),
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

        Ok(MerkleTreePathGadget { path })
    }
}

impl<P, HGadget, ConstraintF> AllocGadget<MerkleTreeTwoPaths<P>, ConstraintF>
    for MerkleTreeTwoPathsGadget<P, HGadget, ConstraintF>
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
        T: Borrow<MerkleTreeTwoPaths<P>>,
    {
        let mut old_path = Vec::new();
        let paths_borrow = value_gen()?;
        let paths = paths_borrow.borrow();
        for (i, &(ref l, ref r)) in paths.old_path.path.iter().enumerate() {
            let l_hash = HGadget::OutputGadget::alloc(
                &mut cs.ns(|| format!("old_path_l_child_{}", i)),
                || Ok(l.clone()),
            )?;
            let r_hash = HGadget::OutputGadget::alloc(
                &mut cs.ns(|| format!("old_path_r_child_{}", i)),
                || Ok(r.clone()),
            )?;
            old_path.push((l_hash, r_hash));
        }
        let mut new_path = Vec::new();
        for (i, &(ref l, ref r)) in paths.new_path.path.iter().enumerate() {
            let l_hash = HGadget::OutputGadget::alloc(
                &mut cs.ns(|| format!("new_path_l_child_{}", i)),
                || Ok(l.clone()),
            )?;
            let r_hash = HGadget::OutputGadget::alloc(
                &mut cs.ns(|| format!("new_path_r_child_{}", i)),
                || Ok(r.clone()),
            )?;
            new_path.push((l_hash, r_hash));
        }
        Ok(MerkleTreeTwoPathsGadget { old_path, new_path })
    }

    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<MerkleTreeTwoPaths<P>>,
    {
        let mut old_path = Vec::new();
        let paths_borrow = value_gen()?;
        let paths = paths_borrow.borrow();
        for (i, &(ref l, ref r)) in paths.old_path.path.iter().enumerate() {
            let l_hash = HGadget::OutputGadget::alloc_input(
                &mut cs.ns(|| format!("old_path_l_child_{}", i)),
                || Ok(l.clone()),
            )?;
            let r_hash = HGadget::OutputGadget::alloc_input(
                &mut cs.ns(|| format!("old_path_r_child_{}", i)),
                || Ok(r.clone()),
            )?;
            old_path.push((l_hash, r_hash));
        }

        let mut new_path = Vec::new();
        for (i, &(ref l, ref r)) in paths.new_path.path.iter().enumerate() {
            let l_hash = HGadget::OutputGadget::alloc_input(
                &mut cs.ns(|| format!("new_path_l_child_{}", i)),
                || Ok(l.clone()),
            )?;
            let r_hash = HGadget::OutputGadget::alloc_input(
                &mut cs.ns(|| format!("new_path_r_child_{}", i)),
                || Ok(r.clone()),
            )?;
            new_path.push((l_hash, r_hash));
        }

        Ok(MerkleTreeTwoPathsGadget { old_path, new_path })
    }
}

#[cfg(test)]
mod test {
    #[cfg(not(feature = "std"))]
    pub(crate) use alloc::collections::HashMap;

    #[cfg(feature = "std")]
    pub(crate) use std::collections::HashMap;

    use crate::{
        crh::{
            pedersen::{constraints::PedersenCRHGadget, PedersenCRH, PedersenWindow},
            FixedLengthCRH, FixedLengthCRHGadget,
        },
        merkle_sparse_tree::*,
    };
    use algebra::jubjub::{Fq, JubJubAffine as JubJub};
    use r1cs_core::ConstraintSystem;
    use rand::SeedableRng;
    use rand_xorshift::XorShiftRng;

    use super::*;
    use r1cs_std::{jubjub::JubJubGadget, test_constraint_system::TestConstraintSystem};

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

    fn generate_merkle_tree(leaves: &HashMap<usize, [u8; 30]>, use_bad_root: bool) -> () {
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);

        let crh_parameters = H::setup(&mut rng).unwrap();
        let tree = JubJubMerkleTree::new(crh_parameters.clone(), leaves).unwrap();
        let root = tree.root();
        let mut satisfied = true;
        for (i, leaf) in leaves.iter() {
            let mut cs = TestConstraintSystem::<Fq>::new();
            let proof = tree.generate_proof(*i, &leaf).unwrap();
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
            let index_g = UInt64::constant((*i) as u64);

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
            cw.check_membership_with_index(
                &mut cs.ns(|| format!("new_witness_check_with_index_{}", i)),
                &crh_parameters,
                &root,
                &leaf_g,
                &index_g,
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
    fn good_root_membership_test() {
        let mut leaves: HashMap<usize, [u8; 30]> = HashMap::new();
        for i in 0..10u8 {
            let input = [i; 30];
            leaves.insert(i as usize, input);
        }
        generate_merkle_tree(&leaves, false);
    }

    #[should_panic]
    #[test]
    fn bad_root_membership_test() {
        let mut leaves: HashMap<usize, [u8; 30]> = HashMap::new();
        for i in 0..10u8 {
            let input = [i; 30];
            leaves.insert(i as usize, input);
        }
        generate_merkle_tree(&leaves, true);
    }

    fn generate_merkle_tree_and_test_update(
        old_leaves: &HashMap<usize, [u8; 2]>,
        new_leaves: &HashMap<usize, [u8; 2]>,
    ) -> () {
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);

        let crh_parameters = H::setup(&mut rng).unwrap();
        let mut tree = JubJubMerkleTree::new(crh_parameters.clone(), old_leaves).unwrap();
        let mut satisfied = true;
        for (i, new_leaf) in new_leaves.iter() {
            let mut cs = TestConstraintSystem::<Fq>::new();

            let old_root = tree.root.unwrap();
            let update_proof = tree.update_and_prove(*i, &new_leaf).unwrap();
            let new_leaf_membership_proof = tree.generate_proof(*i, &new_leaf).unwrap();
            let new_root = tree.root.unwrap();

            assert!(update_proof
                .verify(&crh_parameters, &old_root, &new_root, &new_leaf, *i)
                .unwrap());
            assert!(new_leaf_membership_proof
                .verify_with_index(&crh_parameters, &new_root, &new_leaf, *i)
                .unwrap());

            // Allocate Merkle Tree Root
            let old_root_gadget = <HG as FixedLengthCRHGadget<H, _>>::OutputGadget::alloc(
                &mut cs.ns(|| format!("old_digest_{}", i)),
                || Ok(old_root),
            )
            .unwrap();
            let new_root_gadget = <HG as FixedLengthCRHGadget<H, _>>::OutputGadget::alloc(
                &mut cs.ns(|| format!("new_digest_{}", i)),
                || Ok(new_root),
            )
            .unwrap();

            let constraints_from_digests = cs.num_constraints();
            println!("constraints from digests: {}", constraints_from_digests);

            // Allocate Parameters for CRH
            let crh_parameters = <HG as FixedLengthCRHGadget<H, Fq>>::ParametersGadget::alloc(
                &mut cs.ns(|| format!("new_parameters_{}", i)),
                || Ok(crh_parameters.clone()),
            )
            .unwrap();

            let constraints_from_parameters = cs.num_constraints() - constraints_from_digests;
            println!(
                "constraints from parameters: {}",
                constraints_from_parameters
            );

            // Allocate Leaf
            let leaf_g = UInt8::constant_vec(new_leaf);
            let index_g = UInt64::constant((*i) as u64);

            let constraints_from_leaf =
                cs.num_constraints() - constraints_from_parameters - constraints_from_digests;
            println!("constraints from leaf: {}", constraints_from_leaf);

            // Allocate Merkle Tree Path
            let update_proof_cw = MerkleTreeTwoPathsGadget::<_, HG, _>::alloc(
                &mut cs.ns(|| format!("new_witness_update_{}", i)),
                || Ok(update_proof),
            )
            .unwrap();

            let constraints_from_two_paths = cs.num_constraints()
                - constraints_from_parameters
                - constraints_from_digests
                - constraints_from_leaf;
            println!("constraints from two paths: {}", constraints_from_two_paths);

            let new_leaf_membership_proof_cw = MerkleTreePathGadget::<_, HG, _>::alloc(
                &mut cs.ns(|| format!("new_witness_new_membership_{}", i)),
                || Ok(new_leaf_membership_proof),
            )
            .unwrap();

            let constraints_from_path = cs.num_constraints()
                - constraints_from_parameters
                - constraints_from_digests
                - constraints_from_leaf
                - constraints_from_two_paths;
            println!("constraints from path: {}", constraints_from_path);

            let leaf_g: &[UInt8] = leaf_g.as_slice();
            update_proof_cw
                .check_update(
                    &mut cs.ns(|| format!("update_witness_check_{}", i)),
                    &crh_parameters,
                    &old_root_gadget,
                    &new_root_gadget,
                    &leaf_g,
                    &index_g,
                )
                .unwrap();
            new_leaf_membership_proof_cw
                .check_membership_with_index(
                    &mut cs.ns(|| format!("new_witness_check_with_index_{}", i)),
                    &crh_parameters,
                    &new_root_gadget,
                    &leaf_g,
                    &index_g,
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
                + constraints_from_digests
                + constraints_from_parameters
                + constraints_from_two_paths
                + constraints_from_path;
            println!(
                "number of constraints: {}",
                cs.num_constraints() - setup_constraints
            );
        }

        assert!(satisfied);
    }

    #[test]
    fn good_root_update_test() {
        let mut old_leaves: HashMap<usize, [u8; 2]> = HashMap::new();
        for i in 0..4u8 {
            let input = [i; 2];
            old_leaves.insert(i as usize, input);
        }
        let mut new_leaves: HashMap<usize, [u8; 2]> = HashMap::new();
        for i in 0..8u8 {
            let input = [i + 1; 2];
            new_leaves.insert(i as usize, input);
        }
        generate_merkle_tree_and_test_update(&old_leaves, &new_leaves);
    }
}
