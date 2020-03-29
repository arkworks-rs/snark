use algebra::{Field, PrimeField};
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::prelude::*;

use primitives::{
    crh::FieldBasedHash,
    merkle_tree::field_based_mht::*,
};
use crate::crh::FieldBasedHashGadget;

use std::borrow::Borrow;
use std::marker::PhantomData;

pub struct FieldBasedMerkleTreePathGadget<P, HGadget, ConstraintF>
    where
        P: FieldBasedMerkleTreeConfig,
        P::H: FieldBasedHash<Data = ConstraintF>,
        HGadget: FieldBasedHashGadget<P::H, ConstraintF>,
        ConstraintF: Field,
{
    path: Vec<(HGadget::DataGadget, Boolean)>,
}

impl<P, HGadget, ConstraintF> FieldBasedMerkleTreePathGadget<P, HGadget, ConstraintF>
    where
        P: FieldBasedMerkleTreeConfig,
        P::H: FieldBasedHash<Data = ConstraintF>,
        HGadget: FieldBasedHashGadget<P::H, ConstraintF>,
        ConstraintF: PrimeField,
{
    pub fn check_membership<
        CS: ConstraintSystem<ConstraintF>,
    >(
        &self,
        cs: CS,
        root: &HGadget::DataGadget,
        leaf: &HGadget::DataGadget,
    ) -> Result<(), SynthesisError> {
        self.conditionally_check_membership(cs, root, leaf, &Boolean::Constant(true))
    }

    /// Coherently with the primitive, if `P::HASH_LEAVES` = `true` then we hash the
    /// leaf, otherwise we assume it to be just one FieldGadget element.
    pub fn conditionally_check_membership<
        CS: ConstraintSystem<ConstraintF>,
    >(
        &self,
        mut cs: CS,
        root: &HGadget::DataGadget,
        leaf: &HGadget::DataGadget,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {

        debug_assert!(self.path.len() == P::HEIGHT - 1);

        let mut previous_hash = (*leaf).clone();

        for (i, &(ref sibling_hash, ref direction)) in self.path.iter().enumerate() {

            //Select left hash based on direction
            let lhs = HGadget::DataGadget::conditionally_select(cs.ns(|| format!("Choose left hash {}", i)),
                                                                  direction,
                                                                  &sibling_hash,
                                                                  &previous_hash
            )?;

            //Select right hash based on direction
            let rhs = HGadget::DataGadget::conditionally_select(cs.ns(|| format!("Choose right hash {}", i)),
                                                                  direction,
                                                                  &previous_hash,
                                                                  &sibling_hash
            )?;

            previous_hash = hash_inner_node_gadget::<P::H, HGadget, ConstraintF, _>(
                &mut cs.ns(|| format!("hash_inner_node_{}", i)),
                lhs,
                rhs,
            )?;
        }

        root.conditional_enforce_equal(
            &mut cs.ns(|| "root_is_last"),
            &previous_hash,
            should_enforce,
        )
    }
}

pub struct FieldBasedMerkleTreeGadget<P, HGadget, ConstraintF>
    where
        P: FieldBasedMerkleTreeConfig,
        P::H: FieldBasedHash<Data = ConstraintF>,
        HGadget: FieldBasedHashGadget<P::H, ConstraintF>,
        ConstraintF: PrimeField,
{
    _params:        PhantomData<P>,
    _hash_gadget:   PhantomData<HGadget>,
    _field:         PhantomData<ConstraintF>,
}

impl<P, HGadget, ConstraintF> FieldBasedMerkleTreeGadget<P, HGadget, ConstraintF>
    where
        P: FieldBasedMerkleTreeConfig,
        P::H: FieldBasedHash<Data = ConstraintF>,
        HGadget: FieldBasedHashGadget<P::H, ConstraintF>,
        ConstraintF: PrimeField,
{
    pub fn check_leaves<CS: ConstraintSystem<ConstraintF>>
    (
        cs: CS,
        leaves: &[HGadget::DataGadget],
        root: &HGadget::DataGadget,
    ) -> Result<(), SynthesisError> {
        Self::conditionally_check_leaves(cs, leaves, root, &Boolean::Constant(true))
    }

    /// Starting from all the leaves in the Merkle Tree, reconstructs and enforces
    /// the Merkle Root. NOTE: This works iff Merkle Tree has been created by passing
    /// all leaves (i.e. padding_tree = null).
    pub fn conditionally_check_leaves<CS: ConstraintSystem<ConstraintF>>
    (
        mut cs: CS,
        leaves: &[HGadget::DataGadget],
        root: &HGadget::DataGadget,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        debug_assert!(leaves.len() == 2_usize.pow((P::HEIGHT - 1) as u32));

        let mut prev_level_nodes = leaves.to_vec();
        //Iterate over all levels except the root
        for level in 0..P::HEIGHT-1 {
            let mut curr_level_nodes = vec![];

            //Iterate over all nodes in a level. We assume their number to be even (e.g a power of two)

            for (i, nodes) in prev_level_nodes.chunks(2).enumerate() {
                //Compute parent hash
                let mut children = vec![];
                children.push(nodes[0].clone());
                children.push(nodes[1].clone());
                let parent_hash = HGadget::check_evaluation_gadget(
                    cs.ns(|| format!("hash_children_pair_{}_of_level_{}", i, level)),
                    children.as_slice(),
                )?;
                curr_level_nodes.push(parent_hash);
            }
            prev_level_nodes = curr_level_nodes;
        }
        //At this point, we should have only the root in prev_level_nodes
        //Enforce equality with the root
        debug_assert!(prev_level_nodes.len() == 1);

        //Enforce equality with the root

        root.conditional_enforce_equal(
            &mut cs.ns(|| "root_is_last"),
            &prev_level_nodes[0],
            should_enforce,
        )?;

        Ok(())
    }
}

pub(crate) fn hash_inner_node_gadget<H, HG, ConstraintF, CS>(
    cs: CS,
    left_child: HG::DataGadget,
    right_child: HG::DataGadget,
) -> Result<HG::DataGadget, SynthesisError>
    where
        ConstraintF: Field,
        CS: ConstraintSystem<ConstraintF>,
        H: FieldBasedHash<Data = ConstraintF>,
        HG: FieldBasedHashGadget<H, ConstraintF>,
{
    HG::check_evaluation_gadget(cs, &[left_child, right_child])
}

impl<P, HGadget, ConstraintF> AllocGadget<FieldBasedMerkleTreePath<P>, ConstraintF>
for FieldBasedMerkleTreePathGadget<P, HGadget, ConstraintF>
    where
        P: FieldBasedMerkleTreeConfig,
        P::H: FieldBasedHash<Data = ConstraintF>,
        HGadget: FieldBasedHashGadget<P::H, ConstraintF>,
        ConstraintF: Field,
{
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
        where
            F: FnOnce() -> Result<T, SynthesisError>,
            T: Borrow<FieldBasedMerkleTreePath<P>>,
    {
        let mut path = Vec::new();
        for (i, &(ref sibling, ref d)) in value_gen()?.borrow().path.iter().enumerate() {
            let sibling_hash =
                HGadget::DataGadget::alloc(&mut cs.ns(|| format!("sibling_hash_{}", i)), || {
                    Ok(sibling)
                })?;
            let direction =
                Boolean::alloc(&mut cs.ns(|| format!("direction_bit_{}", i)), || {
                    Ok(d)
                })?;
            path.push((sibling_hash, direction));
        }
        Ok(FieldBasedMerkleTreePathGadget { path })
    }

    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
        where
            F: FnOnce() -> Result<T, SynthesisError>,
            T: Borrow<FieldBasedMerkleTreePath<P>>,
    {
        let mut path = Vec::new();
        for (i, &(ref sibling, ref d)) in value_gen()?.borrow().path.iter().enumerate() {
            let sibling_hash =
                HGadget::DataGadget::alloc_input(&mut cs.ns(|| format!("sibling_hash_{}", i)), || {
                    Ok(sibling)
                })?;
            let direction =
                Boolean::alloc_input(&mut cs.ns(|| format!("direction_bit_{}", i)), || {
                    Ok(d)
                })?;
            path.push((sibling_hash, direction));
        }
        Ok(FieldBasedMerkleTreePathGadget { path })
    }
}

#[cfg(test)]
mod test {
    use primitives::{
        crh::MNT4PoseidonHash,
        merkle_tree::field_based_mht::*,
    };
    use crate::crh::MNT4PoseidonHashGadget;
    use algebra::fields::mnt4753::Fr;
    use r1cs_core::ConstraintSystem;
    use rand::{Rng, SeedableRng};
    use rand_xorshift::XorShiftRng;
    use super::*;
    use r1cs_std::{
        fields::mnt6753::FqGadget,
        test_constraint_system::TestConstraintSystem,
    };

    struct MNT4753FieldBasedMerkleTreeParams;

    impl FieldBasedMerkleTreeConfig for MNT4753FieldBasedMerkleTreeParams {
        const HEIGHT: usize = 6;
        type H = MNT4PoseidonHash;
    }

    type MNT4753FieldBasedMerkleTree = FieldBasedMerkleHashTree<MNT4753FieldBasedMerkleTreeParams>;

    type HG = MNT4PoseidonHashGadget;

    fn check_merkle_paths(leaves: &[Fr], use_bad_root: bool) -> bool {

        let tree = MNT4753FieldBasedMerkleTree::new(leaves).unwrap();
        let root = tree.root();
        let mut satisfied = true;

        //Merkle Path Gadget test
        for (i, leaf) in leaves.iter().enumerate() {
            let mut cs = TestConstraintSystem::<Fr>::new();
            let proof = tree.generate_proof(i, leaf).unwrap();
            assert!(proof.verify(&root, leaf).unwrap());

            // Allocate Merkle Tree Root
            let root = FqGadget::alloc(
                &mut cs.ns(|| format!("new_digest_{}", i)),
                || {
                    if use_bad_root {
                        Ok(Fr::zero())
                    } else {
                        Ok(root)
                    }
                },
            )
                .unwrap();

            // Allocate Leaf
            let leaf_g = FqGadget::alloc(cs.ns(|| "alloc leaf"), || Ok(leaf)).unwrap();

            // Allocate Merkle Tree Path
            let cw = FieldBasedMerkleTreePathGadget::<_, HG, _>::alloc(
                &mut cs.ns(|| format!("new_witness_{}", i)),
                || Ok(proof),
            )
                .unwrap();

            cw.check_membership(
                &mut cs.ns(|| format!("new_witness_check_{}", i)),
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
        }

        satisfied
    }

    fn check_leaves(leaves: &[Fr], use_bad_root: bool) -> bool {

        let tree = MNT4753FieldBasedMerkleTree::new(leaves).unwrap();
        let root = tree.root();

        //Merkle Tree Gadget test
        let mut cs = TestConstraintSystem::<Fr>::new();

        // Allocate Merkle Tree Root
        let root = FqGadget::alloc(
            &mut cs.ns(|| "root_digest_{}"),
            || {
                if use_bad_root {
                    Ok(Fr::zero())
                } else {
                    Ok(root)
                }
            },
        )
            .unwrap();

        //Alloc leaves
        let mut leaves_g = vec![];
        for (i, leaf) in leaves.iter().enumerate() {
            leaves_g.push(FqGadget::alloc(cs.ns(|| format!("alloc leaf_{}", i)), || Ok(leaf)).unwrap());
        }

        //Check MR from leaves
        FieldBasedMerkleTreeGadget::<MNT4753FieldBasedMerkleTreeParams, HG, Fr>::check_leaves(
            &mut cs.ns(|| "check all leaves belong to MT"),
            &leaves_g,
            &root,
        ).unwrap();

        if !cs.is_satisfied() {
            println!(
                "Unsatisfied constraint: {}",
                cs.which_is_unsatisfied().unwrap()
            );
        }

        cs.is_satisfied()
    }

    #[test]
    fn good_root_test() {
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);

        //Test #leaves << 2^HEIGHT
        let mut leaves = Vec::new();
        for _ in 0..4 {
            let f: Fr = rng.gen();
            leaves.push(f);
        }
        assert!(check_merkle_paths(&leaves, false));

        //Test #leaves = 2^HEIGHT - 1
        let mut leaves = Vec::new();
        for _ in 0..16 {
            let f: Fr = rng.gen();
            leaves.push(f);
        }
        assert!(check_merkle_paths(&leaves, false));

        //Test #leaves == 2^HEIGHT
        let mut leaves = Vec::new();
        for _ in 0..32 {
            let f: Fr = rng.gen();
            leaves.push(f);
        }
        assert!(check_merkle_paths(&leaves, false));
        assert!(check_leaves(&leaves, false));
    }

    #[test]
    fn bad_root_test() {
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);

        //Test #leaves << 2^HEIGHT
        let mut leaves = Vec::new();
        for _ in 0..4 {
            let f: Fr = rng.gen();
            leaves.push(f);
        }
        assert!(!check_merkle_paths(&leaves, true));

        //Test #leaves = 2^HEIGHT - 1
        let mut leaves = Vec::new();
        for _ in 0..16 {
            let f: Fr = rng.gen();
            leaves.push(f);
        }
        assert!(!check_merkle_paths(&leaves, true));

        //Test #leaves == 2^HEIGHT
        let mut leaves = Vec::new();
        for _ in 0..32 {
            let f: Fr = rng.gen();
            leaves.push(f);
        }
        assert!(!check_merkle_paths(&leaves, true));
        assert!(!check_leaves(&leaves, true));
    }
}