use algebra::{Field, PrimeField};
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::prelude::*;

use crate::{
    crh::{FieldBasedHash, FieldBasedHashGadget},
    merkle_tree::field_based_mht::*,
};

use std::borrow::Borrow;
use r1cs_std::to_field_gadget_vec::ToConstraintFieldGadget;

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
        leaf: impl ToConstraintFieldGadget<ConstraintF, FieldGadget = HGadget::DataGadget>,
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
        leaf: impl ToConstraintFieldGadget<ConstraintF, FieldGadget = HGadget::DataGadget>,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {

        assert_eq!(self.path.len(), P::HEIGHT - 1);

        let mut previous_hash =
            if P::HASH_LEAVES {
                let leaf_elements = leaf.to_field_gadget_elements()?;
                HGadget::check_evaluation_gadget(
                    cs.ns(|| "hash leaf"),
                    leaf_elements.as_slice(),
                )?
            } else {
                leaf.to_field_gadget_elements().map(|fe_g| fe_g[0].clone())?
            };

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

    use crate::{
        crh::{
            MNT4PoseidonHash, MNT4PoseidonHashGadget
        },
        merkle_tree::field_based_mht::*,
    };
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
        const HASH_LEAVES: bool = true;
        const HEIGHT: usize = 32;
        type H = MNT4PoseidonHash;
    }

    type MNT4753FieldBasedMerkleTree = FieldBasedMerkleHashTree<MNT4753FieldBasedMerkleTreeParams>;

    struct MNT4753FieldBasedMerkleTreeParamsNoHash;

    impl FieldBasedMerkleTreeConfig for MNT4753FieldBasedMerkleTreeParamsNoHash {
        const HASH_LEAVES: bool = false;
        const HEIGHT: usize = 32;
        type H = MNT4PoseidonHash;
    }

    type MNT4753FieldBasedMerkleTreeNoHash = FieldBasedMerkleHashTree<MNT4753FieldBasedMerkleTreeParamsNoHash>;

    type HG = MNT4PoseidonHashGadget;

    fn generate_merkle_tree(leaves: &[Fr], use_bad_root: bool) -> () {

        let tree = MNT4753FieldBasedMerkleTree::new(leaves).unwrap();
        let root = tree.root();
        let mut satisfied = true;
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
                leaf_g,
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
        assert!(satisfied);
    }

    fn generate_merkle_tree_no_hash(leaves: &[Fr], use_bad_root: bool) -> () {

        let tree = MNT4753FieldBasedMerkleTreeNoHash::new(leaves).unwrap();
        let root = tree.root();
        let mut satisfied = true;
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
                leaf_g,
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
        assert!(satisfied);
    }

    #[test]
    fn good_root_test() {
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);
        let mut leaves = Vec::new();
        for _ in 0..4 {
            let f: Fr = rng.gen();
            leaves.push(f);
        }
        generate_merkle_tree(&leaves, false);
        generate_merkle_tree_no_hash(&leaves, false);
    }

    #[should_panic]
    #[test]
    fn bad_root_test() {
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);
        let mut leaves = Vec::new();
        for _ in 0..4 {
            let f: Fr = rng.gen();
            leaves.push(f);
        }
        generate_merkle_tree(&leaves, true);
        generate_merkle_tree_no_hash(&leaves, false);
    }
}