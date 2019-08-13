use algebra::Field;
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::prelude::*;
use r1cs_std::boolean::AllocatedBit;

use crate::crypto_primitives::{mht::HashMembershipProof, CommitmentScheme, FixedLengthCRH};

use crate::gadgets::{commitment::CommitmentGadget, crh::FixedLengthCRHGadget};

use crate::ledger::{CommPath, Digest, LedgerDigest, LedgerWitness};

use std::{borrow::Borrow, marker::PhantomData};

pub trait LCWGadget<C: CommitmentScheme, D: LedgerDigest, CW: LedgerWitness<D>, ConstraintF: Field> {
    type ParametersGadget: AllocGadget<D::Parameters, ConstraintF>;
    type CommitmentGadget: AllocGadget<C::Output, ConstraintF>;
    type DigestGadget: AllocGadget<D, ConstraintF>;
    type WitnessGadget: AllocGadget<CW, ConstraintF>;

    fn check_witness_gadget<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        parameters: &Self::ParametersGadget,
        ledger_state_digest: &Self::DigestGadget,
        commitment: &Self::CommitmentGadget,
        witness: &Self::WitnessGadget,
    ) -> Result<(), SynthesisError> {
        Self::conditionally_check_witness_gadget(
            cs,
            parameters,
            ledger_state_digest,
            commitment,
            witness,
            &Boolean::Constant(true),
        )
    }

    fn conditionally_check_witness_gadget<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        parameters: &Self::ParametersGadget,
        ledger_state_digest: &Self::DigestGadget,
        commitment: &Self::CommitmentGadget,
        witness: &Self::WitnessGadget,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError>;
}

pub struct IdealLedgerGadget<C, H, HGadget, CGadget> {
    #[doc(hidden)]
    _comm_scheme: PhantomData<C>,
    #[doc(hidden)]
    _hash:        PhantomData<H>,
    #[doc(hidden)]
    _comm_gadget: PhantomData<CGadget>,
    #[doc(hidden)]
    _hash_gadget: PhantomData<HGadget>,
}

pub struct CommitmentWitness<
    H: FixedLengthCRH,
    C: CommitmentScheme,
    HGadget: FixedLengthCRHGadget<H, ConstraintF>,
    ConstraintF: Field,
> {
    path:    Vec<(HGadget::OutputGadget, HGadget::OutputGadget)>,
    _crh:    PhantomData<H>,
    _comm:   PhantomData<C>,
    _engine: PhantomData<ConstraintF>,
}

pub struct DigestGadget<H: FixedLengthCRH, HGadget: FixedLengthCRHGadget<H, ConstraintF>, ConstraintF: Field> {
    digest:  HGadget::OutputGadget,
    #[doc(hidden)]
    _crh:    PhantomData<H>,
    #[doc(hidden)]
    _engine: PhantomData<ConstraintF>,
}

impl<C, ConstraintF, H, CGadget, HGadget> LCWGadget<C, Digest<H>, CommPath<H, C::Output>, ConstraintF>
    for IdealLedgerGadget<C, H, HGadget, CGadget>
where
    C: CommitmentScheme,
    C::Output: Eq,
    ConstraintF: Field,
    H: FixedLengthCRH,
    CGadget: CommitmentGadget<C, ConstraintF>,
    HGadget: FixedLengthCRHGadget<H, ConstraintF>,
{
    type ParametersGadget = <HGadget as FixedLengthCRHGadget<H, ConstraintF>>::ParametersGadget;
    type DigestGadget = DigestGadget<H, HGadget, ConstraintF>;

    type CommitmentGadget = <CGadget as CommitmentGadget<C, ConstraintF>>::OutputGadget;
    type WitnessGadget = CommitmentWitness<H, C, HGadget, ConstraintF>;

    /// Given a `leaf` and `path`, check that the `path` is a valid
    /// authentication path for the `leaf` in a Merkle tree.
    /// Note: It is assumed that the root is contained in the `path`.
    fn conditionally_check_witness_gadget<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        parameters: &Self::ParametersGadget,
        root_hash: &Self::DigestGadget,
        commitment: &Self::CommitmentGadget,
        witness: &Self::WitnessGadget,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        assert_eq!(
            witness.path.len(),
            (HashMembershipProof::<H, C::Output>::MAX_HEIGHT - 1) as usize
        );
        // Check that the hash of the given leaf matches the leaf hash in the membership
        // proof.
        let commitment_bits = commitment.to_bytes(&mut cs.ns(|| "commitment_to_bytes"))?;
        let commitment_hash = HGadget::check_evaluation_gadget(
            cs.ns(|| "check_evaluation_gadget"),
            parameters,
            &commitment_bits,
        )?;

        // Check if leaf is one of the bottom-most siblings.
        let leaf_is_left = AllocatedBit::alloc(&mut cs.ns(|| "leaf_is_left"), || {
            Ok(commitment_hash == witness.path[0].0)
        })?
        .into();
        <HGadget::OutputGadget>::conditional_enforce_equal_or(
            &mut cs.ns(|| "check_leaf_is_left"),
            &leaf_is_left,
            &commitment_hash,
            &witness.path[0].0,
            &witness.path[0].1,
            should_enforce,
        )?;

        // Check levels between leaf level and root.
        let mut previous_hash = commitment_hash;
        for (i, &(ref left_hash, ref right_hash)) in witness.path.iter().enumerate() {
            // Check if the previous_hash matches the correct current hash.
            let previous_is_left =
                AllocatedBit::alloc(&mut cs.ns(|| format!("previous_is_left_{}", i)), || {
                    Ok(&previous_hash == left_hash)
                })?
                .into();

            <HGadget::OutputGadget>::conditional_enforce_equal_or(
                &mut cs.ns(|| format!("check_equals_which_{}", i)),
                &previous_is_left,
                &previous_hash,
                left_hash,
                right_hash,
                should_enforce,
            )?;

            previous_hash = hash_inner_node_gadget::<H, HGadget, ConstraintF, _>(
                &mut cs.ns(|| format!("hash_inner_node_{}", i)),
                parameters,
                left_hash,
                right_hash,
            )?;
        }

        root_hash.digest.conditional_enforce_equal(
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

impl<H, HGadget, ConstraintF> AllocGadget<Digest<H>, ConstraintF> for DigestGadget<H, HGadget, ConstraintF>
where
    H: FixedLengthCRH,
    HGadget: FixedLengthCRHGadget<H, ConstraintF>,
    ConstraintF: Field,
{
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<Digest<H>>,
    {
        let digest = HGadget::OutputGadget::alloc(&mut cs.ns(|| "digest"), || {
            Ok(value_gen()?.borrow().0.clone())
        })?;

        Ok(DigestGadget {
            digest,
            _crh: PhantomData,
            _engine: PhantomData,
        })
    }

    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<Digest<H>>,
    {
        let digest = HGadget::OutputGadget::alloc_input(&mut cs.ns(|| "input_digest"), || {
            Ok(value_gen()?.borrow().0.clone())
        })?;
        Ok(DigestGadget {
            digest,
            _crh: PhantomData,
            _engine: PhantomData,
        })
    }
}

impl<H, C, HGadget, ConstraintF> AllocGadget<CommPath<H, C::Output>, ConstraintF>
    for CommitmentWitness<H, C, HGadget, ConstraintF>
where
    H: FixedLengthCRH,
    C: CommitmentScheme,
    HGadget: FixedLengthCRHGadget<H, ConstraintF>,
    ConstraintF: Field,
{
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<CommPath<H, C::Output>>,
    {
        let mut path = Vec::new();
        for (i, &(ref l, ref r)) in value_gen()?.borrow().0.path.iter().enumerate() {
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
        Ok(CommitmentWitness {
            path,
            _crh: PhantomData,
            _comm: PhantomData,
            _engine: PhantomData,
        })
    }

    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<CommPath<H, C::Output>>,
    {
        let mut path = Vec::new();
        for (i, &(ref l, ref r)) in value_gen()?.borrow().0.path.iter().enumerate() {
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

        Ok(CommitmentWitness {
            path,
            _crh: PhantomData,
            _comm: PhantomData,
            _engine: PhantomData,
        })
    }
}

#[cfg(test)]
mod test {
    use std::rc::Rc;

    use crate::crypto_primitives::{
        commitment::{
            pedersen::{PedersenCommitment, PedersenRandomness},
            CommitmentScheme,
        },
        crh::{
            pedersen::{PedersenCRH, PedersenWindow},
            FixedLengthCRH,
        },
        mht::*,
    };
    use algebra::{
        curves::jubjub::JubJubAffine as JubJub,
        fields::jubjub::fr::Fr,
        fields::jubjub::fq::Fq,
        Group
    };
    use rand::{ChaChaRng, Rand, SeedableRng};
    use r1cs_core::ConstraintSystem;

    use super::*;
    use crate::gadgets::{
        commitment::pedersen::PedersenCommitmentGadget,
        crh::{pedersen::PedersenCRHGadget, FixedLengthCRHGadget},
    };
    use r1cs_std::{
        groups::curves::twisted_edwards::jubjub::JubJubGadget,
        test_constraint_system::TestConstraintSystem,
    };

    use crate::ledger::{CommPath, Digest};

    #[derive(Clone)]
    pub(super) struct Window4x256;
    impl PedersenWindow for Window4x256 {
        const WINDOW_SIZE: usize = 4;
        const NUM_WINDOWS: usize = 256;
    }

    type H = PedersenCRH<JubJub, Window4x256>;
    type HG = PedersenCRHGadget<JubJub, Fq, JubJubGadget>;
    type C = PedersenCommitment<JubJub, Window4x256>;
    type CG = PedersenCommitmentGadget<JubJub, Fq, JubJubGadget>;
    type JubJubMHT = MerkleHashTree<H, <C as CommitmentScheme>::Output>;
    type LG = IdealLedgerGadget<C, H, HG, CG>;
    type DG = DigestGadget<H, HG, Fq>;
    type LCWG = CommitmentWitness<H, C, HG, Fq>;

    fn generate_merkle_tree(leaves: &[<C as CommitmentScheme>::Output]) -> () {
        let seed: [u32; 8] = [
            2053759276, 152413135, 1690980041, 4293109333, 2390175708, 686052238, 1844363894,
            1379683288,
        ];
        let mut rng = ChaChaRng::from_seed(&seed);

        let crh_parameters = Rc::new(H::setup(&mut rng).unwrap());
        let tree = JubJubMHT::new(crh_parameters.clone(), &leaves).unwrap();
        let root = tree.root();
        let mut satisfied = true;
        for (i, leaf) in leaves.iter().enumerate() {
            let mut cs = TestConstraintSystem::<Fq>::new();
            let proof = tree.generate_proof(i, &leaf).unwrap();
            assert!(proof.verify(&crh_parameters, &root, &leaf).unwrap());

            let digest = DG::alloc(&mut cs.ns(|| format!("new_digest_{}", i)), || {
                Ok(Digest(root))
            })
            .unwrap();
            let constraints_from_digest = cs.num_constraints();
            println!("constraints from digest: {}", constraints_from_digest);

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

            let comm = <CG as CommitmentGadget<C, Fq>>::OutputGadget::alloc(
                &mut cs.ns(|| format!("new_comm_{}", i)),
                || {
                    let leaf: JubJub = *leaf;
                    Ok(leaf)
                },
            )
            .unwrap();
            let constraints_from_comm =
                cs.num_constraints() - constraints_from_parameters - constraints_from_digest;
            println!("constraints from comm: {}", constraints_from_comm);

            let cw = LCWG::alloc(&mut cs.ns(|| format!("new_witness_{}", i)), || {
                Ok(CommPath(proof))
            })
            .unwrap();
            let constraints_from_path = cs.num_constraints()
                - constraints_from_parameters
                - constraints_from_digest
                - constraints_from_comm;
            println!("constraints from path: {}", constraints_from_path);
            LG::check_witness_gadget(
                &mut cs.ns(|| format!("new_witness_check_{}", i)),
                &crh_parameters,
                &digest,
                &comm,
                &cw,
            )
            .unwrap();
            if !cs.is_satisfied() {
                satisfied = false;
                println!(
                    "Unsatisfied constraint: {}",
                    cs.which_is_unsatisfied().unwrap()
                );
            }
            let setup_constraints = constraints_from_comm
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
    fn mht_gadget_test() {
        let mut leaves = Vec::new();
        let seed: [u32; 8] = [
            2053759276, 152413135, 1690980041, 4293109333, 2390175708, 686052238, 1844363894,
            1379683288,
        ];
        let mut rng = ChaChaRng::from_seed(&seed);
        let comm_parameters = C::setup(&mut rng).unwrap();
        for i in 0..4u8 {
            let r = PedersenRandomness(Fr::rand(&mut rng));
            let input = [i, i, i, i, i, i, i, i];
            leaves.push(C::commit(&comm_parameters, &input, &r).unwrap());
        }
        generate_merkle_tree(&leaves);
        // let mut leaves = Vec::new();
        // for i in 0..100u8 {
        //     leaves.push([i, i, i, i, i, i, i, i]);
        // }
        // generate_merkle_tree(&leaves);
    }

    fn bad_merkle_tree_verify(leaves: &[<C as CommitmentScheme>::Output]) -> () {
        let seed: [u32; 8] = [
            2053759276, 152413135, 1690980041, 4293109333, 2390175708, 686052238, 1844363894,
            1379683288,
        ];
        let mut rng = ChaChaRng::from_seed(&seed);

        let crh_parameters = Rc::new(H::setup(&mut rng).unwrap());
        let tree = JubJubMHT::new(crh_parameters.clone(), &leaves).unwrap();
        let root = tree.root();
        for (i, leaf) in leaves.iter().enumerate() {
            let mut cs = TestConstraintSystem::<Fq>::new();
            let proof = tree.generate_proof(i, &leaf).unwrap();
            assert!(proof.verify(&crh_parameters, &root, &leaf).unwrap());

            let digest = DG::alloc(&mut cs.ns(|| format!("new_digest_{}", i)), || {
                Ok(Digest(JubJub::zero()))
            })
            .unwrap();
            let crh_parameters =
                <HG as FixedLengthCRHGadget<H, Fq>>::ParametersGadget::alloc(
                    &mut cs.ns(|| format!("new_parameters_{}", i)),
                    || Ok(crh_parameters.clone()),
                )
                .unwrap();
            let comm = <CG as CommitmentGadget<C, Fq>>::OutputGadget::alloc(
                &mut cs.ns(|| format!("new_comm_{}", i)),
                || {
                    let leaf = *leaf;
                    Ok(leaf)
                },
            )
            .unwrap();
            let cw = LCWG::alloc(&mut cs.ns(|| format!("new_witness_{}", i)), || {
                Ok(CommPath(proof))
            })
            .unwrap();
            LG::check_witness_gadget(
                &mut cs.ns(|| format!("new_witness_check_{}", i)),
                &crh_parameters,
                &digest,
                &comm,
                &cw,
            )
            .unwrap();
            if !cs.is_satisfied() {
                println!(
                    "Unsatisfied constraints: {}",
                    cs.which_is_unsatisfied().unwrap()
                );
            }
            assert!(cs.is_satisfied());
        }
    }

    #[should_panic]
    #[test]
    fn bad_root_test() {
        let mut leaves = Vec::new();
        let seed: [u32; 8] = [
            2053759276, 152413135, 1690980041, 4293109333, 2390175708, 686052238, 1844363894,
            1379683288,
        ];
        let mut rng = ChaChaRng::from_seed(&seed);
        let comm_parameters = C::setup(&mut rng).unwrap();
        for i in 0..4u8 {
            let r = PedersenRandomness(Fr::rand(&mut rng));
            let input = [i, i, i, i, i, i, i, i];
            leaves.push(C::commit(&comm_parameters, &input, &r).unwrap());
        }
        bad_merkle_tree_verify(&leaves);
    }
}
