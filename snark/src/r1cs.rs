use ark_std::rand::{RngCore, CryptoRng};
use ark_ff::Field;
use ark_relations::r1cs::{R1CS, ConstraintMatrices, Instance, Witness, ConstraintSynthesizer};
use crate::SNARK;

/// A [`SNARKForR1CS`] is a [`SNARK`] for the [`R1CS`] relation.
pub trait SNARKForR1CS<F: Field>: SNARK<R1CS<F>> {
    /// Generate inputs for the SNARK indexer from [`cs`].
    fn indexer_inputs<CS: ConstraintSynthesizer<F>>(cs: &CS) -> ConstraintMatrices<F>;

    /// Generate inputs for the SNARK prover from [`cs`].
    fn prover_inputs<CS: ConstraintSynthesizer<F>>(
        cs: &CS,
    ) -> (
        Instance<F>,
        Witness<F>,
    );

    /// Generate inputs for the SNARK verifier from [`cs`].
    fn verifier_inputs<CS: ConstraintSynthesizer<F>>(cs: &CS) -> Instance<F>;

    /// Generates a proof of satisfaction of the constraint system induced by [`cs`].
    fn prove_with_cs<CS: ConstraintSynthesizer<F>, Rng: RngCore + CryptoRng>(
        pk: &Self::ProvingKey,
        c: CS,
        rng: &mut Rng,
    ) -> Result<Self::Proof, Self::Error> {
        let (instance, witness) = Self::prover_inputs(&c);
        Self::prove(pk, &instance, &witness, rng)
    }


    /// Verify that [`proof`] is a valid proof with respect to [`vk`] and to the
    /// instance induced by [`cs`].
    fn verify_with_cs<CS: ConstraintSynthesizer<F>>(
        vk: &Self::VerifyingKey,
        cs: CS,
        proof: &Self::Proof,
    ) -> Result<bool, Self::Error> {
        let instance = Self::verifier_inputs(&cs);
        Self::verify(vk, &instance, proof)
    }

    /// Checks that [`proof`] is a valid proof of with respect to [`pvk`] and
    /// the instance induced by [`cs`].
    fn verify_with_cs_and_processed_vk<CS: ConstraintSynthesizer<F>>(
        pvk: &Self::ProcessedVerifyingKey,
        cs: CS,
        proof: &Self::Proof,
    ) -> Result<bool, Self::Error> {
        let instance = Self::verifier_inputs(&cs);
        Self::verify_with_processed_vk(pvk, &instance, proof)
    }
}
