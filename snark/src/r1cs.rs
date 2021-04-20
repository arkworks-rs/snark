use ark_relations::r1cs::{R1CS, ConstraintSynthesizer};

/// The basic functionality for a SNARK.
pub trait SNARKForR1CS<F: Field>: SNARK<R1CS<F>> {
    /// Generates a proof of satisfaction of the arithmetic circuit C (specified
    /// as R1CS constraints).
    fn prove<C: ConstraintSynthesizer<F>, Rng: RngCore + CryptoRng>(
        circuit_pk: &Self::ProvingKey,
        witness_generator: C,
        rng: &mut Rng,
    ) -> Result<Self::Proof, Self::Error>;

    /// Checks that `proof` is a valid proof of the satisfaction of circuit
    /// encoded in `circuit_vk`, with respect to the public input `public_input`,
    /// specified as R1CS constraints.
    fn verify<C: ConstraintSynthesizer<F>>(
        circuit_vk: &Self::VerifyingKey,
        public_input_generator: C,
        proof: &Self::Proof,
    ) -> Result<bool, Self::Error> {
        let pvk = Self::process_vk(circuit_vk)?;
        Self::verify_with_processed_vk(&pvk, public_input_generator, proof)
    }

    /// Preprocesses `circuit_vk` to enable faster verification.
    fn process_vk(
        circuit_vk: &Self::VerifyingKey,
    ) -> Result<Self::ProcessedVerifyingKey, Self::Error>;

    /// Checks that `proof` is a valid proof of the satisfaction of circuit
    /// encoded in `circuit_pvk`, with respect to the public input `public_input`,
    /// specified as R1CS constraints.
    fn verify_with_processed_vk<C: ConstraintSynthesizer<F>>(
        circuit_pvk: &Self::ProcessedVerifyingKey,
        public_input_generator: C,
        proof: &Self::Proof,
    ) -> Result<bool, Self::Error>;
}

/// A SNARK with (only) circuit-specific setup.
pub trait CircuitSpecificSetupSNARKForR1CS<F: Field>: CircuitSpecificSetupSNARK<R1CS<F>> + SNARK<R1CS<F>> {
    /// The setup algorithm for circuit-specific SNARKs.
    fn setup<C: ConstraintSynthesizer<F>, Rng: RngCore + CryptoRng>(
        index_generator: C,
        rng: &mut Rng,
    ) -> Result<(Self::ProvingKey, Self::VerifyingKey), Self::Error>;
}

/// A SNARK with universal setup. That is, a SNARK where the trusted setup is
/// circuit-independent.
pub trait UniversalSetupSNARKForR1CS<F: Field>: UniversalSetupSNARK<R1CS<F>> + SNARKForR1CS<F> {
    /// Specifies how to bound the size of public parameters required to
    /// generate the index proving and verification keys for a given
    /// circuit.
    fn universal_setup<Rng: RngCore + CryptoRng>(
        compute_bound: &Self::ComputationBound,
        rng: &mut Rng,
    ) -> Result<Self::PublicParameters, Self::Error>;

    /// Indexes the public parameters according to the circuit `circuit`, and
    /// outputs circuit-specific proving and verification keys.
    fn index<C: ConstraintSynthesizer<F>>(
        pp: &Self::PublicParameters,
        index_generator: C,
        rng: &mut Rng,
    ) -> Result<
        (Self::ProvingKey, Self::VerifyingKey),
        IndexingError<Self::ComputationBound, Self::Error>,
    >;
}
