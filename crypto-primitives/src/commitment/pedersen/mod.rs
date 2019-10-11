use crate::Error;
use algebra::UniformRand;
use algebra::{Field, ToConstraintField};
use algebra::{bytes::ToBytes, groups::Group, BitIterator, FpParameters, PrimeField};

use rand::Rng;
use std::marker::PhantomData;

use super::CommitmentScheme;
use std::io::{Result as IoResult, Write};

pub use crate::crh::pedersen::PedersenWindow;
use crate::crh::{
    pedersen::{PedersenCRH, PedersenParameters as PedersenCRHParameters},
    FixedLengthCRH,
};

#[cfg(feature = "r1cs")]
pub mod constraints;

#[derive(Clone)]
pub struct PedersenParameters<G: Group> {
    pub randomness_generator: Vec<G>,
    pub generators:           Vec<Vec<G>>,
}

pub struct PedersenCommitment<G: Group, W: PedersenWindow> {
    group:  PhantomData<G>,
    window: PhantomData<W>,
}

#[derive(Derivative)]
#[derivative(
    Clone(bound = "G: Group"),
    PartialEq(bound = "G: Group"),
    Debug(bound = "G: Group"),
    Eq(bound = "G: Group"),
    Default(bound = "G: Group")
)]
pub struct PedersenRandomness<G: Group>(pub G::ScalarField);

impl<G: Group> UniformRand for PedersenRandomness<G> {
    #[inline]
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        PedersenRandomness(UniformRand::rand(rng))
    }
}

impl<G: Group> ToBytes for PedersenRandomness<G> {
    fn write<W: Write>(&self, writer: W) -> IoResult<()> {
        self.0.write(writer)
    }
}

impl<G: Group, W: PedersenWindow> CommitmentScheme for PedersenCommitment<G, W> {
    type Parameters = PedersenParameters<G>;
    type Randomness = PedersenRandomness<G>;
    type Output = G;

    fn setup<R: Rng>(rng: &mut R) -> Result<Self::Parameters, Error> {
        let time = start_timer!(|| format!(
            "PedersenCOMM::Setup: {} {}-bit windows; {{0,1}}^{{{}}} -> G",
            W::NUM_WINDOWS,
            W::WINDOW_SIZE,
            W::NUM_WINDOWS * W::WINDOW_SIZE
        ));
        let num_powers = <G::ScalarField as PrimeField>::Params::MODULUS_BITS as usize;
        let randomness_generator = PedersenCRH::<_, W>::generator_powers(num_powers, rng);
        let generators = PedersenCRH::<_, W>::create_generators(rng);
        end_timer!(time);

        Ok(Self::Parameters {
            randomness_generator,
            generators,
        })
    }

    fn commit(
        parameters: &Self::Parameters,
        input: &[u8],
        randomness: &Self::Randomness,
    ) -> Result<Self::Output, Error> {
        let commit_time = start_timer!(|| "PedersenCOMM::Commit");
        // If the input is too long, return an error.
        if input.len() > W::WINDOW_SIZE * W::NUM_WINDOWS {
            panic!("incorrect input length: {:?}", input.len());
        }
        // Pad the input to the necessary length.
        let mut padded_input = Vec::with_capacity(input.len());
        let mut input = input;
        if (input.len() * 8) < W::WINDOW_SIZE * W::NUM_WINDOWS {
            let current_length = input.len();
            padded_input.extend_from_slice(input);
            for _ in current_length..((W::WINDOW_SIZE * W::NUM_WINDOWS) / 8) {
                padded_input.push(0u8);
            }
            input = padded_input.as_slice();
        }
        assert_eq!(parameters.generators.len(), W::NUM_WINDOWS);

        // Invoke Pedersen CRH here, to prevent code duplication.

        let crh_parameters = PedersenCRHParameters {
            generators: parameters.generators.clone(),
        };
        let mut result = PedersenCRH::<_, W>::evaluate(&crh_parameters, &input)?;
        let randomize_time = start_timer!(|| "Randomize");

        // Compute h^r.
        let mut scalar_bits = BitIterator::new(randomness.0.into_repr()).collect::<Vec<_>>();
        scalar_bits.reverse();
        for (bit, power) in scalar_bits
            .into_iter()
            .zip(&parameters.randomness_generator)
        {
            if bit {
                result += power
            }
        }
        end_timer!(randomize_time);
        end_timer!(commit_time);

        Ok(result)
    }
}


impl<ConstraintF: Field, G: Group + ToConstraintField<ConstraintF>> ToConstraintField<ConstraintF> for PedersenParameters<G> {
    #[inline]
    fn to_field_elements(&self) -> Result<Vec<ConstraintF>, Error> {
        Ok(Vec::new())
    }
}
