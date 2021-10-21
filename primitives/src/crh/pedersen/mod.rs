use crate::{bytes_to_bits, CryptoError, Error};
use rand::Rng;
use rayon::prelude::*;
use std::{
    fmt::{Debug, Formatter, Result as FmtResult},
    marker::PhantomData,
};

use crate::crh::FixedLengthCRH;
use algebra::{groups::Group, Field, ToConstraintField};
use serde::{Deserialize, Serialize};

pub trait PedersenWindow: Clone {
    const WINDOW_SIZE: usize;
    const NUM_WINDOWS: usize;
}

#[derive(Clone, Default, Serialize, Deserialize)]
#[serde(bound(deserialize = "G: Group"))]
pub struct PedersenParameters<G: Group> {
    pub generators: Vec<Vec<G>>,
}

pub struct PedersenCRH<G: Group, W: PedersenWindow> {
    group: PhantomData<G>,
    window: PhantomData<W>,
}

impl<G: Group, W: PedersenWindow> PedersenCRH<G, W> {
    pub fn create_generators<R: Rng>(rng: &mut R) -> Vec<Vec<G>> {
        let mut generators_powers = Vec::new();
        for _ in 0..W::NUM_WINDOWS {
            generators_powers.push(Self::generator_powers(W::WINDOW_SIZE, rng));
        }
        generators_powers
    }

    pub fn generator_powers<R: Rng>(num_powers: usize, rng: &mut R) -> Vec<G> {
        let mut cur_gen_powers = Vec::with_capacity(num_powers);
        let mut base = G::rand(rng);
        for _ in 0..num_powers {
            cur_gen_powers.push(base);
            base.double_in_place();
        }
        cur_gen_powers
    }
}

impl<G: Group, W: PedersenWindow> FixedLengthCRH for PedersenCRH<G, W> {
    const INPUT_SIZE_BITS: usize = W::WINDOW_SIZE * W::NUM_WINDOWS;
    type Output = G;
    type Parameters = PedersenParameters<G>;

    fn setup<R: Rng>(rng: &mut R) -> Result<Self::Parameters, Error> {
        let time = start_timer!(|| format!(
            "PedersenCRH::Setup: {} {}-bit windows; {{0,1}}^{{{}}} -> G",
            W::NUM_WINDOWS,
            W::WINDOW_SIZE,
            W::NUM_WINDOWS * W::WINDOW_SIZE
        ));
        let generators = Self::create_generators(rng);
        end_timer!(time);
        Ok(Self::Parameters { generators })
    }

    fn evaluate(parameters: &Self::Parameters, input: &[u8]) -> Result<Self::Output, Error> {
        let eval_time = start_timer!(|| "PedersenCRH::Eval");

        if (input.len() * 8) > W::WINDOW_SIZE * W::NUM_WINDOWS {
            return Err(Box::new(CryptoError::Other(
                format!(
                    "incorrect input length {:?} for window params {:?}x{:?}",
                    input.len(),
                    W::WINDOW_SIZE,
                    W::NUM_WINDOWS
                )
                .to_owned(),
            )));
        }

        let mut padded_input = Vec::with_capacity(input.len());
        let mut input = input;
        // Pad the input if it is not the current length.
        if (input.len() * 8) < W::WINDOW_SIZE * W::NUM_WINDOWS {
            let current_length = input.len();
            padded_input.extend_from_slice(input);
            for _ in current_length..((W::WINDOW_SIZE * W::NUM_WINDOWS) / 8) {
                padded_input.push(0u8);
            }
            input = padded_input.as_slice();
        }

        if parameters.generators.len() != W::NUM_WINDOWS {
            Err(Box::new(CryptoError::Other(
                format!(
                    "Incorrect pp of size {:?}x{:?} for window params {:?}x{:?}",
                    parameters.generators[0].len(),
                    parameters.generators.len(),
                    W::WINDOW_SIZE,
                    W::NUM_WINDOWS
                )
                .to_owned(),
            )))?
        }

        // Compute sum of h_i^{m_i} for all i.
        let result = bytes_to_bits(input)
            .par_chunks(W::WINDOW_SIZE)
            .zip(&parameters.generators)
            .map(|(bits, generator_powers)| {
                let mut encoded = G::zero();
                for (bit, base) in bits.iter().zip(generator_powers.iter()) {
                    if *bit {
                        encoded = encoded + base;
                    }
                }
                encoded
            })
            .reduce(|| G::zero(), |a, b| a + &b);
        end_timer!(eval_time);

        Ok(result)
    }
}

impl<G: Group> Debug for PedersenParameters<G> {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        write!(f, "Pedersen Hash Parameters {{\n")?;
        for (i, g) in self.generators.iter().enumerate() {
            write!(f, "\t  Generator {}: {:?}\n", i, g)?;
        }
        write!(f, "}}\n")
    }
}

impl<G: Group> PedersenParameters<G> {
    pub fn check_consistency(&self) -> bool {
        for (i, p1) in self.generators.iter().enumerate() {
            if p1[0] == G::zero() {
                return false; // infinity generator
            }
            for p2 in self.generators.iter().skip(i + 1) {
                if p1[0] == p2[0] {
                    return false; // duplicate generator
                }
                if p1[0] == p2[0].neg() {
                    return false; // inverse generator
                }
            }
        }
        return true;
    }
}

impl<ConstraintF: Field, G: Group + ToConstraintField<ConstraintF>> ToConstraintField<ConstraintF>
    for PedersenParameters<G>
{
    #[inline]
    fn to_field_elements(&self) -> Result<Vec<ConstraintF>, Error> {
        Ok(Vec::new())
    }
}
