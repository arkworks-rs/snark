use crate::{Error, Vec};
use core::{
    fmt::{Debug, Formatter, Result as FmtResult},
    marker::PhantomData,
};
use rand::Rng;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

use super::pedersen::{bytes_to_bits, PedersenCRH, PedersenWindow};
use crate::crh::FixedLengthCRH;
use algebra_core::{biginteger::BigInteger, fields::PrimeField, groups::Group};
use ff_fft::cfg_chunks;

#[cfg(feature = "r1cs")]
pub mod constraints;

pub const CHUNK_SIZE: usize = 3;

#[derive(Clone, Default)]
pub struct BoweHopwoodPedersenParameters<G: Group> {
    pub generators: Vec<Vec<G>>,
}

pub struct BoweHopwoodPedersenCRH<G: Group, W: PedersenWindow> {
    group:  PhantomData<G>,
    window: PhantomData<W>,
}

impl<G: Group, W: PedersenWindow> BoweHopwoodPedersenCRH<G, W> {
    pub fn create_generators<R: Rng>(rng: &mut R) -> Vec<Vec<G>> {
        let mut generators = Vec::new();
        for _ in 0..W::NUM_WINDOWS {
            let mut generators_for_segment = Vec::new();
            let mut base = G::rand(rng);
            for _ in 0..W::WINDOW_SIZE {
                generators_for_segment.push(base);
                for _ in 0..4 {
                    base.double_in_place();
                }
            }
            generators.push(generators_for_segment);
        }
        generators
    }
}

impl<G: Group, W: PedersenWindow> FixedLengthCRH for BoweHopwoodPedersenCRH<G, W> {
    const INPUT_SIZE_BITS: usize = PedersenCRH::<G, W>::INPUT_SIZE_BITS;
    type Output = G;
    type Parameters = BoweHopwoodPedersenParameters<G>;

    fn setup<R: Rng>(rng: &mut R) -> Result<Self::Parameters, Error> {
        fn calculate_num_chunks_in_segment<F: PrimeField>() -> usize {
            let upper_limit = F::modulus_minus_one_div_two();
            let mut c = 0;
            let mut range = F::BigInt::from(2_u64);
            while range < upper_limit {
                range.muln(4);
                c += 1;
            }

            c
        }

        let maximum_num_chunks_in_segment = calculate_num_chunks_in_segment::<G::ScalarField>();
        if W::WINDOW_SIZE > maximum_num_chunks_in_segment {
            return Err(format!(
                "Bowe-Hopwood hash must have a window size resulting in scalars < (p-1)/2, \
                 maximum segment size is {}",
                maximum_num_chunks_in_segment
            )
            .into());
        }

        let time = start_timer!(|| format!(
            "BoweHopwoodPedersenCRH::Setup: {} segments of {} 3-bit chunks; {{0,1}}^{{{}}} -> G",
            W::NUM_WINDOWS,
            W::WINDOW_SIZE,
            W::WINDOW_SIZE * W::NUM_WINDOWS * CHUNK_SIZE
        ));
        let generators = Self::create_generators(rng);
        end_timer!(time);
        Ok(Self::Parameters { generators })
    }

    fn evaluate(parameters: &Self::Parameters, input: &[u8]) -> Result<Self::Output, Error> {
        let eval_time = start_timer!(|| "BoweHopwoodPedersenCRH::Eval");

        if (input.len() * 8) > W::WINDOW_SIZE * W::NUM_WINDOWS * CHUNK_SIZE {
            panic!(
                "incorrect input length {:?} for window params {:?}x{:?}x{}",
                input.len(),
                W::WINDOW_SIZE,
                W::NUM_WINDOWS,
                CHUNK_SIZE,
            );
        }

        let mut padded_input = Vec::with_capacity(input.len());
        let input = bytes_to_bits(input);
        // Pad the input if it is not the current length.
        padded_input.extend_from_slice(&input);
        if input.len() % CHUNK_SIZE != 0 {
            let current_length = input.len();
            for _ in 0..(CHUNK_SIZE - current_length % CHUNK_SIZE) {
                padded_input.push(false);
            }
        }

        assert_eq!(padded_input.len() % CHUNK_SIZE, 0);

        assert_eq!(
            parameters.generators.len(),
            W::NUM_WINDOWS,
            "Incorrect pp of size {:?} for window params {:?}x{:?}x{}",
            parameters.generators.len(),
            W::WINDOW_SIZE,
            W::NUM_WINDOWS,
            CHUNK_SIZE,
        );
        for generators in parameters.generators.iter() {
            assert_eq!(generators.len(), W::WINDOW_SIZE);
        }
        assert_eq!(CHUNK_SIZE, 3);

        // Compute sum of h_i^{sum of
        // (1-2*c_{i,j,2})*(1+c_{i,j,0}+2*c_{i,j,1})*2^{4*(j-1)} for all j in segment}
        // for all i. Described in section 5.4.1.7 in the Zcash protocol
        // specification.

        let result = cfg_chunks!(padded_input, W::WINDOW_SIZE * CHUNK_SIZE)
            .zip(&parameters.generators)
            .map(|(segment_bits, segment_generators)| {
                cfg_chunks!(segment_bits, CHUNK_SIZE)
                    .zip(segment_generators)
                    .map(|(chunk_bits, generator)| {
                        let mut encoded = generator.clone();
                        if chunk_bits[0] {
                            encoded = encoded + generator;
                        }
                        if chunk_bits[1] {
                            encoded += &generator.double();
                        }
                        if chunk_bits[2] {
                            encoded = encoded.neg();
                        }
                        encoded
                    })
                    .sum::<G>()
            })
            .sum::<G>();

        end_timer!(eval_time);

        Ok(result)
    }
}

impl<G: Group> Debug for BoweHopwoodPedersenParameters<G> {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        write!(f, "Bowe-Hopwood Pedersen Hash Parameters {{\n")?;
        for (i, g) in self.generators.iter().enumerate() {
            write!(f, "\t  Generator {}: {:?}\n", i, g)?;
        }
        write!(f, "}}\n")
    }
}

#[cfg(test)]
mod test {
    use crate::{
        crh::{bowe_hopwood::BoweHopwoodPedersenCRH, pedersen::PedersenWindow},
        FixedLengthCRH,
    };
    use algebra::{jubjub::JubJubProjective, test_rng};

    #[test]
    fn test_simple_bh() {
        #[derive(Clone)]
        struct TestWindow {}
        impl PedersenWindow for TestWindow {
            const WINDOW_SIZE: usize = 63;
            const NUM_WINDOWS: usize = 8;
        }

        let rng = &mut test_rng();
        let params =
            <BoweHopwoodPedersenCRH<JubJubProjective, TestWindow> as FixedLengthCRH>::setup(rng)
                .unwrap();
        <BoweHopwoodPedersenCRH<JubJubProjective, TestWindow> as FixedLengthCRH>::evaluate(
            &params,
            &[1, 2, 3],
        )
        .unwrap();
    }
}
