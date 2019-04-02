use crate::{BigInteger, FpParameters, PairingEngine, PrimeField, ProjectiveCurve};
use rayon::prelude::*;

pub struct FixedBaseMSM;

impl FixedBaseMSM {
    pub fn get_mul_window_size(num_scalars: usize) -> usize {
        if num_scalars < 32 {
            3
        } else {
            (f64::from(num_scalars as u32)).ln().ceil() as usize
        }
    }

    pub fn get_window_table<T: ProjectiveCurve>(
        scalar_size: usize,
        window: usize,
        g: T,
    ) -> Vec<Vec<T>> {
        let in_window = 1 << window;
        let outerc = (scalar_size + window - 1) / window;
        let last_in_window = 1 << (scalar_size - (outerc - 1) * window);

        let mut multiples_of_g = vec![vec![T::zero(); in_window]; outerc];

        let mut gouter = g;
        for outer in 0..outerc {
            let mut ginner = T::zero();
            let cur_in_window = if outer == outerc - 1 {
                last_in_window
            } else {
                in_window
            };
            for inner in 0..cur_in_window {
                multiples_of_g[outer][inner] = ginner;
                ginner = ginner + &gouter;
            }
            for _ in 0..window {
                gouter.double_in_place();
            }
        }
        multiples_of_g
    }

    pub fn windowed_mul<E: PairingEngine, T: ProjectiveCurve>(
        scalar_size: usize,
        window: usize,
        multiples_of_g: &Vec<Vec<T>>,
        scalar: &E::Fr,
    ) -> T {
        let outerc = (scalar_size + window - 1) / window;
        let mut scalar_val = scalar.into_repr().to_bits();
        scalar_val.reverse();

        let mut res = multiples_of_g[0][0];
        for outer in 0..outerc {
            let mut inner = 0usize;
            for i in 0..window {
                if outer * window + i < (<E::Fr as PrimeField>::Params::MODULUS_BITS as usize)
                    && scalar_val[outer * window + i]
                {
                    inner |= 1 << i;
                }
            }
            res = res + &multiples_of_g[outer][inner];
        }
        res
    }

    pub fn batch_mul<E: PairingEngine, T: ProjectiveCurve>(
        scalar_size: usize,
        window: usize,
        table: &Vec<Vec<T>>,
        v: &Vec<E::Fr>,
    ) -> Vec<T> {
        v.par_iter()
            .map(|e| Self::windowed_mul::<E, T>(scalar_size, window, table, e))
            .collect::<Vec<_>>()
    }
}
