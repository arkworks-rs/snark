use crate::{
    cfg_chunks_mut,
    curves::{batch_bucketed_add, BatchGroupArithmeticSlice, BucketPosition, BATCH_SIZE},
    fields::FpParameters,
    AffineCurve, PrimeField, ProjectiveCurve, Vec,
};
use num_traits::identities::Zero;

use core::fmt;

use rand::Rng;
#[cfg(feature = "parallel")]
use {rand::thread_rng, rayon::prelude::*};

#[derive(Debug, Clone)]
pub struct VerificationError;

impl fmt::Display for VerificationError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Verification Error. Not in subgroup")
    }
}

fn verify_points<C: AffineCurve, R: Rng>(
    points: &[C],
    num_buckets: usize,
    _new_security_param: Option<usize>, // Only pass new_security_param if possibly recursing
    rng: &mut R,
) -> Result<(), VerificationError> {
    let n_points = points.len();
    let mut bucket_assign = Vec::with_capacity(points.len());
    for i in 0..n_points {
        bucket_assign.push(BucketPosition {
            bucket: rng.gen_range(0, num_buckets) as u32,
            position: i as u32,
        });
    }
    let _now = timer!();
    let mut buckets = batch_bucketed_add(num_buckets, &mut points.to_vec(), &mut bucket_assign[..]);
    timer_println!(_now, format!("bucketed add({}, {})", num_buckets, n_points));

    // We use the batch_scalar_mul to check the subgroup condition if
    // there are sufficient number of buckets. For SW curves, the number
    // elems for the batch mul to become useful is around 2^24.
    let _now = timer!();
    let verification_failure = if num_buckets >= BATCH_SIZE {
        cfg_chunks_mut!(buckets, BATCH_SIZE).for_each(|e| {
            let length = e.len();
            e[..].batch_scalar_mul_in_place::<<C::ScalarField as PrimeField>::BigInt>(
                &mut vec![C::ScalarField::modulus().into(); length][..],
                4,
            );
        });
        !buckets.iter().all(|&p| p.is_zero())
    } else {
        !buckets
            .iter()
            .all(|&b| b.into_projective().mul(C::ScalarField::modulus()).is_zero())
    };
    timer_println!(_now, "mul by modulus");
    if verification_failure {
        return Err(VerificationError);
    }
    Ok(())
}

fn run_rounds<C: AffineCurve, R: Rng>(
    points: &[C],
    num_buckets: usize,
    num_rounds: usize,
    new_security_param: Option<usize>,
    rng: &mut R,
) -> Result<(), VerificationError> {
    #[cfg(feature = "parallel")]
    {
        if num_rounds > 2 {
            use std::sync::Arc;
            let ref_points = Arc::new(points.to_vec());
            let mut threads = vec![];
            for _ in 0..num_rounds {
                let ref_points_thread = ref_points.clone();
                // We only use std when a multicore environment is available
                threads.push(std::thread::spawn(
                    move || -> Result<(), VerificationError> {
                        let mut rng = &mut thread_rng();
                        verify_points(
                            &ref_points_thread[..],
                            num_buckets,
                            new_security_param,
                            &mut rng,
                        )?;
                        Ok(())
                    },
                ));
            }
            for thread in threads {
                thread.join().unwrap()?;
            }
        } else {
            for _ in 0..num_rounds {
                verify_points(points, num_buckets, new_security_param, rng)?;
            }
        }
    }

    #[cfg(not(feature = "parallel"))]
    for _ in 0..num_rounds {
        verify_points(points, num_buckets, new_security_param, rng)?;
    }

    Ok(())
}

pub fn batch_verify_in_subgroup<C: AffineCurve, R: Rng>(
    points: &[C],
    security_param: usize,
    rng: &mut R,
) -> Result<(), VerificationError> {
    #[cfg(feature = "std")]
    let cost_estimate = (<C::ScalarField as PrimeField>::Params::MODULUS_BITS as f64
        * (0.5 * 7.0 / 6.0 * 0.8 + 1.0 / 5.0))
        .ceil() as usize;
    #[cfg(not(feature = "std"))]
    let cost_estimate = <C::ScalarField as PrimeField>::Params::MODULUS_BITS as usize * 5 / 4;

    let (num_buckets, num_rounds, _) = get_max_bucket(
        security_param,
        points.len(),
        // We estimate the costs of a single scalar multiplication in the batch affine, w-NAF GLV
        // case as 7/6 * 0.5 * n_bits * 0.8 (doubling) + 0.5 * 1/(w + 1) * n_bits
        // (addition) We take into account that doubling in the batch add model is cheaper
        // as it requires less cache use
        cost_estimate,
    );
    run_rounds(points, num_buckets, num_rounds, None, rng)?;
    Ok(())
}

/// We get the greatest power of 2 number of buckets such that we minimise the
/// number of rounds while satisfying the constraint that
/// n_rounds * buckets * next_check_per_elem_cost < n
fn get_max_bucket(
    security_param: usize,
    n_elems: usize,
    next_check_per_elem_cost: usize,
) -> (usize, usize, usize) {
    #[cfg(feature = "std")]
    {
        let mut log2_num_buckets = 1f64;
        let num_rounds = |log2_num_buckets: f64| -> usize {
            (security_param as f64 / log2_num_buckets).ceil() as usize
        };

        while num_rounds(log2_num_buckets)
            * next_check_per_elem_cost
            * (2f64.powf(log2_num_buckets).ceil() as usize)
            < n_elems
            && num_rounds(log2_num_buckets + 0.1) > 1
        {
            log2_num_buckets += 0.1;
        }
        (
            2f64.powf(log2_num_buckets).ceil() as usize, // number of buckets
            num_rounds(log2_num_buckets),                // number of rounds
            log2_num_buckets.ceil() as usize,            // new security param
        )
    }

    #[cfg(not(feature = "std"))]
    {
        let mut log2_num_buckets: u32 = 1;
        let num_rounds = |log2_num_buckets: u32| -> usize {
            (security_param - 1) / (log2_num_buckets as usize) + 1
        };

        while num_rounds(log2_num_buckets)
            * next_check_per_elem_cost
            * (2_i32.pow(log2_num_buckets) as usize)
            < n_elems
            && num_rounds(log2_num_buckets + 1) > 1
        {
            log2_num_buckets += 1;
        }
        (
            2_i32.pow(log2_num_buckets) as usize, // number of buckets
            num_rounds(log2_num_buckets),         // number of rounds
            log2_num_buckets as usize,            // new security param
        )
    }
}
