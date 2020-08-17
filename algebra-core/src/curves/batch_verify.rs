use crate::fields::FpParameters;
use crate::{batch_bucketed_add_split, log2, AffineCurve, PrimeField};
use num_traits::{identities::Zero, Pow};
use rand::thread_rng;
use rand::Rng;
use std::fmt;

#[derive(Debug, Clone)]
pub struct VerificationError;

impl fmt::Display for VerificationError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Verification Error. Not in subgroup")
    }
}

pub fn batch_verify_in_subgroup<C: AffineCurve>(
    points: &[C],
    security_param: usize,
) -> Result<(), VerificationError> {
    let (num_buckets, num_rounds) = get_max_bucket(security_param, points.len());
    let rng = &mut thread_rng();

    for _ in 0..num_rounds {
        let mut bucket_assign = Vec::with_capacity(points.len());
        for _ in 0..points.len() {
            bucket_assign.push(rng.gen_range(0, num_buckets));
        }
        let buckets = batch_bucketed_add_split(num_buckets, points, &bucket_assign[..]);

        if num_buckets <= 3 {
            if !buckets.iter().all(|b| {
                b.mul(<C::ScalarField as PrimeField>::Params::MODULUS) == C::Projective::zero()
            }) {
                return Err(VerificationError);
            }
        } else {
            batch_verify_in_subgroup(&buckets[..], log2(num_buckets) as usize)?;
        }
    }
    Ok(())
}

// We get the greatest power of 2 number of buckets
// such that we minimise the number of rounds
// while satisfying the constraint that number of rounds * buckets * 2 < n

// Number of buckets is always greater than new security param.
// So only need 1 round subsequently
fn get_max_bucket(security_param: usize, n_elems: usize) -> (usize, usize) {
    let mut log2_num_buckets = 1;
    let num_rounds =
        |log2_num_buckets: usize| -> usize { (security_param - 1) / log2_num_buckets + 1 };

    while num_rounds(log2_num_buckets) * 2 * (2.pow(log2_num_buckets) as usize) < n_elems
        && num_rounds(log2_num_buckets) > 1
    {
        log2_num_buckets += 1;
    }
    (
        2.pow(log2_num_buckets) as usize,
        num_rounds(log2_num_buckets),
    )
}
