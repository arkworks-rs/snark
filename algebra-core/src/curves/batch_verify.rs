use crate::fields::FpParameters;
use crate::{
    cfg_chunks_mut,
    curves::{batch_bucketed_add_split, BatchGroupArithmeticSlice},
    AffineCurve, PrimeField, ProjectiveCurve, Vec,
};
use num_traits::{identities::Zero, Pow};

use core::fmt;
#[cfg(feature = "parallel")]
use rand::thread_rng;
use rand::Rng;

const MAX_BUCKETS_FOR_FULL_CHECK: usize = 2;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[derive(Debug, Clone)]
pub struct VerificationError;

impl fmt::Display for VerificationError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Verification Error. Not in subgroup")
    }
}

// Only pass new_security_param if possibly recursing
fn verify_points<C: AffineCurve, R: Rng>(
    points: &[C],
    num_buckets: usize,
    new_security_param: Option<usize>,
    rng: &mut R,
) -> Result<(), VerificationError> {
    let mut bucket_assign = Vec::with_capacity(points.len());
    for _ in 0..points.len() {
        bucket_assign.push(rng.gen_range(0, num_buckets));
    }
    let mut buckets = batch_bucketed_add_split(num_buckets, points, &bucket_assign[..], 12);

    // Check that all the buckets belong to the subgroup, either by calling
    // the batch verify recusively, or by directly checking by multiplying by group order
    // when the number of buckets is small enough
    if num_buckets <= MAX_BUCKETS_FOR_FULL_CHECK || new_security_param == None {
        // We use the batch scalar mul to check the subgroup condition if
        // there are sufficient number of buckets
        let verification_failure = if num_buckets >= 4096 {
            cfg_chunks_mut!(buckets, 4096).for_each(|e| {
                let length = e.len();
                e[..].batch_scalar_mul_in_place::<<C::ScalarField as PrimeField>::BigInt>(
                    &mut vec![C::ScalarField::modulus().into(); length][..],
                    1,
                );
            });
            !buckets.iter().all(|&p| p == C::zero())
        } else {
            !buckets
                .iter()
                .all(|&b| b.mul(C::ScalarField::modulus()) == C::Projective::zero())
        };
        if verification_failure {
            return Err(VerificationError);
        }
    } else {
        // Since !new_security_param.is_none():
        let new_security_param = new_security_param.unwrap();

        // Temporarily commented out until a fix can be found for the recursive version of the test

        // if buckets.len() > 4096 {
        //     batch_verify_in_subgroup_recursive(&buckets[..], new_security_param, rng)?;
        // } else {

        batch_verify_in_subgroup_proj(
            &buckets
                .iter()
                .map(|&p| p.into())
                .collect::<Vec<C::Projective>>()[..],
            new_security_param,
            rng,
        )?;

        // }
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

    #[cfg(not(feature = "parallel"))]
    {
        for _ in 0..num_rounds {
            verify_points(points, num_buckets, new_security_param, rng)?;
        }
    }

    Ok(())
}

pub fn batch_verify_in_subgroup<C: AffineCurve, R: Rng>(
    points: &[C],
    security_param: usize,
    rng: &mut R,
) -> Result<(), VerificationError> {
    let (num_buckets, num_rounds, _) = get_max_bucket(
        security_param,
        points.len(),
        <C::ScalarField as PrimeField>::Params::MODULUS_BITS as usize,
    );
    run_rounds(points, num_buckets, num_rounds, None, rng)?;
    Ok(())
}

/// Temporarily commented out until a fix can be found for the recursive version of the test

// pub fn batch_verify_in_subgroup_recursive<C: AffineCurve, R: Rng>(
//     points: &[C],
//     security_param: usize,
//     rng: &mut R,
// ) -> Result<(), VerificationError> {
//     // we add security for maximum depth, as recursive depth adds additional error to error bound
//     let security_param = security_param + (log2(log2(security_param) as usize) as usize) + 1;
//     let (num_buckets, num_rounds, new_security_param) =
//         get_max_bucket(security_param, points.len(), 2);
//     run_rounds(points, num_buckets, num_rounds, Some(new_security_param), rng)?;
//     Ok(())
// }

pub fn batch_verify_in_subgroup_proj<C: ProjectiveCurve, R: Rng>(
    points: &[C],
    security_param: usize,
    rng: &mut R,
) -> Result<(), VerificationError> {
    let (num_buckets, num_rounds, new_security_param) =
        get_max_bucket(security_param, points.len(), 2);

    for _ in 0..num_rounds {
        let mut bucket_assign = Vec::with_capacity(points.len());
        for _ in 0..points.len() {
            bucket_assign.push(rng.gen_range(0, num_buckets));
        }
        // If our batch size is too small, we do the naive bucket add
        let zero = C::zero();
        let mut buckets = vec![zero; num_buckets];
        for (p, a) in points.iter().zip(bucket_assign) {
            buckets[a].add_assign(p);
        }

        if num_buckets <= MAX_BUCKETS_FOR_FULL_CHECK {
            if !buckets
                .iter()
                .all(|b| b.mul(<C::ScalarField as PrimeField>::Params::MODULUS) == C::zero())
            {
                return Err(VerificationError);
            }
        } else {
            batch_verify_in_subgroup_proj(&buckets[..], new_security_param, rng)?;
        }
    }
    Ok(())
}

// We get the greatest power of 2 number of buckets
// such that we minimise the number of rounds
// while satisfying the constraint that
// number of rounds * buckets * next_check_per_elem_cost < n
fn get_max_bucket(
    security_param: usize,
    n_elems: usize,
    next_check_per_elem_cost: usize,
) -> (usize, usize, usize) {
    let mut log2_num_buckets = 1;
    let num_rounds =
        |log2_num_buckets: usize| -> usize { (security_param - 1) / log2_num_buckets + 1 };

    while num_rounds(log2_num_buckets)
        * next_check_per_elem_cost
        * (2.pow(log2_num_buckets) as usize)
        < n_elems
        && num_rounds(log2_num_buckets) > 1
    {
        log2_num_buckets += 1;
    }
    (
        2.pow(log2_num_buckets) as usize, // number of buckets
        num_rounds(log2_num_buckets),     // number of rounds
        log2_num_buckets,                 // new security param
    )
}
