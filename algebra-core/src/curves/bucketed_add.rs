use crate::{
    curves::{BatchGroupArithmeticSlice, BATCH_SIZE},
    AffineCurve, Vec,
};

#[derive(Copy, Clone, Debug)]
pub struct BucketPosition {
    pub bucket: u32,
    pub position: u32,
}

/// The objective of this function is to identify an addition tree of
/// independent elliptic curve group additions for each bucket, and to batch the
/// independent additions using the batch affine inversion method.

/// The strategy taken is to sort a list of bucket assignments of all the
/// elements (which we can for most intents and purposes, think of as being
/// uniformly random) by bucket, so that indices corresponding to elements that
/// must be added together are physically collocated in memory. Then, in the
/// first round, we proceed to perform independent additions producing
/// intermediate results at the greatest depth for each addition tree (each
/// corresponding to a bucket), and write the result to a new vector. We do so
/// to improve cache locality for future rounds, and take advantage of the
/// CPU-intensive nature of elliptic curve operations along with prfetching to
/// hide the latency of reading from essentially random locations in memory.

/// Subsequently, we perform the additions in place, and the second operands
/// become junk data. Finally, when we only have the buckets left (no more
/// additions left to perform), we copy the result into a destination `res`
/// slice.
#[inline]
pub fn batch_bucketed_add<C: AffineCurve>(
    buckets: usize,
    elems: &[C],
    bucket_positions: &mut [BucketPosition],
) -> Vec<C> {
    assert_eq!(elems.len(), bucket_positions.len());
    assert!(elems.len() > 0);

    let _now = timer!();
    // We sort the bucket positions so that indices of elements assigned
    // to the same bucket are continguous. This way, we can easily identify
    // how to construct the addition tree for that bucket.
    bucket_positions.sort_unstable_by_key(|x| x.bucket);
    timer_println!(_now, "sort");

    let mut len = bucket_positions.len();
    let mut all_ones = true;
    let mut new_len = 0; // len counter
    let mut glob = 0; // global counters
    let mut loc = 1; // local counter
    let mut batch = 0; // batch counter
    let mut instr = Vec::<(u32, u32)>::with_capacity(BATCH_SIZE);
    let mut new_elems = Vec::<C>::with_capacity(elems.len() * 3 / 8);

    let mut scratch_space = Vec::<Option<C>>::with_capacity(BATCH_SIZE / 2);

    let _now = timer!();
    // In the first loop, we copy the results of the first in place addition tree
    // to a local vector, new_elems
    // Subsequently, we perform all the operations in place
    while glob < len {
        let current_bucket = bucket_positions[glob].bucket;
        // We are iterating over elements using a global `glob` counter, and counting
        // how many in a row are being assigned to the same bucket, using the `loc`
        // counter.
        while glob + 1 < len && bucket_positions[glob + 1].bucket == current_bucket {
            glob += 1;
            loc += 1;
        }
        // If the current bucket exceeds buckets, it encodes a noop
        if current_bucket >= buckets as u32 {
            loc = 1;
        } else if loc > 1 {
            // all ones is false if next len is not 1

            // in other words, we have not reached the terminating
            // condition that after the current round of addition
            // there is only one element left in each addition tree

            // This would be the case, if each addition tree had at
            // most 2 elements in the current round.
            if loc > 2 {
                all_ones = false;
            }
            let is_odd = loc % 2 == 1;
            let half = loc / 2;
            // We encode instructions to add adjacent elements
            for i in 0..half {
                instr.push((
                    bucket_positions[glob - (loc - 1) + 2 * i].position,
                    bucket_positions[glob - (loc - 1) + 2 * i + 1].position,
                ));
                // Compactification of buckets
                bucket_positions[new_len + i] = BucketPosition {
                    bucket: current_bucket,
                    position: (new_len + i) as u32,
                };
            }
            // If there are an odd number of elements, the lone element
            // without a partner will be copied over to the `new_elems`
            // vector, a noop which is encoded as !0u32
            if is_odd {
                instr.push((bucket_positions[glob].position, !0u32));
                bucket_positions[new_len + half] = BucketPosition {
                    bucket: current_bucket,
                    position: (new_len + half) as u32,
                };
            }
            // Reset the local_counter and update state

            // We compactify the `bucket_positions` data by shifing left
            // `new_len` is the len of the current compactified vector.

            // We also update the `batch` counter to decide when it is
            // optimal to invoke the batch inversion, i.e. when we have
            // accumulated enough independent additions.
            new_len += half + (loc % 2);
            batch += half;
            loc = 1;

            if batch >= BATCH_SIZE / 2 {
                // We need instructions for copying data in the case
                // of noops. We encode noops/copies as !0u32
                elems[..].batch_add_write(&instr[..], &mut new_elems, &mut scratch_space);

                instr.clear();
                batch = 0;
            }
        } else {
            instr.push((bucket_positions[glob].position, !0u32));
            bucket_positions[new_len] = BucketPosition {
                bucket: current_bucket,
                position: new_len as u32,
            };
            new_len += 1;
        }
        glob += 1;
    }
    if instr.len() > 0 {
        elems[..].batch_add_write(&instr[..], &mut new_elems, &mut scratch_space);
        instr.clear();
    }
    glob = 0;
    batch = 0;
    loc = 1;
    len = new_len;
    new_len = 0;

    // We repeat the above procedure, except, since we are performing the addition
    // trees in place, we do not need to encode noops to force a copy to a new
    // vector.
    while !all_ones {
        all_ones = true;
        while glob < len {
            let current_bucket = bucket_positions[glob].bucket;
            while glob + 1 < len && bucket_positions[glob + 1].bucket == current_bucket {
                glob += 1;
                loc += 1;
            }
            if current_bucket >= buckets as u32 {
                loc = 1;
            } else if loc > 1 {
                // all ones is false if next len is not 1
                if loc != 2 {
                    all_ones = false;
                }
                let is_odd = loc % 2 == 1;
                let half = loc / 2;
                for i in 0..half {
                    instr.push((
                        bucket_positions[glob - (loc - 1) + 2 * i].position,
                        bucket_positions[glob - (loc - 1) + 2 * i + 1].position,
                    ));
                    bucket_positions[new_len + i] = bucket_positions[glob - (loc - 1) + 2 * i];
                }
                if is_odd {
                    bucket_positions[new_len + half] = bucket_positions[glob];
                }
                // Reset the local_counter and update state
                new_len += half + (loc % 2);
                batch += half;
                loc = 1;

                if batch >= BATCH_SIZE / 2 {
                    &mut new_elems[..].batch_add_in_place_same_slice(&instr[..]);
                    instr.clear();
                    batch = 0;
                }
            } else {
                bucket_positions[new_len] = bucket_positions[glob];
                new_len += 1;
            }
            glob += 1;
        }
        if instr.len() > 0 {
            &mut new_elems[..].batch_add_in_place_same_slice(&instr[..]);
            instr.clear();
        }
        glob = 0;
        batch = 0;
        loc = 1;
        len = new_len;
        new_len = 0;
    }
    timer_println!(_now, "addition tree");

    let zero = C::zero();
    let mut res = vec![zero; buckets];

    let _now = timer!();
    for i in 0..len {
        let (pos, buc) = (bucket_positions[i].position, bucket_positions[i].bucket);
        res[buc as usize] = new_elems[pos as usize];
    }
    timer_println!(_now, "reassign");
    res
}
