use crate::{
    cfg_iter_mut,
    curves::{BatchGroupArithmeticSlice, BATCH_SIZE},
    log2, AffineCurve, Vec,
};

#[cfg(feature = "std")]
use {core::cmp::Ordering, std::collections::HashMap, voracious_radix_sort::*};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[derive(Copy, Clone, Debug)]
pub struct BucketPosition {
    pub bucket: u32,
    pub position: u32,
}

impl PartialOrd for BucketPosition {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.bucket.partial_cmp(&other.bucket)
    }
}
impl PartialEq for BucketPosition {
    fn eq(&self, other: &Self) -> bool {
        self.bucket == other.bucket
    }
}

impl Radixable<u32> for BucketPosition {
    type Key = u32;
    #[inline]
    fn key(&self) -> Self::Key {
        self.bucket
    }
}

const RATIO_MULTIPLIER: usize = 2;

#[inline]
#[cfg(feature = "std")]
pub fn batch_bucketed_add_radix<C: AffineCurve>(
    buckets: usize,
    elems: &[C],
    bucket_positions: &mut [BucketPosition],
) -> Vec<C> {
    assert_eq!(elems.len(), bucket_positions.len());
    assert!(elems.len() > 0);

    let mut len = bucket_positions.len();
    let mut all_ones = true;
    let mut new_len = 0; // len counter
    let mut glob = 0; // global counters
    let mut loc = 1; // local counter
    let mut batch = 0; // batch counter
    let mut instr = Vec::<(u32, u32)>::with_capacity(BATCH_SIZE);
    let mut new_elems = Vec::<C>::with_capacity(elems.len() * 3 / 8);

    let mut scratch_space = Vec::<Option<C>>::with_capacity(BATCH_SIZE / 2);

    // In the first loop, we copy the results of the first in place addition tree
    // to a local vector, new_elems
    // Subsequently, we perform all the operations in place

    while glob < len {
        let current_bucket = bucket_positions[glob].bucket;
        while glob + 1 < len && bucket_positions[glob + 1].bucket == current_bucket {
            glob += 1;
            loc += 1;
        }
        if current_bucket >= buckets as u32 {
            loc = 1;
        } else {
            // all ones is false if next len is not 1
            if loc > 2 {
                all_ones = false;
            }
            let is_odd = loc % 2 == 1;
            let half = loc / 2;
            for i in 0..half {
                instr.push((
                    bucket_positions[glob - (loc - 1) + 2 * i].position,
                    bucket_positions[glob - (loc - 1) + 2 * i + 1].position,
                ));
                bucket_positions[new_len + i] = BucketPosition {
                    bucket: current_bucket,
                    position: (new_len + i) as u32,
                };
            }
            if is_odd {
                instr.push((bucket_positions[glob].position, !0u32));
                bucket_positions[new_len + half] = BucketPosition {
                    bucket: current_bucket,
                    position: (new_len + half) as u32,
                };
            }
            // Reset the local_counter and update state
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
    let zero = C::zero();
    let mut res = vec![zero; buckets];

    for i in 0..len {
        let (pos, buc) = (bucket_positions[i].position, bucket_positions[i].bucket);
        res[buc as usize] = new_elems[pos as usize];
    }
    res
}

#[inline]
#[cfg(feature = "std")]
pub fn batch_bucketed_add<C: AffineCurve>(
    buckets: usize,
    elems: &mut [C],
    bucket_assign: &[usize],
) -> Vec<C> {
    let num_split = 2i32.pow(log2(buckets) / 2 + 2) as usize;
    let split_size = (buckets - 1) / num_split + 1;
    let mut bucket_split = vec![Vec::with_capacity(split_size); num_split];

    // Get the inverted index for the positions assigning to each buckets
    for (position, &bucket) in bucket_assign.iter().enumerate() {
        if bucket < buckets {
            bucket_split[bucket / split_size].push((bucket as u32, position as u32));
        }
    }

    let offset = ((elems.len() - 1) / buckets + 1) * RATIO_MULTIPLIER;
    let mut index = vec![0u32; offset * buckets];
    let mut assign_hash = HashMap::<usize, Vec<u32>>::new();

    for split in bucket_split {
        for (bucket, position) in split {
            let bucket = bucket as usize;
            let idx = bucket * offset;
            let n_assignments = index[idx] as usize;
            index[idx] += 1;
            // If we have run out of space for the fixed sized offsets, we add the assignments
            // to a dynamically sized vector stored in a hashmap
            if n_assignments >= offset - 1 {
                let assign_vec = assign_hash
                    .entry(bucket)
                    .or_insert(Vec::with_capacity(offset));
                if n_assignments == offset - 1 {
                    assign_vec.extend_from_slice(&index[idx + 1..idx + offset]);
                }
                assign_vec.push(position);
            } else {
                index[idx + n_assignments + 1] = position;
            }
        }
    }

    // Instructions for indexes for the in place addition tree
    let mut instr: Vec<Vec<(u32, u32)>> = vec![];
    // Find the maximum depth of the addition tree
    let max_depth = index
        .iter()
        .step_by(offset)
        .map(|x| log2(*x as usize))
        .max()
        .unwrap() as usize;

    // Generate in-place addition instructions that implement the addition tree
    // for each bucket from the leaves to the root
    for i in 0..max_depth {
        let mut instr_row = Vec::<(u32, u32)>::with_capacity(buckets);
        for bucket in 0..buckets {
            let idx = bucket * offset;
            let len = index[idx] as usize;

            if len > 1 << (max_depth - i - 1) {
                let new_len = (len - 1) / 2 + 1;
                // We must deal with vector
                if len > offset - 1 {
                    let assign_vec = assign_hash.entry(bucket).or_default();
                    if new_len <= offset - 1 {
                        for j in 0..len / 2 {
                            index[idx + j + 1] = assign_vec[2 * j];
                            instr_row.push((assign_vec[2 * j], assign_vec[2 * j + 1]));
                        }
                        if len % 2 == 1 {
                            index[idx + new_len] = assign_vec[len - 1];
                        }
                        assign_hash.remove(&bucket);
                    } else {
                        for j in 0..len / 2 {
                            assign_vec[j] = assign_vec[2 * j];
                            instr_row.push((assign_vec[2 * j], assign_vec[2 * j + 1]));
                        }
                        if len % 2 == 1 {
                            assign_vec[new_len - 1] = assign_vec[len - 1];
                        }
                    }
                } else {
                    for j in 0..len / 2 {
                        index[idx + j + 1] = index[idx + 2 * j + 1];
                        instr_row.push((index[idx + 2 * j + 1], index[idx + 2 * j + 2]));
                    }
                    if len % 2 == 1 {
                        index[idx + new_len] = index[idx + len];
                    }
                }
                // New length is the ceil of (old_length / 2)
                index[idx] = new_len as u32;
            }
        }
        if instr_row.len() > 0 {
            instr.push(instr_row);
        }
    }

    for instr_row in instr.iter() {
        for instr_chunk in C::get_chunked_instr::<(u32, u32)>(&instr_row[..], BATCH_SIZE).iter() {
            elems[..].batch_add_in_place_same_slice(&instr_chunk[..]);
        }
    }

    let zero = C::zero();
    let mut res = vec![zero; buckets];

    for bucket in 0..buckets {
        if index[offset * bucket] == 1 {
            res[bucket] = elems[index[offset * bucket + 1] as usize];
        } else if index[offset * bucket] > 1 {
            debug_assert!(false, "Did not successfully reduce index");
        }
    }
    res
}

#[cfg(not(feature = "std"))]
pub fn batch_bucketed_add<C: AffineCurve>(
    buckets: usize,
    elems: &mut [C],
    bucket_assign: &[usize],
) -> Vec<C> {
    let num_split = 2i32.pow(log2(buckets) / 2 + 2) as usize;
    let split_size = (buckets - 1) / num_split + 1;
    let ratio = elems.len() / buckets * 2;
    // Get the inverted index for the positions assigning to each bucket
    let mut bucket_split = vec![vec![]; num_split];
    let mut index = vec![Vec::with_capacity(ratio); buckets];

    for (position, &bucket) in bucket_assign.iter().enumerate() {
        // Check the bucket assignment is valid
        if bucket < buckets {
            // index[bucket].push(position);
            bucket_split[bucket / split_size].push((bucket, position));
        }
    }

    for split in bucket_split {
        for (bucket, position) in split {
            index[bucket].push(position as u32);
        }
    }

    // Instructions for indexes for the in place addition tree
    let mut instr: Vec<Vec<(u32, u32)>> = vec![];
    // Find the maximum depth of the addition tree
    let max_depth = index.iter()
        // log_2
        .map(|x| log2(x.len()))
        .max().unwrap();

    // Generate in-place addition instructions that implement the addition tree
    // for each bucket from the leaves to the root
    for i in 0..max_depth {
        let mut instr_row = Vec::<(u32, u32)>::with_capacity(buckets);
        for to_add in index.iter_mut() {
            if to_add.len() > 1 << (max_depth - i - 1) {
                let mut new_to_add = vec![];
                for j in 0..(to_add.len() / 2) {
                    new_to_add.push(to_add[2 * j]);
                    instr_row.push((to_add[2 * j], to_add[2 * j + 1]));
                }
                if to_add.len() % 2 == 1 {
                    new_to_add.push(*to_add.last().unwrap());
                }
                *to_add = new_to_add;
            }
        }
        instr.push(instr_row);
    }

    for instr_row in instr.iter() {
        for instr in C::get_chunked_instr::<(u32, u32)>(&instr_row[..], BATCH_SIZE).iter() {
            elems[..].batch_add_in_place_same_slice(&instr[..]);
        }
    }

    let zero = C::zero();
    let mut res = vec![zero; buckets];

    for (i, to_add) in index.iter().enumerate() {
        if to_add.len() == 1 {
            res[i] = elems[to_add[0] as usize];
        } else if to_add.len() > 1 {
            debug_assert!(false, "Did not successfully reduce to_add");
        }
    }
    res
}

// We make the batch bucket add cache-oblivious by splitting the problem
// into sub problems recursively
pub fn batch_bucketed_add_split<C: AffineCurve>(
    buckets: usize,
    elems: &[C],
    bucket_assign: &[usize],
    target_n_buckets_hint: usize,
) -> Vec<C> {
    let split_size = if buckets > 1 << target_n_buckets_hint {
        1 << target_n_buckets_hint
    } else {
        buckets
    };
    let num_split = (buckets - 1) / split_size + 1;
    let mut elem_split = vec![vec![]; num_split];
    let mut bucket_split = vec![vec![]; num_split];

    let split_window = 1 << 5;
    let split_split = (num_split - 1) / split_window + 1;

    for i in 0..split_split {
        for (position, &bucket) in bucket_assign.iter().enumerate() {
            let split_index = bucket / split_size;
            // Check the bucket assignment is valid
            if bucket < buckets
                && split_index >= i * split_window
                && split_index < (i + 1) * split_window
            {
                bucket_split[split_index].push(bucket % split_size);
                elem_split[split_index].push(elems[position]);
            }
        }
    }

    let res = cfg_iter_mut!(elem_split)
        .zip(cfg_iter_mut!(bucket_split))
        .filter(|(e, _)| e.len() > 0)
        .map(|(elems, buckets)| batch_bucketed_add(split_size, &mut elems[..], &buckets[..]))
        .flatten()
        .collect();
    res
}
