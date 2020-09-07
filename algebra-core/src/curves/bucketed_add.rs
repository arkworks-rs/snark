use crate::{
    cfg_iter_mut,
    curves::{BatchGroupArithmeticSlice, BATCH_SIZE},
    log2, AffineCurve, Vec,
};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[cfg(feature = "std")]
use std::collections::HashMap;

const RATIO_MULTIPLIER: usize = 2;

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
