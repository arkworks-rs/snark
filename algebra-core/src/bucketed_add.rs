use crate::{cfg_iter_mut, curves::BatchGroupArithmeticSlice, AffineCurve};

// #[cfg(feature = "parallel")]
// use rayon::prelude::*;

const BATCH_ADD_SIZE: usize = 4096;

// We make the batch bucket add cache-oblivious by splitting the problem
// into sub problems recursively
pub fn batch_bucketed_add_split<C: AffineCurve>(
    buckets: usize,
    elems: &[C],
    bucket_assign: &[usize],
    bucket_size: usize,
) -> Vec<C> {
    let split_size = if buckets >= 1 << 26 {
        1 << 16
    } else {
        1 << bucket_size
    };
    let num_split = (buckets - 1) / split_size + 1;
    // println!("{}, {}", split_size, num_split);
    let mut elem_split = vec![vec![]; num_split];
    let mut bucket_split = vec![vec![]; num_split];

    let now = std::time::Instant::now();

    for (position, &bucket) in bucket_assign.iter().enumerate() {
        bucket_split[bucket / split_size].push(bucket % split_size);
        elem_split[bucket / split_size].push(elems[position]);
    }

    // println!(
    //     "\nAssign bucket and elem split: {:?}",
    //     now.elapsed().as_micros()
    // );

    let now = std::time::Instant::now();

    let res = if split_size < 1 << (bucket_size + 1) {
        cfg_iter_mut!(elem_split)
            .zip(cfg_iter_mut!(bucket_split))
            .map(|(elems, bucket)| batch_bucketed_add(split_size, &mut elems[..], &bucket[..]))
            .flatten()
            .collect()
    } else {
        // println!("CALLING RECURSIVE");
        elem_split
            .iter()
            .zip(bucket_split.iter())
            .map(|(elems, bucket)| {
                batch_bucketed_add_split(split_size, &elems[..], &bucket[..], bucket_size)
            })
            .flatten()
            .collect()
    };

    // println!("Bucketed add: {:?}", now.elapsed().as_micros());
    res
}

pub fn batch_bucketed_add<C: AffineCurve>(
    buckets: usize,
    elems: &mut [C],
    bucket_assign: &[usize],
) -> Vec<C> {
    let num_split = if buckets >= 1 << 14 { 4096 } else { 1 };
    let split_size = buckets / num_split;
    let ratio = elems.len() / buckets * 2;
    // Get the inverted index for the positions assigning to each bucket
    let now = std::time::Instant::now();
    let mut bucket_split = vec![vec![]; num_split];
    let mut index = vec![Vec::with_capacity(ratio); buckets];

    // We use two levels of assignments to help with cache locality.
    for (position, &bucket) in bucket_assign.iter().enumerate() {
        bucket_split[bucket / split_size].push((bucket, position));
    }

    for split in bucket_split {
        for (bucket, position) in split {
            index[bucket].push(position);
        }
    }
    // println!("\nGenerate Inverted Index: {:?}", now.elapsed().as_micros());

    // Instructions for indexes for the in place addition tree
    let mut instr: Vec<Vec<(usize, usize)>> = vec![];
    // Find the maximum depth of the addition tree
    let max_depth = index.iter()
        // log_2
        .map(|x| crate::log2(x.len()))
        .max().unwrap();

    let now = std::time::Instant::now();
    // Generate in-place addition instructions that implement the addition tree
    // for each bucket from the leaves to the root
    for i in 0..max_depth {
        let mut instr_row = Vec::<(usize, usize)>::with_capacity(buckets);
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
    // println!("Generate Instr: {:?}", now.elapsed().as_micros());

    let now = std::time::Instant::now();
    // let mut elems_mut_1 = elems.to_vec();

    for instr_row in instr.iter() {
        for instr in C::get_chunked_instr::<(usize, usize)>(&instr_row[..], BATCH_ADD_SIZE).iter() {
            elems[..].batch_add_in_place_same_slice(&instr[..]);
        }
    }
    // println!("Batch add in place: {:?}", now.elapsed().as_micros());

    let now = std::time::Instant::now();
    let zero = C::zero();
    let mut res = vec![zero; buckets];

    for (i, to_add) in index.iter().enumerate() {
        if to_add.len() > 1 {
            panic!("Did not successfully reduce to_add");
        } else if to_add.len() == 1 {
            res[i] = elems[to_add[0]];
        }
    }

    // println!("Reassign: {:?}", now.elapsed().as_micros());
    res
}
