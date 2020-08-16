use crate::{AffineCurve, curves::BatchGroupArithmeticSlice};
use std::cmp::Ordering;
use voracious_radix_sort::*;

const BATCH_ADD_SIZE: usize = 4096;

#[derive(Copy, Clone, Debug)]
struct ReverseIndex {
    pos: usize,
    bucket: u64,
}

impl Radixable<u64> for ReverseIndex {
    type Key = u64;
    #[inline]
    fn key(&self) -> Self::Key {
        self.bucket
    }
}

impl PartialOrd for ReverseIndex {
    fn partial_cmp(&self, other: &ReverseIndex) -> Option<Ordering> {
        self.bucket.partial_cmp(&other.bucket)
    }
}

impl PartialEq for ReverseIndex {
    fn eq(&self, other: &Self) -> bool {
        self.bucket == other.bucket
    }
}

pub fn batch_bucketed_add<C: AffineCurve>(
    buckets: usize,
    elems: &[C],
    bucket_assign: &[usize],
) -> Vec<C> {

    let now = std::time::Instant::now();
    // let mut index = vec![Vec::with_capacity(8); buckets];
    // for (position, &bucket) in bucket_assign.iter().enumerate() {
    //     index[bucket].push(position);
    // }
    // Instead of the above, we do a radix sort by bucket value instead, and store offsets

    let mut index = vec![Vec::with_capacity(8); buckets];
    let mut to_sort = bucket_assign.iter()
        .enumerate()
        .map(|(pos, bucket)| ReverseIndex{ pos, bucket: *bucket as u64 })
        .collect::<Vec<ReverseIndex>>();
    to_sort.voracious_stable_sort();
    to_sort.iter().for_each(|x| index[x.bucket as usize].push(x.pos));

    println!("Generate Index: {:?}", now.elapsed().as_micros());

    // Instructions for indexes for the in place addition tree
    let mut instr: Vec<Vec<(usize, usize)>> = vec![];
    // Find the maximum depth of the addition tree
    let max_depth = index.iter()
        // log_2
        .map(|x| crate::log2(x.len()))
        .max().unwrap();

    let now = std::time::Instant::now();
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
    println!("Generate Instr: {:?}", now.elapsed().as_micros());

    let now = std::time::Instant::now();
    let mut elems_mut_1 = elems.to_vec();

    for instr_row in instr.iter() {
        for chunk in instr_row.chunks(BATCH_ADD_SIZE) {
            elems_mut_1[..].batch_add_in_place_same_slice(chunk.to_vec());
        }
    }
    println!("Batch add in place: {:?}", now.elapsed().as_micros());


    let now = std::time::Instant::now();
    let zero = C::zero();
    let mut res = vec![zero; buckets];


    for (i, to_add) in index.iter().enumerate() {
        if to_add.len() > 1 {
            panic!("Did not successfully reduce to_add");
        } else if to_add.len() == 1 {
            res[i] = elems_mut_1[to_add[0]];
        }
    }

    println!("Reassign: {:?}", now.elapsed().as_micros());
    res
}
