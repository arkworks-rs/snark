
pub mod bw6_761_g1_bucket_add_kernel {
    use accel::*;
    use rayon::prelude::*;

    use algebra::{BigInteger, FpParameters, Zero};
    use algebra_core::{curves::{ProjectiveCurve, AffineCurve}, fields::PrimeField};

    #[kernel_mod]
    pub mod batch_add_write {
        pub unsafe fn batch_add_write(

        )
    }

    pub fn batch_add_in_place_same_slice(

    )

    pub fn run_kernel(
        buckets: usize,
        elems: &[G1Affine],
        bucket_positions: &mut [BucketPosition],
    ) -> Vec<G1Affine> {
        run_kernel_inner::<G1Affine>(buckets, elems, bucket_positions)
    }

    pub fn run_kernel_inner<C: AffineCurve>(
        buckets: usize,
        elems: DeviceMemory<C>,
        bucket_positions: &mut [BucketPosition],
        ctx: &Context,
    ) -> Vec<C> {
        assert_eq!(elems.len(), bucket_positions.len());
        assert!(elems.len() > 0);

        const BATCH_SIZE: usize = (elems.len() - 1) / 16 + 1;

        let _now = timer!();
        dlsd_radixsort(bucket_positions, 8);
        timer_println!(_now, "radixsort");

        let mut len = bucket_positions.len();
        let mut all_ones = true;
        let mut new_len = 0; // len counter
        let mut glob = 0; // global counters
        let mut loc = 1; // local counter
        let mut batch = 0; // batch counter
        let mut instr = DeviceMemory::<(u32, u32)>::zeros(BATCH_SIZE + 1024);
        let mut new_elems = Vec::<C>::with_capacity(elems.len() * 3 / 8);

        let mut scratch_space = Vec::<Option<C>>::with_capacity(BATCH_SIZE / 2);

        let _now = timer!();
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
            } else if loc > 1 {
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
                    batch_add_write_kernel::batch_add_write(&elems[..], &instr[..], &mut new_elems, &mut scratch_space);

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
            batch_add_write_kernel::batch_add_write(&elems[..], &instr[..], &mut new_elems, &mut scratch_space);
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
                        batch_add_in_place_same_slice_kernel::batch_add_in_place_same_slice(
                            &mut new_elems[..],
                            &instr[..]
                        );
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
                batch_add_in_place_same_slice_kernel::batch_add_in_place_same_slice(
                    &mut new_elems[..],
                    &instr[..]
                );
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
}
