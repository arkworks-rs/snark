use super::*;
use crate::{
    crh::{
        poseidon::batched_crh::PoseidonBatchHash,
        BatchFieldBasedHash
    },
    UpdatableBatchFieldBasedHash
};
use std::marker::PhantomData;

// Updatable Batch Poseidon Hash: each update() call will add a new batch for which computing
// the hash. The hashes for each batch are computed in parallel whenever the number of batch
// reaches cpu_loads size (a variable used to evenly distribute the work among the available
// cpus, given by hashes_per_core * number_of_cores, where hashes_per_core is a variable set by
// the user and defaults to 1. The user may set the value of hashes_per_core as desired depending
// on its system and use case), and the 'outputs' vector will be extended with the new results.
// The finalize() function will process the eventual remaining batches in 'pending' Vec and will
// return the result along with all the previous outputs. Since all the outputs are kept in memory
// for consistency with the specifications of BatchPoseidonHash, but they may become unneeded,
// the clear() function will free the memory from the output Vec.
pub struct UpdatableBatchPoseidonHash<F: PrimeField + MulShort, P: PoseidonParameters<Fr = F>>{
    outputs: Vec<F>,
    pending: Vec<F>,
    cpu_load: usize,
    _parameters: PhantomData<P>,
}

impl<F, P> UpdatableBatchPoseidonHash<F, P>
    where
        F: PrimeField + MulShort ,
        P: PoseidonParameters<Fr = F>,
{
    pub fn new(hashes_per_core: Option<usize>) -> Self {
        let cpu_load = rayon::current_num_threads() * hashes_per_core.unwrap_or(1);
        Self {
            outputs: Vec::with_capacity((cpu_load * P::R)/2),
            pending: Vec::with_capacity(cpu_load * P::R),
            cpu_load,
            _parameters: PhantomData,
        }
    }

    #[inline]
    fn apply_permutation(&mut self) {
        let output = PoseidonBatchHash::<F, P>::batch_evaluate(
            self.pending.as_slice()
        ).unwrap();
        self.outputs.extend_from_slice(output.as_slice());
    }

    #[inline]
    fn _finalize(&self) -> Vec<F> {
        let mut outputs = self.outputs.clone();
        let output = PoseidonBatchHash::<F, P>::batch_evaluate(
            self.pending.as_slice()
        ).unwrap();
        outputs.extend_from_slice(output.as_slice());
        outputs
    }
}

impl<F, P> UpdatableBatchFieldBasedHash for UpdatableBatchPoseidonHash<F, P>
    where
        F: PrimeField + MulShort,
        P: PoseidonParameters<Fr = F>,
{
    type Data = F;
    type Parameters = P;

    #[inline]
    fn update(&mut self, input: &[Self::Data]) -> &mut Self {
        assert_eq!(input.len(), P::R);
        self.pending.extend_from_slice(input);

        if self.pending.len() == self.cpu_load * P::R {
            self.apply_permutation();
            self.pending.clear();
        }
        self
    }

    fn finalize(&self) -> Vec<Self::Data> {
        if !self.pending.is_empty() {
            self._finalize()
        } else {
            self.outputs.clone()
        }
    }

    fn clear(&mut self) {
        self.outputs.clear();
    }
}

#[cfg(test)]
mod test {
    use algebra::{fields::{mnt4753::Fr as MNT4Fr, mnt6753::Fr as MNT6Fr}, UniformRand};
    use rand::SeedableRng;
    use rand_xorshift::XorShiftRng;
    use super::*;

    type MNT4PoseidonBatchHash = PoseidonBatchHash<MNT4Fr, MNT4753PoseidonParameters>;
    type MNT6PoseidonBatchHash = PoseidonBatchHash<MNT6Fr, MNT6753PoseidonParameters>;

    type UpdatableMNT4PoseidonBatchHash = UpdatableBatchPoseidonHash<MNT4Fr, MNT4753PoseidonParameters>;
    type UpdatableMNT6PoseidonBatchHash = UpdatableBatchPoseidonHash<MNT6Fr, MNT6753PoseidonParameters>;

    #[test]
    fn test_updatable_batch_poseidon_hash_mnt4(){
        let samples = 100;
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        // Test consistency with PoseidonHash and UpdatablePoseidonHash
        for i in 1..=samples {
            let input = vec![MNT4Fr::rand(&mut rng); MNT4753PoseidonParameters::R * i];
            let hash_output = MNT4PoseidonBatchHash::batch_evaluate(input.as_slice()).unwrap();
            let updatable_hash_output = {
                let mut updatable_hash = UpdatableMNT4PoseidonBatchHash::new(None);
                for input in input.chunks(MNT4753PoseidonParameters::R) {
                    updatable_hash.update(input);
                }
                updatable_hash.finalize()
            };
            assert_eq!(hash_output, updatable_hash_output, "Hashes output with {} inputs must be equal", i);
        }

        // Test finalize() holding the state and allowing updates in between different calls to it
        let input = vec![MNT4Fr::rand(&mut rng); 4];
        let h_out = MNT4PoseidonBatchHash::batch_evaluate(input.as_slice()).unwrap();

        let mut uh =  UpdatableMNT4PoseidonBatchHash::new(None);
        uh.update(&input[0..2]);
        uh.finalize();
        uh.update(&input[2..4]);
        assert_eq!(h_out, uh.finalize());

        //Test finalize() being idempotent
        assert_eq!(h_out, uh.finalize());
    }

    #[test]
    fn test_updatable_poseidon_hash_mnt6(){
        let samples = 100;
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        // Test consistency with PoseidonHash and UpdatablePoseidonHash
        for i in 1..=samples {
            let input = vec![MNT6Fr::rand(&mut rng); MNT6753PoseidonParameters::R * i];
            let hash_output = MNT6PoseidonBatchHash::batch_evaluate(input.as_slice()).unwrap();
            let updatable_hash_output = {
                let mut updatable_hash = UpdatableMNT6PoseidonBatchHash::new(None);
                for input in input.chunks(MNT6753PoseidonParameters::R) {
                    updatable_hash.update(input);
                }
                updatable_hash.finalize()
            };
            assert_eq!(hash_output, updatable_hash_output, "Hashes output with {} inputs must be equal", i);
        }

        // Test finalize() holding the state and allowing updates in between different calls to it
        let input = vec![MNT6Fr::rand(&mut rng); 4];
        let h_out = MNT6PoseidonBatchHash::batch_evaluate(input.as_slice()).unwrap();

        let mut uh =  UpdatableMNT6PoseidonBatchHash::new(None);
        uh.update(&input[0..2]);
        uh.finalize();
        uh.update(&input[2..4]);
        assert_eq!(h_out, uh.finalize());

        //Test finalize() being idempotent
        assert_eq!(h_out, uh.finalize());
    }
}