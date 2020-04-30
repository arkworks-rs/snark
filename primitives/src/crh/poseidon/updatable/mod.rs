use super::*;
use crate::UpdatableFieldBasedHash;
use std::marker::PhantomData;

pub mod batched;

pub struct UpdatablePoseidonHash<F: PrimeField + MulShort, P: PoseidonParameters<Fr = F>>{
    state: Vec<F>,
    pending: Vec<F>,
    _parameters: PhantomData<P>,
}

impl<F, P> UpdatablePoseidonHash<F, P>
    where
        F: PrimeField + MulShort ,
        P: PoseidonParameters<Fr = F>,
{
    pub fn new(personalization: Option<Vec<F>>) -> Self {
        let mut state = Vec::with_capacity(P::T);
        for i in 0..P::T {
            state.push(P::AFTER_ZERO_PERM[i]);
        }
        Self {
            state,
            pending: Vec::with_capacity(P::R),
            _parameters: PhantomData,
        }
    }

    fn apply_permutation(&mut self) {
        for (input, state) in self.pending.iter().zip(self.state.iter_mut()) {
            *state += input;
        }
        self.state[P::R] += &P::C2;
        PoseidonHash::<F, P>::poseidon_perm(&mut self.state);
    }

    fn _finalize(&self) -> F {
        let mut state = self.state.clone();
        for (input, s) in self.pending.iter().zip(state.iter_mut()) {
            *s += input;
        }
        state[P::R] += &P::C2;
        PoseidonHash::<F, P>::poseidon_perm(&mut state);
        state[0]
    }

}

impl<F, P> UpdatableFieldBasedHash for UpdatablePoseidonHash<F, P>
    where
        F: PrimeField + MulShort,
        P: PoseidonParameters<Fr = F>,
{
    type Data = F;
    type Parameters = P;

    fn update(&mut self, input: &Self::Data) -> &mut Self {
        self.pending.push(input.clone());
        if self.pending.len() == P::R {
            self.apply_permutation();
            self.pending.clear();
        }
        self
    }

    fn finalize(&self) -> Self::Data {
        if !self.pending.is_empty() {
            self._finalize()
        } else {
            self.state[0]
        }
    }
}

#[cfg(test)]
mod test {
    use algebra::{
        fields::{mnt4753::Fr as MNT4Fr, mnt6753::Fr as MNT6Fr},
        UniformRand,
    };
    use rand::SeedableRng;
    use rand_xorshift::XorShiftRng;
    use super::*;

    type UpdatableMNT4PoseidonHash = UpdatablePoseidonHash<MNT4Fr, MNT4753PoseidonParameters>;
    type UpdatableMNT6PoseidonHash = UpdatablePoseidonHash<MNT6Fr, MNT6753PoseidonParameters>;

    #[test]
    fn test_updatable_poseidon_hash_mnt4(){
        let samples = 100;
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        for i in 1..=samples {
            let input = vec![MNT4Fr::rand(&mut rng); i];
            let hash_output = MNT4PoseidonHash::evaluate(input.as_slice()).unwrap();
            let updatable_hash_output = {
                let mut updatable_hash = UpdatableMNT4PoseidonHash::new(None);
                for input in input.iter() {
                    updatable_hash.update(input);
                }
                updatable_hash.finalize()
            };
            assert_eq!(hash_output, updatable_hash_output, "Hashes output with {} inputs must be equal", i);
        }
    }

    #[test]
    fn test_updatable_poseidon_hash_mnt6(){
        let samples = 100;
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        for i in 1..=samples {
            let input = vec![MNT6Fr::rand(&mut rng); i];
            let hash_output = MNT6PoseidonHash::evaluate(input.as_slice()).unwrap();
            let updatable_hash_output = {
                let mut updatable_hash = UpdatableMNT6PoseidonHash::new(None);
                for input in input.iter() {
                    updatable_hash.update(input);
                }
                updatable_hash.finalize()
            };
            assert_eq!(hash_output, updatable_hash_output, "Hashes output with {} inputs must be equal", i);
        }
    }
}