use super::*;
use std::marker::PhantomData;
use crate::UpdatableFieldBasedHash;

pub mod batched;
pub use batched::*;

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
        let mut instance = Self {
            state,
            pending: Vec::with_capacity(P::R),
            _parameters: PhantomData,
        };

        // If personalization Vec is not multiple of the rate, we pad it with zero field elements.
        // This will allow eventually to precompute the constants of the initial state. This
        // is exactly as doing H(personalization, padding, ...). NOTE: this way of personalizing
        // the hash is not mentioned in https://eprint.iacr.org/2019/458.pdf
        if personalization.is_some(){
            let personalization = personalization.unwrap();
            let padding = personalization.len() % P::R;

            for p in personalization.into_iter(){
                instance.update(p);
            }

            for _ in 0..padding {
                instance.update(F::zero());
            }
            assert_eq!(instance.pending.len(), 0);
        }
        instance
    }

    #[inline]
    fn apply_permutation(&mut self) {
        for (input, state) in self.pending.iter().zip(self.state.iter_mut()) {
            *state += input;
        }
        self.state[P::R] += &P::C2;
        PoseidonHash::<F, P>::poseidon_perm(&mut self.state);
    }

    #[inline]
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

    fn update(&mut self, input: Self::Data) -> &mut Self {
        self.pending.push(input);
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
    use algebra::{fields::{mnt4753::Fr as MNT4Fr, mnt6753::Fr as MNT6Fr}, UniformRand, Field};
    use rand::SeedableRng;
    use rand_xorshift::XorShiftRng;
    use super::*;

    type UpdatableMNT4PoseidonHash = UpdatablePoseidonHash<MNT4Fr, MNT4753PoseidonParameters>;
    type UpdatableMNT6PoseidonHash = UpdatablePoseidonHash<MNT6Fr, MNT6753PoseidonParameters>;

    #[test]
    fn test_updatable_poseidon_hash_mnt4(){
        let samples = 100;
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        // Test consistency with PoseidonHash and UpdatablePoseidonHash
        for i in 1..=samples {
            let input = vec![MNT4Fr::rand(&mut rng); i];
            let hash_output = MNT4PoseidonHash::evaluate(input.as_slice()).unwrap();
            let updatable_hash_output = {
                let mut updatable_hash = UpdatableMNT4PoseidonHash::new(None);
                for input in input.into_iter() {
                    updatable_hash.update(input);
                }
                updatable_hash.finalize()
            };
            assert_eq!(hash_output, updatable_hash_output, "Hashes output with {} inputs must be equal", i);
        }

        // Test finalize() holding the state and allowing updates in between different calls to it
        let input = vec![MNT4Fr::rand(&mut rng); 2];
        let h_out = MNT4PoseidonHash::evaluate(input.as_slice()).unwrap();

        let mut uh =  UpdatableMNT4PoseidonHash::new(None);
        uh.update(input[0]);
        uh.finalize();
        uh.update(input[1]);
        assert_eq!(h_out, uh.finalize());

        //Test finalize() being idempotent
        assert_eq!(h_out, uh.finalize());

        // Test initializing UpdatablePoseidonHash with personalization is the same as concatenating
        // to PoseidonHash input the personalization and the padding.
        let input = MNT4Fr::rand(&mut rng);
        let samples = 10;
        for i in 1..=samples{
            let personalization = vec![MNT4Fr::rand(&mut rng); i];

            let mut hash_input = personalization.clone();
            let padding = vec![MNT4Fr::zero(); personalization.len() % MNT4753PoseidonParameters::R];
            hash_input.extend_from_slice(padding.as_slice());
            hash_input.push(input.clone());
            let hash_output = MNT4PoseidonHash::evaluate(hash_input.as_slice()).unwrap();

            let updatable_hash_output = {
                let mut updatable_hash = UpdatableMNT4PoseidonHash::new(Some(personalization));
                updatable_hash.update(input);
                updatable_hash.finalize()
            };
            assert_eq!(hash_output,
                       updatable_hash_output,
                       "Hashes output with {} elements of personalization must be equal", i
            );
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
                for input in input.into_iter() {
                    updatable_hash.update(input);
                }
                updatable_hash.finalize()
            };
            assert_eq!(hash_output, updatable_hash_output, "Hashes output with {} inputs must be equal", i);
        }

        // Test finalize() holding the state and allowing updates in between different calls to it
        let input = vec![MNT6Fr::rand(&mut rng); 2];
        let h_out = MNT6PoseidonHash::evaluate(input.as_slice()).unwrap();

        let mut uh =  UpdatableMNT6PoseidonHash::new(None);
        uh.update(input[0]);
        uh.finalize();
        uh.update(input[1]);
        assert_eq!(h_out, uh.finalize());

        //Test finalize() being idempotent
        assert_eq!(h_out, uh.finalize());

        // Test initializing UpdatablePoseidonHash with personalization is the same as concatenating
        // to PoseidonHash input the personalization and the padding.
        let input = MNT6Fr::rand(&mut rng);
        let samples = 10;
        for i in 1..=samples{
            let personalization = vec![MNT6Fr::rand(&mut rng); i];

            let mut hash_input = personalization.clone();
            let padding = vec![MNT6Fr::zero(); personalization.len() % MNT6753PoseidonParameters::R];
            hash_input.extend_from_slice(padding.as_slice());
            hash_input.push(input.clone());
            let hash_output = MNT6PoseidonHash::evaluate(hash_input.as_slice()).unwrap();

            let updatable_hash_output = {
                let mut updatable_hash = UpdatableMNT6PoseidonHash::new(Some(personalization));
                updatable_hash.update(input);
                updatable_hash.finalize()
            };
            assert_eq!(hash_output,
                       updatable_hash_output,
                       "Hashes output with {} elements of personalization must be equal", i
            );
        }
    }
}