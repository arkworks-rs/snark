extern crate rand;
extern crate rayon;

use algebra::{PrimeField, MulShort};

use std::marker::PhantomData;

use crate::{
    crh::{
        FieldBasedHash,
        FieldBasedHashParameters,
    }, CryptoError, Error
};

pub mod batched_crh;

pub mod parameters;
pub use self::parameters::*;

#[derive(Derivative)]
#[derivative(
    Clone(bound = ""),
    Debug(bound = ""),
)]
pub struct PoseidonHash<F: PrimeField, P: PoseidonParameters<Fr = F>>{
    state: Vec<F>,
    pending: Vec<F>,
    input_size: Option<usize>,
    updates_ctr: usize,
    mod_rate: bool,
    _parameters: PhantomData<P>,
}

pub trait PoseidonParameters: 'static + FieldBasedHashParameters + Clone {

    const T: usize;  // Number of S-Boxes
    const R_F:i32;   // Number of full rounds
    const R_P:i32;   // Number of partial rounds
    const ZERO:Self::Fr;   // The zero element in the field
    const AFTER_ZERO_PERM: &'static[Self::Fr]; // State vector after a zero permutation
    const ROUND_CST: &'static[Self::Fr];  // Array of round constants
    const MDS_CST: &'static[Self::Fr];  // The MDS matrix
    const MDS_CST_SHORT: &'static[Self::Fr];  // The MDS matrix for fast matrix multiplication

}

// Function that does the scalar multiplication
// It uses Montgomery multiplication
// Constants are defined such that the result is x * t * 2^n mod M,
// that is the Montgomery representation of the operand x * t mod M, and t is the 64-bit constant
#[allow(dead_code)]
#[inline]
pub fn scalar_mul<F: PrimeField + MulShort<F, Output = F>, P: PoseidonParameters<Fr=F>> (res: &mut F, state: &mut[F], mut start_idx_cst: usize) {

    state.iter().for_each(|x| {
        let elem = x.mul(&P::MDS_CST[start_idx_cst]);
        start_idx_cst += 1;
        *res += &elem;
    });
}

// Function that does the mix matrix
#[allow(dead_code)]
#[inline]
pub fn matrix_mix<F: PrimeField + MulShort<F, Output = F>, P: PoseidonParameters<Fr=F>>  (state: &mut Vec<F>) {

    // the new state where the result will be stored initialized to zero elements
    let mut new_state = vec![F::zero(); P::T];

    let mut idx_cst = 0;
    for i in 0..P::T {
        scalar_mul::<F,P>(&mut new_state[i], state, idx_cst);
        idx_cst += P::T;
    }
    *state = new_state;
}

// Function that does the scalar multiplication
// It uses a partial Montgomery multiplication defined as PM(x, t) = x * t * 2^-64 mod M
// t is a 64-bit matrix constant. In the algorithm, the constants are represented in
// partial Montgomery representation, i.e. t * 2^64 mod M
#[inline]
pub fn scalar_mul_fast<F: PrimeField + MulShort<F, Output = F>, P: PoseidonParameters<Fr=F>> (res: &mut F, state: &mut[F], mut start_idx_cst: usize) {
    state.iter().for_each(|&x| {
        let elem = P::MDS_CST_SHORT[start_idx_cst].mul_short(x);
        start_idx_cst += 1;
        *res += &elem;
    });
}

// Function that does the mix matrix with fast algorithm
#[inline]
pub fn matrix_mix_short<F: PrimeField + MulShort<F, Output = F>, P: PoseidonParameters<Fr=F>> (state: &mut Vec<F>) {

    // the new state where the result will be stored initialized to zero elements
    let mut new_state = vec![F::zero(); P::T];

    let mut idx_cst = 0;
    for i in 0..P::T {
        scalar_mul_fast::<F,P>(&mut new_state[i], state, idx_cst);
        idx_cst += P::T;
    }
    *state = new_state;
}

impl<F: PrimeField + MulShort<F, Output = F>, P: PoseidonParameters<Fr=F>> PoseidonHash<F, P> {

    fn _init(constant_size: Option<usize>, mod_rate: bool, personalization: Option<&[F]>) -> Self
    {
        let mut state = Vec::with_capacity(P::T);
        for i in 0..P::T {
            state.push(P::AFTER_ZERO_PERM[i]);
        }


        let mut instance = Self {
            state,
            pending: Vec::with_capacity(P::R),
            input_size: constant_size,
            updates_ctr: 0,
            mod_rate,
            _parameters: PhantomData,
        };

        // If personalization Vec is not multiple of the rate, we pad it with zero field elements.
        // This will allow eventually to precompute the constants of the initial state. This
        // is exactly as doing H(personalization, padding, ...). NOTE: this way of personalizing
        // the hash is not mentioned in https://eprint.iacr.org/2019/458.pdf
        if personalization.is_some(){
            // Use a support variable-length non mod rate instance
            let mut personalization_instance = Self::init_variable_length(false, None);
            let personalization = personalization.unwrap();

            // Apply personalization
            for &p in personalization.into_iter(){
                personalization_instance.update(p);
            }

            // Apply padding (according to the variable length input strategy)
            personalization_instance.update(F::one());
            for _ in personalization_instance.pending.len()..P::R {
                personalization_instance.update(F::zero());
            }
            assert_eq!(personalization_instance.pending.len(), 0);

            // Set the new initial state
            instance.state = personalization_instance.state;
        }
        instance
    }

    #[inline]
    fn apply_permutation_if_needed(&mut self) {
        if self.pending.len() == P::R {
            for (input, state) in self.pending.iter().zip(self.state.iter_mut()) {
                *state += input;
            }
            Self::poseidon_perm(&mut self.state);
            self.pending.clear();
        }
    }

    #[inline]
    fn get_hash(mut state: Vec<F>, inputs: Vec<F>) -> F {
        for (input, s) in inputs.iter().zip(state.iter_mut()) {
            *s += input;
        }
        Self::poseidon_perm(&mut state);
        state[0]
    }

    #[inline]
    // Padding strategy is described in https://eprint.iacr.org/2019/458.pdf(Section 4.2)
    fn pad_and_finalize(&self) -> F
    {
        // Constant input length instance
        if self.input_size.is_some() {

            // The constant size is modulus rate, so we already have the hash output in state[0]
            // as a permutation is applied each time pending reaches P::R length
            if self.pending.is_empty() {
                self.state[0].clone()
            }

            // Pending is not empty: pad with 0s up to rate then compute the hash,
            else {
                Self::get_hash(self.state.clone(), self.pending.clone())
            }
        }

        // Variable input length instance
        else {

            // The input is of variable length, but always modulus rate: result is already available
            // in state[0] as a permutation is applied each time pending reaches P::R length
            if self.mod_rate {
                self.state[0].clone()
            }

            // The input is of variable length, but not modulus rate: we always need to apply
            // padding. Pad with a single 1 and then 0s up to rate. Compute hash.
            else {
                let mut pending = self.pending.clone(); // Can also be empty if the input happens to be mod rate
                pending.push(F::one());
                Self::get_hash(self.state.clone(), pending)
            }
        }
    }

    /// Inversion S-Box only !
    pub(crate) fn poseidon_perm (state: &mut Vec<F>) {

        // index that goes over the round constants
        let mut round_cst_idx = 0;

        // First full rounds
        for _i in 0..P::R_F {

            // Add the round constants to the state vector
            for d in state.iter_mut() {
                let rc = P::ROUND_CST[round_cst_idx];
                *d += &rc;
                round_cst_idx += 1;
            }

            // Apply the S-BOX to each of the elements of the state vector
            // Use Montgomery simultaneous inversion
            let mut w: Vec<P::Fr> = Vec::new();
            let mut accum_prod = P::Fr::one();

            w.push(accum_prod);

            // Calculate the intermediate partial products
            for d in state.iter() {
                accum_prod = accum_prod * &d;
                w.push(accum_prod);
            }

            if accum_prod == P::Fr::zero() {
                // At least one of the S-Boxes is zero
                // Calculate inverses individually
                for d in state.iter_mut() {
                    // The S-BOX is an inversion function
                    if *d != P::Fr::zero() {
                        *d = (*d).inverse().unwrap();
                    }
                }
            } else {
                let mut w_bar = accum_prod.inverse().unwrap();

                // Extract the individual inversions
                let mut idx: i64 = w.len() as i64 - P::R as i64;
                for d in state.iter_mut().rev() {
                    let tmp = d.clone();
                    *d = w_bar * &w[idx as usize];
                    w_bar = w_bar * &tmp;
                    idx -= 1;
                }
            }

            // Perform the matrix mix
            matrix_mix_short::<F, P>(state);

        }

        // Partial rounds
        for _i in 0..P::R_P {

            // Add the round constants to the state vector
            for d in state.iter_mut() {
                let rc = P::ROUND_CST[round_cst_idx];
                *d += &rc;
                round_cst_idx += 1;
            }

            // Apply S-BOX only to the first element of the state vector
            if state[0]!= P::Fr::zero() {
                state[0] = state[0].inverse().unwrap();
            }

            // Apply the matrix mix
            matrix_mix_short::<F, P>(state);
        }

        // Second full rounds
        for _i in 0..P::R_F {

            // Add the round constants
            for d in state.iter_mut() {
                let rc = P::ROUND_CST[round_cst_idx];
                *d += &rc;
                round_cst_idx += 1;
            }

            // Apply the S-BOX to each of the elements of the state vector
            // Use Montgomery simultaneous inversion
            let mut w: Vec<P::Fr> = Vec::new();
            let mut accum_prod = P::Fr::one();

            w.push(accum_prod);

            // Calculate the intermediate partial products
            for d in state.iter() {
                accum_prod = accum_prod * &d;
                w.push(accum_prod);
            }

            if accum_prod == P::Fr::zero() {
                // At least one of the S-Boxes is zero
                // Calculate inverses individually
                for d in state.iter_mut() {
                    // The S-BOX is an inversion function
                    if *d != P::Fr::zero() {
                        *d = (*d).inverse().unwrap();
                    }
                }
            } else {
                let mut w_bar = accum_prod.inverse().unwrap();

                // Extract the individual inversions
                let mut idx: i64 = w.len() as i64 - P::R as i64;
                for d in state.iter_mut().rev() {
                    let tmp = d.clone();
                    *d = w_bar * &w[idx as usize];
                    w_bar = w_bar * &tmp;
                    idx -= 1;
                }
            }

            // Perform the matrix mix
            matrix_mix_short::<F, P>(state);
        }
    }
}

impl<F, P> FieldBasedHash for PoseidonHash<F, P>
    where
        F: PrimeField + MulShort<F, Output = F>,
        P: PoseidonParameters<Fr = F>,
{
    type Data = F;
    type Parameters = P;

    fn init_constant_length(input_size: usize, personalization: Option<&[Self::Data]>) -> Self {
        Self::_init(
            Some(input_size),
            input_size % P::R == 0, // Not taken into account, can be any
            personalization
        )
    }

    fn init_variable_length(mod_rate: bool, personalization: Option<&[Self::Data]>) -> Self {
        Self::_init(None, mod_rate, personalization)
    }

    fn update(&mut self, input: Self::Data) -> &mut Self {
        self.pending.push(input);
        self.updates_ctr += 1;
        self.apply_permutation_if_needed();
        self
    }

    fn finalize(&self) -> Result<Self::Data, Error> {

        let error_condition =

            // Constant input length instance, but the size of the input is different from the declared one
            (self.input_size.is_some() && self.updates_ctr != self.input_size.unwrap())

            ||

            // Variable modulus rate input length instance, but the size of the input is not modulus rate
            (self.input_size.is_none() && self.mod_rate && self.updates_ctr % P::R != 0);

        // If one of the conditions above is true, we must throw an error
        if error_condition
        {
            Err(Box::new(CryptoError::HashingError("attempt to finalize with an input of invalid size".to_owned())))
        }

        // Otherwise pad if needed (according to the Self instance type) and return the hash output
        else
        {
            Ok(self.pad_and_finalize())
        }
    }

    fn reset(&mut self, personalization: Option<&[Self::Data]>) -> &mut Self {
        let new_instance = Self::_init(self.input_size, self.mod_rate, personalization);
        *self = new_instance;
        self
    }
}

#[cfg(test)]
mod test {
    use algebra::{Field, PrimeField, MulShort};
    use crate::crh::{
        FieldBasedHash,
        test::{
            constant_length_field_based_hash_test, variable_length_field_based_hash_test
        }
    };
    use crate::{FieldBasedHashParameters, PoseidonParameters, PoseidonHash};

    fn generate_inputs<F: PrimeField>(num: usize) -> Vec<F>{
        let mut inputs = Vec::with_capacity(num);
        for i in 1..=num {
            let input = F::from(i as u32);
            inputs.push(input);
        }
        inputs
    }

    fn poseidon_permutation_regression_test<F: PrimeField + MulShort<F, Output = F>, P: PoseidonParameters<Fr = F>>(
        start_states: Vec<Vec<F>>,
        end_states:   Vec<Vec<F>>,
    )
    {
        // Regression test
        start_states.into_iter().zip(end_states).enumerate().for_each(|(i, (mut start_state, end_state))| {
            PoseidonHash::<F, P>::poseidon_perm(&mut start_state);
            assert_eq!(
                start_state,
                end_state,
                "Incorrect end state {}:\n Expected\n{:?}\n, Found\n {:?}\n", i, start_state, end_state);
        });
    }

    fn test_routine<F: PrimeField, H: FieldBasedHash<Data = F>>(
        num_samples:  usize,
    )
    {
        let rate = <H::Parameters as FieldBasedHashParameters>::R;
        for i in 0..num_samples {

            let ins = generate_inputs::<F>(i + 1);

            // Constant length
            {
                let mut digest = H::init_constant_length(i + 1, None);

                constant_length_field_based_hash_test::<H>(
                    &mut digest,
                    ins.clone()
                );
            }

            // Variable length
            {
                let mod_rate = (i + 1) % rate == 0;
                let mut digest = H::init_variable_length(mod_rate, None);

                variable_length_field_based_hash_test::<H>(
                    &mut digest,
                    ins.clone(),
                    mod_rate
                );

                // Test also case in which mod_rate is false but the input happens to be mod rate
                if mod_rate {
                    let mut digest = H::init_variable_length(!mod_rate, None);
                    variable_length_field_based_hash_test::<H>(
                        &mut digest,
                        ins,
                        !mod_rate
                    );
                }
            }
        }
    }

    #[cfg(feature = "mnt4_753")]
    #[test]
    fn test_poseidon_hash_mnt4() {
        use algebra::{
            biginteger::BigInteger768,
            fields::mnt4753::Fr as MNT4753Fr
        };
        use crate::crh::poseidon::parameters::mnt4753::{
            MNT4PoseidonHash, MNT4753PoseidonParameters
        };

        // Test vectors are computed via the script in ./parameters/scripts/permutation_mnt4fr.sage
        let start_states = vec![
            vec![MNT4753Fr::zero(); 3],
            vec![
                MNT4753Fr::new(BigInteger768([0xf770348fbe4e29b6,0xfefd6b30dfb52494,0xec61827e5cf9425,0xc6288db72079112c,0xd70e11f75c351bac,0x2e4657caf8648c8e,0x7f9f3a94358aa2f7,0xee7f886bb42e8eab,0xe5ae5d4ec1b0796f,0xd056464cb38777c6,0xf3d7cd676c74ae38,0x120d49a741c34,])),
                MNT4753Fr::new(BigInteger768([0x96de60f9741f78b7,0xa98cc9495bb4615e,0xc4b3aeadfd321c2c,0x40e4b75eb8fe1116,0x1396ee290297e819,0x9762744e4cfded19,0xbedcef99b43ee15a,0x8b84865c31d378a0,0xf5468754aa4a4c4e,0xfd715c8245c2e124,0x31cb5bb04a339986,0xdaf306180aed,])),
                MNT4753Fr::new(BigInteger768([0x7e874134d509e406,0x729d013268020212,0x8b362dd530097799,0xae5054da3ad04250,0xe2e7413bd0fcbe5f,0xad08673f2f925bee,0xfb93f0ee8900d97e,0x2c1d037343b00151,0xd3dac3f2b1139f55,0x154e788ae1aca4cc,0x663269814fb52d57,0x676d9c4d8329,])),
            ],
            vec![
                MNT4753Fr::new(BigInteger768([0xa26b0bc72724d615,0x729202dca25403d4,0x1b2ff6dc78c46b5e,0xed529329c88557ec,0xa7264c3cd1f1ca2d,0xa9f0e2b1e57c800f,0x2322b96082d360ec,0x138d00037c082f1c,0x6c25792c21edce0a,0x75723fc00d8d1bc3,0xf60868fea31de240,0x14e224d41e354,])),
                MNT4753Fr::new(BigInteger768([0x21c229d68cde6f3f,0xf96b852ba3677e55,0x815b51e9b5e329c2,0xedec4ec2b77a9d36,0x44e0217411a0dea9,0x724a35de8cbd3141,0x8008cb5f0106c484,0x921855777c1c9cd3,0xd87d5b5babb7c9ab,0x603fc082a06ed8c4,0xe589b5a1adea946e,0x129d1f84a0c66,])),
                MNT4753Fr::new(BigInteger768([0x80794339ccdf973f,0x8f537759fc1b1aca,0x7997a170b362d649,0x7b1cddf6db6ca199,0x6b25316a81753330,0xa143d6d50bd07ebf,0x4d65e4fd6f8587d6,0x572c858cf606bd90,0x245465ba33e044b1,0x86f9aaa423b9390,0x8ee2bbed6bda13a6,0x7fa83fcd7a59,])),
            ],
            vec![
                MNT4753Fr::new(BigInteger768([0x275345cd3949fba9,0xaa492ccf37b80d9,0xdd9c6b17371c879a,0x846303d5b851d739,0x8d2b1b900c8c2227,0x780824b721514171,0xe08b4ffffb8a4f71,0xc69a0eb1b3f3ad,0x409578a5de88b1df,0xef2b552006465afb,0x2539560ecdf8147,0x134fe3e183dcd,])),
                MNT4753Fr::new(BigInteger768([0xf7f3c59f70e5b72a,0xec1ae7ed077f2d99,0xbbf075b432e1a2d8,0xf32012c620b8cd09,0x81e964a2687b8654,0x43082373cc23c4f6,0x494428fd5d2b9d5,0xed89d49a5f32ca1a,0x8d2c7f6937d4bc08,0x8aa8316d21567c0c,0x5e2c9cde56f4c802,0x6422f65bc889,])),
                MNT4753Fr::new(BigInteger768([0x44238a7e541cdf0,0xc09a1bda2e310a6d,0xef2001005bbaf873,0x1fd97ee19fea97eb,0xce43458dee7839cd,0x735d8cff80565348,0xca740dd90f883e06,0x8825f23c63c39a44,0xe80c50eb3548e408,0xddc815aae7e6a432,0x519048208b84f07f,0x50d352305dca,])),
            ],
            vec![
                MNT4753Fr::new(BigInteger768([0x911b5559a3eeb52d,0x482afb0b1b566e49,0x3983c4efc4fb37da,0x3288b81e77372d01,0xc69bd18751793c34,0x103f732ca150f840,0xbe72b866f7fd8512,0x19f4e9f908c9d1bf,0xb7976427cfc0fe4e,0xc9f43b7c2ad54601,0x3f2eb373787a291,0x9d3dd62a7475,])),
                MNT4753Fr::new(BigInteger768([0x799693496d2180d4,0x9c8364f338a500b7,0x37a57ca5674e1252,0x2c19b0502325bead,0x32b30a126f41f5ac,0x8bcd51ff52cedf29,0x9e04cb66d8d16160,0x59e8aaadbc99fab6,0xbd046f342e99d386,0x4488dd3ce29590aa,0xdcc2bb0149b02eaa,0x1543d162aa244,])),
                MNT4753Fr::new(BigInteger768([0xbb41e5acd82643f9,0x4042aec0d83f7624,0x2c14ed2f563bb21e,0x9cee7ec494eb57e9,0x41eec6c2b0056ac2,0xd1ea7cfa30f223ef,0xf148c377c2fba415,0xb3b56ee96972c9cb,0x82c3e44086911217,0x9ef750feb5842cc6,0x9f33c28feb810dc0,0x727b9f80e6df,])),
            ],
        ];

        let end_states = vec![
            vec![
                MNT4753Fr::new(BigInteger768([0x4f54c026da6ed8f0,0x12700bf5ad94f6c9,0x23a3fa62e9c042c1,0x2394c785581c75e7,0x839626f16bd60d08,0xb29828eef68c9bd4,0xd1479004b0f71d2,0x9d1a0dffdd1e7b00,0x9f1df2af9215e68c,0xc562186972253d2e,0xf6b8c66a6f3999b0,0xa040e4e0ff92,])),
                MNT4753Fr::new(BigInteger768([0xb0258a782c08064,0x6a04841f8be4990a,0xda027778a67d713b,0xb88be63be3cac9b4,0xde929c2510a321e5,0xc0d9dd704886213e,0xfbe0efc728d44f11,0x77c8d6422b5eb368,0x2827d5af4fe0fbad,0xb90c8793bc2a9d21,0xf9ce1fdde5140214,0x15a64a6345311,])),
                MNT4753Fr::new(BigInteger768([0xde9731dd4ad29db3,0x86caaccf88b402a1,0xe5e77eee08fca8a2,0x1dd9e752e50aad07,0x2d0f73cfb9508a83,0xb2b6ab08f14d96eb,0x224833c17d87490d,0x4e7d2e13141aaa55,0x1796b61e1cc3563,0xdbeb6f5ed60179f,0xb0633f07c680eda2,0x601b999b7143,])),
            ],
            vec![
                MNT4753Fr::new(BigInteger768([0xe749d7517ebe099b,0xc6abeacc602cf0bf,0x958f4b91f3c3b22d,0x9a295b36c4a6ea9e,0xd3925085d5ae2179,0xf23a8b4284968652,0x8018232a8a8fd30b,0x34533842150d4c6a,0xf0531c8f2f4a3dd4,0xeaab2b7956c6e7cb,0x9fc2b52eb516b457,0x7e2c759defce,])),
                MNT4753Fr::new(BigInteger768([0xfc5dab1dedb49656,0x78deb85913893c98,0x6088942fdbff357e,0xb3c15f514de46072,0x5dc205c3ccd4df39,0x591d9320bec689a6,0x99a7765caae47a86,0x2fcfe60a560fa3ed,0x43e2f302b5852456,0x5b4087eaa01f39c6,0xcc7db3f671985b7d,0x1272366ae322b,])),
                MNT4753Fr::new(BigInteger768([0xc23a10d72a73058e,0x7125f89599d62e8e,0x944ffd3948d3b453,0xc1513ee7ef29c1d2,0xdf1ddf8a25a2233,0x193c0cac56b49055,0xcb23ffde25ea2bd6,0x6d4a4ad2f3e415af,0x7da1b50b3731057,0x30f2f41a6746bd09,0x2a3cfda1f9885424,0xe6f1af34a223,])),
            ],
            vec![
                MNT4753Fr::new(BigInteger768([0xbfcb18d74e65c563,0x722359395bfeb077,0xb8e0b7abddb9a694,0xc830a386c2854b6b,0x53d7c0704e145ce,0xbe91d2a17d6f8874,0x2b49e38e1b99292a,0xc7e2cb48c2be1151,0xa5e54b3a714aad54,0xf634e385fe3d9b90,0x66f9a11a59535867,0x1425351d064a2,])),
                MNT4753Fr::new(BigInteger768([0x4a28ff3c4fecbb8d,0x60a639f0a2a002d9,0x5149d27ed99128c1,0x6dacfe4ce235b503,0xf21ef2fe6f344e69,0xbac70a5d64a033de,0x54f1cb89e291c8e6,0x2548230a2b8eeb67,0x763440a89ffdc8de,0x3ac6435a7c2b7922,0xacb97881f998663d,0x8ae31b1e760f,])),
                MNT4753Fr::new(BigInteger768([0x9dfe82b5a7baefa5,0x14bff3144e3c4f00,0xcbb47c1db66e74c4,0x8c3d330245b24464,0x3be7110fcc0f2674,0xb4a9281c6d349356,0xa4894a010cef488c,0x2abe0a21b8a83ca7,0xf9e9d807e418b54,0x439e4046be879838,0x3204e13287f737d5,0x3098a5738444,])),
            ],
            vec![
                MNT4753Fr::new(BigInteger768([0x470bac44ae262597,0x37c75eb3f00758fb,0xae77bbd563b5fac6,0xa22469cb36563eb5,0x4db9a5ea229af500,0xf6848cf2a64ad4a5,0x3a4611a0ed9e6243,0xf63fb5b6489325dd,0x1a9c90dd1544863f,0xdab1cb220fdf73d4,0xb9ec40309591932b,0x141777a73c602,])),
                MNT4753Fr::new(BigInteger768([0xedab7a7bd3a0061b,0x32d0ba278e569bec,0x83a9e0f317060812,0x29acd35e4d33cdb6,0x3f13496b623a9cde,0xa565606e05e4a5d,0xba87579c189af741,0x45bcb5fbad648a4e,0x32e8658135401638,0xbc853abb54e732b5,0xc37855ec443e12d3,0x1ad1ff8f54ad6,])),
                MNT4753Fr::new(BigInteger768([0xaba94817dccf0311,0x601cdff2f1e54d9e,0x6a0d8ab8a097a5b6,0x51d8c83d12239512,0x92f9ef537fc921e8,0x688b9fe86605c3ae,0x250ebdd755ad043c,0x29d412ee38a1e765,0xb31f5447678264b4,0x7d053f0ea44d854b,0x5d83d881795db690,0x397b9db5b588,])),
            ],
            vec![
                MNT4753Fr::new(BigInteger768([0xf0afca787979dcae,0x42fbae09a94724f3,0xce13b6f47a98712e,0x68faa457e317c516,0x7f77afa6123189da,0xf24b93d153626436,0xa40c88d389b68cfd,0x9b032ff8170c5c10,0xb90fa1c19b5affe3,0xc6cb43fb1342f46b,0x73a8195215425b8a,0x16cfda5a32fef,])),
                MNT4753Fr::new(BigInteger768([0xd864f5bc7dbdbe12,0xd316f0a8460332b6,0xada86ced0ff99e99,0x80860702b69fbf79,0xe4a85e8c6fe21f02,0xdc253a82c99e4359,0x538ca29cb25f1740,0xb4b3b0c1728477d2,0x2ae092fa5a67319a,0xf11e69b6ea6e795b,0xbd153a2d52cd7fe1,0x172ce347450d4,])),
                MNT4753Fr::new(BigInteger768([0x16d7536835c3972f,0x6e1897915f2ecc3e,0xa12771652da6c8b8,0xaf97a5aaa35b7313,0xae2a361cddc23c31,0xefc41bde8666d6dc,0x6cdd6c01057a661,0x7235dca1f39f8bc6,0x6332b45ab259d,0x851fb01167d8a74a,0x1c840faa9ad5c9b7,0xfe4f5c82b740,])),
            ]
        ];

        poseidon_permutation_regression_test::<MNT4753Fr, MNT4753PoseidonParameters>(start_states, end_states);
        test_routine::<MNT4753Fr, MNT4PoseidonHash>(3)
    }

    #[cfg(feature = "mnt6_753")]
    #[test]
    fn test_poseidon_hash_mnt6() {
        use algebra::{
            biginteger::BigInteger768,
            fields::mnt6753::Fr as MNT6753Fr
        };
        use crate::crh::poseidon::parameters::mnt6753::{
            MNT6PoseidonHash, MNT6753PoseidonParameters,
        };

        // Test vectors are computed via the script in ./parameters/scripts/permutation_mnt6fr.sage
        let start_states = vec![
            vec![MNT6753Fr::zero(); 3],
            vec![
                MNT6753Fr::new(BigInteger768([0x2045f548c283a386,0x9f90d7b623ef9965,0x634e1e0bcd6ce5f1,0xed09fb1cd92e2f48,0xa4b92193ab3a4c,0xc38d823f5e556d81,0x93e8a09384f1d5f0,0xa463757a137a2127,0xc948555766dabe44,0x3246e78f29a70bfe,0x21ebc006f85e213,0x18c2c2170055e,])),
                MNT6753Fr::new(BigInteger768([0x5abb4a33f5026781,0xa3510b40fb1bd1e7,0xce8ae77f3e0e9a1d,0xd1375569096b196a,0x107156721a5241bd,0x82b75d6eb65ccdc,0x9f6a6933bbe7d8ad,0x9335a61a85fe8998,0x5179ec766656404c,0x8052414d46077e77,0xb77841abce4c69c,0x10e71d39ef7ee,])),
                MNT6753Fr::new(BigInteger768([0xf76a1c08fa236882,0x25e1b757eb33ed43,0x1f63d4997a13c8b1,0xe23eae7ea2605b4b,0xe8c20feb190f9dd,0xa63856368a5c24f9,0x114eaf0c94cc670b,0xe858d17f6da22272,0x9b5443cadda8156a,0xfe92bd2a3eefc8b3,0x2c8a4defc4a4ff9,0x19cc15d056674,])),
            ],
            vec![
                MNT6753Fr::new(BigInteger768([0x1b940903c57e8e7f,0xbd38cde2e8e16008,0xe18d1abcfe05990a,0x8e86b1ca3a0ee1f5,0x33a31929417f05f9,0x170be227265f62bd,0x29e22c2b9864352a,0x901db3c41b27245e,0xc3bc6e6cfce69e3c,0x498f01eea65c0215,0xbf86a87e3005b3db,0x90f488bd8e09,])),
                MNT6753Fr::new(BigInteger768([0xb2d9ad48cbb812ba,0xc53cb754a7a02d89,0x89f52c6630ad8f86,0xe623c68f3610652f,0x198f83682c814e5d,0xfb78854e850e95fb,0x46e398cb56c27f78,0x81e60dab3991f035,0x3babbc1fe35f4f30,0x8056c683be44ffab,0x167af8aceb070f00,0x1a2572baaf46d,])),
                MNT6753Fr::new(BigInteger768([0x6242acf3bfbe2c6e,0x7afcb4878b2fcab1,0xccdee01e7839e6ff,0x8ebef555a3fcaeb9,0xa627b970cb4d56d2,0xb672bd365dab0d61,0x71f74eef13dab0fd,0x5a138a0bd718f4c3,0x7d08a2cf2ef0747c,0x8a0cdeefcdfded66,0xfe18f6573bbabadb,0x12c02029e0030,])),
            ],
            vec![
                MNT6753Fr::new(BigInteger768([0xf2ca60b5bb454d9f,0xb4ae3ba59e4a711,0x62154368b888061c,0x6214f711b35b4f9,0x5dd4d44dc9d4f0ad,0x4304e1c271f64602,0x80d4e3b0e1025ae3,0x5316732d6accc44d,0x24fc5d7d7bba464e,0x12d10c9485d208a1,0xca6df371c62a8872,0x86ce9f608bae,])),
                MNT6753Fr::new(BigInteger768([0xcdf0f7492613b504,0x455faa0e541fa1e6,0xb77242df6b8a68be,0x3b5435160d723cb6,0x77b8914a813586bf,0xc17dabd68e491d26,0xa85720ce2231af9d,0xd19e81cea5d64c41,0x56c90bfdb43ce182,0x9ff4ff3aba6a9a01,0x8875906afee26342,0x16a993a8df862,])),
                MNT6753Fr::new(BigInteger768([0xad98e2452d8be558,0xed19ce15ee0069d3,0xf889b49a8ad1016e,0x42760a3cbfb291b7,0x3d94e422b333dc5d,0xc27cbbac2884c097,0x851fd495c84543e9,0xf9b100c34675f211,0x11eae122f8ff1706,0xf3eecc4f60743020,0x38fc6ca1e5d1b4a7,0xffa8124e7034,])),
            ],
            vec![
                MNT6753Fr::new(BigInteger768([0x376743561f755f41,0xf0a8830457e9879b,0xa134b300b8f2d67b,0x1806324852aa9feb,0xdb98705dbf283859,0x565bca638d85ee57,0x1c6b6fe1fe752b0,0xd9eb23176d5d6110,0x5c5e86b5774422e2,0xd6fdf4c92ea236a1,0xeb2a915f44b72fa3,0x195c5f80dbf29,])),
                MNT6753Fr::new(BigInteger768([0x4c674244dfb68ecc,0x24a104856852ac3f,0x83e10d6c10dd3f4f,0xe99fe1f0d8553d3c,0x2d371b923253d5c0,0x14594932de89a19e,0xfd4589d2f8e53f17,0xe2ba2c7b929a53b3,0x3891f35b974a36ec,0xf17f8749ca140c09,0x6be74c21301f7c9e,0x13de4e1311a04,])),
                MNT6753Fr::new(BigInteger768([0xc366ce203caca4b0,0xe1d195b5bf3af54e,0x24b93c34bd0043ee,0x91559c070b29c53a,0xe866e46830168ff8,0xaeeda2129518cab7,0x37f8bb28ae15d7f3,0x5811fb22acd02c55,0xce7d805057f58acc,0x3a80df0b2af5f4fd,0x4dc7c29c8f6bed72,0xe511723afdb9,])),
            ]
        ];

        let end_states = vec![
            vec![
                MNT6753Fr::new(BigInteger768([0xef99f18ca1164fb0,0x1bf161755d689806,0x83ee017c500c6964,0x8abab822f92200c0,0x4b64884b9cc7eef9,0x53d4a2f13e17017c,0x551b8da2668dad8a,0x9939a48a0191c96c,0x2e1d80ef403671a0,0xb037bb60fbeb0212,0x6a22eba60581eb12,0x6ec196c9026d,])),
                MNT6753Fr::new(BigInteger768([0x18c4207483ba0f2f,0x6c50abc8aca74de3,0x7c1acfd6686351c,0xf367937c1356e91f,0xcdbf0447592ec1,0xe13763baac982387,0x2e1f904290e7045f,0xb6ffbcccd73c1092,0xfae22550de44cf2c,0x14c26231e52c7eae,0x471836049049f3b7,0xdc46826797ae,])),
                MNT6753Fr::new(BigInteger768([0x2ee4a96e4cda5f6f,0x7442a7b7f51fdbfc,0x23d03839ab7d811,0x1f873a8c0ddfd7a4,0x872f14e24612551a,0xd43181c852d5f78b,0xb2ff35a74130d2cd,0xd64aaa80f389157,0xb954953b8d35d74,0x37aba7a7212e96c,0xcce2fff62e11a3d4,0xfb3f9157120d,])),
            ],
            vec![
                MNT6753Fr::new(BigInteger768([0x626e4d0e6e3e1936,0x7c99da459f8385d0,0xbd84a2fb934889a6,0xff40b1979118e180,0x76cb8b37a32cce54,0x6c389f3f88157389,0xb9f0135ec3d92cc2,0xfd6a928e603a79be,0x5472af35b978d0a6,0x109995c9831f98c2,0x976c556bfe34da5a,0xf838693b701,])),
                MNT6753Fr::new(BigInteger768([0x58fb485fd781fcc6,0xd92a60427ce67147,0x2cca412943d39ade,0xc55d3362bac1743,0xcb8dcfa4ae0fcda1,0x25bde06b8f99facd,0x2d30b30add5faa3e,0xbe0ebdda1ba7458d,0x296f6010c1db1c7b,0x506364ec0031a00e,0x24c13847d3fe6ab7,0xea0c23423f1a,])),
                MNT6753Fr::new(BigInteger768([0xc36816e6dafa2f57,0x554255a6e34a34d4,0x29f17ff72b3c5695,0xae97815a3cc05077,0x64a0824e4b9b1aae,0x267cf597a9a556ef,0x8d8c67fc33757cbc,0xad2db4d1a3c73012,0xf3fcee4d169de439,0xfc4632cd5cb31baf,0xe1420a2c4e68de6,0x1bd34ad51cd02,])),
            ],
            vec![
                MNT6753Fr::new(BigInteger768([0x160dacef01b59531,0x313dd55066586bd8,0xdcc16189ec00c953,0xcc44967095828982,0x1066ee6f582ba5ea,0x3d879be40c078337,0xb9cb0ef83e1b4a51,0xc9b91de1e758c41,0xe578ceb8440e2bb8,0x3d6f2d210d4278df,0x2bab83b243a3335a,0x1afd20a9dbdc7,])),
                MNT6753Fr::new(BigInteger768([0x3a7ee60628dc201d,0xae1dcd081da757a,0xde1625ce6e93bc19,0xfb1a64dd14c0ae77,0x1bb5eba30eb2f202,0xdf064e762ce2f903,0x9abc764fb4c55d03,0x6db04d43d811c05d,0x87d85ec650763745,0x1bdcd095b0e1ada2,0x8681985565baa005,0x154d78a914323,])),
                MNT6753Fr::new(BigInteger768([0x101437542e4c39d4,0xcbdcf8d57d75fdd2,0x40996ed826c3b401,0xe492943442e0833b,0xf088ed10c7619f8c,0xb8e27256e0a69172,0x7112494180a5924,0x58d0e045a50972e9,0x4285049c582ed300,0xba0daceb8ab6d3c0,0x5ebb479b97c4c24d,0x820fdfe15d33,])),
            ],
            vec![
                MNT6753Fr::new(BigInteger768([0x645f79445d3423f1,0x699de15f996c470c,0x3740c3b7e7818751,0xac5c029dba988fd2,0x7342c873ecef9aee,0x4ff8cedd8fa15877,0xa9f8d05cc0c37cdb,0x6342d403e9995fcc,0xcd1206bec9b26855,0x9c7d8a00045eb24d,0x9c63e4f9f6757a65,0x1b358d82afeb4,])),
                MNT6753Fr::new(BigInteger768([0x5c47dc04494f4bd2,0x9c673cd9289d41af,0x162259acba9d8d18,0x62cad4f296328097,0x8aaf9e1700b7c75d,0x55e78bf0544350b2,0x4f68ebcc4892c902,0xdab2889f96fa7b5b,0x2a03de10d75b9f18,0x1ea1e16fc08e4df6,0x6acecbff7d2f538,0x9435d0a83b56,])),
                MNT6753Fr::new(BigInteger768([0x57c48852f8169d69,0x770318c8f24e3ac0,0xa0305f4306f0fbf4,0xf24a6cdad69062c1,0x193310c1c542ab5e,0x34b6461663f4fe2a,0xe7a085a783023999,0xb5ce7b9c96faf8e0,0x7552f4cfa41a306a,0x2f174937af08a752,0x1a0cef0caa379120,0xaf994027adab,])),
            ],
            vec![
                MNT6753Fr::new(BigInteger768([0x9720001bc9352497,0x26db06d9f4127454,0x9cce839d50eab099,0xba25501620cf63a9,0x795125f6eb018f87,0x694e8cec73b544f8,0xdb77a066d8a2cdd5,0x7aabd5789a9eafe3,0x178cc6b3542ceaa6,0xa6ac0cd365b9c275,0x122759efe8da9356,0x8e1dde78adb9,])),
                MNT6753Fr::new(BigInteger768([0xa9c2b63431ec99e7,0xb05d41809af7e5dc,0x2cbd97c762aecb7,0x4d41c4687b6d4477,0x8381b288c0dbf80,0x50d30f6e9cd8073e,0xbd5d9a24ab8be9f5,0x53f6ff54d29bfaf6,0xdfcf47396745930f,0xf9624d429b121957,0x2eff2dd22352fa1c,0x8062baa0e970,])),
                MNT6753Fr::new(BigInteger768([0x686af5fafbfbf6ea,0x1e1c039393b53fbf,0x395bda15104e42d7,0x86bd133dc0ecd7de,0xe6edda60379dd98,0xa4b50608cd0cbda3,0x71914eaa21572,0x716fc727079df56d,0x92d198f1997ebcb0,0x2bc460bbd690afcc,0xed78f65c0b4e499e,0x2bfad26243bd,])),
            ]
        ];

        poseidon_permutation_regression_test::<MNT6753Fr, MNT6753PoseidonParameters>(start_states, end_states);
        test_routine::<MNT6753Fr, MNT6PoseidonHash>(3)
    }
}