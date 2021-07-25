extern crate rand;
extern crate rayon;

use algebra::{Field, PrimeField};

use std::{
    ops::Mul,
    marker::PhantomData
};

use crate::{
    crh::{
        FieldBasedHash,
        FieldBasedHashParameters,
        SBox,
    }, CryptoError, Error
};

pub mod batched_crh;

pub mod parameters;
pub use self::parameters::*;

pub mod sbox;
pub use self::sbox::*;

pub trait PoseidonParameters: 'static + FieldBasedHashParameters + Clone {
    const T: usize;  // Number of S-Boxes
    const R_F:i32;   // Number of full rounds
    const R_P:i32;   // Number of partial rounds
    const ZERO:Self::Fr;   // The zero element in the field
    const AFTER_ZERO_PERM: &'static[Self::Fr]; // State vector after a zero permutation
    const ROUND_CST: &'static[Self::Fr];  // Array of round constants
    const MDS_CST: &'static[Self::Fr];  // The MDS matrix

    /// Add round constants to `state` starting from `start_idx_cst`, modifying `state` in place.
    #[inline]
    fn add_round_constants(
        state: &mut [<Self as FieldBasedHashParameters>::Fr],
        start_idx_cst: &mut usize
    )
    {
        for d in state.iter_mut() {
            let rc = Self::ROUND_CST[*start_idx_cst];
            *d += &rc;
            *start_idx_cst += 1;
        }
    }

    /// Perform scalar multiplication between vectors `res` and `state`,
    /// modifying `res` in place.
    #[inline]
    fn dot_product(
        res: &mut <Self as FieldBasedHashParameters>::Fr,
        state: &mut [<Self as FieldBasedHashParameters>::Fr],
        mut start_idx_cst: usize
    ) {
        state.iter().for_each(|x| {
            let elem = x.mul(&Self::MDS_CST[start_idx_cst]);
            start_idx_cst += 1;
            *res += &elem;
        });
    }

    /// Perform matrix mix on `state`, modifying `state` in place.
    #[inline]
    fn matrix_mix(state: &mut Vec<<Self as FieldBasedHashParameters>::Fr>)
    {
        // the new state where the result will be stored initialized to zero elements
        let mut new_state = vec![<Self as FieldBasedHashParameters>::Fr::zero(); Self::T];

        let mut idx_cst = 0;
        for i in 0..Self::T {
            Self::dot_product(&mut new_state[i], state, idx_cst);
            idx_cst += Self::T;
        }
        *state = new_state;
    }
}

pub trait PoseidonShortParameters: PoseidonParameters {
    /// MDS matrix supporting short Montgomery multiplication with respect to the short
    /// Montgomery constant R_2=2^64
    const MDS_CST_SHORT: &'static [Self::Fr];
}

#[derive(Derivative)]
#[derivative(
Clone(bound = ""),
Debug(bound = ""),
)]
pub struct PoseidonHash<F: PrimeField, P: PoseidonParameters<Fr = F>, SB: SBox<Field = F, Parameters = P>>{
    state: Vec<F>,
    pending: Vec<F>,
    input_size: Option<usize>,
    updates_ctr: usize,
    mod_rate: bool,
    _parameters: PhantomData<P>,
    _sbox: PhantomData<SB>,
}

impl<F, P, SB> PoseidonHash<F, P, SB>
    where
        F: PrimeField,
        P: PoseidonParameters<Fr = F>,
        SB: SBox<Field = F, Parameters = P>,
{
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
            _sbox: PhantomData,
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
            while personalization_instance.pending.len() != 0 {
                personalization_instance.update(F::zero());
            }
            debug_assert_eq!(personalization_instance.pending.len(), 0);

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

    pub(crate) fn poseidon_perm (state: &mut Vec<F>) {

        // index that goes over the round constants
        let round_cst_idx = &mut 0;

        // First full rounds
        for _i in 0..P::R_F {
            // Add the round constants to the state vector
            P::add_round_constants(state, round_cst_idx);

            // Apply the S-BOX to each of the elements of the state vector
            SB::apply_full(state);

            // Perform matrix mix
            P::matrix_mix(state);
        }

        // Partial rounds
        for _i in 0..P::R_P {
            // Add the round constants to the state vector
            P::add_round_constants(state, round_cst_idx);

            // Apply S-BOX only to the first element of the state vector
            SB::apply_partial(state);

            // Perform matrix mix
            P::matrix_mix(state);
        }

        // Second full rounds
        for _i in 0..P::R_F {
            // Add the round constants
            P::add_round_constants(state, round_cst_idx);

            // Apply the S-BOX to each of the elements of the state vector
            SB::apply_full(state);

            // Perform matrix mix
            P::matrix_mix(state);
        }
    }
}

impl<F, P, SB> FieldBasedHash for PoseidonHash<F, P, SB>
    where
        F: PrimeField,
        P: PoseidonParameters<Fr = F>,
        SB: SBox<Field = F, Parameters = P>,
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
    use algebra::{Field, PrimeField};
    use crate::crh::{
        FieldBasedHash, SBox,
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

    fn poseidon_permutation_regression_test<F: PrimeField, P: PoseidonParameters<Fr = F>, SB: SBox<Field = F, Parameters = P>>(
        start_states: Vec<Vec<F>>,
        end_states:   Vec<Vec<F>>,
    )
    {
        // Regression test
        start_states.into_iter().zip(end_states).enumerate().for_each(|(i, (mut start_state, end_state))| {
            PoseidonHash::<F, P, SB>::poseidon_perm(&mut start_state);
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
            MNT4PoseidonHash, MNT4753PoseidonParameters, MNT4InversePoseidonSBox
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

        poseidon_permutation_regression_test::<MNT4753Fr, MNT4753PoseidonParameters, MNT4InversePoseidonSBox>(
            start_states, end_states
        );
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
            MNT6PoseidonHash, MNT6753PoseidonParameters, MNT6InversePoseidonSBox,
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

        poseidon_permutation_regression_test::<MNT6753Fr, MNT6753PoseidonParameters, MNT6InversePoseidonSBox>(
            start_states, end_states
        );
        test_routine::<MNT6753Fr, MNT6PoseidonHash>(3)
    }

    #[cfg(feature = "bn_382")]
    #[test]
    fn test_poseidon_hash_bn382_fr() {
        use algebra::{
            biginteger::BigInteger384,
            fields::bn_382::Fr as BN382Fr
        };
        use crate::crh::poseidon::parameters::bn382::{
            BN382FrPoseidonHash, BN382FrPoseidonParameters, BN382FrQuinticSbox
        };

        // Test vectors are computed via the script in ./parameters/scripts/permutation_bn382.sage
        let start_states = vec![
            vec![BN382Fr::zero(); 3],
            vec![
                BN382Fr::new(BigInteger384([0x3c2fc28fee546f2f,0x46673e3f762a05a9,0xf3a4196fb7f077b0,0xea452bd940906dd0,0x61a33d4ae39ee3a1,0xb0fe1409f6b6ad,])),
                BN382Fr::new(BigInteger384([0x4ca6094864916d57,0x2ca4974e3b8e4d6,0x4a05915e47dd8a18,0x5e6888ec4e9811ed,0xb1ddb601c4144f40,0x45c5e2dccf92992,])),
                BN382Fr::new(BigInteger384([0x20e98c5412e8a53c,0xb2df907a45237f4,0x89db0df005eb52fb,0xc77948ae1a2a2cda,0xf5ddb01fdc5f2ca4,0x17cd7c819448cb46,])),
            ],
            vec![
                BN382Fr::new(BigInteger384([0x541a0ca7bfcc881f,0xf88b8f238697be3c,0x36e61e96d2fb8d14,0x1a3edaa7cbaee4cb,0x55a2ae58ee66a979,0x100171f764d62113,])),
                BN382Fr::new(BigInteger384([0x4abbd93002288653,0x37e17a329d1fa261,0xcd880c8eaf7a18b9,0xb0c2cd616408d2cf,0x5e938101f5333493,0x22a361e49171b56c,])),
                BN382Fr::new(BigInteger384([0x75efbb3b47ed610d,0x872b59023b1582f,0x154f1c9f55385a05,0x130ecac1483ed87c,0xc9c4f03d0a0e838,0x11985516d2a7f963,])),
            ],
            vec![
                BN382Fr::new(BigInteger384([0x5676a2aa9827db5c,0x42dffaf55931d898,0x4df6a1acb8359ba6,0xb6c57235a1057d95,0x8c80b33063239cec,0xc7219289e5b6fbe,])),
                BN382Fr::new(BigInteger384([0xe17739d0259851fc,0x8cf3a336d885e861,0x5147bc1f93978a33,0x371ff88b2aaa0b59,0xc4fac8e7e213807e,0x49b925c3136f71b,])),
                BN382Fr::new(BigInteger384([0xe2fdff72d46fcb0a,0x4067d514d9cb9ecf,0xf51b3c5b3bc11d00,0xeed7a12d7d42ee4c,0xc8bb6a1b0a079aa7,0xd047e537eb9ac58,])),
            ],
            vec![
                BN382Fr::new(BigInteger384([0xcd0f5560277aad4d,0xffe03011802d3fd1,0xf74446bb0aa3e8e2,0x5e1f3daa54d09f36,0x459daf13600a2960,0x1cd498d82eb74a2d,])),
                BN382Fr::new(BigInteger384([0xc1ca68ef9f0d7346,0x78bfeb6a95ea63e5,0xce164dee9dba93a,0x60f8dbaa8634a63a,0xfcfb923ab4911528,0x93128aeb82dbf04,])),
                BN382Fr::new(BigInteger384([0x16c5a3b0f84a1808,0x1bb720c4473aa741,0xe3dd83f67121d1fb,0x31dc7f9ff20507b8,0xc86761e6ec443333,0x6f67c54083f05db,])),
            ]
        ];

        let end_states = vec![
            vec![
                BN382Fr::new(BigInteger384([0x3600ae9dea9ba41a,0x17a35e3bedfc2e1,0x7ee93e40052b3867,0xc555ef28f09e84e9,0x1ef349664ad402cf,0x1c49706c59f09b25,])),
                BN382Fr::new(BigInteger384([0xbb6f865c755b9100,0x6f6ccbea5f0c5847,0x4cfd3606c21c2573,0x3512ec3dc6889f67,0xc7981de6b0710b5f,0x109fe23f817aa0cf,])),
                BN382Fr::new(BigInteger384([0x39de13d041934215,0x5370089a3da7c4fe,0x512952ce97e48c03,0xe1c26f50f4c9c4c1,0x1f008942e907b93e,0x1910b7b5453ff08f,])),
            ],
            vec![
                BN382Fr::new(BigInteger384([0xf3b93ceda4f3a5c,0x5dcd6b6bc043fd10,0x8d383811267393b4,0x66f48dee2d1b12df,0xbdb9d022d8ef1832,0x3d7e58786b39ef4,])),
                BN382Fr::new(BigInteger384([0x44aa122585436d31,0x28935d91839eef2b,0xda2ba836d955d3fe,0x200274d572c207a8,0x68ea32c32bf9e76c,0x1b6e87d7d7bd71b6,])),
                BN382Fr::new(BigInteger384([0x65ba9efee2204115,0x81b822106a189c40,0x72b7d6e504e281b4,0xa51d8ac7dd820df0,0x1ea1f1cb92430cbc,0x23a85bdeb2d2dd16,])),
            ],
            vec![
                BN382Fr::new(BigInteger384([0x47c7598c44ff8d16,0x9a2a4de7e4caa199,0xa64228ccfb671b,0xe507c52bab4c227c,0xa03bae146874c577,0x142abb97131a15ce,])),
                BN382Fr::new(BigInteger384([0x6e5c0a1b6c74884d,0xf5bb78ce31dc03be,0xe12a8aea2fdbfb1c,0x27806b8e798e5047,0xdb908a200b3040d9,0xe722e2590de5b3d,])),
                BN382Fr::new(BigInteger384([0xa3c9528966e64486,0x475589fea46633f1,0xd74899c26b7cc411,0x1771d0995b78fb5d,0xf4e48a25c61e9202,0x13751c53efdbf754,])),
            ],
            vec![
                BN382Fr::new(BigInteger384([0xbaa9e75cb23bcf05,0x35d727f254ae75d7,0xacb20d326450e2b8,0x177c73eda4c84fdb,0x51f291a5f9dd6033,0x8788cee947e9501,])),
                BN382Fr::new(BigInteger384([0xf1b326ebc984ec0,0x866f44f24cf07054,0x5f070db622ccd3da,0xceb0f26208090d9e,0xdd7bd626dbb1d31e,0xa8a45f03c973521,])),
                BN382Fr::new(BigInteger384([0x32e4799fc1db07b1,0xbbdcbf7c6b9e2f24,0xf7cbd541b37e4650,0xd8143503afc7320a,0x75a91583524c9a16,0x1c9f9295f8bce898,])),
            ],
            vec![
                BN382Fr::new(BigInteger384([0xd5db324446615bf8,0x96f94dd6887732e0,0x56020c6319093a3a,0x5ef153e7bc15f69b,0x1b87643733a4b798,0x16787d5e34111ed,])),
                BN382Fr::new(BigInteger384([0x4558ed95b354fe81,0x31ba491852c4023,0x98af2996db40ba92,0xb4c3ac53e548ec3d,0x96c9e81d713719ea,0x1eefdfa3b6b479ae,])),
                BN382Fr::new(BigInteger384([0x3340788274f54c1f,0xb2d040485d2fd9d6,0xd7df55b13440dbf3,0x856bf5fc77c7f48b,0x48cf9764e0e67a05,0x1816ef21b6373a7,])),
            ]
        ];

        poseidon_permutation_regression_test::<BN382Fr, BN382FrPoseidonParameters, BN382FrQuinticSbox>(
            start_states, end_states
        );
        test_routine::<BN382Fr, BN382FrPoseidonHash>(3)
    }

    #[cfg(feature = "bn_382")]
    #[test]
    fn test_poseidon_hash_bn382_fq() {
        use algebra::{
            biginteger::BigInteger384,
            fields::bn_382::Fq as BN382Fq
        };
        use crate::crh::poseidon::parameters::bn382_dual::{
            BN382FqPoseidonHash, BN382FqPoseidonParameters, BN382FqQuinticSbox
        };

        // Test vectors are computed via the script in ./parameters/scripts/permutation_bn382dual.sage
        let start_states = vec![
            vec![BN382Fq::zero(); 3],
            vec![
                BN382Fq::new(BigInteger384([0x239d004c236ddb,0x88d83e760e8bd5bb,0x2ca0f68190713e45,0x8f6a964f924c8fff,0x62a854d505daa3f3,0x295e6179129332c,])),
                BN382Fq::new(BigInteger384([0xa4a0f24f69849fbf,0x751e1bb2f93df901,0x6955afa141342da,0x3a242cea266d1ac2,0x4f838810d428645,0x397b9821248dd08,])),
                BN382Fq::new(BigInteger384([0x5985d03eb267a372,0x6491f79810a21027,0xe65805fff01b641a,0x3aa8f9b916f74025,0x7ed27d962144ab7f,0x17f25f1815f2512c,])),
            ],
            vec![
                BN382Fq::new(BigInteger384([0x86ead0985648077a,0x7a50ef9f2086cc9d,0x69c612dbec57975e,0x8647aacd9ab88959,0x3a5fabf8692b8d12,0xbdb03daf76eb57,])),
                BN382Fq::new(BigInteger384([0x4409ea78e288db0a,0xc14e0a759b5fd26b,0x1ae0285264db243b,0xf2be0cf31a448a05,0xd103243aef14ada3,0x1189adbc498d1570,])),
                BN382Fq::new(BigInteger384([0xfa5e0c518b29c440,0x28cbb1257edeb8a6,0x7a120a8c0658b3b5,0x13040f12fb2249f8,0xb71143b9ada3922c,0x1ee9611738dbe1b3,])),
            ],
            vec![
                BN382Fq::new(BigInteger384([0x3bbae40afacfabc1,0x518f05b12a86d30,0xa7c6c267a8c546f3,0xdff2338e035d8c38,0x45cad929932db574,0x179803640786a069,])),
                BN382Fq::new(BigInteger384([0xe4c488029c73ab3d,0x9cbea6f936421688,0xa733a951138f8904,0x9566d6bc3392168,0xe102fe13109c07ae,0x1e4c4733f9c926f1,])),
                BN382Fq::new(BigInteger384([0xbeeabdfd33d7d4d4,0x258d58e0edf24637,0x644767bec95dd149,0x780c156441e1c292,0xb0b849ce82fd90a2,0xb189d134bfa9ced,])),
            ],
            vec![
                BN382Fq::new(BigInteger384([0x1d9b42d2a2f73ea8,0xb9b3cf9c1e9aea41,0xa3c8780de2c255f2,0xff9617a521bc6a15,0x3dfe0e09411bbce1,0x1872aac1dea2aba8,])),
                BN382Fq::new(BigInteger384([0x166383182fda3435,0x3125ac12879ae7e6,0x425286423e9432b,0x796686a5176807f8,0x826f8b280eb7669c,0x37172d9cb2e8efd,])),
                BN382Fq::new(BigInteger384([0x2fc9a35fae8c69f3,0x8dca72688a8fa1c4,0xe7a690c67ed759d6,0xcde98c6072dd8eb4,0xa4bd01fd0dbe1bcd,0xf556423e114180e,])),
            ]
        ];

        let end_states = vec![
            vec![
                BN382Fq::new(BigInteger384([0x27dcc9c1f001c02d,0x7fc9de4b5ab915ed,0x7c6832557c4a410d,0x320b95a8fa27bf32,0xe5c89c9c09bd67e5,0x65748e22de4f8c5,])),
                BN382Fq::new(BigInteger384([0x7cdb27778c5d6796,0xad588ee542be3389,0x68e926bfdd6398ec,0xe432240624573240,0x2766c91ade70f83f,0x170646120652b37c,])),
                BN382Fq::new(BigInteger384([0xcada65af3ba4e9c4,0x7e4561e9933627cd,0x8cb8757ddb2e0730,0x610ecc5beda633e0,0x984de49537e8c3ec,0x1349deb07a8f6f52,]))
            ],
            vec![
                BN382Fq::new(BigInteger384([0xcfd422c316b20422,0xf15801f500d95821,0x360f5beb123f7d4e,0xfc13f1eabfe897f0,0xc70e46eea3b47d2c,0x14eb20b8f8cc25e5,])),
                BN382Fq::new(BigInteger384([0x18fb3a5f70545729,0xadc0d9cd0b986c7b,0xc0f502215de819a9,0x21bff5966fdde339,0xc39b173777b1f86b,0x1e01840238fce37a,])),
                BN382Fq::new(BigInteger384([0x70fd0a437704dfb5,0xc0afdaef11a41929,0x8a3d1c5e46648541,0x97c16c79daeb557d,0xd18b01c167ec00e6,0x10d02b9f59132a1d,])),
            ],
            vec![
                BN382Fq::new(BigInteger384([0x1035143aba9695b9,0xf532c66887edbfcd,0xa6bd2998470d554f,0x831687ccd8a703ff,0xb75bed9a7ae1bab5,0x8b4c6d206c82fb8,])),
                BN382Fq::new(BigInteger384([0x7e3d0019dd9387ab,0x746b6db1b8c19f4b,0x2964ec70d389adf6,0x8333f2f4045ebb5f,0x31832aff0cd42bc1,0x16572d68fc8031d5,])),
                BN382Fq::new(BigInteger384([0x208b12c54d10bf3b,0xce12a04a890b4859,0x24fc1c25be961547,0xf8e6e4ee5cf48107,0x43a590c19365296e,0x58b7ff26592e23f,])),
            ],
            vec![
                BN382Fq::new(BigInteger384([0x4f52c2933ef585f2,0x93b9868fb78ca000,0x390b415d3dda671c,0x7376e52933a4470,0x6f4cb578d987419,0xc539440279dc102,])),
                BN382Fq::new(BigInteger384([0x3b8ed76f186a092f,0xfc7e9b70f3a206d0,0xa3bbb0c1436c65a2,0xfe0aeae213ba4473,0x9d8ff7b60fe2b888,0x35cb00af8ae79df,])),
                BN382Fq::new(BigInteger384([0x1e27e68b262adee9,0xd4b7220a4be055ae,0x4ac1d5ab2530b8b,0x34e9beab4c8c6260,0xa37fff7e0bb5c229,0xa75e8ec286abe8e,])),
            ],
            vec![
                BN382Fq::new(BigInteger384([0x9b735d2a3353402c,0xd4547e70eb8130fa,0x2438c5a8bed96075,0x32fdf7691a26f030,0xa1f649648c34ed64,0x22a1ead2ba837f97,])),
                BN382Fq::new(BigInteger384([0x39b0c7a9271496c,0xcfec5f805bdb5e00,0xa9aead920a13442d,0xf8c824e2dedc3993,0x81b407a948baa360,0x205d8c200fb40967,])),
                BN382Fq::new(BigInteger384([0xf9cc3c9cf970f38c,0xaf92136db468bbb9,0xb1c839b8e1eb9561,0xf92e59ecbe79cc84,0x34c857f5954e45f8,0x8344e8ada34f5d1,])),
            ]
        ];

        poseidon_permutation_regression_test::<BN382Fq, BN382FqPoseidonParameters, BN382FqQuinticSbox>(
            start_states, end_states
        );
        test_routine::<BN382Fq, BN382FqPoseidonHash>(3)
    }

    #[cfg(feature = "tweedle")]
    #[test]
    fn test_poseidon_hash_tweedle_fr() {
        use algebra::{
            biginteger::BigInteger256,
            fields::tweedle::Fr as TweedleFr
        };
        use crate::crh::poseidon::parameters::tweedle_dee::{
            TweedleFrPoseidonHash, TweedleFrPoseidonParameters, TweedleFrQuinticSbox
        };

        // Test vectors are computed via the script in ./parameters/scripts/permutation_deefr.sage
        let start_states = vec![
            vec![TweedleFr::zero(); 3],
            vec![
                TweedleFr::new(BigInteger256([0x2d9ced12b8448fa3,0xe47617895bcb1def,0xdb309341af8fc9bc,0x3518ed3d596d9b3d,])),
                TweedleFr::new(BigInteger256([0x2f00b53bfb408372,0x6de08091d9994983,0x30787444ac8639a3,0x18b1a8fe589e66ad,])),
                TweedleFr::new(BigInteger256([0xbbff40a91825c30d,0xa82ca4dd45ed43cd,0x3ce8daf6c9c21029,0x10c0f7735f33aa7a,])),
            ],
            vec![
                TweedleFr::new(BigInteger256([0x5f37a0bd77589e1f,0x5473621f06e318b0,0x134c69d294364fc2,0x17ce475fc0918e98,])),
                TweedleFr::new(BigInteger256([0xf997aedfd435a00c,0xff8244711a05ace4,0x111f3729665dfce3,0x12e06c5d75a20f44,])),
                TweedleFr::new(BigInteger256([0x4fe219488f716f3b,0x47994803d7aa1b4b,0x83c0b9401250e3df,0xc55e3e5129040af,])),
            ],
            vec![
                TweedleFr::new(BigInteger256([0x1c88b7f17d83e522,0x63bbb3d972a8a79,0x3cd3b269e9148e61,0x107064754c2219f6,])),
                TweedleFr::new(BigInteger256([0xd98347c19ef61123,0x8c2f919a2ce03104,0x19a6ebeb17c8d50b,0x211359dab98e662b,])),
                TweedleFr::new(BigInteger256([0x6fca9aeca36a6a90,0x9a5901d4db4cb38b,0xb7a625b6fa9c1d25,0x1c0c5a9e4863c446,])),
            ],
            vec![
                TweedleFr::new(BigInteger256([0x52cc4aa39d8838b8,0x412ba25c63120ebb,0x667515874f0074d6,0x1d2f166897ea99e,])),
                TweedleFr::new(BigInteger256([0x466265a678233c51,0xd6b41807e24ee39f,0xee5874453e9c291c,0x1b0bbd1b8e79ea9d,])),
                TweedleFr::new(BigInteger256([0x49d2b1885d136bf6,0xfebba4a8e8c0595b,0xa5b4ca600f485e66,0x27c2b78d22e855c0,])),
            ],
        ];

        let end_states = vec![
            vec![
                TweedleFr::new(BigInteger256([0x85614442a60ac11a,0x55a43ca8180d2e08,0x43f61ff197080ac4,0x19d87eb89a42aaf1,])),
                TweedleFr::new(BigInteger256([0xa2f6b5a9a16d3790,0xc947563b131a126c,0x52c19607bb4b6640,0xc4604a460df1c57,])),
                TweedleFr::new(BigInteger256([0x7d8f3c1679a9cbe2,0xb09fdc38ee15fe77,0x810720bf23be8578,0x2ab876d1a0abfa95,])),
            ],
            vec![
                TweedleFr::new(BigInteger256([0xc4a37b8664180077,0xd8390d652933725e,0xaafa5d29eb656edb,0x296682761320f48c,])),
                TweedleFr::new(BigInteger256([0x2fffbed47e729020,0x6d243b1d399f42dd,0x2bcea2d0461856d7,0x2fc6f9c7c62a5088,])),
                TweedleFr::new(BigInteger256([0x8b617097039cbf5f,0xc3e9594e65f53809,0x96f163d2a6e08e55,0x1283bbfbfafe0185,])),
            ],
            vec![
                TweedleFr::new(BigInteger256([0xb0e21925172f0ba3,0x22bb8d3720914af7,0x31ee2b9a26424619,0x2184d5590df49e25,])),
                TweedleFr::new(BigInteger256([0x4f525fe270112fb8,0x59d975c2bc66f456,0x1740475c80005233,0x3f44acd2d334fee9,])),
                TweedleFr::new(BigInteger256([0xda02921fa73b4778,0xb9b7c2742272dbeb,0xb3491dacb990965c,0x3cffd4206f4264e,])),
            ],
            vec![
                TweedleFr::new(BigInteger256([0x9a5d804c8f8980d7,0x60f4ba8f01fccce4,0x95428b68f3a9eba3,0x3108ed7e0636e1e7,])),
                TweedleFr::new(BigInteger256([0xf5e24f59c7e404d7,0xf4a10531d95222b1,0xb55cfa77a621836f,0x15f7c485bf9b2bf1,])),
                TweedleFr::new(BigInteger256([0xf65bd157052e1b45,0x180aa5b7e51b8a46,0xe451d510b5cf9dae,0x7cdd9f00493bc73,])),
            ],
            vec![
                TweedleFr::new(BigInteger256([0x7c080f4b62e78aab,0xc6294e279a622677,0xcabd73efb2584d6d,0x10186a71cc08159e,])),
                TweedleFr::new(BigInteger256([0xdb3d4f4a63e1324d,0x6705ae25ff9b471f,0xccae1d131341f589,0x1b31cd963165eccc,])),
                TweedleFr::new(BigInteger256([0x9860019e6edc3f2f,0x14ca7a30bb1a5c36,0xf4e9f4abe3f7ef0c,0x143d7bf07e7f54c7,])),
            ],
        ];

        poseidon_permutation_regression_test::<TweedleFr, TweedleFrPoseidonParameters, TweedleFrQuinticSbox>(
            start_states, end_states
        );
        test_routine::<TweedleFr, TweedleFrPoseidonHash>(3)
    }

    #[cfg(feature = "tweedle")]
    #[test]
    fn test_poseidon_hash_tweedle_fq() {
        use algebra::{
            biginteger::BigInteger256,
            fields::tweedle::Fq as TweedleFq
        };
        use crate::crh::poseidon::parameters::tweedle_dum::{
            TweedleFqPoseidonHash, TweedleFqPoseidonParameters, TweedleFqQuinticSbox
        };

        // Test vectors are computed via the script in ./parameters/scripts/permutation_dumfr.sage
        let start_states = vec![
            vec![
                TweedleFq::zero(); 3
            ],
            vec![
                TweedleFq::new(BigInteger256([0x530261dfc524611d,0xebde2e5e0c454577,0x31c9a2fd3288dbd8,0x22faf97cf0bfa8ed,])),
                TweedleFq::new(BigInteger256([0x25f47e32d936f0c0,0x9c88b0ffb8d56acc,0x3c1a4050825c76ac,0xf81aaaddfb679df,])),
                TweedleFq::new(BigInteger256([0x129cb322f4812820,0x5b218d2750d9cc33,0x5baa3f8af95e185b,0xf5713c92c9b59a5,])),
            ],
            vec![
                TweedleFq::new(BigInteger256([0x8c70fb5700e28179,0x58d04dff4aeb7baa,0x7d229f69585bbc4c,0x1a53f352bbb741f,])),
                TweedleFq::new(BigInteger256([0x983971f4bc40e955,0xf9c4aa245dc69370,0xc90afb10e865d7fa,0x25c68f3eda91e782,])),
                TweedleFq::new(BigInteger256([0x553902e820896d7e,0xea7238f532c5b890,0x66c31bc5cacadbb5,0x11fbf51d7acd7811,])),
            ],
            vec![
                TweedleFq::new(BigInteger256([0x8c5101f47ede0f2b,0xdde609c8ee90d5e9,0xf53611e4c9658d0b,0x9b8ad64dd287d37,])),
                TweedleFq::new(BigInteger256([0xe79daeebc658d0a,0x3019b7ed8cae3dd8,0xe4966f5f01879f27,0x2f1328f79025e70c,])),
                TweedleFq::new(BigInteger256([0x49ad0534394806ae,0x6ab073974f741a93,0x3e043b146513dfe5,0x29b158cd24e843e4,])),
            ],
            vec![
                TweedleFq::new(BigInteger256([0x3a410990938e76ed,0x4bd4f247c6c2215b,0xe815c6d61abfe6f9,0x94daa5bcfb9eb6f,])),
                TweedleFq::new(BigInteger256([0x3787fbb0c8dcfe1a,0xf67406e5daf43fae,0x7a5fc8f335f28767,0x18ff0f241943eec8,])),
                TweedleFq::new(BigInteger256([0xc72a940881085fd6,0x7096ba03e87353af,0x32decb002f5a4e83,0x492cc5ac858b06a,])),
            ],
        ];

        let end_states = vec![
            vec![
                TweedleFq::new(BigInteger256([0x46ef7b471f039f54,0x7516283cc67869f2,0x561a6334ba7a39f1,0x293842a1538ac01b,])),
                TweedleFq::new(BigInteger256([0x6f10ff3b97995e3b,0x7650f70901d51a88,0x9f13555ea4caf2eb,0x14ed7f5560a0a1e1,])),
                TweedleFq::new(BigInteger256([0x815126351fe00f44,0x921a5f3ad5a6e83c,0x5f614c0b1bdaf5f7,0x7733c69a8892f0e,])),
            ],
            vec![
                TweedleFq::new(BigInteger256([0xf39ca6429f499eb1,0x69657c642b509baa,0xbb0a2f6bb3a44a7b,0x1b0f054ee6b06ee5,])),
                TweedleFq::new(BigInteger256([0x9eab499dc61a7d92,0x457d1a9027e66bd4,0x74f80311cef652a5,0x2f0dc832cc821ed,])),
                TweedleFq::new(BigInteger256([0xe5949837b34cdd97,0x2fdd08e41ac8e36f,0xbfcb6768fbb981d,0x1521b70d21fc43fb,])),
            ],
            vec![
                TweedleFq::new(BigInteger256([0x21fb36a475c20033,0x6a938adf93ceda77,0xa05bc36806e89296,0x1cd7a0d468136dd3,])),
                TweedleFq::new(BigInteger256([0x6295c60c77022ca5,0x440a39652987ef94,0xbe9a8f921e81b656,0x3ade3ff16b820c56,])),
                TweedleFq::new(BigInteger256([0x62f4df55b1158a3d,0x6787fff1b51e08ed,0x47b46cd1709e9d30,0x3c4bbad805b5838c,])),
            ],
            vec![
                TweedleFq::new(BigInteger256([0xf0b39ffa74b62183,0x9c87a4fea04e092a,0xe7ef4462efcf6492,0x1495692d563b0275,])),
                TweedleFq::new(BigInteger256([0x1758eeffd0793b03,0x37e1f13b2b104aa,0x71c181dd5d62c9d,0x3448bf7ebad19d00,])),
                TweedleFq::new(BigInteger256([0x63feeddf9fd791f,0xcf11513a74efebf6,0xc046e6ff5b45f4af,0x13a773bcdaabf9b1,])),
            ],
            vec![
                TweedleFq::new(BigInteger256([0x6f2ad1eed8b08a65,0x23e051559fea114f,0x6e9855acf367f614,0x1f6ff3e5034d9adb,])),
                TweedleFq::new(BigInteger256([0xc76c27513034009f,0xf08aae84a5bdaf00,0xb4614eed8e6839d5,0x18b4587f29cdb052,])),
                TweedleFq::new(BigInteger256([0xa5a9c19386d171db,0x57321c0b6d91fa65,0xaa19cb2f60d37e5b,0x12a05d4caaa7d0ca,])),
            ],
        ];

        poseidon_permutation_regression_test::<TweedleFq, TweedleFqPoseidonParameters, TweedleFqQuinticSbox>(
            start_states, end_states
        );
        test_routine::<TweedleFq, TweedleFqPoseidonHash>(3)
    }
}