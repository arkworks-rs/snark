extern crate rand;
extern crate rayon;

use algebra::PrimeField;

use std::marker::PhantomData;

use crate::crh::{
    FieldBasedHash,
    FieldBasedHashParameters,
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
    const C2:Self::Fr;     // The constant to add in the position corresponding to the capacity
    const AFTER_ZERO_PERM: &'static[Self::Fr]; // State vector after a zero permutation
    const ROUND_CST: &'static[Self::Fr];  // Array of round constants
    const MDS_CST: &'static[Self::Fr];  // The MDS matrix
}

#[derive(Debug)]
pub struct PoseidonHash<F: PrimeField, P: PoseidonParameters<Fr = F>, SB: PoseidonSBox<P>>{
    state: Vec<F>,
    pending: Vec<F>,
    _parameters: PhantomData<P>,
    _sbox: PhantomData<SB>,
}

impl<F, P, SB> PoseidonHash<F, P, SB>
    where
        F: PrimeField,
        P: PoseidonParameters<Fr = F>,
        SB: PoseidonSBox<P>,
{
    #[inline]
    fn apply_permutation(&mut self) {
        for (input, state) in self.pending.iter().zip(self.state.iter_mut()) {
            *state += input;
        }
        self.state[P::R] += &P::C2;
        Self::poseidon_perm(&mut self.state);
    }

    #[inline]
    fn _finalize(&self) -> F {
        let mut state = self.state.clone();
        for (input, s) in self.pending.iter().zip(state.iter_mut()) {
            *s += input;
        }
        state[P::R] += &P::C2;
        Self::poseidon_perm(&mut state);
        state[0]
    }

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
            SB::apply_full(state, false)
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
            SB::apply_partial(state);
        }

        // Second full rounds
        // Process only to R_F - 1 iterations. The last iteration does not contain a matrix mix
        for _i in 0..(P::R_F-1) {

            // Add the round constants
            for d in state.iter_mut() {
                let rc = P::ROUND_CST[round_cst_idx];
                *d += &rc;
                round_cst_idx += 1;
            }

            // Apply the S-BOX to each of the elements of the state vector
            SB::apply_full(state, false);
        }

        // Last full round does not perform the matrix_mix
        // Add the round constants
        for d in state.iter_mut() {
            let rc = P::ROUND_CST[round_cst_idx];
            *d += &rc;
            round_cst_idx += 1;
        }

        // Apply the S-BOX to each of the elements of the state vector
        SB::apply_full(state, true);
    }
}

impl<F, P, SB> FieldBasedHash for PoseidonHash<F, P, SB>
    where
        F: PrimeField,
        P: PoseidonParameters<Fr = F>,
        SB: PoseidonSBox<P>,
{
    type Data = F;
    type Parameters = P;

    fn init(personalization: Option<&[Self::Data]>) -> Self {
        let mut state = Vec::with_capacity(P::T);
        for i in 0..P::T {
            state.push(P::AFTER_ZERO_PERM[i]);
        }
        let mut instance = Self {
            state,
            pending: Vec::with_capacity(P::R),
            _parameters: PhantomData,
            _sbox: PhantomData,
        };

        // If personalization Vec is not multiple of the rate, we pad it with zero field elements.
        // This will allow eventually to precompute the constants of the initial state. This
        // is exactly as doing H(personalization, padding, ...). NOTE: this way of personalizing
        // the hash is not mentioned in https://eprint.iacr.org/2019/458.pdf
        if personalization.is_some(){
            let personalization = personalization.unwrap();

            for &p in personalization.into_iter(){
                instance.update(p);
            }

            let padding = if personalization.len() % P::R != 0 {
                P::R - ( personalization.len() % P::R )
            } else {
                0
            };

            for _ in 0..padding {
                instance.update(F::zero());
            }
            assert_eq!(instance.pending.len(), 0);
        }
        instance
    }

    // Note: `Field` implements the `Copy` trait, therefore invoking this function won't
    // cause a moving of ownership for `input`, but just a copy. Another copy is
    // performed below in `self.pending.push(input);`
    // We can reduce this to one copy by passing a reference to `input`, but from an
    // interface point of view this is not logically correct: someone calling this
    // functions will likely not use the `input` anymore in most of the cases
    // (in the other cases he can just clone it).
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

    fn reset(&mut self, personalization: Option<&[Self::Data]>) -> &mut Self {
        let new_instance = Self::init(personalization);
        *self = new_instance;
        self
    }
}

#[cfg(test)]
mod test {
    use algebra::{
        fields::{
            mnt4753::Fr as MNT4753Fr,
            mnt6753::Fr as MNT6753Fr,
        },
    };
    use std::str::FromStr;
    use algebra::biginteger::BigInteger768;
    use crate::crh::{
        poseidon::parameters::{
            mnt4753::MNT4PoseidonHash,
            mnt6753::MNT6PoseidonHash,
        },
        FieldBasedHash,
    };

    #[test]
    fn test_poseidon_hash_mnt4() {

        // Regression test
        let expected_output = MNT4753Fr::new(BigInteger768([120759599714708995, 15132412086599307425, 1270378153255747692, 3280164418217209635, 5680179791594071572, 2475152338055275001, 9455820118751334058, 6363436228419696186, 3538976751580678769, 14987158621073838958, 10703097083485496843, 48481977539350]));
        let mut poseidon_digest = MNT4PoseidonHash::init(None);
        let input = [MNT4753Fr::from_str("1").unwrap(), MNT4753Fr::from_str("2").unwrap()];

        let output = poseidon_digest
            .update(input[0])
            .update(input[1])
            .finalize();
        assert_eq!(output, expected_output, "Outputs do not match for MNT4753.");

        // Test finalize() holding the state and allowing updates in between different calls to it
        poseidon_digest
            .reset(None)
            .update(input[0].clone());
        poseidon_digest.finalize();
        poseidon_digest.update(input[1].clone());
        assert_eq!(output, poseidon_digest.finalize());

        //Test finalize() being idempotent
        assert_eq!(output, poseidon_digest.finalize());
    }

    #[test]
    fn test_poseidon_hash_mnt4_single_element() {
        let expected_output = MNT4753Fr::new(BigInteger768([10133114337753187244, 13011129467758174047, 14520750556687040981, 911508844858788085, 1859877757310385382, 9602832310351473622, 8300303689130833769, 981323167857397563, 5760566649679562093, 8644351468476031499, 10679665778836668809, 404482168782668]));
        let mut poseidon_digest = MNT4PoseidonHash::init(None);
        poseidon_digest.update(MNT4753Fr::from_str("1").unwrap());
        let output = poseidon_digest.finalize();
        assert_eq!(output, expected_output, "Outputs do not match for MNT4753.");
    }

    #[test]
    fn test_poseidon_hash_mnt4_three_element() {
        let expected_output = MNT4753Fr::new(BigInteger768([5991160601160569512, 9804741598782512164, 8257389273544061943, 15170134696519047397, 9908596892162673198, 7815454566429677811, 9000639780203615183, 8443915450757188195, 1987926952117715938, 17724141978374492147, 13890449093436164383, 191068391234529]));
        let mut poseidon_digest = MNT4PoseidonHash::init(None);

        for i in 1..=3{
            poseidon_digest.update(MNT4753Fr::from(i as u64));
        }

        let output = poseidon_digest.finalize();
        assert_eq!(output, expected_output, "Outputs do not match for MNT4753.");
    }

    #[test]
    fn test_poseidon_hash_mnt6() {
        let expected_output = MNT6753Fr::new(BigInteger768([8195238283171732026, 13694263410588344527, 1885103367289967816, 17142467091011072910, 13844754763865913168, 14332001103319040991, 8911700442280604823, 6452872831806760781, 17467681867740706391, 5384727593134901588, 2343350281633109128, 244405261698305]));
        let mut poseidon_digest = MNT6PoseidonHash::init(None);
        let input = [MNT6753Fr::from_str("1").unwrap(), MNT6753Fr::from_str("2").unwrap()];

        let output = poseidon_digest
            .update(input[0])
            .update(input[1])
            .finalize();
        assert_eq!(output, expected_output, "Outputs do not match for MNT6753.");

        // Test finalize() holding the state and allowing updates in between different calls to it
        poseidon_digest
            .reset(None)
            .update(input[0].clone());
        poseidon_digest.finalize();
        poseidon_digest.update(input[1].clone());
        assert_eq!(output, poseidon_digest.finalize());

        //Test finalize() being idempotent
        assert_eq!(output, poseidon_digest.finalize());
    }

    #[test]
    fn test_poseidon_hash_mnt6_single_element() {
        let expected_output = MNT6753Fr::new(BigInteger768([9820480440897423048, 13953114361017832007, 6124683910518350026, 12198883805142820977, 16542063359667049427, 16554395404701520536, 6092728884107650560, 1511127385771028618, 14755502041894115317, 9806346309586473535, 5880260960930089738, 191119811429922]));
        let mut poseidon_digest = MNT6PoseidonHash::init(None);
        let input = MNT6753Fr::from_str("1").unwrap();
        poseidon_digest.update(input);
        let output = poseidon_digest.finalize();
        assert_eq!(output, expected_output, "Outputs do not match for MNT6753.");
    }

    #[test]
    fn test_poseidon_hash_mnt6_three_element() {
        let expected_output = MNT6753Fr::new(BigInteger768([13800884891843937189, 3814452749758584714, 14612220153016028606, 15886322817426727111, 12444362646204085653, 5214641378156871899, 4248022398370599899, 5982332416470364372, 3842784910369906888, 11445718704595887413, 5723531295320926061, 101830932453997]));
        let mut poseidon_digest = MNT6PoseidonHash::init(None);

        for i in 1..=3{
            let input = MNT6753Fr::from(i as u64);
            poseidon_digest.update(input);
        }

        let output = poseidon_digest.finalize();
        assert_eq!(output, expected_output, "Outputs do not match for MNT6753.");
    }

    use algebra::{
        fields::bn_382::{
            Fr as BN382Fr, Fq as BN382Fq,
        },
        biginteger::BigInteger384,
    };
    use crate::crh::parameters::bn382::*;

    #[test]
    fn test_poseidon_hash_bn382_fr() {
        let expected_output = BN382Fr::new(BigInteger384([5374955110091081208, 9708994766202121080, 14988884941712225891, 5210165913215347951, 13114182334648522197, 392522167697949297]));

        let mut digest = BN382FrPoseidonHash::init(None);
        digest.update(BN382Fr::from_str("1").unwrap());
        digest.update(BN382Fr::from_str("2").unwrap());
        let output = digest.finalize();

        assert_eq!(output, expected_output, "Outputs do not match for BN382Fr");
    }

    #[test]
    fn test_poseidon_hash_bn382_fq() {
        let expected_output = BN382Fq::new(BigInteger384([10704305393280846886, 13510271104066299406, 8759721062701909552, 14597420682011858322, 7770486455870140465, 1389855295932765543]));

        let mut digest = BN382FqPoseidonHash::init(None);
        digest.update(BN382Fq::from_str("1").unwrap());
        digest.update(BN382Fq::from_str("2").unwrap());
        let output = digest.finalize();

        assert_eq!(output, expected_output, "Outputs do not match for BN382Fq");
    }

}