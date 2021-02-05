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
use crate::{AlgebraicSponge, SpongeMode};

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

#[derive(Derivative)]
#[derivative(
Clone(bound = ""),
Debug(bound = "")
)]
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
    pub fn new(state: Vec<F>, pending: Vec<F>) -> Self {
        Self {
            state,
            pending,
            _parameters: PhantomData,
            _sbox: PhantomData
        }
    }

    #[inline]
    fn apply_permutation(&mut self, add_c2: bool) {
        for (input, state) in self.pending.iter().zip(self.state.iter_mut()) {
            *state += input;
        }
        if add_c2 { self.state[P::R] += &P::C2 };
        Self::poseidon_perm(&mut self.state);
    }

    #[inline]
    fn _finalize(&self, state: &mut Vec<F>, add_c2: bool) -> F {
        for (input, s) in self.pending.iter().zip(state.iter_mut()) {
            *s += input;
        }
        if add_c2 { state[P::R] += &P::C2; }
        Self::poseidon_perm(state);
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
        assert_eq!(P::T - P::R, 1, "The assumption that the capacity is one field element is not satisfied.");

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
            self.apply_permutation(true);
            self.pending.clear();
        }
        self
    }

    fn finalize(&self) -> Self::Data {
        if !self.pending.is_empty() {
            self._finalize(&mut self.state.clone(), true)
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

#[derive(Derivative)]
#[derivative(
Clone(bound = ""),
Debug(bound = "")
)]
pub struct PoseidonSponge<F: PrimeField, P: PoseidonParameters<Fr = F>, SB: PoseidonSBox<P>> {
    pub(crate) mode:   SpongeMode,
    pub(crate) digest: PoseidonHash<F, P, SB>,
}

impl<F, P, SB> PoseidonSponge<F, P, SB>
    where
        F: PrimeField,
        P: PoseidonParameters<Fr = F>,
        SB: PoseidonSBox<P>,
{
    pub fn new(mode: SpongeMode, digest: PoseidonHash<F, P, SB>) -> Self {
        Self {
            mode,
            digest
        }
    }

    fn clear_pending_and_apply_permutation(&mut self) {
        self.digest.apply_permutation(false);
        self.digest.pending.clear();
    }

    pub fn get_pending(&self) -> &[F] {
        &self.digest.pending
    }
}

impl<F, P, SB> AlgebraicSponge<F> for PoseidonSponge<F, P, SB>
    where
        F: PrimeField,
        P: PoseidonParameters<Fr = F>,
        SB: PoseidonSBox<P>,
{
    fn init() -> Self {
        let digest = PoseidonHash::<F, P, SB>::init(None);
        let mode = SpongeMode::Absorbing;
        Self { mode, digest }
    }

    fn get_state(&self) -> &[F] {
        self.digest.state.as_slice()
    }

    fn set_state(&mut self, state: Vec<F>) {
        assert_eq!(state.len(), P::T);
        self.digest.state = state;
    }

    fn get_mode(&self) -> &SpongeMode {
        &self.mode
    }

    fn set_mode(&mut self, mode: SpongeMode) {
        self.mode = mode;
    }

    fn absorb(&mut self, elems: Vec<F>) {

        if elems.len() > 0 {
            match self.mode {
                // If we were absorbing keep doing it
                SpongeMode::Absorbing => {
                    elems.into_iter().for_each(|f| {
                        self.digest.pending.push(f);
                        if self.digest.pending.len() == P::R {
                            // Apply a permutation when we reach rate field elements
                            self.clear_pending_and_apply_permutation();
                        }
                    })
                },

                // If we were squeezing, change the state into absorbing
                SpongeMode::Squeezing => {
                    self.mode = SpongeMode::Absorbing;
                    self.absorb(elems);
                }
            }
        }
    }

    fn squeeze(&mut self, num: usize) -> Vec<F> {
        let mut outputs = Vec::with_capacity(num);

        if num > 0 {
            match self.mode {
                SpongeMode::Absorbing => {
                    // If pending is empty and we were in absorbing, it means that a Poseidon
                    // permutation was applied just before calling squeeze(), (unless you absorbed
                    // nothing, but that is handled) therefore it's wasted to apply another
                    // permutation, and we can directly add state[0] to the outputs
                    if self.digest.pending.len() == 0 {
                        outputs.push(self.digest.state[0].clone());
                    }

                    // If pending is not empty and we were absorbing, then we need to add the
                    // pending elements to the state and then apply a permutation
                    else {
                        self.clear_pending_and_apply_permutation();
                        outputs.push(self.digest.state[0].clone());
                    }
                    self.mode = SpongeMode::Squeezing;
                    outputs.append(&mut self.squeeze(num - 1));
                },

                // If we were squeezing, then squeeze the required number of field elements
                SpongeMode::Squeezing => {
                    for _ in 0..num {
                        PoseidonHash::<F, P, SB>::poseidon_perm(&mut self.digest.state);
                        outputs.push(self.digest.state[0].clone());
                    }
                }
            }
        }
        outputs
    }
}

impl<F, P, SB> From<Vec<F>> for PoseidonSponge<F, P, SB>
    where
        F: PrimeField,
        P: PoseidonParameters<Fr = F>,
        SB: PoseidonSBox<P>,
{
    fn from(other: Vec<F>) -> Self {
        assert_eq!(other.len(), P::T);
        let digest = PoseidonHash::<F, P, SB>{
            state: other,
            pending: Vec::with_capacity(P::R),
            _parameters: PhantomData,
            _sbox: PhantomData
        };
        let mode = SpongeMode::Absorbing;
        Self { mode, digest }
    }
}

#[cfg(test)]
mod test {
    use algebra::{fields::{
        mnt4753::Fr as MNT4753Fr,
        mnt6753::Fr as MNT6753Fr,
    }, PrimeField};
    use algebra::biginteger::BigInteger768;
    use crate::crh::{
        poseidon::parameters::{
            mnt4753::{
                MNT4PoseidonHash, MNT4PoseidonSponge,
            },
            mnt6753::{
                MNT6PoseidonHash, MNT6PoseidonSponge,
            },
        },
        test::{field_based_hash_test, algebraic_sponge_test}
    };

    fn generate_inputs<F: PrimeField>(num: usize) -> Vec<F>{
        let mut inputs = Vec::with_capacity(num);
        for i in 1..=num {
            let input = F::from(i as u32);
            inputs.push(input);
        }
        inputs
    }

    #[cfg(feature = "mnt4_753")]
    #[test]
    fn test_poseidon_hash_mnt4() {
        field_based_hash_test::<MNT4PoseidonHash>(
            None,
            generate_inputs(1),
            MNT4753Fr::new(BigInteger768([10133114337753187244, 13011129467758174047, 14520750556687040981, 911508844858788085, 1859877757310385382, 9602832310351473622, 8300303689130833769, 981323167857397563, 5760566649679562093, 8644351468476031499, 10679665778836668809, 404482168782668]))
        );

        field_based_hash_test::<MNT4PoseidonHash>(
            None,
            generate_inputs(2),
            MNT4753Fr::new(BigInteger768([120759599714708995, 15132412086599307425, 1270378153255747692, 3280164418217209635, 5680179791594071572, 2475152338055275001, 9455820118751334058, 6363436228419696186, 3538976751580678769, 14987158621073838958, 10703097083485496843, 48481977539350]))
        );

        field_based_hash_test::<MNT4PoseidonHash>(
            None,
            generate_inputs(3),
            MNT4753Fr::new(BigInteger768([5991160601160569512, 9804741598782512164, 8257389273544061943, 15170134696519047397, 9908596892162673198, 7815454566429677811, 9000639780203615183, 8443915450757188195, 1987926952117715938, 17724141978374492147, 13890449093436164383, 191068391234529]))

        );

        algebraic_sponge_test::<MNT4PoseidonSponge, _>(
            generate_inputs(5),
            MNT4753Fr::new(BigInteger768([5523198498380909748, 115671896985227974, 1569791264643314974, 10995686465166995133, 13403013599916971011, 9712036026598684290, 2998254759663594264, 8111306964576313791, 14787788173217046374, 5019183223370964031, 2046072629858084037, 254417771852919])),
        );
    }

    #[cfg(feature = "mnt6_753")]
    #[test]
    fn test_poseidon_hash_mnt6() {
        let expected_output = MNT6753Fr::new(BigInteger768([9820480440897423048, 13953114361017832007, 6124683910518350026, 12198883805142820977, 16542063359667049427, 16554395404701520536, 6092728884107650560, 1511127385771028618, 14755502041894115317, 9806346309586473535, 5880260960930089738, 191119811429922]));
        field_based_hash_test::<MNT6PoseidonHash>(
            None,
            generate_inputs(1),
            expected_output
        );

        let expected_output = MNT6753Fr::new(BigInteger768([8195238283171732026, 13694263410588344527, 1885103367289967816, 17142467091011072910, 13844754763865913168, 14332001103319040991, 8911700442280604823, 6452872831806760781, 17467681867740706391, 5384727593134901588, 2343350281633109128, 244405261698305]));
        field_based_hash_test::<MNT6PoseidonHash>(
            None,
            generate_inputs(2),
            expected_output
        );

        let expected_output = MNT6753Fr::new(BigInteger768([13800884891843937189, 3814452749758584714, 14612220153016028606, 15886322817426727111, 12444362646204085653, 5214641378156871899, 4248022398370599899, 5982332416470364372, 3842784910369906888, 11445718704595887413, 5723531295320926061, 101830932453997]));
        field_based_hash_test::<MNT6PoseidonHash>(
            None,
            generate_inputs(3),
            expected_output
        );

        algebraic_sponge_test::<MNT6PoseidonSponge, _>(
            generate_inputs(5),
            MNT6753Fr::new(BigInteger768([4415375432859735514, 15277127492869068266, 6420528934628268057, 5636828761846368316, 15914884428340991861, 12211035015422291435, 9014434969167954921, 15632055196340174537, 1936740514626417030, 2037074531769378557, 3262429322356753926, 159992276980186])),
        );
    }

    use algebra::{
        fields::bn_382::{
            Fr as BN382Fr, Fq as BN382Fq,
        },
        biginteger::BigInteger384,
    };
    use crate::crh::parameters::bn382::*;

    #[cfg(feature = "bn_382")]
    #[test]
    fn test_poseidon_hash_bn382_fr() {
        let expected_output = BN382Fr::new(BigInteger384([5374955110091081208, 9708994766202121080, 14988884941712225891, 5210165913215347951, 13114182334648522197, 392522167697949297]));
        field_based_hash_test::<BN382FrPoseidonHash>(
            None,
            generate_inputs(2),
            expected_output
        );

        algebraic_sponge_test::<BN382FrPoseidonSponge, _>(
            generate_inputs(5),
            BN382Fr::new(BigInteger384([10936988494830768508, 13847454081411578749, 11142634849204675481, 15630701206177827142, 9268091518839142138, 1696903081923115093])),
        );
    }

    #[cfg(feature = "bn_382")]
    #[test]
    fn test_poseidon_hash_bn382_fq() {
        let expected_output = BN382Fq::new(BigInteger384([10704305393280846886, 13510271104066299406, 8759721062701909552, 14597420682011858322, 7770486455870140465, 1389855295932765543]));
        field_based_hash_test::<BN382FqPoseidonHash>(
            None,
            generate_inputs(2),
            expected_output
        );

        algebraic_sponge_test::<BN382FqPoseidonSponge, _>(
            generate_inputs(5),
            BN382Fq::new(BigInteger384([1720359977283232686, 13054139547712885489, 5187034200847661218, 4192055198669470320, 9342360683435217824, 2331312681920757379])),
        );
    }

    use algebra::{
        fields::tweedle::{
            fr::Fr as tweedleFr, fq::Fq as tweedleFq,
        },
        biginteger::BigInteger256,
    };
    use crate::crh::parameters::tweedle::*;

    #[cfg(feature = "tweedle")]
    #[test]
    fn test_poseidon_hash_tweedle_fr() {
        let expected_output = tweedleFr::new(BigInteger256([10853721058648678319, 16685221779982166148, 7657542961996224896, 317411048550701368]));
        field_based_hash_test::<TweedleFrPoseidonHash>(
            None,
            generate_inputs(2),
            expected_output
        );

        algebraic_sponge_test::<TweedleFrPoseidonSponge, _>(
            generate_inputs(5),
            tweedleFr::new(BigInteger256([7978538512485120357, 2855094910189988323, 8391520218117106983, 4530816604245346005]))
        );
    }

    #[cfg(feature = "tweedle")]
    #[test]
    fn test_poseidon_hash_tweedle_fq() {
        let expected_output = tweedleFq::new(BigInteger256([9400878458790897114, 15068972336691613232, 13673927707991766433, 1032370839092625161]));
        field_based_hash_test::<TweedleFqPoseidonHash>(
            None,
            generate_inputs(2),
            expected_output
        );
        algebraic_sponge_test::<TweedleFqPoseidonSponge, _>(
            generate_inputs(5),
            tweedleFq::new(BigInteger256([6129396713371884441, 13029129558030886129, 9912298868009899566, 385156517461505756]))
        );
    }
}