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
use crate::AlgebraicSponge;

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

impl<F, P, SB> AlgebraicSponge<F> for PoseidonHash<F, P, SB>
    where
        F: PrimeField,
        P: PoseidonParameters<Fr = F>,
        SB: PoseidonSBox<P>,
{
    fn new() -> Self {
        Self::init(None)
    }

    fn absorb(&mut self, elems: Vec<F>) {
        elems.into_iter().for_each(|f| {
            self.pending.push(f);
            if self.pending.len() == P::R {
                self.apply_permutation(false);
                self.pending.clear();
            }
        })
    }

    fn squeeze(&self, num: usize) -> Vec<F> {
        let mut state = self.state.clone();
        for (input, s) in self.pending.iter().zip(state.iter_mut()) {
            *s += input;
        }
        let mut output = Vec::with_capacity(num);
        for _ in 0..num {
            Self::poseidon_perm(&mut state);
            output.push(state[0].clone());
        }
        output
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
            mnt4753::MNT4PoseidonHash,
            mnt6753::MNT6PoseidonHash,
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

        algebraic_sponge_test::<MNT4PoseidonHash, _>(
            generate_inputs(5),
            vec![
                MNT4753Fr::new(BigInteger768([5523198498380909748, 115671896985227974, 1569791264643314974, 10995686465166995133, 13403013599916971011, 9712036026598684290, 2998254759663594264, 8111306964576313791, 14787788173217046374, 5019183223370964031, 2046072629858084037, 254417771852919])),
                MNT4753Fr::new(BigInteger768([16127953875295592987, 11894863219027726126, 17701631661641187911, 14870672795897948467, 12477792521878046393, 7933118744616181152, 16504185330957910903, 4440754255912752398, 2147815826017548215, 15556058995384616380, 7316328326052908756, 488587201854075])),
                MNT4753Fr::new(BigInteger768([17093615125500557324, 5853740299004434536, 14185517698974334849, 15081129818867394374, 4393775450585412195, 15317350057972385858, 602096567370343479, 171645566249469757, 213587553721144028, 16770187529782220214, 10656253922486106039, 183907258171817])),
                MNT4753Fr::new(BigInteger768([12712297803466521417, 5718281485803573501, 13963200307320285636, 13279387328772122744, 11400771971860800684, 18129978631072081884, 16152730030736233432, 4443584430842639200, 10168547976596639262, 8852760702638429316, 6245616001490582246, 290914933918621])),
                MNT4753Fr::new(BigInteger768([10000480287694253354, 5873338108557543580, 14250396374254092226, 15940814885061001617, 7624733801679967429, 10913153568558026344, 10086885331318146225, 5264350148957133820, 9338271530149930761, 15261863060094932371, 15830626088040155644, 456905335294658]))
            ]
        );
    }

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

        algebraic_sponge_test::<MNT6PoseidonHash, _>(
            generate_inputs(5),
            vec![
                MNT6753Fr::new(BigInteger768([4415375432859735514, 15277127492869068266, 6420528934628268057, 5636828761846368316, 15914884428340991861, 12211035015422291435, 9014434969167954921, 15632055196340174537, 1936740514626417030, 2037074531769378557, 3262429322356753926, 159992276980186])),
                MNT6753Fr::new(BigInteger768([893755365716394547, 9594068362757650017, 10451094618402356015, 6337497625726267689, 4160278191670220893, 6010491150125800070, 828252792235881440, 2673899143833542881, 13863607618619424503, 10731000306034932792, 12767597931406600360, 152088456622435])),
                MNT6753Fr::new(BigInteger768([2814201016530476081, 8511720540093155210, 4852589771558571678, 13110129496171929816, 4735221329350745796, 12663409230803838109, 11568212364226125976, 1693070022169342918, 2888127900540921339, 4204338349274483900, 15311382501005232626, 438220753471979])),
                MNT6753Fr::new(BigInteger768([1775158153332815483, 2026417705400082199, 11123169942942274855, 14330277722886619530, 7836191979464539572, 4865965751577341882, 16998150036867812011, 2739207342391276882, 4247579097230352822, 6464967486406868296, 10032824663714337053, 463362661909267])),
                MNT6753Fr::new(BigInteger768([17595647052173133652, 5444801609492079542, 17111761943624358533, 9849391399095026693, 263233265735065436, 1281866738651078101, 2333944581006551597, 13500573177506341637, 10203063908384511927, 15054379420590094900, 13913006388276371623, 451783078041284]))
            ]
        );
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
        field_based_hash_test::<BN382FrPoseidonHash>(
            None,
            generate_inputs(2),
            expected_output
        );

        algebraic_sponge_test::<BN382FrPoseidonHash, _>(
            generate_inputs(5),
            vec![
                BN382Fr::new(BigInteger384([10936988494830768508, 13847454081411578749, 11142634849204675481, 15630701206177827142, 9268091518839142138, 1696903081923115093])),
                BN382Fr::new(BigInteger384([696581192457039113, 7735026456746748023, 2858087404435942973, 8551026292503659964, 10649071471792522052, 1622107112684056397])),
                BN382Fr::new(BigInteger384([13847543644949650964, 4311663852800204305, 4992676739644232980, 4027896082677996619, 12253082643656956806, 534636992911167858])),
                BN382Fr::new(BigInteger384([13260537455622274637, 17371511872105869645, 1623345855391475943, 1872797973978834237, 7609952746869694608, 847112821092614604])),
                BN382Fr::new(BigInteger384([2836595273011902376, 5748666575747773318, 13745164138914401597, 11231420449817912461, 3298276789500282573, 1513701281573365844]))
            ]
        );
    }

    #[test]
    fn test_poseidon_hash_bn382_fq() {
        let expected_output = BN382Fq::new(BigInteger384([10704305393280846886, 13510271104066299406, 8759721062701909552, 14597420682011858322, 7770486455870140465, 1389855295932765543]));
        field_based_hash_test::<BN382FqPoseidonHash>(
            None,
            generate_inputs(2),
            expected_output
        );

        algebraic_sponge_test::<BN382FqPoseidonHash, _>(
            generate_inputs(5),
            vec![
                BN382Fq::new(BigInteger384([1720359977283232686, 13054139547712885489, 5187034200847661218, 4192055198669470320, 9342360683435217824, 2331312681920757379])),
                BN382Fq::new(BigInteger384([4268335443493036543, 5752606093574636763, 7217511812669325408, 5182020843125294140, 5400661559267345372, 2261580635543305540])),
                BN382Fq::new(BigInteger384([1487056021032845247, 4971413697289841114, 829985713515651434, 17668516547215117574, 6183480290866190120, 2494711740844925696])),
                BN382Fq::new(BigInteger384([12440731074295045789, 13776548474636848476, 9611723257742392848, 14113646477500207966, 2479383220929757679, 2462018425115041259])),
                BN382Fq::new(BigInteger384([13447602317280700149, 7333172335365879673, 3401105980077026956, 7168529825556409731, 18013175927097546368, 369466729188096152])),
            ]
        );
    }
}