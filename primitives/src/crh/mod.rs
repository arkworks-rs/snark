use algebra::bytes::ToBytes;
use rand::Rng;
use std::hash::Hash;

pub mod bowe_hopwood;
pub mod injective_map;
pub mod pedersen;

use crate::Error;


pub trait FixedLengthCRH {
    const INPUT_SIZE_BITS: usize;
    type Output: ToBytes + Clone + Eq + Hash + Default;
    type Parameters: Clone + Default;

    fn setup<R: Rng>(r: &mut R) -> Result<Self::Parameters, Error>;
    fn evaluate(parameters: &Self::Parameters, input: &[u8]) -> Result<Self::Output, Error>;
}

//Temporary mock of Poseidon interfaces

use algebra::Field;

pub trait FieldBasedHashParameters: Sized + Clone{
    type Fr: Field;
}

pub trait FieldBasedHash {
    type Data: Field;
    type Parameters: FieldBasedHashParameters<Fr = Self::Data>;

    fn evaluate(input: &[Self::Data]) -> Result<Self::Data, Error>;
}

pub trait PoseidonParameters: FieldBasedHashParameters {
    //Constants here
}

use algebra::fields::mnt4753::Fr as MNT4753Fr;
use algebra::fields::mnt6753::Fr as MNT6753Fr;
use std::marker::PhantomData;

#[derive(Derivative)]
#[derivative(Clone)]
pub struct MNT4HashParameters;

impl FieldBasedHashParameters for MNT4HashParameters{
    type Fr = MNT4753Fr;
}

impl PoseidonParameters for MNT4HashParameters {}

#[derive(Derivative)]
#[derivative(Clone)]
pub struct MNT6HashParameters;

impl FieldBasedHashParameters for MNT6HashParameters{
    type Fr = MNT6753Fr;
}

impl PoseidonParameters for MNT6HashParameters {}

pub struct PoseidonHash<F: Field, P: PoseidonParameters<Fr = F>>{
    _field:      PhantomData<F>,
    _parameters: PhantomData<P>,
}

impl<F: Field, P: PoseidonParameters<Fr = F>> FieldBasedHash for PoseidonHash<F, P>{
    type Data = F;
    type Parameters = P;

    fn evaluate(input: &[F]) -> Result<F, Error> {
        //Dummy impl, just for test
        let mut res = F::zero();
        input.iter().for_each(|f| res += f);
        Ok(res)
    }
}

pub type MNT4PoseidonHash = PoseidonHash<MNT4753Fr, MNT4HashParameters>;
pub type MNT6PoseidonHash = PoseidonHash<MNT6753Fr, MNT6HashParameters>;
