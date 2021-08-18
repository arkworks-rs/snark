//! The scalar field of the SECG secp256k1, a 256 bit prime field.
use crate::{
    biginteger::BigInteger320 as BigInteger,
    fields::{Fp320, Fp320Parameters, FpParameters},
};

pub struct FrParameters;

pub type Fr = Fp320<FrParameters>;

// As the prime modulus has excatly 256 bit we need to oversize by another
// 64  bit register, yielding Fp320.
impl Fp320Parameters for FrParameters {}
impl FpParameters for FrParameters {
    type BigInt = BigInteger;

    /// p = 115792089237316195423570985008687907852837564279074904382605163141518161494337
    const MODULUS: BigInteger = BigInteger([
        0xbfd25e8cd0364141,
        0xbaaedce6af48a03b,
        0xfffffffffffffffe,
        0xffffffffffffffff,
        0x0,
    ]);

    const MODULUS_BITS: u32 = 256;

    const REPR_SHAVE_BITS: u32 = 64;

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    const TWO_ADICITY: u32 = 6;

    /// (p - 1)/2 = 
    /// = 57896044618658097711785492504343953926418782139537452191302581570759080747168
    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0xdfe92f46681b20a0,
        0x5d576e7357a4501d,
        0xffffffffffffffff,
        0x7fffffffffffffff,
        0x0,
    ]);
    
    /// T = (p - 1)/2^TWO_ADICITY = 
    /// = 1809251394333065553493296640760748560200586941860545380978205674086221273349
    const T: BigInteger = BigInteger([
        0xeeff497a3340d905,
        0xfaeabb739abd2280,
        0xffffffffffffffff,
        0x3ffffffffffffff,
        0x0,
    ]);

    /// (T - 1) / 2 = 
    /// = 904632097166532776746648320380374280100293470930272690489102837043110636674
    const T_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0x777fa4bd19a06c82,
        0xfd755db9cd5e9140,
        0xffffffffffffffff,
        0x1ffffffffffffff,
        0x0,
    ]);

    /// R = 2^320 mod q = 
    /// = 7976748203231275684456616952498544216114824026705293737984
    const R: BigInteger = BigInteger([
        0x0000000000000000,
        0x402da1732fc9bebf,
        0x4551231950b75fc4,
        0x0000000000000001,
        0x0,
    ]);

    /// R2 = (2^320)^2 mod q
    const R2: BigInteger = BigInteger([
        0x1e004f504dfd7f79,
        0x8fcf59774a052ea,
        0x27c4120fc94e1653,
        0x3c1a6191e5702644,
        0x0,
    ]);

    /// INV = -q^{-1} (mod 2^64)
    const INV: u64 = 5408259542528602431;

    /// GENERATOR = 5 (Montgomery rep.)
    const GENERATOR: BigInteger = BigInteger([
        0x0000000000000000,
        0x40e4273feef0b9bb,
        0x5a95af7e9394ded5,
        0x0000000000000006,
        0x0,
    ]);

    /// ROOT_OF_UNITY = GENERATOR^T (Montgomery rep.)
    const ROOT_OF_UNITY: BigInteger = BigInteger([
        0xd8fdbe6bab1827c2,
        0x55c9ffadb706d40d,
        0xdae7a417e6a410dd,
        0xd3d93fa528627010,
        0x0,
    ]);

}