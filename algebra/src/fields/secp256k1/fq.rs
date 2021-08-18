//! The base field of the SECG secp256k1, a 256 bit prime field. 
use crate::{
    biginteger::BigInteger320 as BigInteger,
    fields::{FpParameters, Fp320, Fp320Parameters},
};

pub type Fq = Fp320<FqParameters>;

pub struct FqParameters;

// As the prime modulus has excatly 256 bit we need to oversize by another
// 64  bit register, yielding Fp320.
impl Fp320Parameters for FqParameters {}
impl FpParameters for FqParameters {
    type BigInt = BigInteger;

    /// p = 2^256 - 2^32 - 2^9 - 2^8 - 2^7 - 2^6 - 2^4 - 1 =
    ///  = 115792089237316195423570985008687907853269984665640564039457584007908834671663
    const MODULUS: BigInteger = BigInteger([
        0xfffffffefffffc2f,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0x0,
    ]);

    const MODULUS_BITS: u32 = 256;

    const REPR_SHAVE_BITS: u32 = 64;

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    const TWO_ADICITY: u32 = 1;

    /// (p-1)/2 = 
    ///  = 57896044618658097711785492504343953926634992332820282019728792003954417335831
    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0xffffffff7ffffe17,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0x7fffffffffffffff,
        0x0,
    ]);

    /// T = (p-1)/2^twoadicity = MODULUS_MINUS_ONE_DIV_TWO
    const T: BigInteger = BigInteger([
        0xffffffff7ffffe17,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0x7fffffffffffffff,
        0x0,
    ]);

    /// (T - 1) / 2 = 
    /// = 28948022309329048855892746252171976963317496166410141009864396001977208667915
    const T_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0xffffffffbfffff0b,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0x3fffffffffffffff,
        0x0,
    ]);

    /// Montgomery constant R = 2^320 mod p = 
    ///  = 79228180536733297607775879168
    const R: BigInteger = BigInteger([
        0x0000000000000000,
        0x00000001000003d1,
        0x0000000000000000,
        0x0000000000000000,
        0x0,
    ]);

    /// R2 = (2^320)^2 mod p = 18446752466076602529
    const R2: BigInteger = BigInteger([
        0x0000000000000000,
        0x0000000000000000,
        0x000007a2000e90a1,
        0x0000000000000001,
        0x0,
    ]);

    /// INV = -p^{-1} (mod 2^64)
    const INV: u64 = 15580212934572586289;

    /// GENERATOR = 5 (Montgomery rep.)
    const GENERATOR: BigInteger = BigInteger([
        0x0000000000000000,
        0x0000000500001315,
        0x0000000000000000,
        0x0000000000000000,
        0x0,
    ]);

    /// ROOT_OF_UNITY = GENERATOR^T (Montgomery rep.) 
    /// = -1 mod p
    const ROOT_OF_UNITY: BigInteger = BigInteger([
        0xfffffffefffffc2f,
        0xfffffffefffffc2e,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0x0,
    ]);

}