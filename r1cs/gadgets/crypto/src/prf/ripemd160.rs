//! Circuits for the [RIPEMD160] hash function
//!
//! [RIPEMD160]: https://tools.ietf.org/html/rfc6234

use r1cs_std::boolean::Boolean;
use r1cs_std::eq::MultiEq;
use r1cs_std::uint32::UInt32;
use r1cs_std::uint8::UInt8;
use r1cs_core::{ConstraintSystem, SynthesisError};
use algebra::PrimeField;

use crate::sha256::{sha256_ch_boolean, triop};

/// Outputs K[round_idx] and K'[round_idx]
fn get_round_constants(round_idx: usize) -> (UInt32, UInt32) {
    let (k, k_prime): (u32, u32) = match round_idx {
        round_idx if round_idx <= 15 =>                    (0x00000000, 0x50a28be6),
        round_idx if round_idx >= 16 && round_idx <= 31 => (0x5a827999, 0x5c4dd124),
        round_idx if round_idx >= 32 && round_idx <= 47 => (0x6ed9eba1, 0x6d703ef3),
        round_idx if round_idx >= 48 && round_idx <= 63 => (0x8f1bbcdc, 0x7a6d76e9),
        round_idx if round_idx >= 64 && round_idx <= 79 => (0xa953fd4e, 0x00000000),
        _ => unreachable!()
    };
    (UInt32::constant(k), UInt32::constant(k_prime))
}

#[allow(clippy::unreadable_literal)]
const IV: [u32; 5] = [
    0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476, 0xc3d2e1f0
];

/// Amount for first rotate left
const S: [usize; 80] = [
    11, 14, 15, 12, 5, 8, 7, 9, 11, 13, 14, 15, 6, 7, 9, 8,
    7, 6, 8, 13, 11, 9, 7, 15, 7, 12, 15, 9, 11, 7, 13, 12,
    11, 13, 6, 7, 14, 9, 13, 15, 14, 8, 13, 6, 5, 12, 7, 5,
    11, 12, 14, 15, 14, 15, 9, 8, 9, 14, 5, 6, 8, 6, 5, 12,
    9, 15, 5, 11, 6, 8, 13, 12, 5, 12, 13, 14, 11, 8, 5, 6
];

/// Amount for second rotate left
const S_PRIME: [usize; 80] = [
    8, 9, 9, 11, 13, 15, 15, 5, 7, 7, 8, 11, 14, 14, 12, 6,
    9, 13, 15, 7, 12, 8, 9, 11, 7, 7, 12, 7, 6, 15, 13, 11,
    9, 7, 15, 11, 8, 6, 6, 14, 12, 13, 5, 14, 13, 13, 7, 5,
    15, 5, 8, 11, 14, 14, 6, 14, 6, 9, 12, 9, 12, 5, 15, 8,
    8, 5, 12, 9, 12, 5, 14, 6, 8, 13, 6, 5, 15, 13, 11, 11
];

/// First message word selection
const R: [usize; 80] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
    7, 4, 13, 1, 10, 6, 15, 3, 12, 0, 9, 5, 2, 14, 11, 8,
    3, 10, 14, 4, 9, 15, 8, 1, 2, 7, 0, 6, 13, 11, 5, 12,
    1, 9, 11, 10, 0, 8, 12, 4, 13, 3, 7, 15, 14, 5, 6, 2,
    4, 0, 5, 9, 7, 12, 2, 10, 14, 1, 3, 8, 11, 6, 15, 13
];

/// Second message word selection
const R_PRIME: [usize; 80] = [
    5, 14, 7, 0, 9, 2, 11, 4, 13, 6, 15, 8, 1, 10, 3, 12,
    6, 11, 3, 7, 0, 13, 5, 10, 14, 15, 8, 12, 4, 9, 1, 2,
    15, 5, 1, 3, 7, 14, 6, 9, 11, 8, 12, 2, 10, 0, 4, 13,
    8, 6, 4, 1, 3, 11, 15, 0, 5, 12, 2, 13, 9, 7, 10, 14,
    12, 15, 10, 4, 1, 5, 8, 7, 6, 2, 13, 14, 0, 3, 9, 11
];

/// Apply the proper round function, selected with round_idx, to inputs a, b and c
fn apply_round_function<ConstraintF, CS>(
    cs:        CS,
    a:         &UInt32,
    b:         &UInt32,
    c:         &UInt32,
    round_idx: usize,
) -> Result<UInt32, SynthesisError>
    where
        ConstraintF: PrimeField,
        CS: ConstraintSystem<ConstraintF>
{
    // Select and return the proper round function according to round idx
    let (tri_fn, circ_fn): (fn(u32, u32, u32) -> u32, fn(&mut CS, usize, &Boolean, &Boolean, &Boolean) -> Result<Boolean, SynthesisError>) = match round_idx {
        // f(j, a, b, c) = a XOR b XOR c (0 <= j <= 15)
        round_idx if round_idx <= 15 => {
            (
                |a, b, c| a ^ b ^ c,
                |cs, i, a, b, c| {
                    let result = Boolean::xor(cs.ns(|| format!("A XOR B {}", i)), &a, &b)?;
                    Boolean::xor(cs.ns(|| format!("A XOR B XOR C {}", i)), &result,&c)
                }
            )
        },
        // f(j, a, b, c) = (a AND b) OR (NOT(a) AND c)  (16 <= j <= 31)
        // Note: It's the same as the sha256_ch function e.g. (a AND b) XOR (NOT(a) AND c)
        // since the two logical expressions have the same truth table
        round_idx if round_idx >= 16 && round_idx <= 31 => {
            (
                |a, b, c| (a & b) | ((!a) & c),
                |cs, i, a, b, c| sha256_ch_boolean(cs.ns(|| format!("ch {}", i)), a, b, c)
            )
        },
        // f(j, a, b, c) = (a OR NOT(b)) XOR c (32 <= j <= 47)
        round_idx if round_idx >= 32 && round_idx <= 47 => {
            (
                |a, b, c| (a | (!b)) ^ c,
                |cs, i, a, b, c| {
                    let t = Boolean::or(cs.ns(|| format!("A OR NOT B {}", i)), &a, &b.not())?;
                    Boolean::xor(cs.ns(|| format!("(A OR NOT B) XOR C {}", i)), &t, &c)
                }
            )
        },
        // f(j, a, b, c) = (c AND a) OR (NOT(c) AND b)  (48 <= j <= 63)
        // Note: It's the same as the second round function, but with permuted variables
        round_idx if round_idx >= 48 && round_idx <= 63 => {
            (
                |a, b, c| (c & a) | ((!c) & b),
                |cs, i, a, b, c| sha256_ch_boolean(cs.ns(|| format!("permuted ch {}", i)), c, a, b)
            )
        },
        // f(j, a, b, c) = (b OR NOT(c)) XOR a (64 <= j <= 79)
        // Note: It's the same as the third round function but with permuted variables
        round_idx if round_idx >= 64 && round_idx <= 79 => {
            (
                |a, b, c| (b | (!c)) ^ a,
                |cs, i, a, b, c| {
                    let t = Boolean::or(cs.ns(|| format!("B OR NOT C {}", i)), &b, &c.not())?;
                    Boolean::xor(cs.ns(|| format!("(B OR NOT C) XOR a {}", i)), &t, &a)
                }
            )
        },
        _ => unreachable!()
    };

    triop(cs, a, b, c, tri_fn, circ_fn)
}

pub fn ripemd160_block_no_padding<ConstraintF, CS>(
    mut cs: CS,
    input: &[Boolean],
) -> Result<Vec<Boolean>, SynthesisError>
    where
        ConstraintF: PrimeField,
        CS: ConstraintSystem<ConstraintF>,
{
    assert_eq!(input.len(), 512);

    Ok(
        ripemd160_compression_function(&mut cs, &input, &get_ripemd160_iv())?
            .into_iter()
            .flat_map(|e| e.to_bits_le())
            .collect(),
    )
}

pub fn ripemd160<ConstraintF, CS>(mut cs: CS, input: &[Boolean]) -> Result<Vec<Boolean>, SynthesisError>
    where
        ConstraintF: PrimeField,
        CS: ConstraintSystem<ConstraintF>,
{
    assert!(input.len() % 8 == 0);

    let mut padded = input.to_vec();
    let plen = padded.len() as u64;
    // append 0x80
    padded.append(&mut UInt8::constant(1).into_bits_be());

    // append K '0' bits, where K is the minimum number >= 0 such that L + 1 + K + 64 is a multiple of 512
    while (padded.len() + 64) % 512 != 0 {
        padded.push(Boolean::constant(false));
    }
    // append L as a 64-bit little-endian integer, making the total post-processed length a multiple of 512 bits
    for b in (0..64).map(|i| (plen >> i) & 1 == 1) {
        padded.push(Boolean::constant(b));
    }
    assert!(padded.len() % 512 == 0);

    let mut cur = get_ripemd160_iv();
    for (i, block) in padded.chunks(512).enumerate() {
        cur = ripemd160_compression_function(cs.ns(|| format!("block {}", i)), block, &cur)?;
    }

    Ok(cur.into_iter().flat_map(|e| e.to_bits_le()).collect())
}

fn get_ripemd160_iv() -> Vec<UInt32> {
    IV.iter().map(|&v| UInt32::constant(v)).collect()
}

fn ripemd160_compression_function<ConstraintF, CS>(
    cs: CS,
    input: &[Boolean],
    current_hash_value: &[UInt32],
) -> Result<Vec<UInt32>, SynthesisError>
    where
        ConstraintF: PrimeField,
        CS: ConstraintSystem<ConstraintF>,
{
    assert_eq!(input.len(), 512);
    assert_eq!(current_hash_value.len(), 5);

    let mut a = current_hash_value[0].clone();
    let mut b = current_hash_value[1].clone();
    let mut c = current_hash_value[2].clone();
    let mut d = current_hash_value[3].clone();
    let mut e = current_hash_value[4].clone();
    let mut a_prime = a.clone();
    let mut b_prime = b.clone();
    let mut c_prime = c.clone();
    let mut d_prime = d.clone();
    let mut e_prime = e.clone();

    let x = input
        .chunks(32)
        .map(|e| UInt32::from_bits_le(e))
        .collect::<Vec<_>>();

    let mut cs = MultiEq::new(cs);

    for i in 0..80 {

        let cs = &mut cs.ns(|| format!("compression round {}", i));
        let mut t = {
            let f = apply_round_function(cs.ns(|| format!("first round function {}", i)), &b, &c, &d, i)?;
            let selected_input_word = x[R[i]].clone();
            let result = UInt32::addmany(
                cs.ns(|| format!("fist rotl(a + f + x + k) {}", i)),
                &[a, f, selected_input_word, get_round_constants(i).0]
            )?.rotl(S[i]);
            UInt32::addmany(
                cs.ns(|| format!("compute first T {}", i)),
                &[result, e.clone()]
            )
        }?;

        a = e;
        e = d;
        d = c.rotl(10);
        c = b;
        b = t;

        t = {
            let f = apply_round_function(
                cs.ns(|| format!("second round function {}", i)),
                &b_prime,
                &c_prime,
                &d_prime,
                79 - i
            )?;
            let selected_input_word = x[R_PRIME[i]].clone();
            let result = UInt32::addmany(
                cs.ns(|| format!("second rotl(a + f + x + k) {}", i)),
                &[a_prime, f, selected_input_word, get_round_constants(i).1]
            )?.rotl(S_PRIME[i]);
            UInt32::addmany(
                cs.ns(|| format!("compute second T {}", i)),
                &[result, e_prime.clone()]
            )
        }?;

        a_prime = e_prime;
        e_prime = d_prime;
        d_prime = c_prime.rotl(10);
        c_prime = b_prime;
        b_prime = t;
    }

    // Final addition round

    let h0 = UInt32::addmany(
        cs.ns(|| "compute final t"),
        &[current_hash_value[1].clone(), c, d_prime]
    )?;

    let h1 = UInt32::addmany(
        cs.ns(|| "new h1"),
        &[current_hash_value[2].clone(), d, e_prime],
    )?;

    let h2 = UInt32::addmany(
        cs.ns(|| "new h2"),
        &[current_hash_value[3].clone(), e, a_prime],
    )?;

    let h3 = UInt32::addmany(
        cs.ns(|| "new h3"),
        &[current_hash_value[4].clone(), a, b_prime],
    )?;

    let h4 = UInt32::addmany(
        cs.ns(|| "new h4"),
        &[current_hash_value[0].clone(), b, c_prime],
    )?;

    Ok(vec![h0, h1, h2, h3, h4])
}


#[cfg(test)]
mod test {
    use super::*;
    use r1cs_std::{alloc::AllocGadget, boolean::AllocatedBit, test_constraint_system::TestConstraintSystem};
    use algebra::fields::bls12_381::Fr;

    use rand::{RngCore, SeedableRng};
    use rand_xorshift::XorShiftRng;

    #[test]
    fn test_full_block() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x3d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        let iv = get_ripemd160_iv();

        let mut cs = TestConstraintSystem::<Fr>::new();
        let input_bits: Vec<_> = (0..512)
            .map(|i| {
                Boolean::from(
                    AllocatedBit::alloc(
                        cs.ns(|| format!("input bit {}", i)),
                        || Ok(rng.next_u32() % 2 != 0),
                    )
                        .unwrap(),
                )
            })
            .collect();

        ripemd160_compression_function(cs.ns(|| "ripemd160"), &input_bits, &iv).unwrap();

        assert!(cs.is_satisfied());
        assert_eq!(cs.num_constraints() - 512, 18797);
    }

     #[test]
     fn native_test() {
         use ripemd160::{Digest, Ripemd160};

         let mut rng = XorShiftRng::from_seed([
             0x59, 0x62, 0xbe, 0x3d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
             0xbc, 0xe5,
         ]);

         for input_len in (0..32).chain((32..256).filter(|a| a % 8 == 0)) {
             let mut h = Ripemd160::new();
             let data: Vec<u8> = (0..input_len).map(|_| rng.next_u32() as u8).collect();
             h.update(&data);
             let hash_result = h.finalize();

             let mut cs = TestConstraintSystem::<Fr>::new();
             let mut input_bits = vec![];

             for (byte_i, input_byte) in data.into_iter().enumerate() {
                 for bit_i in 0..8 {
                     let cs = cs.ns(|| format!("input bit {} {}", byte_i, bit_i));

                     input_bits.push(
                         AllocatedBit::alloc(cs, || Ok((input_byte >> bit_i) & 1u8 == 1u8))
                             .unwrap()
                             .into(),
                     );
                 }
             }

             let r = ripemd160(&mut cs, &input_bits).unwrap();

             assert!(cs.is_satisfied());

             let mut s = hash_result
                 .iter()
                 .flat_map(|&byte| (0..8).map(move |i| (byte >> i) & 1u8 == 1u8));

             for b in r {
                 match b {
                     Boolean::Is(b) => {
                         assert!(s.next().unwrap() == b.get_value().unwrap());
                     }
                     Boolean::Not(b) => {
                         assert!(s.next().unwrap() != b.get_value().unwrap());
                     }
                     Boolean::Constant(b) => {
                         assert!(input_len == 0);
                         assert!(s.next().unwrap() == b);
                     }
                 }
             }
         }
     }

    #[test]
    fn compare_against_test_vectors() {
        let test_inputs = [
            "",
            "a",
            "abc",
            "message digest",
            "abcdefghijklmnopqrstuvwxyz",
            "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq",
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789",
            &(0..8).map(|_| "1234567890").collect::<String>(),
            "The quick brown fox jumps over the lazy dog",
            "The quick brown fox jumps over the lazy cog",
        ];

        let test_outputs = [
            "9c1185a5c5e9fc54612808977ee8f548b2258d31",
            "0bdc9d2d256b3ee9daae347be6f4dc835a467ffe",
            "8eb208f7e05d987a9b044a8e98c6b087f15a0bfc",
            "5d0689ef49d2fae572b881b123a85ffa21595f36",
            "f71c27109c692c1b56bbdceb5b9d2865b3708dbc",
            "12a053384a9c0c88e405a06c27dcf49ada62eb2b",
            "b0e20b6e3116640286ed3a87a5713079b21f5189",
            "9b752e45573d4b39f4dbd3323cab82bf63326bfb",
            "37f332f68db77bd9d7edd4969571ad671cf9dd3b",
            "132072df690933835eb8b6ad0b77e7b6f14acad7",
        ];

        for (test_input, test_output) in test_inputs.iter().zip(test_outputs.iter()) {
            let mut cs = TestConstraintSystem::<Fr>::new();
            let mut input_bits = vec![];

            for (byte_i, input_byte) in test_input.as_bytes().into_iter().enumerate() {
                for bit_i in 0..8 {
                    let cs = cs.ns(|| format!("input bit {} {}", byte_i, bit_i));

                    input_bits.push(
                        AllocatedBit::alloc(cs, || Ok((input_byte >> bit_i) & 1u8 == 1u8))
                            .unwrap()
                            .into(),
                    );
                }
            }

            let r = ripemd160(&mut cs, &input_bits).unwrap();

            assert!(cs.is_satisfied());

            let expected_output = hex::decode(test_output).unwrap();
            let mut s = expected_output
                .iter()
                .flat_map(|&byte| (0..8).map(move |i| (byte >> i) & 1u8 == 1u8));

            for b in r {
                match b {
                    Boolean::Is(b) => {
                        assert!(s.next().unwrap() == b.get_value().unwrap());
                    }
                    Boolean::Not(b) => {
                        assert!(s.next().unwrap() != b.get_value().unwrap());
                    }
                    Boolean::Constant(b) => {
                        assert!(input_bits.len() == 0);
                        assert!(s.next().unwrap() == b);
                    }
                }
            }
        }
    }
}