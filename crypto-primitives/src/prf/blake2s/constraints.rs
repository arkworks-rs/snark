use algebra_core::PrimeField;
use r1cs_core::{ConstraintSystemRef, Namespace, SynthesisError};

use crate::{prf::PRFGadget, Vec};
use r1cs_std::prelude::*;

use core::borrow::Borrow;

// 2.1.  Parameters
// The following table summarizes various parameters and their ranges:
//               | BLAKE2b          | BLAKE2s          |
// --------------+------------------+------------------+
// Bits in word  | w = 64           | w = 32           |
// Rounds in F   | r = 12           | r = 10           |
// Block bytes   | bb = 128         | bb = 64          |
// Hash bytes    | 1 <= nn <= 64    | 1 <= nn <= 32    |
// Key bytes     | 0 <= kk <= 64    | 0 <= kk <= 32    |
// Input bytes   | 0 <= ll < 2**128 | 0 <= ll < 2**64  |
// --------------+------------------+------------------+
// G Rotation    | (R1, R2, R3, R4) | (R1, R2, R3, R4) |
// constants =   | (32, 24, 16, 63) | (16, 12,  8,  7) |
// --------------+------------------+------------------+
//

const R1: usize = 16;
const R2: usize = 12;
const R3: usize = 8;
const R4: usize = 7;

// Round     |  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 |
// ----------+-------------------------------------------------+
// SIGMA[0]  |  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 |
// SIGMA[1]  | 14 10  4  8  9 15 13  6  1 12  0  2 11  7  5  3 |
// SIGMA[2]  | 11  8 12  0  5  2 15 13 10 14  3  6  7  1  9  4 |
// SIGMA[3]  |  7  9  3  1 13 12 11 14  2  6  5 10  4  0 15  8 |
// SIGMA[4]  |  9  0  5  7  2  4 10 15 14  1 11 12  6  8  3 13 |
// SIGMA[5]  |  2 12  6 10  0 11  8  3  4 13  7  5 15 14  1  9 |
// SIGMA[6]  | 12  5  1 15 14 13  4 10  0  7  6  3  9  2  8 11 |
// SIGMA[7]  | 13 11  7 14 12  1  3  9  5  0 15  4  8  6  2 10 |
// SIGMA[8]  |  6 15 14  9 11  3  0  8 12  2 13  7  1  4 10  5 |
// SIGMA[9]  | 10  2  8  4  7  6  1  5 15 11  9 14  3 12 13  0 |
// ----------+-------------------------------------------------+
//

const SIGMA: [[usize; 16]; 10] = [
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
    [14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3],
    [11, 8, 12, 0, 5, 2, 15, 13, 10, 14, 3, 6, 7, 1, 9, 4],
    [7, 9, 3, 1, 13, 12, 11, 14, 2, 6, 5, 10, 4, 0, 15, 8],
    [9, 0, 5, 7, 2, 4, 10, 15, 14, 1, 11, 12, 6, 8, 3, 13],
    [2, 12, 6, 10, 0, 11, 8, 3, 4, 13, 7, 5, 15, 14, 1, 9],
    [12, 5, 1, 15, 14, 13, 4, 10, 0, 7, 6, 3, 9, 2, 8, 11],
    [13, 11, 7, 14, 12, 1, 3, 9, 5, 0, 15, 4, 8, 6, 2, 10],
    [6, 15, 14, 9, 11, 3, 0, 8, 12, 2, 13, 7, 1, 4, 10, 5],
    [10, 2, 8, 4, 7, 6, 1, 5, 15, 11, 9, 14, 3, 12, 13, 0],
];

// 3.1.  Mixing Function G
// The G primitive function mixes two input words, "x" and "y", into
// four words indexed by "a", "b", "c", and "d" in the working vector
// v[0..15].  The full modified vector is returned.  The rotation
// constants (R1, R2, R3, R4) are given in Section 2.1.
// FUNCTION G( v[0..15], a, b, c, d, x, y )
// |
// |   v[a] := (v[a] + v[b] + x) mod 2**w
// |   v[d] := (v[d] ^ v[a]) >>> R1
// |   v[c] := (v[c] + v[d])     mod 2**w
// |   v[b] := (v[b] ^ v[c]) >>> R2
// |   v[a] := (v[a] + v[b] + y) mod 2**w
// |   v[d] := (v[d] ^ v[a]) >>> R3
// |   v[c] := (v[c] + v[d])     mod 2**w
// |   v[b] := (v[b] ^ v[c]) >>> R4
// |
// |   RETURN v[0..15]
// |
// END FUNCTION.
//

fn mixing_g<ConstraintF: PrimeField>(
    v: &mut [UInt32<ConstraintF>],
    a: usize,
    b: usize,
    c: usize,
    d: usize,
    x: &UInt32<ConstraintF>,
    y: &UInt32<ConstraintF>,
) -> Result<(), SynthesisError> {
    v[a] = UInt32::addmany(&[v[a].clone(), v[b].clone(), x.clone()])?;
    v[d] = v[d].xor(&v[a])?.rotr(R1);
    v[c] = UInt32::addmany(&[v[c].clone(), v[d].clone()])?;
    v[b] = v[b].xor(&v[c])?.rotr(R2);
    v[a] = UInt32::addmany(&[v[a].clone(), v[b].clone(), y.clone()])?;
    v[d] = v[d].xor(&v[a])?.rotr(R3);
    v[c] = UInt32::addmany(&[v[c].clone(), v[d].clone()])?;
    v[b] = v[b].xor(&v[c])?.rotr(R4);

    Ok(())
}

// 3.2.  Compression Function F
// Compression function F takes as an argument the state vector "h",
// message block vector "m" (last block is padded with zeros to full
// block size, if required), 2w-bit offset counter "t", and final block
// indicator flag "f".  Local vector v[0..15] is used in processing.  F
// returns a new state vector.  The number of rounds, "r", is 12 for
// BLAKE2b and 10 for BLAKE2s.  Rounds are numbered from 0 to r - 1.
// FUNCTION F( h[0..7], m[0..15], t, f )
// |
// |      // Initialize local work vector v[0..15]
// |      v[0..7] := h[0..7]              // First half from state.
// |      v[8..15] := IV[0..7]            // Second half from IV.
// |
// |      v[12] := v[12] ^ (t mod 2**w)   // Low word of the offset.
// |      v[13] := v[13] ^ (t >> w)       // High word.
// |
// |      IF f = TRUE THEN                // last block flag?
// |      |   v[14] := v[14] ^ 0xFF..FF   // Invert all bits.
// |      END IF.
// |
// |      // Cryptographic mixing
// |      FOR i = 0 TO r - 1 DO           // Ten or twelve rounds.
// |      |
// |      |   // Message word selection permutation for this round.
// |      |   s[0..15] := SIGMA[i mod 10][0..15]
// |      |
// |      |   v := G( v, 0, 4,  8, 12, m[s[ 0]], m[s[ 1]] )
// |      |   v := G( v, 1, 5,  9, 13, m[s[ 2]], m[s[ 3]] )
// |      |   v := G( v, 2, 6, 10, 14, m[s[ 4]], m[s[ 5]] )
// |      |   v := G( v, 3, 7, 11, 15, m[s[ 6]], m[s[ 7]] )
// |      |
// |      |   v := G( v, 0, 5, 10, 15, m[s[ 8]], m[s[ 9]] )
// |      |   v := G( v, 1, 6, 11, 12, m[s[10]], m[s[11]] )
// |      |   v := G( v, 2, 7,  8, 13, m[s[12]], m[s[13]] )
// |      |   v := G( v, 3, 4,  9, 14, m[s[14]], m[s[15]] )
// |      |
// |      END FOR
// |
// |      FOR i = 0 TO 7 DO               // XOR the two halves.
// |      |   h[i] := h[i] ^ v[i] ^ v[i + 8]
// |      END FOR.
// |
// |      RETURN h[0..7]                  // New state.
// |
// END FUNCTION.
//

fn blake2s_compression<ConstraintF: PrimeField>(
    h: &mut [UInt32<ConstraintF>],
    m: &[UInt32<ConstraintF>],
    t: u64,
    f: bool,
) -> Result<(), SynthesisError> {
    assert_eq!(h.len(), 8);
    assert_eq!(m.len(), 16);

    // static const uint32_t blake2s_iv[8] =
    // {
    // 0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A,
    // 0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19
    // };
    //

    let mut v = Vec::with_capacity(16);
    v.extend_from_slice(h);
    v.push(UInt32::constant(0x6A09E667));
    v.push(UInt32::constant(0xBB67AE85));
    v.push(UInt32::constant(0x3C6EF372));
    v.push(UInt32::constant(0xA54FF53A));
    v.push(UInt32::constant(0x510E527F));
    v.push(UInt32::constant(0x9B05688C));
    v.push(UInt32::constant(0x1F83D9AB));
    v.push(UInt32::constant(0x5BE0CD19));

    assert_eq!(v.len(), 16);

    v[12] = v[12].xor(&UInt32::constant(t as u32))?;
    v[13] = v[13].xor(&UInt32::constant((t >> 32) as u32))?;

    if f {
        v[14] = v[14].xor(&UInt32::constant(u32::max_value()))?;
    }

    for i in 0..10 {
        let s = SIGMA[i % 10];

        mixing_g(&mut v, 0, 4, 8, 12, &m[s[0]], &m[s[1]])?;
        mixing_g(&mut v, 1, 5, 9, 13, &m[s[2]], &m[s[3]])?;
        mixing_g(&mut v, 2, 6, 10, 14, &m[s[4]], &m[s[5]])?;
        mixing_g(&mut v, 3, 7, 11, 15, &m[s[6]], &m[s[7]])?;
        mixing_g(&mut v, 0, 5, 10, 15, &m[s[8]], &m[s[9]])?;
        mixing_g(&mut v, 1, 6, 11, 12, &m[s[10]], &m[s[11]])?;
        mixing_g(&mut v, 2, 7, 8, 13, &m[s[12]], &m[s[13]])?;
        mixing_g(&mut v, 3, 4, 9, 14, &m[s[14]], &m[s[15]])?;
    }

    for i in 0..8 {
        h[i] = h[i].xor(&v[i])?;
        h[i] = h[i].xor(&v[i + 8])?;
    }

    Ok(())
}

// FUNCTION BLAKE2( d[0..dd-1], ll, kk, nn )
// |
// |     h[0..7] := IV[0..7]          // Initialization Vector.
// |
// |     // Parameter block p[0]
// |     h[0] := h[0] ^ 0x01010000 ^ (kk << 8) ^ nn
// |
// |     // Process padded key and data blocks
// |     IF dd > 1 THEN
// |     |       FOR i = 0 TO dd - 2 DO
// |     |       |       h := F( h, d[i], (i + 1) * bb, FALSE )
// |     |       END FOR.
// |     END IF.
// |
// |     // Final block.
// |     IF kk = 0 THEN
// |     |       h := F( h, d[dd - 1], ll, TRUE )
// |     ELSE
// |     |       h := F( h, d[dd - 1], ll + bb, TRUE )
// |     END IF.
// |
// |     RETURN first "nn" bytes from little-endian word array h[].
// |
// END FUNCTION.
//

pub fn evaluate_blake2s<ConstraintF: PrimeField>(
    input: &[Boolean<ConstraintF>],
) -> Result<Vec<UInt32<ConstraintF>>, SynthesisError> {
    assert!(input.len() % 8 == 0);
    let mut parameters = [0; 8];
    parameters[0] = 0x01010000 ^ 32;
    evaluate_blake2s_with_parameters(input, &parameters)
}

pub fn evaluate_blake2s_with_parameters<F: PrimeField>(
    input: &[Boolean<F>],
    parameters: &[u32; 8],
) -> Result<Vec<UInt32<F>>, SynthesisError> {
    assert!(input.len() % 8 == 0);

    let mut h = Vec::with_capacity(8);
    h.push(UInt32::constant(0x6A09E667).xor(&UInt32::constant(parameters[0]))?);
    h.push(UInt32::constant(0xBB67AE85).xor(&UInt32::constant(parameters[1]))?);
    h.push(UInt32::constant(0x3C6EF372).xor(&UInt32::constant(parameters[2]))?);
    h.push(UInt32::constant(0xA54FF53A).xor(&UInt32::constant(parameters[3]))?);
    h.push(UInt32::constant(0x510E527F).xor(&UInt32::constant(parameters[4]))?);
    h.push(UInt32::constant(0x9B05688C).xor(&UInt32::constant(parameters[5]))?);
    h.push(UInt32::constant(0x1F83D9AB).xor(&UInt32::constant(parameters[6]))?);
    h.push(UInt32::constant(0x5BE0CD19).xor(&UInt32::constant(parameters[7]))?);

    let mut blocks: Vec<Vec<UInt32<F>>> = vec![];

    for block in input.chunks(512) {
        let mut this_block = Vec::with_capacity(16);
        for word in block.chunks(32) {
            let mut tmp = word.to_vec();
            while tmp.len() < 32 {
                tmp.push(Boolean::constant(false));
            }
            this_block.push(UInt32::from_bits_le(&tmp));
        }
        while this_block.len() < 16 {
            this_block.push(UInt32::constant(0));
        }
        blocks.push(this_block);
    }

    if blocks.is_empty() {
        blocks.push((0..16).map(|_| UInt32::constant(0)).collect());
    }

    for (i, block) in blocks[0..blocks.len() - 1].iter().enumerate() {
        blake2s_compression(&mut h, block, ((i as u64) + 1) * 64, false)?;
    }

    blake2s_compression(
        &mut h,
        &blocks[blocks.len() - 1],
        (input.len() / 8) as u64,
        true,
    )?;

    Ok(h)
}

use crate::prf::Blake2s;

pub struct Blake2sGadget;
#[derive(Clone, Debug)]
pub struct OutputVar<ConstraintF: PrimeField>(pub Vec<UInt8<ConstraintF>>);

impl<ConstraintF: PrimeField> EqGadget<ConstraintF> for OutputVar<ConstraintF> {
    #[tracing::instrument(target = "r1cs")]
    fn is_eq(&self, other: &Self) -> Result<Boolean<ConstraintF>, SynthesisError> {
        self.0.is_eq(&other.0)
    }

    /// If `should_enforce == true`, enforce that `self` and `other` are equal; else,
    /// enforce a vacuously true statement.
    #[tracing::instrument(target = "r1cs")]
    fn conditional_enforce_equal(
        &self,
        other: &Self,
        should_enforce: &Boolean<ConstraintF>,
    ) -> Result<(), SynthesisError> {
        self.0.conditional_enforce_equal(&other.0, should_enforce)
    }

    /// If `should_enforce == true`, enforce that `self` and `other` are not equal; else,
    /// enforce a vacuously true statement.
    #[tracing::instrument(target = "r1cs")]
    fn conditional_enforce_not_equal(
        &self,
        other: &Self,
        should_enforce: &Boolean<ConstraintF>,
    ) -> Result<(), SynthesisError> {
        self.0
            .as_slice()
            .conditional_enforce_not_equal(other.0.as_slice(), should_enforce)
    }
}

impl<ConstraintF: PrimeField> ToBytesGadget<ConstraintF> for OutputVar<ConstraintF> {
    #[inline]
    fn to_bytes(&self) -> Result<Vec<UInt8<ConstraintF>>, SynthesisError> {
        Ok(self.0.clone())
    }
}

impl<ConstraintF: PrimeField> AllocVar<[u8; 32], ConstraintF> for OutputVar<ConstraintF> {
    #[tracing::instrument(target = "r1cs", skip(cs, f))]
    fn new_variable<T: Borrow<[u8; 32]>>(
        cs: impl Into<Namespace<ConstraintF>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let bytes = f().map(|b| *b.borrow()).unwrap_or([0u8; 32]);
        match mode {
            AllocationMode::Constant => Ok(Self(UInt8::constant_vec(&bytes))),
            AllocationMode::Input => UInt8::new_input_vec(cs, &bytes).map(Self),
            AllocationMode::Witness => UInt8::new_witness_vec(cs, &bytes).map(Self),
        }
    }
}

impl<F: PrimeField> R1CSVar<F> for OutputVar<F> {
    type Value = [u8; 32];

    fn cs(&self) -> ConstraintSystemRef<F> {
        self.0.cs()
    }

    fn value(&self) -> Result<Self::Value, SynthesisError> {
        let mut value = [0u8; 32];
        for (val_i, self_i) in value.iter_mut().zip(&self.0) {
            *val_i = self_i.value()?;
        }
        Ok(value)
    }
}

impl<F: PrimeField> PRFGadget<Blake2s, F> for Blake2sGadget {
    type OutputVar = OutputVar<F>;

    #[tracing::instrument(target = "r1cs", skip(cs))]
    fn new_seed(cs: impl Into<Namespace<F>>, seed: &[u8; 32]) -> Vec<UInt8<F>> {
        let ns = cs.into();
        let cs = ns.cs();
        UInt8::new_witness_vec(r1cs_core::ns!(cs, "New Blake2s seed"), seed).unwrap()
    }

    #[tracing::instrument(target = "r1cs", skip(seed, input))]
    fn evaluate(seed: &[UInt8<F>], input: &[UInt8<F>]) -> Result<Self::OutputVar, SynthesisError> {
        assert_eq!(seed.len(), 32);
        let input: Vec<_> = seed
            .iter()
            .chain(input)
            .flat_map(|b| b.to_bits_le().unwrap())
            .collect();
        let result: Vec<_> = evaluate_blake2s(&input)?
            .into_iter()
            .flat_map(|int| int.to_bytes().unwrap())
            .collect();
        Ok(OutputVar(result))
    }
}

#[cfg(test)]
mod test {
    use algebra::ed_on_bls12_381::Fq as Fr;
    use rand::{Rng, SeedableRng};
    use rand_xorshift::XorShiftRng;

    use crate::prf::blake2s::{constraints::evaluate_blake2s, Blake2s as B2SPRF};
    use blake2::VarBlake2s;
    use r1cs_core::ConstraintSystem;

    use super::Blake2sGadget;
    use r1cs_std::prelude::*;

    #[test]
    fn test_blake2s_constraints() {
        let cs = ConstraintSystem::<Fr>::new_ref();
        let input_bits: Vec<_> = (0..512)
            .map(|_| Boolean::new_witness(r1cs_core::ns!(cs, "input bit"), || Ok(true)).unwrap())
            .collect();
        evaluate_blake2s(&input_bits).unwrap();
        assert!(cs.is_satisfied().unwrap());
        assert_eq!(cs.num_constraints(), 21792);
    }

    #[test]
    fn test_blake2s_prf() {
        use crate::prf::{PRFGadget, PRF};
        use rand::Rng;

        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
        let cs = ConstraintSystem::<Fr>::new_ref();

        let mut seed = [0u8; 32];
        rng.fill(&mut seed);

        let mut input = [0u8; 32];
        rng.fill(&mut input);

        let seed_var = Blake2sGadget::new_seed(cs.clone(), &seed);
        let input_var =
            UInt8::new_witness_vec(r1cs_core::ns!(cs, "declare_input"), &input).unwrap();
        let out = B2SPRF::evaluate(&seed, &input).unwrap();
        let actual_out_var = <Blake2sGadget as PRFGadget<_, Fr>>::OutputVar::new_witness(
            r1cs_core::ns!(cs, "declare_output"),
            || Ok(out),
        )
        .unwrap();

        let output_var = Blake2sGadget::evaluate(&seed_var, &input_var).unwrap();
        output_var.enforce_equal(&actual_out_var).unwrap();

        if !cs.is_satisfied().unwrap() {
            println!(
                "which is unsatisfied: {:?}",
                cs.which_is_unsatisfied().unwrap()
            );
        }
        assert!(cs.is_satisfied().unwrap());
    }

    #[test]
    fn test_blake2s_precomp_constraints() {
        // Test that 512 fixed leading bits (constants)
        // doesn't result in more constraints.

        let cs = ConstraintSystem::<Fr>::new_ref();
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
        let input_bits: Vec<_> = (0..512)
            .map(|_| Boolean::constant(rng.gen()))
            .chain((0..512).map(|_| {
                Boolean::new_witness(r1cs_core::ns!(cs, "input bit"), || Ok(true)).unwrap()
            }))
            .collect();
        evaluate_blake2s(&input_bits).unwrap();
        assert!(cs.is_satisfied().unwrap());
        assert_eq!(cs.num_constraints(), 21792);
    }

    #[test]
    fn test_blake2s_constant_constraints() {
        let cs = ConstraintSystem::<Fr>::new_ref();
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
        let input_bits: Vec<_> = (0..512)
            .map(|_| Boolean::<Fr>::constant(rng.gen()))
            .collect();
        evaluate_blake2s(&input_bits).unwrap();
        assert_eq!(cs.num_constraints(), 0);
    }

    #[test]
    fn test_blake2s() {
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        for input_len in (0..32).chain((32..256).filter(|a| a % 8 == 0)) {
            use digest::*;
            let mut h = VarBlake2s::new_keyed(&[], 32);

            let data: Vec<u8> = (0..input_len).map(|_| rng.gen()).collect();

            h.input(&data);

            let mut hash_result = Vec::with_capacity(h.output_size());
            h.variable_result(|res| hash_result.extend_from_slice(res));

            let cs = ConstraintSystem::<Fr>::new_ref();

            let mut input_bits = vec![];

            for input_byte in data.into_iter() {
                for bit_i in 0..8 {
                    let cs = r1cs_core::ns!(cs, "input bit");

                    input_bits.push(
                        Boolean::new_witness(cs, || Ok((input_byte >> bit_i) & 1u8 == 1u8))
                            .unwrap(),
                    );
                }
            }

            let r = evaluate_blake2s(&input_bits).unwrap();

            assert!(cs.is_satisfied().unwrap());

            let mut s = hash_result
                .iter()
                .flat_map(|&byte| (0..8).map(move |i| (byte >> i) & 1u8 == 1u8));

            for chunk in r {
                for b in chunk.to_bits_le() {
                    match b {
                        Boolean::Is(b) => {
                            assert!(s.next().unwrap() == b.value().unwrap());
                        }
                        Boolean::Not(b) => {
                            assert!(s.next().unwrap() != b.value().unwrap());
                        }
                        Boolean::Constant(b) => {
                            assert!(input_len == 0);
                            assert!(s.next().unwrap() == b);
                        }
                    }
                }
            }
        }
    }
}
