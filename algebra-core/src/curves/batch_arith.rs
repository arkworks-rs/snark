use crate::{biginteger::BigInteger, AffineCurve, Field, Vec};
use core::ops::Neg;
use num_traits::Zero;

/// We use a batch size that is big enough to amortise the cost of the actual inversion
/// close to zero while not straining the CPU cache by generating and fetching from
/// large w-NAF tables and slices [G]
pub const BATCH_SIZE: usize = 4096;
/// 0 == Identity; 1 == Neg; 2 == GLV; 3 == GLV + Neg
pub const ENDO_CODING_BITS: usize = 2;

#[inline(always)]
pub fn decode_endo_from_u32(index_code: u32) -> (usize, u8) {
    (
        index_code as usize >> ENDO_CODING_BITS,
        index_code as u8 % 4,
    )
}

pub trait BatchGroupArithmetic
where
    Self: Sized + Clone + Copy + Zero + Neg<Output = Self>,
{
    type BBaseField: Field;

    /*
    We use the w-NAF method, achieving point density of approximately 1/(w + 1)
    and requiring storage of only 2^(w - 1).
    Refer to e.g. Improved Techniques for Fast Exponentiation, Section 4
    Bodo M¨oller 2002. https://www.bmoeller.de/pdf/fastexp-icisc2002.pdf
    */

    /// Computes [[p, 3 * p, ..., (2^w - 1) * p], ..., [q, 3* q,  ..., ]]
    /// We need to manipulate the offsets when using the table
    fn batch_wnaf_tables(bases: &[Self], w: usize) -> Vec<Self> {
        let half_size = 1 << (w - 1);
        let batch_size = bases.len();

        let zero = Self::zero();
        let mut tables = vec![zero; half_size * batch_size];

        let mut a_2 = bases.to_vec();
        let mut tmp = bases.to_vec();

        let instr = (0..batch_size).map(|x| x as u32).collect::<Vec<_>>();
        Self::batch_double_in_place(&mut a_2, &instr[..], None);

        for i in 0..half_size {
            if i != 0 {
                let instr = (0..batch_size)
                    .map(|x| (x as u32, x as u32))
                    .collect::<Vec<_>>();
                Self::batch_add_in_place(&mut tmp, &mut a_2.to_vec()[..], &instr[..]);
            }

            for (elem_id, &p) in tmp.iter().enumerate() {
                tables[elem_id * half_size + i] = p.clone();
            }
        }
        tables
    }

    /// Computes the vectorised version of the wnaf integer recoding
    /// Optionally takes a slice of booleans which indicate whether that
    /// scalar is negative. If so, it negates the recoding.
    /// Mutates scalars in place
    fn batch_wnaf_opcode_recoding<BigInt: BigInteger>(
        scalars: &mut [BigInt],
        w: usize,
        negate: Option<&[bool]>,
    ) -> Vec<Vec<Option<i16>>> {
        assert!(w > 0);
        let batch_size = scalars.len();
        let window_size: i16 = 1 << (w + 1);
        let half_window_size: i16 = 1 << w;

        let mut op_code_vectorised =
            Vec::<Vec<Option<i16>>>::with_capacity(scalars[0].as_ref().len() * 64);

        let mut all_none = false;

        match negate {
            None => {
                while !all_none {
                    let mut opcode_row = Vec::with_capacity(batch_size);
                    for s in scalars.iter_mut() {
                        if s.is_zero() {
                            opcode_row.push(None);
                        } else {
                            let op = if s.is_odd() {
                                let mut z: i16 = (s.as_ref()[0] % (1 << (w + 1))) as i16;

                                if z < half_window_size {
                                    s.sub_noborrow(&BigInt::from(z as u64));
                                } else {
                                    z = z - window_size;
                                    s.add_nocarry(&BigInt::from((-z) as u64));
                                }
                                z
                            } else {
                                0
                            };
                            opcode_row.push(Some(op));
                            s.div2();
                        }
                    }
                    all_none = opcode_row.iter().all(|x| x.is_none());
                    if !all_none {
                        op_code_vectorised.push(opcode_row);
                    }
                }
            }
            Some(bools) => {
                while !all_none {
                    let mut opcode_row = Vec::with_capacity(batch_size);
                    for (s, neg) in scalars.iter_mut().zip(bools) {
                        if s.is_zero() {
                            opcode_row.push(None);
                        } else {
                            let op = if s.is_odd() {
                                let mut z: i16 = (s.as_ref()[0] % (1 << (w + 1))) as i16;
                                if z < half_window_size {
                                    s.sub_noborrow(&BigInt::from(z as u64));
                                } else {
                                    z = z - window_size;
                                    s.add_nocarry(&BigInt::from((-z) as u64));
                                }
                                if *neg {
                                    -z
                                } else {
                                    z
                                }
                            } else {
                                0
                            };
                            opcode_row.push(Some(op));
                            s.div2();
                        }
                    }
                    all_none = opcode_row.iter().all(|x| x.is_none());
                    if !all_none {
                        op_code_vectorised.push(opcode_row);
                    }
                }
            }
        }
        op_code_vectorised
    }

    /*
    We define a series of batched primitive EC ops, each of which is most suitable
    to a given scenario.

    We encode the indexes as u32s to save on fetch latency via better cacheing. The
    principle we are applying is that the len of the batch ops should never exceed
    about 2^20, and the table size would never exceed 2^10, so 32 bits will always
    be enough
    */

    /// Mutates bases to be doubled in place
    /// Accepts optional scratch space which might help by reducing the
    /// number of heap allocations for the Vector-based scratch_space
    fn batch_double_in_place(
        bases: &mut [Self],
        index: &[u32],
        scratch_space: Option<&mut Vec<Self::BBaseField>>,
    );

    /// Mutates bases in place and stores result in the first operand.
    /// The element corresponding to the second operand becomes junk data.
    fn batch_add_in_place_same_slice(bases: &mut [Self], index: &[(u32, u32)]);

    /// Mutates bases in place and stores result in bases.
    /// The elements in other become junk data.
    fn batch_add_in_place(bases: &mut [Self], other: &mut [Self], index: &[(u32, u32)]);

    /// Adds elements in bases with elements in other (for instance, a table), utilising
    /// a scratch space to store intermediate results.
    fn batch_add_in_place_read_only(
        _bases: &mut [Self],
        _other: &[Self],
        _index: &[(u32, u32)],
        _scratch_space: &mut Vec<Self>,
    ) {
        unimplemented!()
    }

    fn batch_add_write(
        _lookup: &[Self],
        _index: &[(u32, u32)],
        _new_elems: &mut Vec<Self>,
        _scratch_space: &mut Vec<Option<Self>>,
    ) {
        unimplemented!()
    }

    /// Performs a batch scalar multiplication using the w-NAF encoding
    /// utilising the primitive batched ops
    fn batch_scalar_mul_in_place<BigInt: BigInteger>(
        mut bases: &mut [Self],
        scalars: &mut [BigInt],
        w: usize,
    ) {
        let opcode_vectorised = Self::batch_wnaf_opcode_recoding::<BigInt>(scalars, w, None);
        let tables = Self::batch_wnaf_tables(bases, w);
        let half_size = 1 << (w - 1);

        // Set all points to 0;
        let zero = Self::zero();
        for p in bases.iter_mut() {
            *p = zero;
        }

        for opcode_row in opcode_vectorised.iter().rev() {
            let index_double: Vec<_> = opcode_row
                .iter()
                .enumerate()
                .filter(|x| x.1.is_some())
                .map(|x| x.0 as u32)
                .collect();

            Self::batch_double_in_place(&mut bases, &index_double[..], None);

            let mut add_ops: Vec<Self> = opcode_row
                .iter()
                .enumerate()
                .filter(|(_, op)| op.is_some() && op.unwrap() != 0)
                .map(|(i, op)| {
                    let idx = op.unwrap();
                    if idx > 0 {
                        tables[i * half_size + (idx as usize) / 2].clone()
                    } else {
                        tables[i * half_size + (-idx as usize) / 2].clone().neg()
                    }
                })
                .collect();

            let index_add: Vec<_> = opcode_row
                .iter()
                .enumerate()
                .filter(|(_, op)| op.is_some() && op.unwrap() != 0)
                .map(|x| x.0)
                .enumerate()
                .map(|(x, y)| (y as u32, x as u32))
                .collect();

            Self::batch_add_in_place(&mut bases, &mut add_ops[..], &index_add[..]);
        }
    }

    /// Chunks vectorised instructions into a size that does not require
    /// storing a lot of intermediate state
    fn get_chunked_instr<T: Clone>(instr: &[T], batch_size: usize) -> Vec<Vec<T>> {
        let mut res = Vec::new();

        let rem = instr.chunks_exact(batch_size).remainder();
        let mut chunks = instr.chunks_exact(batch_size).peekable();

        if chunks.len() == 0 {
            res.push(rem.to_vec());
        }

        while let Some(chunk) = chunks.next() {
            let chunk = if chunks.peek().is_none() {
                [chunk, rem].concat()
            } else {
                chunk.to_vec()
            };
            res.push(chunk);
        }
        res
    }
}

/// We make the syntax for performing batch ops on slices cleaner
/// by defining a corresponding trait and impl for [G] rather than on G
pub trait BatchGroupArithmeticSlice<G: AffineCurve> {
    fn batch_double_in_place(&mut self, index: &[u32]);

    fn batch_add_in_place_same_slice(&mut self, index: &[(u32, u32)]);

    fn batch_add_in_place(&mut self, other: &mut Self, index: &[(u32, u32)]);

    fn batch_add_write(
        &self,
        index: &[(u32, u32)],
        new_elems: &mut Vec<G>,
        scratch_space: &mut Vec<Option<G>>,
    );

    fn batch_scalar_mul_in_place<BigInt: BigInteger>(&mut self, scalars: &mut [BigInt], w: usize);
}

impl<G: AffineCurve> BatchGroupArithmeticSlice<G> for [G] {
    fn batch_double_in_place(&mut self, index: &[u32]) {
        G::batch_double_in_place(self, index, None);
    }

    fn batch_add_in_place_same_slice(&mut self, index: &[(u32, u32)]) {
        G::batch_add_in_place_same_slice(self, index);
    }

    fn batch_add_in_place(&mut self, other: &mut Self, index: &[(u32, u32)]) {
        G::batch_add_in_place(self, other, index);
    }

    fn batch_add_write(
        &self,
        index: &[(u32, u32)],
        new_elems: &mut Vec<G>,
        scratch_space: &mut Vec<Option<G>>,
    ) {
        G::batch_add_write(self, index, new_elems, scratch_space);
    }

    fn batch_scalar_mul_in_place<BigInt: BigInteger>(&mut self, scalars: &mut [BigInt], w: usize) {
        G::batch_scalar_mul_in_place(self, scalars, w);
    }
}
