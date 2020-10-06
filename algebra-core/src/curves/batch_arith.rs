use crate::{biginteger::BigInteger, AffineCurve, Field, Vec};
use core::ops::Neg;
use either::Either;
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

        let mut two_a = bases.to_vec();
        let instr = (0..batch_size).map(|x| x as u32).collect::<Vec<_>>();
        Self::batch_double_in_place(&mut two_a, &instr[..], None);

        let mut tables = Vec::<Self>::with_capacity(half_size * batch_size);
        tables.extend_from_slice(bases);
        let mut scratch_space = Vec::<Option<Self>>::with_capacity((batch_size - 1) / 2 + 1);

        for i in 1..half_size {
            let instr = (0..batch_size)
                .map(|x| (((i - 1) * batch_size + x) as u32, x as u32))
                .collect::<Vec<_>>();
            Self::batch_add_write_read_self(
                &two_a[..],
                &instr[..],
                &mut tables,
                &mut scratch_space,
            );
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

        if negate.is_some() {
            assert_eq!(scalars.len(), negate.unwrap().len()); // precompute bounds check
        }

        let f = false;
        while !all_none {
            let iter = match negate {
                None => Either::Left(core::iter::repeat(&f).take(batch_size)),
                Some(bools) => Either::Right(bools.iter()),
            };
            let mut opcode_row = Vec::with_capacity(batch_size);
            for (s, &neg) in scalars.iter_mut().zip(iter) {
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
                        if neg {
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
        bases: &mut [Self],
        other: &[Self],
        index: &[(u32, u32)],
        scratch_space: &mut Vec<Self>,
    );

    /// Lookups up group elements according to index, and either adds and writes or simply
    /// writes them to new_elems, using scratch space to store intermediate values. Scratch
    /// space is always cleared after use.
    fn batch_add_write(
        lookup: &[Self],
        index: &[(u32, u32)],
        new_elems: &mut Vec<Self>,
        scratch_space: &mut Vec<Option<Self>>,
    );

    /// Similar to batch_add_write, only that the lookup for the first operand is performed
    /// in new_elems rather than lookup
    fn batch_add_write_read_self(
        lookup: &[Self],
        index: &[(u32, u32)],
        new_elems: &mut Vec<Self>,
        scratch_space: &mut Vec<Option<Self>>,
    );

    /// Performs a batch scalar multiplication using the w-NAF encoding
    /// utilising the primitive batched ops
    fn batch_scalar_mul_in_place<BigInt: BigInteger>(
        mut bases: &mut [Self],
        scalars: &mut [BigInt],
        w: usize,
    ) {
        let batch_size = bases.len();
        let opcode_vectorised = Self::batch_wnaf_opcode_recoding::<BigInt>(scalars, w, None);
        let tables = Self::batch_wnaf_tables(bases, w);

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
                        tables[(idx as usize) / 2 * batch_size + i].clone()
                    } else {
                        tables[(-idx as usize) / 2 * batch_size + i].clone().neg()
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
