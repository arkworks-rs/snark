use crate::{AffineCurve, biginteger::BigInteger};
use num_traits::Zero;
use core::ops::Neg;

pub trait BatchGroupArithmetic
where
    Self: Sized + Clone + Copy + Zero + Neg<Output = Self>,
{
    // This function consumes the scalars
    // We can make this more generic in the future to use other than u16.

    // TODO: Generalise to A != 0
    // Computes [-p, p, -3p, 3p, ..., -2^wp, 2^wp]
    fn batch_wnaf_tables(bases: &[Self], w: usize) -> Vec<Vec<Self>> {
        let half_size = 1 << w;
        let batch_size = bases.len();

        let mut tables = vec![Vec::<Self>::with_capacity(half_size); batch_size];

        let mut a_2 = bases.to_vec();
        let mut tmp = bases.to_vec();

        let instr = (0..batch_size).collect::<Vec<usize>>();
        Self::batch_double_in_place(&mut a_2, &instr[..]);

        for i in 0..half_size {
            if i != 0 {
                let instr = (0..batch_size)
                    .map(|x| (x, x))
                    .collect::<Vec<(usize, usize)>>();
                Self::batch_add_in_place(&mut tmp, &mut a_2.to_vec()[..], &instr[..]);
            }

            for (table, p) in tables.iter_mut().zip(&tmp) {
                table.push(p.clone());
            }
        }
        tables
    }

    // This function mutates the scalars in place
    // We can make this more generic in the future to use other than i16.
    fn batch_wnaf_opcode_recoding<BigInt: BigInteger + AsRef<[u64]>>(
        scalars: &mut [BigInt],
        w: usize,
    ) -> Vec<Vec<Option<i16>>> {
        assert!(w > 0);
        let batch_size = scalars.len();
        let window_size: i16 = 1 << (w + 1);
        let half_window_size: i16 = 1 << w;

        let mut op_code_vectorised =
            Vec::<Vec<Option<i16>>>::with_capacity(scalars[0].as_ref().len() * 64);

        let mut all_none = false;
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
        op_code_vectorised
    }

    // This function consumes the second op as it mutates it in place
    // to prevent memory allocation
    fn batch_double_in_place(bases: &mut [Self], index: &[usize]);

    fn batch_add_in_place_same_slice(bases: &mut [Self], index: &[(usize, usize)]);

    fn batch_add_in_place(bases: &mut [Self], other: &mut [Self], index: &[(usize, usize)]);

    fn batch_scalar_mul_in_place<BigInt: BigInteger>(
        mut bases: &mut [Self],
        scalars: &mut [BigInt],
        w: usize,
    ) {
        let opcode_vectorised = Self::batch_wnaf_opcode_recoding::<BigInt>(scalars, w);
        let tables = Self::batch_wnaf_tables(bases, w);

        // Set all points to 0;
        let zero = Self::zero();
        for p in bases.iter_mut() {
            *p = zero;
        }

        for opcode_row in opcode_vectorised.iter().rev() {
            let index_double: Vec<usize> = opcode_row
                .iter()
                .enumerate()
                .filter(|x| x.1.is_some())
                .map(|x| x.0)
                .collect();

            Self::batch_double_in_place(&mut bases, &index_double[..]);

            let mut add_ops: Vec<Self> = tables
                .iter()
                .zip(opcode_row)
                .filter(|(_, op)| op.is_some() && op.unwrap() != 0)
                .map(|(t, op)| {
                    let idx = op.unwrap();
                    if idx > 0 {
                        t[(idx as usize) / 2].clone()
                    } else {
                        t[((-idx) as usize) / 2].clone().neg()
                    }
                })
                .collect();

            let index_add: Vec<(usize, usize)> = opcode_row
                .iter()
                .enumerate()
                .filter(|(_, op)| op.is_some() && op.unwrap() != 0)
                .map(|x| x.0)
                .enumerate()
                .map(|(x, y)| (y, x))
                .collect();

            Self::batch_add_in_place(&mut bases, &mut add_ops[..], &index_add[..]);
        }
    }

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

// We make the syntax cleaner by defining corresponding trait and impl for [G]
pub trait BatchGroupArithmeticSlice<G: AffineCurve> {
    fn batch_wnaf_tables(&self, w: usize) -> Vec<Vec<G>>;

    fn batch_wnaf_opcode_recoding<BigInt: BigInteger + AsRef<[u64]>>(
        scalars: &mut [BigInt],
        w: usize,
    ) -> Vec<Vec<Option<i16>>>;

    fn batch_double_in_place(&mut self, index: &[usize]);

    fn batch_add_in_place_same_slice(&mut self, index: &[(usize, usize)]);

    fn batch_add_in_place(&mut self, other: &mut Self, index: &[(usize, usize)]);

    fn batch_scalar_mul_in_place<BigInt: BigInteger>(&mut self, scalars: &mut [BigInt], w: usize);
}

impl<G: AffineCurve> BatchGroupArithmeticSlice<G> for [G] {
    fn batch_wnaf_tables(&self, w: usize) -> Vec<Vec<G>> {
        G::batch_wnaf_tables(self, w)
    }

    fn batch_wnaf_opcode_recoding<BigInt: BigInteger + AsRef<[u64]>>(
        scalars: &mut [BigInt],
        w: usize,
    ) -> Vec<Vec<Option<i16>>> {
        G::batch_wnaf_opcode_recoding::<BigInt>(scalars, w)
    }

    fn batch_double_in_place(&mut self, index: &[usize]) {
        G::batch_double_in_place(self, index);
    }

    fn batch_add_in_place_same_slice(&mut self, index: &[(usize, usize)]) {
        G::batch_add_in_place_same_slice(self, index);
    }

    fn batch_add_in_place(&mut self, other: &mut Self, index: &[(usize, usize)]) {
        G::batch_add_in_place(self, other, index);
    }

    fn batch_scalar_mul_in_place<BigInt: BigInteger>(&mut self, scalars: &mut [BigInt], w: usize) {
        G::batch_scalar_mul_in_place(self, scalars, w);
    }
}
