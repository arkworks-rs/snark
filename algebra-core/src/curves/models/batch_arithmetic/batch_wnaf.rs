use crate::batch_group_ops::*;
use zexe_algebra_core::curves::short_weierstrass_jacobian::GroupAffine;
use zexe_algebra_core::biginteger::BigInteger;

// Since w-NAF is constant time, we can proceed in lockstep.
// sources?

// Should this be computed
// Will a table lookup be thwart large batch sizes?
// Need to find suitable batch size to amortise inversion

// We should bench the wnaf_tables using the generic add_or_double_in_place
// and the custom one that aims to reduce memory access.

// Computes [-p, p, -3p, 3p, ..., -2^wp, 2^wp]
pub fn batch_wnaf_tables<P>(w: usize, a: Vec<GroupAffine<P>>) -> Vec<Vec<G>>{
    let half_size = 1 << (w - 1);
    let batch_size = a.len();

    let mut tables = vec![Vec::with_capacity(half_size << 1); batch_size];

    let a_2 = batch_double_in_place_with_edge_cases(&mut a.copy());
    let tmp = a.copy();

    for (p, table) in tmp.zip(tables) { table.append(p); }
    for i in 1..half_size {
        batch_add_in_place_with_edge_cases(&mut tmp, a_2.copy());
        for (p, table) in tmp.zip(tables) {
            table.append(p.neg());
            table.append(p);
        }
    }
    tables
}

// This function consumes the scalars
// We can make this more generic in the future to use other than u16.
pub fn batch_wnaf_opcode_recoding<BigInt: BigInteger>(
    mut scalars: Vec<BigInt>,
    w: usize
) -> Vec<Vec<Option<u16>>> {
    let batch_size = scalars.len();
    let window_size: u16 = 1 << (w + 1);
    let half_window_size: u16 = 1 << w;

    let op_code_vectorised = Vec::<Vec<Option<u16>>>::with_capacity(scalars[0].len() * 64);

    let all_none = false;
    while !all_none {
        let mut opcode_row = Vec::with_capacity(batch_size);

        for s in scalars {
            if s.is_zero() {
                opcode_row.push(None);
            } else {
                let op = if s.is_odd() {
                    let mut z = (s.0[0] % window_size) as u16;

                    if z < half_window_size {
                        s.sub_noborrow(&BigInteger::from(z as u64));
                    } else {
                        let tmp = window_size - z as i16;
                        s.add_nocarry(&BigInteger::from(tmp as u64));
                        z = tmp - 1; // z = 0, 2, ..., 2^(w+1) - 2
                    }
                    z
                } else {
                    window_size // We encode 0s to be 2^(w+1)
                };
                opcode_row.push(Some(op));
                s.div2();
            }
        }

        all_none = opcode_row.all(|x| x.is_none());
        if !all_none {
            op_code_vectorised.push(opcode_row);
        } else {
            break;
        }
    }

    op_code_vectorised
}
