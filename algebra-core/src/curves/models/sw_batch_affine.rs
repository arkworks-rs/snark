#[macro_export]
macro_rules! impl_sw_batch_affine {
    ($GroupAffine: ident) => {
        #[cfg(feature = "prefetch")]
        macro_rules! prefetch_slice {
            ($slice_1: ident, $slice_2: ident, $prefetch_iter: ident) => {
                if let Some((idp_1, idp_2)) = $prefetch_iter.next() {
                    prefetch::<Self>(&mut $slice_1[*idp_1 as usize]);
                    prefetch::<Self>(&mut $slice_2[*idp_2 as usize]);
                }
            };

            ($slice_1: ident, $prefetch_iter: ident) => {
                if let Some((idp_1, _)) = $prefetch_iter.next() {
                    prefetch::<Self>(&mut $slice_1[*idp_1 as usize]);
                }
            };
        }

        #[cfg(feature = "prefetch")]
        macro_rules! prefetch_slice_endo {
            ($slice_1: ident, $slice_2: ident, $prefetch_iter: ident) => {
                if let Some((idp_1, idp_2)) = $prefetch_iter.next() {
                    let (idp_2, _) = decode_endo_from_u32(*idp_2);
                    prefetch::<Self>(&mut $slice_1[*idp_1 as usize]);
                    prefetch::<Self>(&$slice_2[idp_2]);
                }
            };
        }

        #[cfg(feature = "prefetch")]
        macro_rules! prefetch_slice_write {
            ($slice_1: ident, $slice_2: ident, $prefetch_iter: ident) => {
                if let Some((idp_1, idp_2)) = $prefetch_iter.next() {
                    prefetch::<Self>(&$slice_1[*idp_1 as usize]);
                    if *idp_2 != !0u32 {
                        prefetch::<Self>(&$slice_2[*idp_2 as usize]);
                    }
                }
            };
        }

        macro_rules! batch_add_loop_1 {
            ($a: ident, $b: ident, $half: ident, $inversion_tmp: ident) => {
                if $a.is_zero() || $b.is_zero() {
                    ();
                } else if $a.x == $b.x {
                    $half = match $half {
                        None => P::BaseField::one().double().inverse(),
                        _ => $half,
                    };
                    let h = $half.unwrap();

                    // Double
                    // In our model, we consider self additions rare.
                    // So we consider it inconsequential to make them more expensive
                    // This costs 1 modular mul more than a standard squaring,
                    // and one amortised inversion
                    if $a.y == $b.y {
                        let x_sq = $b.x.square();
                        $b.x -= &$b.y; // x - y
                        $a.x = $b.y.double(); // denominator = 2y
                        $a.y = x_sq.double() + &x_sq + &P::COEFF_A; // numerator = 3x^2 + a
                        $b.y -= &(h * &$a.y); // y - (3x^2 + $a./2
                        $a.y *= &$inversion_tmp; // (3x^2 + a) * tmp
                        $inversion_tmp *= &$a.x; // update tmp
                    } else {
                        // No inversions take place if either operand is zero
                        $a.infinity = true;
                        $b.infinity = true;
                    }
                } else {
                    // We can recover x1 + x2 from this. Note this is never 0.
                    $a.x -= &$b.x; // denominator = x1 - x2
                    $a.y -= &$b.y; // numerator = y1 - y2
                    $a.y *= &$inversion_tmp; // (y1 - y2)*tmp
                    $inversion_tmp *= &$a.x // update tmp
                }
            };
        }

        macro_rules! batch_add_loop_2 {
            ($a: ident, $b: ident, $inversion_tmp: ident) => {
                if $a.is_zero() {
                    *$a = $b;
                } else if !$b.is_zero() {
                    let lambda = $a.y * &$inversion_tmp;
                    $inversion_tmp *= &$a.x; // Remove the top layer of the denominator

                    // x3 = l^2 - x1 - x2 or for squaring: 2y + l^2 + 2x - 2y = l^2 - 2x
                    $a.x += &$b.x.double();
                    $a.x = lambda.square() - &$a.x;
                    // y3 = l*(x2 - x3) - y2 or
                    // for squaring: (3x^2 + a)/2y(x - y - x3) - (y - (3x^2 + a)/2) = l*(x - x3) - y
                    $a.y = lambda * &($b.x - &$a.x) - &$b.y;
                }
            };
        }

        impl<P: Parameters> BatchGroupArithmetic for $GroupAffine<P> {
            type BBaseField = P::BaseField;
            /// This implementation of batch group ops takes particular
            /// care to make most use of points fetched from memory to prevent reallocations

            /// It is inspired by Aztec's approach:
            /// https://github.com/AztecProtocol/barretenberg/blob/
            /// c358fee3259a949da830f9867df49dc18768fa26/barretenberg/
            /// src/aztec/ecc/curves/bn254/scalar_multiplication/scalar_multiplication.cpp

            // We require extra scratch space, and since we want to prevent allocation/deallocation overhead
            // we pass it externally for when this function is called many times
            #[inline]
            fn batch_double_in_place(
                bases: &mut [Self],
                index: &[u32],
                scratch_space: Option<&mut Vec<Self::BBaseField>>,
            ) {
                let mut inversion_tmp = P::BaseField::one();

                let mut _scratch_space_inner = if scratch_space.is_none() {
                    Vec::with_capacity(index.len())
                } else {
                    vec![]
                };
                let scratch_space = match scratch_space {
                    Some(vec) => vec,
                    None => &mut _scratch_space_inner,
                };

                debug_assert!(scratch_space.len() == 0);

                #[cfg(feature = "prefetch")]
                let mut prefetch_iter = index.iter();
                #[cfg(feature = "prefetch")]
                prefetch_iter.next();

                for idx in index.iter() {
                    // Prefetch next group into cache
                    #[cfg(feature = "prefetch")]
                    if let Some(idp) = prefetch_iter.next() {
                        prefetch::<Self>(&mut bases[*idp as usize]);
                    }
                    let mut a = &mut bases[*idx as usize];
                    if !a.is_zero() {
                        if a.y.is_zero() {
                            a.infinity = true;
                        } else {
                            let x_sq = a.x.square();
                            let x_sq_3 = x_sq.double() + &x_sq + &P::COEFF_A; // numerator = 3x^2 + a
                            scratch_space.push(x_sq_3 * &inversion_tmp); // (3x^2 + a) * tmp
                            inversion_tmp *= &a.y.double(); // update tmp
                        }
                    }
                }

                inversion_tmp = inversion_tmp.inverse().unwrap(); // this is always in Fp*

                #[cfg(feature = "prefetch")]
                let mut prefetch_iter = index.iter().rev();
                #[cfg(feature = "prefetch")]
                prefetch_iter.next();

                for idx in index.iter().rev() {
                    #[cfg(feature = "prefetch")]
                    if let Some(idp) = prefetch_iter.next() {
                        prefetch::<Self>(&mut bases[*idp as usize]);
                    }
                    let mut a = &mut bases[*idx as usize];
                    if !a.is_zero() {
                        let z = scratch_space.pop().unwrap();
                        #[cfg(feature = "prefetch")]
                        if let Some(e) = scratch_space.last() {
                            prefetch::<P::BaseField>(e);
                        }
                        let lambda = z * &inversion_tmp;
                        inversion_tmp *= &a.y.double(); // Remove the top layer of the denominator

                        // x3 = l^2 + 2x
                        let x3 = &(lambda.square() - &a.x.double());
                        // y3 = l*(x - x3) - y
                        a.y = lambda * &(a.x - x3) - &a.y;
                        a.x = *x3;
                    }
                }

                debug_assert!(scratch_space.len() == 0);

                // We reset the vector
                // Clearing is really unnecessary, but we can do it anyway
                scratch_space.clear();
            }

            #[inline]
            fn batch_add_in_place(bases: &mut [Self], other: &mut [Self], index: &[(u32, u32)]) {
                let mut inversion_tmp = P::BaseField::one();
                let mut half = None;

                #[cfg(feature = "prefetch")]
                let mut prefetch_iter = index.iter();
                #[cfg(feature = "prefetch")]
                prefetch_iter.next();

                // We run two loops over the data separated by an inversion
                for (idx, idy) in index.iter() {
                    #[cfg(feature = "prefetch")]
                    prefetch_slice!(bases, other, prefetch_iter);

                    let (mut a, mut b) = (&mut bases[*idx as usize], &mut other[*idy as usize]);
                    batch_add_loop_1!(a, b, half, inversion_tmp);
                }

                inversion_tmp = inversion_tmp.inverse().unwrap(); // this is always in Fp*

                #[cfg(feature = "prefetch")]
                let mut prefetch_iter = index.iter().rev();
                #[cfg(feature = "prefetch")]
                prefetch_iter.next();

                for (idx, idy) in index.iter().rev() {
                    #[cfg(feature = "prefetch")]
                    prefetch_slice!(bases, other, prefetch_iter);
                    let (mut a, b) = (&mut bases[*idx as usize], other[*idy as usize]);
                    batch_add_loop_2!(a, b, inversion_tmp)
                }
            }

            #[inline]
            fn batch_add_in_place_same_slice(bases: &mut [Self], index: &[(u32, u32)]) {
                let mut inversion_tmp = P::BaseField::one();
                let mut half = None;

                #[cfg(feature = "prefetch")]
                let mut prefetch_iter = index.iter();
                #[cfg(feature = "prefetch")]
                {
                    prefetch_iter.next();
                }

                // We run two loops over the data separated by an inversion
                for (idx, idy) in index.iter() {
                    #[cfg(feature = "prefetch")]
                    prefetch_slice!(bases, bases, prefetch_iter);
                    let (mut a, mut b) = if idx < idy {
                        let (x, y) = bases.split_at_mut(*idy as usize);
                        (&mut x[*idx as usize], &mut y[0])
                    } else {
                        let (x, y) = bases.split_at_mut(*idx as usize);
                        (&mut y[0], &mut x[*idy as usize])
                    };
                    batch_add_loop_1!(a, b, half, inversion_tmp);
                }

                inversion_tmp = inversion_tmp.inverse().unwrap(); // this is always in Fp*

                #[cfg(feature = "prefetch")]
                let mut prefetch_iter = index.iter().rev();
                #[cfg(feature = "prefetch")]
                {
                    prefetch_iter.next();
                }

                for (idx, idy) in index.iter().rev() {
                    #[cfg(feature = "prefetch")]
                    prefetch_slice!(bases, bases, prefetch_iter);
                    let (mut a, b) = if idx < idy {
                        let (x, y) = bases.split_at_mut(*idy as usize);
                        (&mut x[*idx as usize], y[0])
                    } else {
                        let (x, y) = bases.split_at_mut(*idx as usize);
                        (&mut y[0], x[*idy as usize])
                    };
                    batch_add_loop_2!(a, b, inversion_tmp);
                }
            }

            #[inline]
            fn batch_add_in_place_read_only(
                bases: &mut [Self],
                other: &[Self],
                index: &[(u32, u32)],
                scratch_space: &mut Vec<Self>,
            ) {
                let mut inversion_tmp = P::BaseField::one();
                let mut half = None;

                #[cfg(feature = "prefetch")]
                let mut prefetch_iter = index.iter();
                #[cfg(feature = "prefetch")]
                prefetch_iter.next();

                // We run two loops over the data separated by an inversion
                for (idx, idy) in index.iter() {
                    let (idy, endomorphism) = decode_endo_from_u32(*idy);
                    #[cfg(feature = "prefetch")]
                    prefetch_slice_endo!(bases, other, prefetch_iter);

                    let mut a = &mut bases[*idx as usize];

                    // Apply endomorphisms according to encoding
                    let mut b = if endomorphism % 2 == 1 {
                        other[idy].neg()
                    } else {
                        other[idy]
                    };

                    batch_add_loop_1!(a, b, half, inversion_tmp);
                    scratch_space.push(b);
                }

                inversion_tmp = inversion_tmp.inverse().unwrap(); // this is always in Fp*

                #[cfg(feature = "prefetch")]
                let mut prefetch_iter = index.iter().rev();
                #[cfg(feature = "prefetch")]
                prefetch_iter.next();

                for (idx, _) in index.iter().rev() {
                    #[cfg(feature = "prefetch")]
                    {
                        prefetch_slice!(bases, prefetch_iter);
                        let len = scratch_space.len();
                        if len > 0 {
                            prefetch::<Self>(&mut scratch_space[len - 1]);
                        }
                    }
                    let (mut a, b) = (&mut bases[*idx as usize], scratch_space.pop().unwrap());
                    batch_add_loop_2!(a, b, inversion_tmp);
                }
            }

            fn batch_add_write(
                lookup: &[Self],
                index: &[(u32, u32)],
                new_elems: &mut Vec<Self>,
                scratch_space: &mut Vec<Option<Self>>,
            ) {
                let mut inversion_tmp = P::BaseField::one();
                let mut half = None;

                #[cfg(feature = "prefetch")]
                let mut prefetch_iter = index.iter();
                #[cfg(feature = "prefetch")]
                {
                    prefetch_iter.next();
                }

                // We run two loops over the data separated by an inversion
                for (idx, idy) in index.iter() {
                    #[cfg(feature = "prefetch")]
                    prefetch_slice_write!(lookup, lookup, prefetch_iter);

                    if *idy == !0u32 {
                        new_elems.push(lookup[*idx as usize]);
                        scratch_space.push(None);
                    } else {
                        let (mut a, mut b) = (lookup[*idx as usize], lookup[*idy as usize]);
                        batch_add_loop_1!(a, b, half, inversion_tmp);
                        new_elems.push(a);
                        scratch_space.push(Some(b));
                    }
                }

                inversion_tmp = inversion_tmp.inverse().unwrap(); // this is always in Fp*

                for (a, op_b) in new_elems.iter_mut().rev().zip(scratch_space.iter().rev()) {
                    match op_b {
                        Some(b) => {
                            let b_ = *b;
                            batch_add_loop_2!(a, b_, inversion_tmp);
                        }
                        None => (),
                    };
                }
                scratch_space.clear();
            }

            fn batch_add_write_read_self(
                lookup: &[Self],
                index: &[(u32, u32)],
                new_elems: &mut Vec<Self>,
                scratch_space: &mut Vec<Option<Self>>,
            ) {
                let mut inversion_tmp = P::BaseField::one();
                let mut half = None;

                #[cfg(feature = "prefetch")]
                let mut prefetch_iter = index.iter();
                #[cfg(feature = "prefetch")]
                prefetch_iter.next();

                // We run two loops over the data separated by an inversion
                for (idx, idy) in index.iter() {
                    #[cfg(feature = "prefetch")]
                    prefetch_slice_write!(new_elems, lookup, prefetch_iter);

                    if *idy == !0u32 {
                        new_elems.push(lookup[*idx as usize]);
                        scratch_space.push(None);
                    } else {
                        let (mut a, mut b) = (new_elems[*idx as usize], lookup[*idy as usize]);
                        batch_add_loop_1!(a, b, half, inversion_tmp);
                        new_elems.push(a);
                        scratch_space.push(Some(b));
                    }
                }

                inversion_tmp = inversion_tmp.inverse().unwrap(); // this is always in Fp*

                for (a, op_b) in new_elems.iter_mut().rev().zip(scratch_space.iter().rev()) {
                    match op_b {
                        Some(b) => {
                            let b_ = *b;
                            batch_add_loop_2!(a, b_, inversion_tmp);
                        }
                        None => (),
                    };
                }
                scratch_space.clear();
            }
        }
    };
}
