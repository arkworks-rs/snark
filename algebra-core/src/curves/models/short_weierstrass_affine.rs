#[macro_export]
macro_rules! specialise_affine_to_proj {
    ($GroupProjective: ident) => {
        #[cfg(feature = "prefetch")]
        use crate::prefetch;
        use crate::{
            biginteger::BigInteger,
            curves::batch_arith::{decode_endo_from_usize, ENDO_CODING_BITS},
        };

        #[derive(Derivative)]
        #[derivative(
            Copy(bound = "P: Parameters"),
            Clone(bound = "P: Parameters"),
            PartialEq(bound = "P: Parameters"),
            Eq(bound = "P: Parameters"),
            Debug(bound = "P: Parameters"),
            Hash(bound = "P: Parameters")
        )]

        pub struct GroupAffine<P: Parameters> {
            pub x: P::BaseField,
            pub y: P::BaseField,
            pub infinity: bool,
            #[derivative(Debug = "ignore")]
            _params: PhantomData<P>,
        }

        impl<P: Parameters> AffineCurve for GroupAffine<P> {
            const COFACTOR: &'static [u64] = P::COFACTOR;
            type BaseField = P::BaseField;
            type ScalarField = P::ScalarField;
            type Projective = $GroupProjective<P>;

            fn prime_subgroup_generator() -> Self {
                Self::new(
                    P::AFFINE_GENERATOR_COEFFS.0,
                    P::AFFINE_GENERATOR_COEFFS.1,
                    false,
                )
            }

            fn from_random_bytes(bytes: &[u8]) -> Option<Self> {
                P::BaseField::from_random_bytes_with_flags(bytes).and_then(|(x, flags)| {
                    let infinity_flag_mask = SWFlags::Infinity.u8_bitmask();
                    let positive_flag_mask = SWFlags::PositiveY.u8_bitmask();
                    // if x is valid and is zero and only the infinity flag is set, then parse this
                    // point as infinity. For all other choices, get the original point.
                    if x.is_zero() && flags == infinity_flag_mask {
                        Some(Self::zero())
                    } else {
                        let is_positive = flags & positive_flag_mask != 0;
                        Self::get_point_from_x(x, is_positive)
                    }
                })
            }

            fn mul<S: Into<<Self::ScalarField as PrimeField>::BigInt>>(
                &self,
                by: S,
            ) -> Self::Projective {
                if P::GLV {
                    let w = 3;
                    let mut res = Self::Projective::zero();
                    let self_proj = self.into_projective();
                    impl_glv_mul!(Self::Projective, P, w, self_proj, res, by);
                    res
                } else {
                    let bits = BitIterator::new(by.into());
                    self.mul_bits(bits)
                }
            }

            #[inline]
            fn mul_by_cofactor_to_projective(&self) -> Self::Projective {
                self.scale_by_cofactor()
            }

            fn mul_by_cofactor_inv(&self) -> Self {
                self.mul(P::COFACTOR_INV).into()
            }
        }

        #[cfg(feature = "prefetch")]
        macro_rules! prefetch_slice {
            ($slice_1: ident, $slice_2: ident, $prefetch_iter: ident) => {
                if let Some((idp_1, idp_2)) = $prefetch_iter.next() {
                    prefetch::<Self>(&mut $slice_1[*idp_1]);
                    prefetch::<Self>(&mut $slice_2[*idp_2]);
                }
            };

            ($slice_1: ident, $prefetch_iter: ident) => {
                if let Some((idp_1, _)) = $prefetch_iter.next() {
                    prefetch::<Self>(&mut $slice_1[*idp_1]);
                }
            };
        }

        #[cfg(feature = "prefetch")]
        macro_rules! prefetch_slice_endo {
            ($slice_1: ident, $slice_2: ident, $prefetch_iter: ident) => {
                if let Some((idp_1, idp_2)) = $prefetch_iter.next() {
                    let (idp_2, _) = decode_endo_from_usize(*idp_2);
                    prefetch::<Self>(&mut $slice_1[*idp_1]);
                    prefetch::<Self>(&$slice_2[idp_2]);
                }
            };
        }

        impl<P: Parameters> BatchGroupArithmetic for GroupAffine<P> {
            type BBaseField = P::BaseField;
            // This implementation of batch group ops takes particular
            // care to make most use of points fetched from memory to prevent reallocations
            // It is adapted from Aztec's code.

            // https://github.com/AztecProtocol/barretenberg/blob/standardplonkjson/barretenberg/src/
            // aztec/ecc/curves/bn254/scalar_multiplication/scalar_multiplication.cpp

            // We require extra scratch space, and since we want to prevent allocation/deallocation overhead
            // we pass it externally for when this function is called many times
            #[inline]
            fn batch_double_in_place(
                bases: &mut [Self],
                index: &[usize],
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
                        prefetch::<Self>(&mut bases[*idp]);
                    }
                    let mut a = &mut bases[*idx];
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
                let mut scratch_space_counter = (0..scratch_space.len()).rev();
                #[cfg(feature = "prefetch")]
                {
                    prefetch_iter.next();
                    scratch_space_counter.next();
                }

                for idx in index.iter().rev() {
                    #[cfg(feature = "prefetch")]
                    if let Some(idp) = prefetch_iter.next() {
                        prefetch::<Self>(&mut bases[*idp]);
                    }
                    let mut a = &mut bases[*idx];
                    if !a.is_zero() {
                        #[cfg(feature = "prefetch")]
                        if let Some(idp) = scratch_space_counter.next() {
                            prefetch::<P::BaseField>(&mut scratch_space[idp]);
                        }
                        let z = scratch_space.pop().unwrap();
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
            fn batch_add_in_place(
                bases: &mut [Self],
                other: &mut [Self],
                index: &[(usize, usize)],
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
                    prefetch_slice!(bases, other, prefetch_iter);

                    let (mut a, mut b) = (&mut bases[*idx], &mut other[*idy]);
                    if a.is_zero() || b.is_zero() {
                        continue;
                    } else if a.x == b.x {
                        half = match half {
                            None => P::BaseField::one().double().inverse(),
                            _ => half,
                        };
                        let h = half.unwrap();

                        // Double
                        // In our model, we consider self additions rare.
                        // So we consider it inconsequential to make them more expensive
                        // This costs 1 modular mul more than a standard squaring,
                        // and one amortised inversion
                        if a.y == b.y {
                            let x_sq = b.x.square();
                            b.x -= &b.y; // x - y
                            a.x = b.y.double(); // denominator = 2y
                            a.y = x_sq.double() + &x_sq + &P::COEFF_A; // numerator = 3x^2 + a
                            b.y -= &(h * &a.y); // y - (3x^2 + a)/2
                            a.y *= &inversion_tmp; // (3x^2 + a) * tmp
                            inversion_tmp *= &a.x; // update tmp
                        } else {
                            // No inversions take place if either operand is zero
                            a.infinity = true;
                            b.infinity = true;
                        }
                    } else {
                        // We can recover x1 + x2 from this. Note this is never 0.
                        a.x -= &b.x; // denominator = x1 - x2
                        a.y -= &b.y; // numerator = y1 - y2
                        a.y *= &inversion_tmp; // (y1 - y2)*tmp
                        inversion_tmp *= &a.x // update tmp
                    }
                }

                inversion_tmp = inversion_tmp.inverse().unwrap(); // this is always in Fp*

                #[cfg(feature = "prefetch")]
                let mut prefetch_iter = index.iter().rev();
                #[cfg(feature = "prefetch")]
                prefetch_iter.next();

                for (idx, idy) in index.iter().rev() {
                    #[cfg(feature = "prefetch")]
                    prefetch_slice!(bases, other, prefetch_iter);
                    let (mut a, b) = (&mut bases[*idx], other[*idy]);

                    if a.is_zero() {
                        *a = b;
                    } else if !b.is_zero() {
                        let lambda = a.y * &inversion_tmp;
                        inversion_tmp *= &a.x; // Remove the top layer of the denominator

                        // x3 = l^2 - x1 - x2 or for squaring: 2y + l^2 + 2x - 2y = l^2 - 2x
                        a.x += &b.x.double();
                        a.x = lambda.square() - &a.x;
                        // y3 = l*(x2 - x3) - y2 or
                        // for squaring: (3x^2 + a)/2y(x - y - x3) - (y - (3x^2 + a)/2) = l*(x - x3) - y
                        a.y = lambda * &(b.x - &a.x) - &b.y;
                    }
                }
            }

            #[inline]
            fn batch_add_in_place_same_slice(bases: &mut [Self], index: &[(usize, usize)]) {
                let mut inversion_tmp = P::BaseField::one();
                let mut half = None;

                #[cfg(feature = "prefetch")]
                let mut prefetch_iter = index.iter();
                #[cfg(feature = "prefetch")]
                {
                    prefetch_iter.next();
                    prefetch_iter.next();
                }

                // We run two loops over the data separated by an inversion
                for (idx, idy) in index.iter() {
                    #[cfg(feature = "prefetch")]
                    prefetch_slice!(bases, bases, prefetch_iter);
                    let (mut a, mut b) = if idx < idy {
                        let (x, y) = bases.split_at_mut(*idy);
                        (&mut x[*idx], &mut y[0])
                    } else {
                        let (x, y) = bases.split_at_mut(*idx);
                        (&mut y[0], &mut x[*idy])
                    };
                    if a.is_zero() || b.is_zero() {
                        continue;
                    } else if a.x == b.x {
                        half = match half {
                            None => P::BaseField::one().double().inverse(),
                            _ => half,
                        };
                        let h = half.unwrap();

                        // Double
                        // In our model, we consider self additions rare.
                        // So we consider it inconsequential to make them more expensive
                        // This costs 1 modular mul more than a standard squaring,
                        // and one amortised inversion
                        if a.y == b.y {
                            let x_sq = b.x.square();
                            b.x -= &b.y; // x - y
                            a.x = b.y.double(); // denominator = 2y
                            a.y = x_sq.double() + &x_sq + &P::COEFF_A; // numerator = 3x^2 + a
                            b.y -= &(h * &a.y); // y - (3x^2 + a)/2
                            a.y *= &inversion_tmp; // (3x^2 + a) * tmp
                            inversion_tmp *= &a.x; // update tmp
                        } else {
                            // No inversions take place if either operand is zero
                            a.infinity = true;
                            b.infinity = true;
                        }
                    } else {
                        // We can recover x1 + x2 from this. Note this is never 0.
                        a.x -= &b.x; // denominator = x1 - x2
                        a.y -= &b.y; // numerator = y1 - y2
                        a.y *= &inversion_tmp; // (y1 - y2)*tmp
                        inversion_tmp *= &a.x // update tmp
                    }
                }

                inversion_tmp = inversion_tmp.inverse().unwrap(); // this is always in Fp*

                #[cfg(feature = "prefetch")]
                let mut prefetch_iter = index.iter().rev();
                #[cfg(feature = "prefetch")]
                {
                    prefetch_iter.next();
                    prefetch_iter.next();
                }

                for (idx, idy) in index.iter().rev() {
                    #[cfg(feature = "prefetch")]
                    prefetch_slice!(bases, bases, prefetch_iter);
                    let (mut a, b) = if idx < idy {
                        let (x, y) = bases.split_at_mut(*idy);
                        (&mut x[*idx], y[0])
                    } else {
                        let (x, y) = bases.split_at_mut(*idx);
                        (&mut y[0], x[*idy])
                    };
                    if a.is_zero() {
                        *a = b;
                    } else if !b.is_zero() {
                        let lambda = a.y * &inversion_tmp;
                        inversion_tmp *= &a.x; // Remove the top layer of the denominator

                        // x3 = l^2 - x1 - x2 or for squaring: 2y + l^2 + 2x - 2y = l^2 - 2x
                        a.x += &b.x.double();
                        a.x = lambda.square() - &a.x;
                        // y3 = l*(x2 - x3) - y2 or
                        // for squaring: (3x^2 + a)/2y(x - y - x3) - (y - (3x^2 + a)/2) = l*(x - x3) - y
                        a.y = lambda * &(b.x - &a.x) - &b.y;
                    }
                }
            }

            #[inline]
            fn batch_add_in_place_read_only(
                bases: &mut [Self],
                other: &[Self],
                index: &[(usize, usize)],
                scratch_space: Option<&mut Vec<Self>>,
            ) {
                let mut inversion_tmp = P::BaseField::one();
                let mut half = None;

                let mut _scratch_space_inner = if scratch_space.is_none() {
                    Vec::<Self>::with_capacity(index.len())
                } else {
                    vec![]
                };
                let scratch_space = match scratch_space {
                    Some(vec) => vec,
                    None => &mut _scratch_space_inner,
                };

                #[cfg(feature = "prefetch")]
                let mut prefetch_iter = index.iter();
                #[cfg(feature = "prefetch")]
                prefetch_iter.next();

                // We run two loops over the data separated by an inversion
                for (idx, idy) in index.iter() {
                    #[cfg(feature = "prefetch")]
                    prefetch_slice_endo!(bases, other, prefetch_iter);
                    let (idy, endomorphism) = decode_endo_from_usize(*idy);

                    let mut a = &mut bases[*idx];

                    // Apply endomorphisms according to encoding
                    let mut b = if endomorphism % 2 == 1 {
                        other[idy].neg()
                    } else {
                        other[idy]
                    };

                    if P::GLV {
                        if endomorphism >> 1 == 1 {
                            P::glv_endomorphism_in_place(&mut b.x);
                        }
                    }

                    if a.is_zero() || b.is_zero() {
                        scratch_space.push(b);
                        continue;
                    } else if a.x == b.x {
                        half = match half {
                            None => P::BaseField::one().double().inverse(),
                            _ => half,
                        };
                        let h = half.unwrap();

                        // Double
                        // In our model, we consider self additions rare.
                        // So we consider it inconsequential to make them more expensive
                        // This costs 1 modular mul more than a standard squaring,
                        // and one amortised inversion
                        if a.y == b.y {
                            let x_sq = b.x.square();
                            b.x -= &b.y; // x - y
                            a.x = b.y.double(); // denominator = 2y
                            a.y = x_sq.double() + &x_sq + &P::COEFF_A; // numerator = 3x^2 + a
                            b.y -= &(h * &a.y); // y - (3x^2 + a)/2
                            a.y *= &inversion_tmp; // (3x^2 + a) * tmp
                            inversion_tmp *= &a.x; // update tmp
                        } else {
                            // No inversions take place if either operand is zero
                            a.infinity = true;
                            b.infinity = true;
                        }
                    } else {
                        // We can recover x1 + x2 from this. Note this is never 0.
                        a.x -= &b.x; // denominator = x1 - x2
                        a.y -= &b.y; // numerator = y1 - y2
                        a.y *= &inversion_tmp; // (y1 - y2)*tmp
                        inversion_tmp *= &a.x // update tmp
                    }
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
                    let (mut a, b) = (&mut bases[*idx], scratch_space.pop().unwrap());

                    if a.is_zero() {
                        *a = b;
                    } else if !b.is_zero() {
                        let lambda = a.y * &inversion_tmp;
                        inversion_tmp *= &a.x; // Remove the top layer of the denominator

                        // x3 = l^2 - x1 - x2 or for squaring: 2y + l^2 + 2x - 2y = l^2 - 2x
                        a.x += &b.x.double();
                        a.x = lambda.square() - &a.x;
                        // y3 = l*(x2 - x3) - y2 or
                        // for squaring: (3x^2 + a)/2y(x - y - x3) - (y - (3x^2 + a)/2) = l*(x - x3) - y
                        a.y = lambda * &(b.x - &a.x) - &b.y;
                    }
                }
            }

            fn batch_scalar_mul_in_place<BigInt: BigInteger>(
                mut bases: &mut [Self],
                scalars: &mut [BigInt],
                w: usize,
            ) {
                debug_assert!(bases.len() == scalars.len());
                if P::GLV {
                    let mut scratch_space = Vec::<Self::BBaseField>::with_capacity(bases.len());
                    let mut scratch_space_group = Vec::<Self>::with_capacity(bases.len() / w);
                    use itertools::{EitherOrBoth::*, Itertools};
                    let k_vec: Vec<_> = scalars
                        .iter()
                        .map(|k| {
                            P::glv_scalar_decomposition(
                                <P::ScalarField as PrimeField>::BigInt::from_slice(k.as_ref()),
                            )
                        })
                        .collect();

                    let mut k1_scalars: Vec<_> = k_vec.iter().map(|x| (x.0).1).collect();
                    let k1_negates: Vec<_> = k_vec.iter().map(|x| (x.0).0).collect();
                    let mut k2_scalars: Vec<_> = k_vec.iter().map(|x| (x.1).1).collect();
                    let k2_negates: Vec<_> = k_vec.iter().map(|x| (x.1).0).collect();

                    let opcode_vectorised_k1 = Self::batch_wnaf_opcode_recoding(
                        &mut k1_scalars[..],
                        w,
                        Some(k1_negates.as_slice()),
                    );
                    let opcode_vectorised_k2 = Self::batch_wnaf_opcode_recoding(
                        &mut k2_scalars[..],
                        w,
                        Some(k2_negates.as_slice()),
                    );

                    let tables = Self::batch_wnaf_tables(bases, w);
                    let half_size = 1 << w;
                    let batch_size = bases.len();

                    // Set all points to 0;
                    let zero = Self::zero();
                    for p in bases.iter_mut() {
                        *p = zero;
                    }
                    let noop_vec = vec![None; batch_size];

                    for (opcode_row_k1, opcode_row_k2) in opcode_vectorised_k1
                        .iter()
                        .zip_longest(opcode_vectorised_k2.iter())
                        .map(|x| match x {
                            Both(a, b) => (a, b),
                            Left(a) => (a, &noop_vec),
                            Right(b) => (&noop_vec, b),
                        })
                        .rev()
                    {
                        let index_double: Vec<usize> = opcode_row_k1
                            .iter()
                            .zip(opcode_row_k2.iter())
                            .enumerate()
                            .filter(|x| (x.1).0.is_some() || (x.1).1.is_some())
                            .map(|x| x.0)
                            .collect();

                        Self::batch_double_in_place(
                            &mut bases,
                            &index_double[..],
                            Some(&mut scratch_space),
                        );

                        let index_add_k1: Vec<(usize, usize)> = opcode_row_k1
                            .iter()
                            .enumerate()
                            .filter(|(_, op)| op.is_some() && op.unwrap() != 0)
                            .map(|(i, op)| {
                                let idx = op.unwrap();
                                if idx > 0 {
                                    (i, (i * half_size + (idx as usize) / 2) << ENDO_CODING_BITS)
                                } else {
                                    (
                                        i,
                                        ((i * half_size + (-idx as usize) / 2) << ENDO_CODING_BITS)
                                            + 1,
                                    )
                                }
                            })
                            .collect();

                        Self::batch_add_in_place_read_only(
                            &mut bases,
                            &tables[..],
                            &index_add_k1[..],
                            Some(&mut scratch_space_group),
                        );

                        let index_add_k2: Vec<(usize, usize)> = opcode_row_k2
                            .iter()
                            .enumerate()
                            .filter(|(_, op)| op.is_some() && op.unwrap() != 0)
                            .map(|(i, op)| {
                                let idx = op.unwrap();
                                if idx > 0 {
                                    (
                                        i,
                                        ((i * half_size + (idx as usize) / 2) << ENDO_CODING_BITS)
                                            + 2,
                                    )
                                } else {
                                    (
                                        i,
                                        ((i * half_size + (-idx as usize) / 2) << ENDO_CODING_BITS)
                                            + 3,
                                    )
                                }
                            })
                            .collect();

                        Self::batch_add_in_place_read_only(
                            &mut bases,
                            &tables[..],
                            &index_add_k2[..],
                            Some(&mut scratch_space_group),
                        );
                    }
                } else {
                    let mut scratch_space = Vec::<Self::BBaseField>::with_capacity(bases.len());
                    let opcode_vectorised =
                        Self::batch_wnaf_opcode_recoding::<BigInt>(scalars, w, None);
                    let tables = Self::batch_wnaf_tables(bases, w);
                    let half_size = 1 << w;

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

                        Self::batch_double_in_place(
                            &mut bases,
                            &index_double[..],
                            Some(&mut scratch_space),
                        );

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
            }
        }

        impl<P: Parameters> GroupAffine<P> {
            pub fn new(x: P::BaseField, y: P::BaseField, infinity: bool) -> Self {
                Self {
                    x,
                    y,
                    infinity,
                    _params: PhantomData,
                }
            }

            pub fn scale_by_cofactor(&self) -> <Self as AffineCurve>::Projective {
                self.mul_bits(BitIterator::new(P::COFACTOR))
            }

            pub(crate) fn mul_bits<S: AsRef<[u64]>>(
                &self,
                bits: BitIterator<S>,
            ) -> <Self as AffineCurve>::Projective {
                let mut res = <Self as AffineCurve>::Projective::zero();
                for i in bits {
                    res.double_in_place();
                    if i {
                        res.add_assign_mixed(&self)
                    }
                }
                res
            }

            /// Attempts to construct an affine point given an x-coordinate. The
            /// point is not guaranteed to be in the prime order subgroup.
            ///
            /// If and only if `greatest` is set will the lexicographically
            /// largest y-coordinate be selected.
            #[allow(dead_code)]
            pub fn get_point_from_x(x: P::BaseField, greatest: bool) -> Option<Self> {
                // Compute x^3 + ax + b
                let x3b = P::add_b(&((x.square() * &x) + &P::mul_by_a(&x)));

                x3b.sqrt().map(|y| {
                    let negy = -y;

                    let y = if (y < negy) ^ greatest { y } else { negy };
                    Self::new(x, y, false)
                })
            }

            /// Checks that the current point is on the elliptic curve.
            pub fn is_on_curve(&self) -> bool {
                if self.is_zero() {
                    true
                } else {
                    // Check that the point is on the curve
                    let y2 = self.y.square();
                    let x3b = P::add_b(&((self.x.square() * &self.x) + &P::mul_by_a(&self.x)));
                    y2 == x3b
                }
            }

            /// Checks that the current point is in the prime order subgroup given
            /// the point on the curve.
            pub fn is_in_correct_subgroup_assuming_on_curve(&self) -> bool {
                self.mul_bits(BitIterator::new(P::ScalarField::characteristic()))
                    .is_zero()
            }
        }

        impl<P: Parameters> Display for GroupAffine<P> {
            fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
                if self.infinity {
                    write!(f, "GroupAffine(Infinity)")
                } else {
                    write!(f, "GroupAffine(x={}, y={})", self.x, self.y)
                }
            }
        }

        impl<P: Parameters> Zero for GroupAffine<P> {
            fn zero() -> Self {
                Self::new(P::BaseField::zero(), P::BaseField::one(), true)
            }

            fn is_zero(&self) -> bool {
                self.infinity
            }
        }

        impl<P: Parameters> Add<Self> for GroupAffine<P> {
            type Output = Self;
            fn add(self, other: Self) -> Self {
                let mut copy = self;
                copy += &other;
                copy
            }
        }

        impl<'a, P: Parameters> AddAssign<&'a Self> for GroupAffine<P> {
            fn add_assign(&mut self, other: &'a Self) {
                let mut s_proj = <Self as AffineCurve>::Projective::from(*self);
                s_proj.add_assign_mixed(other);
                *self = s_proj.into();
            }
        }

        impl<P: Parameters> Neg for GroupAffine<P> {
            type Output = Self;

            #[inline]
            fn neg(self) -> Self {
                if !self.is_zero() {
                    Self::new(self.x, -self.y, false)
                } else {
                    self
                }
            }
        }

        impl<P: Parameters> ToBytes for GroupAffine<P> {
            #[inline]
            fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
                self.x.write(&mut writer)?;
                self.y.write(&mut writer)?;
                self.infinity.write(writer)
            }
        }

        impl<P: Parameters> FromBytes for GroupAffine<P> {
            #[inline]
            fn read<R: Read>(mut reader: R) -> IoResult<Self> {
                let x = P::BaseField::read(&mut reader)?;
                let y = P::BaseField::read(&mut reader)?;
                let infinity = bool::read(reader)?;
                Ok(Self::new(x, y, infinity))
            }
        }

        impl<P: Parameters> Default for GroupAffine<P> {
            #[inline]
            fn default() -> Self {
                Self::zero()
            }
        }

        impl_sw_curve_serializer!(Parameters);
    };
}

/// Implements GLV mul for a single element with a wNAF table
#[macro_export]
macro_rules! impl_glv_mul {
    ($Projective: ty, $P: ident, $w: ident, $self_proj: ident, $res: ident, $by: ident) => {
        // In the future, make this a GLV parameter entry
        let wnaf_recoding =
            |s: &mut <Self::ScalarField as PrimeField>::BigInt, is_neg: bool| -> Vec<i16> {
                let window_size: i16 = 1 << ($w + 1);
                let half_window_size: i16 = 1 << $w;

                let mut recoding = Vec::<i16>::with_capacity(s.num_bits() as usize / ($w + 1));

                while !s.is_zero() {
                    let op = if s.is_odd() {
                        let mut z: i16 = (s.as_ref()[0] % (1 << ($w + 1))) as i16;

                        if z < half_window_size {
                            s.sub_noborrow(&(z as u64).into());
                        } else {
                            z = z - window_size;
                            s.add_nocarry(&((-z) as u64).into());
                        }
                        if is_neg {
                            -z
                        } else {
                            z
                        }
                    } else {
                        0
                    };
                    recoding.push(op);
                    s.div2();
                }
                recoding
            };

        let ((k1_neg, mut k1), (k2_neg, mut k2)) = $P::glv_scalar_decomposition($by.into());
        let mut wnaf_table_k1 = Vec::<$Projective>::with_capacity(1 << $w);
        let double = $self_proj.double();
        wnaf_table_k1.push($self_proj);
        for _ in 1..(1 << $w) {
            wnaf_table_k1.push(*wnaf_table_k1.last().unwrap() + &double);
        }
        let mut wnaf_table_k2 = wnaf_table_k1.clone();
        wnaf_table_k2
            .iter_mut()
            .for_each(|p| $P::glv_endomorphism_in_place(&mut p.x));

        let k1_ops = wnaf_recoding(&mut k1, k1_neg);
        let k2_ops = wnaf_recoding(&mut k2, k2_neg);

        if k1_ops.len() > k2_ops.len() {
            for &op in k1_ops[k2_ops.len()..].iter().rev() {
                $res.double_in_place();
                if op > 0 {
                    $res += &wnaf_table_k1[(op as usize) / 2];
                } else if op < 0 {
                    $res += &wnaf_table_k1[(-op as usize) / 2].neg();
                }
            }
        } else {
            for &op in k2_ops[k1_ops.len()..].iter().rev() {
                $res.double_in_place();
                if op > 0 {
                    $res += &wnaf_table_k2[(op as usize) / 2];
                } else if op < 0 {
                    $res += &wnaf_table_k2[(-op as usize) / 2].neg();
                }
            }
        }
        for (&op1, &op2) in k1_ops.iter().zip(k2_ops.iter()).rev() {
            $res.double_in_place();
            if op1 > 0 {
                $res += &wnaf_table_k1[(op1 as usize) / 2];
            } else if op1 < 0 {
                $res += &wnaf_table_k1[(-op1 as usize) / 2].neg();
            }
            if op2 > 0 {
                $res += &wnaf_table_k2[(op2 as usize) / 2];
            } else if op2 < 0 {
                $res += &wnaf_table_k2[(-op2 as usize) / 2].neg();
            }
        }
    };
}
