#[macro_export]
macro_rules! specialise_affine_to_proj {
    ($GroupProjective: ident) => {
        use crate::{biginteger::BigInteger, fields::FpParameters};

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

            fn mul<S: Into<<Self::ScalarField as PrimeField>::BigInt>>(&self, by: S) -> Self::Projective {
                let bits = BitIterator::new(by.into());
                self.mul_bits(bits)
            }

            #[inline]
            fn mul_by_cofactor_to_projective(&self) -> Self::Projective {
                self.scale_by_cofactor()
            }

            fn mul_by_cofactor_inv(&self) -> Self {
                self.mul(P::COFACTOR_INV).into()
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

        impl<'a, P: Parameters + FpParameters> BatchArithmetic<'a, GroupAffine<P>> for [GroupAffine<P>] {
            // Computes [-p, p, -3p, 3p, ..., -2^wp, 2^wp]
            fn batch_wnaf_tables(&self, w: usize) -> Vec<Vec<GroupAffine<P>>> {
                let half_size = 1 << w;
                let batch_size = self.len();

                let mut tables: vec![Vec::<GroupAffine<P>>::with_capacity(half_size << 1); batch_size];

                let mut a_2 = vec![];
                a_2.clone_from_slice(&self[..]);
                let mut tmp = vec![];
                tmp.clone_from_slice(&self[..]);

                a_2.batch_double_in_place_with_edge_cases(|x| x.iter_mut());

                for i in 0..half_size {
                    if i != 0 {
                        tmp.batch_add_in_place_with_edge_cases(
                            tmp.iter_mut().zip(a_2.iter_mut())
                        );
                    }

                    for (&mut table, p) in tables.iter_mut().zip(tmp) {
                        table.push(p.neg());
                        table.push(p);
                    }
                }
                // deref coercion
                // let res: &[Self] = &tables;
                // *res
                tables
            }

            // This function consumes the scalars
            // We can make this more generic in the future to use other than u16.
            fn batch_wnaf_opcode_recoding<BigInt: BigInteger + AsRef<[u64]>>(
                scalars: &mut [BigInt],
                w: usize
            ) -> Vec<Vec<Option<u16>>> {
                let batch_size = scalars.len();
                let window_size: u16 = 1 << (w + 1);
                let half_window_size: u16 = 1 << w;

                let op_code_vectorised = Vec::<Vec<Option<u16>>>::with_capacity(scalars[0].as_ref().len() * 64);

                let all_none = false;
                while !all_none {
                    let mut opcode_row = Vec::with_capacity(batch_size);

                    for s in scalars {
                        if s.is_zero() {
                            opcode_row.push(None);
                        } else {
                            let op = if s.is_odd() {
                                let mut z: u16 = (s.as_ref()[0] as u16) % window_size;

                                if z < half_window_size {
                                    s.sub_noborrow(&BigInt::from(z as u64));
                                } else {
                                    let tmp = window_size - z;
                                    s.add_nocarry(&BigInt::from(tmp as u64));
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

                    let all_none = opcode_row.iter().all(|x| x.is_none());
                    if !all_none {
                        op_code_vectorised.push(opcode_row);
                    } else {
                        break;
                    }
                }
                op_code_vectorised
            }

            // This implementation of batch group ops takes particular
            // care to make most use of points fetched from memory
            // And to reuse memory to prevent reallocations
            // It is directly adapted from Aztec's code.

            // https://github.com/AztecProtocol/barretenberg/blob/standardplonkjson/barretenberg/src/
            // aztec/ecc/curves/bn254/scalar_multiplication/scalar_multiplication.cpp

            fn batch_double_in_place_with_edge_cases<I>(&mut self, f: F) -> ()
            where
                F: FnMut(&mut Self) -> I,
                I: Iterator<Item = (&'a mut GroupAffine<P>, &'a mut GroupAffine<P>)> + DoubleEndedIterator
            {
                let mut inversion_tmp = P::BaseField::one();
                let mut scratch_space = Vec::new(); // with_capacity? How to get size?
                // We run two loops over the data separated by an inversion
                for a in f(self) {
                    if !a.is_zero() {
                        if a.y.is_zero() {
                            a.infinity = true;
                        } else {
                            let x_sq = a.x.square();
                            let x_sq_3 = *x_sq.double_in_place() + &x_sq; // numerator = 3x^2
                            scratch_space.push(x_sq_3 * &inversion_tmp); // 3x^2 * tmp
                            inversion_tmp *= &a.y.double(); // update tmp
                        }
                    }
                }

                inversion_tmp = inversion_tmp.inverse().unwrap(); // this is always in Fp*

                for a in f(self).rev() {
                    if !a.is_zero() {
                        let lambda = scratch_space.pop().unwrap() * &inversion_tmp;
                        inversion_tmp *= &a.x; // Remove the top layer of the denominator

                        // x3 = l^2 + 2x
                        let x3 = &(lambda.square() + &a.x.double());
                        // y3 = l*(x - x3) - y
                        a.y = lambda * &(a.x - x3) - &a.y;
                        a.x = *x3;
                    }
                }
            }

            // // May not be secure...
            // fn batch_double_in_place<'a, I>(op_iter: I) -> ()
            // where
            //     I: Iterator<Item = &'a Self>,
            // {
            //     let mut inversion_tmp = P::BaseField::one();
            //     let mut scratch_space = Vec::with_capacity(op_iter.size_hint().0);
            //     // We run two loops over the data separated by an inversion
            //     for &a in op_iter {
            //         let x_sq = a.x.square();
            //         let x_sq_3 = x_sq.double_in_place() + x_sq; // numerator = 3x^2
            //         scratch_space.push(x_sq_3 * inversion_tmp); // 3x^2 * tmp
            //         inversion_tmp *= a.x.double(); // update tmp
            //     }
            //
            //     inversion_tmp = inversion_tmp.inverse().unwrap(); // this is always in Fp*
            //
            //     for &a in op_iter.rev() {
            //         let lambda = scratch_space.pop() * inversion_tmp;
            //         inversion_tmp *= a.x; // Remove the top layer of the denominator
            //
            //         // x3 = l^2 + 2x
            //         let x3 = lambda.square_in_place() + a.x.double();
            //         // y3 = l*(x - x3) - y
            //         a.y = lambda * (a.x - x3) - a.y;
            //         a.x = x3;
            //     }
            // }

            // This implementation takes particular care to make most use of points fetched from memory
            // And to reuse memory to prevent reallocations
            // It is directly adapted from Aztec's code.

            // https://github.com/AztecProtocol/barretenberg/blob/standardplonkjson/barretenberg/src/
            // aztec/ecc/curves/bn254/scalar_multiplication/scalar_multiplication.cpp

            // This function consumes the second op as it mutates it in place
            // to prevent memory allocation

            fn batch_add_in_place_with_edge_cases<I>(&mut self, &mut other: Self, f: F) -> ()
            where
                F: FnMut(&mut Self, &mut Self) -> I,
                I: Iterator<Item = (&'a mut GroupAffine<P>, &'a GroupAffine<P>)> + DoubleEndedIterator
            {
                let mut inversion_tmp = P::BaseField::one();
                // let half = P::BaseField::from_repr(P::MODULUS_MINUS_ONE_DIV_TWO) + P::BaseField::one(); // (p + 1)/2 * 2 = 1
                // We run two loops over the data separated by an inversion
                for (a, b) in f(self) {
                    if a.is_zero() || b.is_zero() {
                        continue;
                    } else if a.x == b.x {
                        // double.
                        // In our model, we consider self additions rare.
                        // So we consider it inconsequential to make them more expensive
                        // This costs 1 modular mul more than a standard squaring
                        if a.y == b.y {
                            let x_sq = b.x.square();
                            b.x -= &b.y; // x - y
                            a.x = b.y.double(); // denominator = 2y
                            a.y = *x_sq.double_in_place() + &x_sq; // numerator = 3x^2
                            // b.y -= half * &a.y; // y - 3x^2/2
                            a.y *= &inversion_tmp; // 3x^2 * tmp
                            inversion_tmp *= &a.x; // update tmp
                        } else {
                            // No inversions take place if either operand is zero
                            a.infinity = true;
                            b.infinity = true;
                        }
                    } else {
                        a.x -= &b.x; // denominator = x1 - x2. We can recover x1 + x2 from this. Note this is never 0.
                        a.y -= &b.y; // numerator = y1 - y2
                        a.y *= &inversion_tmp; // (y1 - y2)*tmp
                        inversion_tmp *= &a.x // update tmp
                    }
                }

                inversion_tmp = inversion_tmp.inverse().unwrap(); // this is always in Fp*

                for (a, b) in f(self).rev() {
                    if a.is_zero() {
                        a = b;
                    } else if !b.is_zero() {
                        let lambda = a.y * &inversion_tmp;
                        inversion_tmp *= &a.x; // Remove the top layer of the denominator

                        // x3 = l^2 + x1 + x2 or for squaring: 2y + l^2 + 2x - 2y = l^2 + 2x
                        a.x += &(lambda.square() + &b.x.double());
                        // y3 = l*(x2 - x3) - y2 or for squaring: 3x^2/2y(x - y - x3) - (y - 3x^2/2) = l*(x - x3) - y
                        a.y = lambda * &(b.x - &a.x) - &b.y;
                    }
                }
            }

            // // This function consumes b_vec as it mutates it in place
            // // to prevent memory allocations
            // fn batch_add_in_place<'a, I>(op_iter: I)
            // where
            //     I: Iterator<Item = (&'a Self, Self)>,
            // {
            //     let mut inversion_tmp = P::BaseField::one();
            //     // We run two loops over the data separated by an inversion
            //     // let mut scratch_space = Vec::<AffineElement>::with_capacity(a_vec.len());
            //     for (&a, b) in op_iter {
            //         a.x -= b.x; // denominator = x1 - x2. We can recover x1 + x2 from this. Note this is never 0.
            //         a.y -= b.y; // numerator = y1 - y2
            //         a.y *= &inversion_tmp; // (y1 - y2)*tmp
            //         inversion_tmp *= a.x // update tmp
            //     }
            //
            //     inversion_tmp = &inversion_tmp.inverse().unwrap(); // this is always in Fp*
            //
            //     for (&a, b) in op_iter.rev() {
            //         let lambda = a.y * inversion_tmp;
            //         inversion_tmp *= &a.x; // Remove the top layer of the denominator
            //         a.x += lambda.square_in_place() + b.x.double();   // x3 = l^2 + x1 + x2
            //         a.y = lambda * (b.x - a.x) - b.y;     // y3 = l*(x2 - x3) - y2
            //     }
            // }

            fn batch_scalar_mul_in_place<BigInt: BigInteger>(
                &mut self,
                w: usize,
                scalars: &mut [BigInt],
            ) {
                let no_op: u16 = 1 << (w + 1); // noop is encoded as window_size
                let opcode_vectorised = Self::batch_wnaf_opcode_recoding::<BigInt>(scalars, w);
                let tables = self.batch_wnaf_tables(w);

                for opcode_row in opcode_vectorised.iter().rev() {

                    let double_iter = |points| {
                        points.iter_mut().zip(opcode_row)
                            .filter(|(p, op)| op.is_some())
                            .map(|x| x.0)
                    }

                    self.batch_double_in_place_with_edge_cases(double_iter);

                    // Copying to this vector might be really stupid...
                    let mut op2: Vec<GroupAffine<P>> = Vec::with_capacity(self.len() / w);

                    let add_iter = |points| {
                        points.iter_mut().zip(opcode_row).zip(tables.iter())
                            .filter(|((_, op), _)| op.is_some() && op.unwrap() != no_op)
                            .map(|((p, op), t)| {
                                op2.push(t[op.unwrap() as usize].clone());
                                (p, op2.last_mut().unwrap())
                            })
                    }

                    self.batch_add_in_place_with_edge_cases(add_iter);
                }
            }

            // fn batch_scalar_mul_in_place_glv<BigInt: BigInteger>(
            //     w: usize,
            //     points: &mut Vec<Self>,
            //     scalars: &mut Vec<BigInt>,
            // ) {
            //     assert_eq!(points.len(), scalars.len());
            //     let batch_size = points.len();
            //     let mut k1 = scalars;
            //     // let (mut k1, mut k2) = Self::batch_glv_decomposition(scalars);
            //
            //     // let p2 = points.map(|p| p.glv_endomorphism());
            //     Self::batch_scalar_mul_in_place::<BigInt>(w, points, &mut k1);
            //     // Self::batch_scalar_mul_in_place(w, p2, k2);
            //     // Self::batch_add_in_place_with_edge_cases(points, p2);
            // }
        }

        impl_sw_curve_serializer!(Parameters);
    }
}
