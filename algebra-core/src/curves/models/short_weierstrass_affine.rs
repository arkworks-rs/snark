#[macro_export]
macro_rules! specialise_affine_to_proj {
    ($GroupProjective: ident) => {
        #[cfg(feature = "prefetch")]
        use crate::prefetch;
        use crate::{
            biginteger::BigInteger,
            curves::batch_arith::{decode_endo_from_u32, ENDO_CODING_BITS},
        };

        #[derive(Derivative)]
        #[derivative(
            Copy(bound = "P: SWModelParameters"),
            Clone(bound = "P: SWModelParameters"),
            PartialEq(bound = "P: SWModelParameters"),
            Eq(bound = "P: SWModelParameters"),
            Debug(bound = "P: SWModelParameters"),
            Hash(bound = "P: SWModelParameters")
        )]
        #[repr(C)]
        pub struct GroupAffine<P: SWModelParameters> {
            pub infinity: bool,
            pub x: P::BaseField,
            pub y: P::BaseField,
            #[derivative(Debug = "ignore")]
            _params: PhantomData<P>,
        }

        impl<P: SWModelParameters> GroupAffine<P> {
            #[inline(always)]
            pub fn has_glv() -> bool {
                P::has_glv()
            }

            #[inline(always)]
            pub fn glv_endomorphism_in_place(elem: &mut <Self as AffineCurve>::BaseField) {
                P::glv_endomorphism_in_place(elem);
            }

            #[inline]
            pub fn glv_scalar_decomposition(
                k: <<Self as AffineCurve>::ScalarField as PrimeField>::BigInt,
            ) -> (
                (
                    bool,
                    <<Self as AffineCurve>::ScalarField as PrimeField>::BigInt,
                ),
                (
                    bool,
                    <<Self as AffineCurve>::ScalarField as PrimeField>::BigInt,
                ),
            ) {
                P::glv_scalar_decomposition(k)
            }
        }

        impl<P: SWModelParameters> AffineCurve for GroupAffine<P> {
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
                if P::has_glv() {
                    let w = 4;
                    let mut res = Self::Projective::zero();
                    let self_proj = self.into_projective();
                    impl_glv_mul!(Self::Projective, P, w, self_proj, res, by);
                    res
                } else {
                    let bits = BitIteratorBE::new(by.into());
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

        impl<P: SWModelParameters> GroupAffine<P> {
            pub fn new(x: P::BaseField, y: P::BaseField, infinity: bool) -> Self {
                Self {
                    x,
                    y,
                    infinity,
                    _params: PhantomData,
                }
            }

            pub fn scale_by_cofactor(&self) -> <Self as AffineCurve>::Projective {
                self.mul_bits(BitIteratorBE::new(P::COFACTOR))
            }

            pub(crate) fn mul_bits<S: AsRef<[u64]>>(
                &self,
                bits: BitIteratorBE<S>,
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
                self.mul_bits(BitIteratorBE::new(P::ScalarField::characteristic()))
                    .is_zero()
            }
        }

        impl<P: SWModelParameters> Display for GroupAffine<P> {
            fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
                if self.infinity {
                    write!(f, "GroupAffine(Infinity)")
                } else {
                    write!(f, "GroupAffine(x={}, y={})", self.x, self.y)
                }
            }
        }

        impl<P: SWModelParameters> Zero for GroupAffine<P> {
            fn zero() -> Self {
                Self::new(P::BaseField::zero(), P::BaseField::one(), true)
            }

            fn is_zero(&self) -> bool {
                self.infinity
            }
        }

        impl<P: SWModelParameters> Add<Self> for GroupAffine<P> {
            type Output = Self;
            fn add(self, other: Self) -> Self {
                let mut copy = self;
                copy += &other;
                copy
            }
        }

        impl<'a, P: SWModelParameters> AddAssign<&'a Self> for GroupAffine<P> {
            fn add_assign(&mut self, other: &'a Self) {
                let mut s_proj = <Self as AffineCurve>::Projective::from(*self);
                s_proj.add_assign_mixed(other);
                *self = s_proj.into();
            }
        }

        impl<P: SWModelParameters> Neg for GroupAffine<P> {
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

        impl_sw_batch_affine!(GroupAffine);

        impl<P: SWModelParameters> ToBytes for GroupAffine<P> {
            #[inline]
            fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
                self.x.write(&mut writer)?;
                self.y.write(&mut writer)?;
                self.infinity.write(writer)
            }
        }

        impl<P: SWModelParameters> FromBytes for GroupAffine<P> {
            #[inline]
            fn read<R: Read>(mut reader: R) -> IoResult<Self> {
                let x = P::BaseField::read(&mut reader)?;
                let y = P::BaseField::read(&mut reader)?;
                let infinity = bool::read(reader)?;
                Ok(Self::new(x, y, infinity))
            }
        }

        impl<P: SWModelParameters> Default for GroupAffine<P> {
            #[inline]
            fn default() -> Self {
                Self::zero()
            }
        }

        impl_sw_curve_serializer!(SWModelParameters);
    };
}

/// Implements GLV mul for a single element with a wNAF tables
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
        for _ in 1..(1 << ($w - 1)) {
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
