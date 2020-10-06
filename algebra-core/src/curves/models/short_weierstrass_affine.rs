#[macro_export]
macro_rules! specialise_affine_to_proj {
    ($GroupProjective: ident) => {
        use crate::batch_arith::decode_endo_from_u32;
        #[cfg(feature = "prefetch")]
        use crate::prefetch;

        #[derive(Derivative)]
        #[derivative(
            Copy(bound = "P: Parameters"),
            Clone(bound = "P: Parameters"),
            PartialEq(bound = "P: Parameters"),
            Eq(bound = "P: Parameters"),
            Debug(bound = "P: Parameters"),
            Hash(bound = "P: Parameters")
        )]
        #[repr(C)]
        pub struct GroupAffine<P: Parameters> {
            pub infinity: bool,
            pub x: P::BaseField,
            pub y: P::BaseField,
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
                let bits = BitIteratorBE::new(by.into());
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

        impl_sw_batch_affine!(GroupAffine);

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
