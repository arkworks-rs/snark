macro_rules! impl_Fp {
    ($Fp:ident, $FpParameters:ident, $BigInteger:ident, $BigIntegerType:ty, $limbs:expr) => {
        pub trait $FpParameters: FpParameters<BigInt = $BigIntegerType> {}

        #[derive(Derivative)]
        #[derivative(
            Default(bound = ""),
            Hash(bound = ""),
            Clone(bound = ""),
            Copy(bound = ""),
            Debug(bound = ""),
            PartialEq(bound = ""),
            Eq(bound = "")
        )]
        pub struct $Fp<P>(
            pub $BigIntegerType,
            #[derivative(Debug = "ignore")]
            #[doc(hidden)]
            pub PhantomData<P>,
        );

        impl<P> $Fp<P> {
            #[inline]
            pub const fn new(element: $BigIntegerType) -> Self {
                Self(element, PhantomData)
            }
        }

        impl<P: $FpParameters> $Fp<P> {
            #[inline]
            pub(crate) fn is_valid(&self) -> bool {
                self.0 < P::MODULUS
            }

            #[inline]
            fn reduce(&mut self) {
                if !self.is_valid() {
                    self.0.sub_noborrow(&P::MODULUS);
                }
            }
        }

        impl<P: $FpParameters> Zero for $Fp<P> {
            #[inline]
            fn zero() -> Self {
                $Fp::<P>($BigInteger::from(0), PhantomData)
            }

            #[inline]
            fn is_zero(&self) -> bool {
                self.0.is_zero()
            }
        }

        impl<P: $FpParameters> One for $Fp<P> {
            #[inline]
            fn one() -> Self {
                $Fp::<P>(P::R, PhantomData)
            }

            #[inline]
            fn is_one(&self) -> bool {
                self.0 == P::R
            }
        }

        impl<P: $FpParameters> Field for $Fp<P> {
            type BasePrimeField = Self;

            #[inline]
            fn double(&self) -> Self {
                let mut temp = *self;
                temp.double_in_place();
                temp
            }

            #[inline]
            fn double_in_place(&mut self) -> &mut Self {
                // This cannot exceed the backing capacity.
                self.0.mul2();
                // However, it may need to be reduced.
                self.reduce();
                self
            }

            #[inline]
            fn characteristic<'a>() -> &'a [u64] {
                P::MODULUS.as_ref()
            }

            #[inline]
            fn from_random_bytes_with_flags(bytes: &[u8]) -> Option<(Self, u8)> {
                let mut result_bytes = [0u8; $limbs * 8];
                for (result_byte, in_byte) in result_bytes.iter_mut().zip(bytes.iter()) {
                    *result_byte = *in_byte;
                }

                let mask: u64 = 0xffffffffffffffff >> P::REPR_SHAVE_BITS;
                // the flags will be at the same byte with the lowest shaven bits or the one after
                let flags_byte_position: usize = 7 - P::REPR_SHAVE_BITS as usize / 8;
                let flags_mask: u8 = ((1 << P::REPR_SHAVE_BITS % 8) - 1) << (8 - P::REPR_SHAVE_BITS % 8);
                // take the last 8 bytes and pass the mask
                let last_bytes = &mut result_bytes[($limbs - 1) * 8..];
                let mut flags: u8 = 0;
                for (i, (b, m)) in last_bytes.iter_mut().zip(&mask.to_le_bytes()).enumerate() {
                    if i == flags_byte_position {
                        flags = *b & flags_mask
                    }
                    *b &= m;
                }

                Self::deserialize(&mut &result_bytes[..]).ok().map(|f| (f, flags))
            }

            #[inline]
            fn square(&self) -> Self {
                let mut temp = self.clone();
                temp.square_in_place();
                temp
            }

            impl_field_square_in_place!($limbs);

            #[inline]
            fn inverse(&self) -> Option<Self> {
                if self.is_zero() {
                    None
                } else {
                    // Guajardo Kumar Paar Pelzl
                    // Efficient Software-Implementation of Finite Fields with Applications to
                    // Cryptography
                    // Algorithm 16 (BEA for Inversion in Fp)

                    let one = $BigInteger::from(1);

                    let mut u = self.0;
                    let mut v = P::MODULUS;
                    let mut b = $Fp::<P>(P::R2, PhantomData); // Avoids unnecessary reduction step.
                    let mut c = Self::zero();

                    while u != one && v != one {
                        while u.is_even() {
                            u.div2();

                            if b.0.is_even() {
                                b.0.div2();
                            } else {
                                b.0.add_nocarry(&P::MODULUS);
                                b.0.div2();
                            }
                        }

                        while v.is_even() {
                            v.div2();

                            if c.0.is_even() {
                                c.0.div2();
                            } else {
                                c.0.add_nocarry(&P::MODULUS);
                                c.0.div2();
                            }
                        }

                        if v < u {
                            u.sub_noborrow(&v);
                            b.sub_assign(&c);
                        } else {
                            v.sub_noborrow(&u);
                            c.sub_assign(&b);
                        }
                    }

                    if u == one {
                        Some(b)
                    } else {
                        Some(c)
                    }
                }
            }

            fn inverse_in_place(&mut self) -> Option<&mut Self> {
                if let Some(inverse) = self.inverse() {
                    *self = inverse;
                    Some(self)
                } else {
                    None
                }
            }

            #[inline]
            fn frobenius_map(&mut self, _: usize) {
                // No-op: No effect in a prime field.
            }
        }

        impl<P: $FpParameters> PrimeField for $Fp<P> {
            type Params = P;
            type BigInt = $BigIntegerType;

            #[inline]
            fn from_repr(r: $BigIntegerType) -> Option<Self> {
                let mut r = $Fp(r, PhantomData);
                if r.is_zero() {
                    Some(r)
                } else if r.is_valid() {
                    r *= &$Fp(P::R2, PhantomData);
                    Some(r)
                } else {
                    None
                }
            }

            impl_field_into_repr!($limbs, $BigIntegerType);
        }

        impl<P: $FpParameters> FftField for $Fp<P> {
            type FftParams = P;

            #[inline]
            fn two_adic_root_of_unity() -> Self {
                $Fp::<P>(P::TWO_ADIC_ROOT_OF_UNITY, PhantomData)
            }

            #[inline]
            fn large_subgroup_root_of_unity() -> Option<Self> {
                Some($Fp::<P>(P::LARGE_SUBGROUP_ROOT_OF_UNITY?, PhantomData))
            }

            #[inline]
            fn multiplicative_generator() -> Self {
                $Fp::<P>(P::GENERATOR, PhantomData)
            }
        }

        impl<P: $FpParameters> SquareRootField for $Fp<P> {
            #[inline]
            fn legendre(&self) -> LegendreSymbol {
                use crate::fields::LegendreSymbol::*;

                // s = self^((MODULUS - 1) // 2)
                let s = self.pow(P::MODULUS_MINUS_ONE_DIV_TWO);
                if s.is_zero() {
                    Zero
                } else if s.is_one() {
                    QuadraticResidue
                } else {
                    QuadraticNonResidue
                }
            }

            #[inline]
            fn sqrt(&self) -> Option<Self> {
                sqrt_impl!(Self, P, self)
            }

            fn sqrt_in_place(&mut self) -> Option<&mut Self> {
                (*self).sqrt().map(|sqrt| {
                    *self = sqrt;
                    self
                })
            }
        }

        impl<P: $FpParameters> Ord for $Fp<P> {
            #[inline(always)]
            fn cmp(&self, other: &Self) -> Ordering {
                self.into_repr().cmp(&other.into_repr())
            }
        }

        impl<P: $FpParameters> PartialOrd for $Fp<P> {
            #[inline(always)]
            fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
                Some(self.cmp(other))
            }
        }

        impl_prime_field_from_int!($Fp, u128, $FpParameters);
        impl_prime_field_from_int!($Fp, u64, $FpParameters);
        impl_prime_field_from_int!($Fp, u32, $FpParameters);
        impl_prime_field_from_int!($Fp, u16, $FpParameters);
        impl_prime_field_from_int!($Fp, u8, $FpParameters);
        impl_prime_field_from_int!($Fp, bool, $FpParameters);

        impl_prime_field_standard_sample!($Fp, $FpParameters);

        impl<P: $FpParameters> ToBytes for $Fp<P> {
            #[inline]
            fn write<W: Write>(&self, writer: W) -> IoResult<()> {
                self.into_repr().write(writer)
            }
        }

        impl<P: $FpParameters> FromBytes for $Fp<P> {
            #[inline]
            fn read<R: Read>(reader: R) -> IoResult<Self> {
                $BigInteger::read(reader).and_then(|b|
                    match $Fp::from_repr(b) {
                        Some(f) => Ok(f),
                        None => Err(crate::error("FromBytes::read failed")),
                    })
            }
        }

        impl<P: $FpParameters> FromStr for $Fp<P> {
            type Err = ();

            /// Interpret a string of numbers as a (congruent) prime field element.
            /// Does not accept unnecessary leading zeroes or a blank string.
            fn from_str(s: &str) -> Result<Self, Self::Err> {
                if s.is_empty() {
                    return Err(());
                }

                if s == "0" {
                    return Ok(Self::zero());
                }

                let mut res = Self::zero();

                let ten = Self::from(<Self as PrimeField>::BigInt::from(10));

                let mut first_digit = true;

                for c in s.chars() {
                    match c.to_digit(10) {
                        Some(c) => {
                            if first_digit {
                                if c == 0 {
                                    return Err(());
                                }

                                first_digit = false;
                            }

                            res.mul_assign(&ten);
                            let digit = Self::from(u64::from(c));
                            res.add_assign(&digit);
                        },
                        None => {
                            return Err(());
                        },
                    }
                }
                if !res.is_valid() {
                    Err(())
                } else {
                    Ok(res)
                }
            }
        }

        impl<P: $FpParameters> Display for $Fp<P> {
            #[inline]
            fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
                write!(f, stringify!($Fp"({})"), self.into_repr())
            }
        }

        impl<P: $FpParameters> Neg for $Fp<P> {
            type Output = Self;
            #[inline]
            #[must_use]
            fn neg(self) -> Self {
                if !self.is_zero() {
                    let mut tmp = P::MODULUS.clone();
                    tmp.sub_noborrow(&self.0);
                    $Fp::<P>(tmp, PhantomData)
                } else {
                    self
                }
            }
        }

        impl<'a, P: $FpParameters> Add<&'a $Fp<P>> for $Fp<P> {
            type Output = Self;

            #[inline]
            fn add(self, other: &Self) -> Self {
                let mut result = self.clone();
                result.add_assign(other);
                result
            }
        }

        impl<'a, P: $FpParameters> Sub<&'a $Fp<P>> for $Fp<P> {
            type Output = Self;

            #[inline]
            fn sub(self, other: &Self) -> Self {
                let mut result = self.clone();
                result.sub_assign(other);
                result
            }
        }

        impl<'a, P: $FpParameters> Mul<&'a $Fp<P>> for $Fp<P> {
            type Output = Self;

            #[inline]
            fn mul(self, other: &Self) -> Self {
                let mut result = self.clone();
                result.mul_assign(other);
                result
            }
        }

        impl<'a, P: $FpParameters> Div<&'a $Fp<P>> for $Fp<P> {
            type Output = Self;

            #[inline]
            fn div(self, other: &Self) -> Self {
                let mut result = self.clone();
                result.mul_assign(&other.inverse().unwrap());
                result
            }
        }

        impl_additive_ops_from_ref!($Fp, $FpParameters);
        impl_multiplicative_ops_from_ref!($Fp, $FpParameters);

        impl<'a, P: $FpParameters> AddAssign<&'a Self> for $Fp<P> {
            #[inline]
            fn add_assign(&mut self, other: &Self) {
                // This cannot exceed the backing capacity.
                self.0.add_nocarry(&other.0);
                // However, it may need to be reduced
                self.reduce();
            }
        }

        impl<'a, P: $FpParameters> SubAssign<&'a Self> for $Fp<P> {
            #[inline]
            fn sub_assign(&mut self, other: &Self) {
                // If `other` is larger than `self`, add the modulus to self first.
                if other.0 > self.0 {
                    self.0.add_nocarry(&P::MODULUS);
                }
                self.0.sub_noborrow(&other.0);
            }
        }

        impl<'a, P: $FpParameters> MulAssign<&'a Self> for $Fp<P> {
            impl_field_mul_assign!($limbs);
        }

        impl<'a, P: $FpParameters> DivAssign<&'a Self> for $Fp<P> {
            #[inline]
            fn div_assign(&mut self, other: &Self) {
                self.mul_assign(&other.inverse().unwrap());
            }
        }
    }
}
