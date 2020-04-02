macro_rules! impl_Fp {
    ($Fp:ident, $FpParameters:ident, $limbs:expr) => {
        pub trait $FpParameters: FpParameters<BigInt = BigInteger> {}

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
            pub BigInteger,
            #[derivative(Debug = "ignore")]
            #[doc(hidden)]
            pub PhantomData<P>,
        );

        impl<P> $Fp<P> {
            #[inline]
            pub const fn new(element: BigInteger) -> Self {
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
                $Fp::<P>(BigInteger::from(0), PhantomData)
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

                    let one = BigInteger::from(1);

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
            type BigInt = BigInteger;

            #[inline]
            fn from_repr(r: BigInteger) -> Self {
                let mut r = $Fp(r, PhantomData);
                if r.is_valid() {
                    r.mul_assign(&$Fp(P::R2, PhantomData));
                    r
                } else {
                    Self::zero()
                }
            }

            impl_field_into_repr!($limbs);

            #[inline]
            fn from_random_bytes(bytes: &[u8]) -> Option<Self> {
                let mut result_bytes = vec![0u8; $limbs * 8];
                for (result_byte, in_byte) in result_bytes.iter_mut().zip(bytes.iter()) {
                    *result_byte = *in_byte;
                }
                BigInteger::read(result_bytes.as_slice())
                    .ok()
                    .and_then(|mut res| {
                        res.as_mut()[$limbs-1] &= 0xffffffffffffffff >> P::REPR_SHAVE_BITS;
                        let result = Self::new(res);
                        if result.is_valid() {
                            Some(result)
                        } else {
                            None
                        }
                    })
            }

            #[inline]
            fn multiplicative_generator() -> Self {
                $Fp::<P>(P::GENERATOR, PhantomData)
            }

            #[inline]
            fn root_of_unity() -> Self {
                $Fp::<P>(P::ROOT_OF_UNITY, PhantomData)
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
                BigInteger::read(reader).map($Fp::from_repr)
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

                let ten = Self::from_repr(<Self as PrimeField>::BigInt::from(10));

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
                            res.add_assign(&Self::from_repr(<Self as PrimeField>::BigInt::from(
                                u64::from(c),
                            )));
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

macro_rules! impl_field_mul_assign {
    ($limbs:expr) => {
        #[inline]
        #[unroll_for_loops]
        fn mul_assign(&mut self, other: &Self) {
            let mut r = [0u64; $limbs];
            let mut carry1 = 0u64;
            let mut carry2 = 0u64;

            for i in 0..$limbs {
                r[0] = fa::mac(r[0], (self.0).0[0], (other.0).0[i], &mut carry1);
                let k = r[0].wrapping_mul(P::INV);
                fa::mac_discard(r[0], k, P::MODULUS.0[0], &mut carry2);
                for j in 1..$limbs {
                    r[j] = fa::mac_with_carry(r[j], (self.0).0[j], (other.0).0[i], &mut carry1);
                    r[j-1] = fa::mac_with_carry(r[j], k, P::MODULUS.0[j], &mut carry2);
                }
                r[$limbs-1] = carry1 + carry2;
            }
            for i in 0..$limbs {
                (self.0).0[i] = r[i];
            }
            self.reduce();
        }
    }
}

macro_rules! impl_field_into_repr {
    ($limbs:expr) => {
        #[inline]
        #[unroll_for_loops]
        fn into_repr(&self) -> BigInteger {
            let mut tmp = *self;
            let mut r = (tmp.0).0;
            // Montgomery Reduction
            for i in 0..$limbs {
                let k = r[i].wrapping_mul(P::INV);
                let mut carry = 0;

                fa::mac_with_carry(r[i], k, P::MODULUS.0[0], &mut carry);
                for j in 1..$limbs {
                    r[(j+i)%$limbs] = fa::mac_with_carry(r[(j+i)%$limbs], k, P::MODULUS.0[j], &mut carry);
                }
                r[($limbs+i)%$limbs] = carry;
            }
            for i in 0..$limbs {
                (tmp.0).0[i] = r[i];
            }
            tmp.0
        }
    }
}

macro_rules! impl_field_square_in_place {
    ($limbs: expr) => {
        #[inline]
        #[unroll_for_loops]
        fn square_in_place(&mut self) -> &mut Self {
            let mut r = [0u64; $limbs*2];

            for i in 0..$limbs-1 {
                let mut carry = 0;
                for j in (i+1)..$limbs {
                    r[i+j] = fa::mac_with_carry(r[i+j], (self.0).0[i], (self.0).0[j], &mut carry);
                }
                r[$limbs+i] = carry;
            }

            r[$limbs*2-1] = r[$limbs*2-2] >> 63;
            for i in 0..$limbs*2-3 {
                r[$limbs*2-2-i] = (r[$limbs*2-2-i] << 1) | (r[$limbs*2-3-i] >> 63)
            }
            r[1] = r[1] << 1;

            let mut carry = 0;
            for i in 0..$limbs {
                r[2*i] = fa::mac_with_carry(r[2*i], (self.0).0[i], (self.0).0[i], &mut carry);
                r[2*i+1] = fa::adc(r[2*i+1], 0, &mut carry);
            }
            // Montgomery reduction
            let mut _carry2 = 0;
            for i in 0..$limbs {
                let k = r[i].wrapping_mul(P::INV);
                let mut carry = 0;
                fa::mac_with_carry(r[i], k, P::MODULUS.0[0], &mut carry);
                for j in 1..$limbs {
                    r[j+i] = fa::mac_with_carry(r[j+i], k, P::MODULUS.0[j], &mut carry);
                }
                r[$limbs+i] = fa::adc(r[$limbs+i], _carry2, &mut carry);
                _carry2 = carry;
            }

            for i in 0..$limbs {
                (self.0).0[i] = r[$limbs+i];
            }
            self.reduce();
            self
        }
    }
}

macro_rules! impl_field_bigint_conv {
    ($field: ident, $bigint: ident, $params: ident) => {
        impl<P: $params> Into<$bigint> for $field<P> {
            fn into(self) -> $bigint {
                self.into_repr()
            }
        }

        impl<P: $params> From<$bigint> for $field<P> {
            fn from(int: $bigint) -> Self {
                Self::from_repr(int)
            }
        }
    };
}

macro_rules! impl_prime_field_standard_sample {
    ($field: ident, $params: ident) => {
        impl<P: $params> rand::distributions::Distribution<$field<P>>
            for rand::distributions::Standard
        {
            #[inline]
            fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> $field<P> {
                loop {
                    let mut tmp = $field(rng.sample(rand::distributions::Standard), PhantomData);
                    // Mask away the unused bits at the beginning.
                    tmp.0
                        .as_mut()
                        .last_mut()
                        .map(|val| *val &= core::u64::MAX >> P::REPR_SHAVE_BITS);

                    if tmp.is_valid() {
                        return tmp;
                    }
                }
            }
        }
    };
}

macro_rules! impl_prime_field_from_int {
    ($field: ident, u128, $params: ident) => {
        impl<P: $params> From<u128> for $field<P> {
            fn from(other: u128) -> Self {
                let upper = (other >> 64) as u64;
                let lower = ((other << 64) >> 64) as u64;
                let mut default_int = P::BigInt::default();
                default_int.0[0] = lower;
                default_int.0[1] = upper;
                Self::from_repr(default_int)
            }
        }
    };
    ($field: ident, $int: ident, $params: ident) => {
        impl<P: $params> From<$int> for $field<P> {
            fn from(other: $int) -> Self {
                Self::from_repr(P::BigInt::from(u64::from(other)))
            }
        }
    };
}

macro_rules! sqrt_impl {
    ($Self:ident, $P:tt, $self:expr) => {{
        use crate::fields::LegendreSymbol::*;
        // https://eprint.iacr.org/2012/685.pdf (page 12, algorithm 5)
        // Actually this is just normal Tonelli-Shanks; since `P::Generator`
        // is a quadratic non-residue, `P::ROOT_OF_UNITY = P::GENERATOR ^ t`
        // is also a quadratic non-residue (since `t` is odd).
        match $self.legendre() {
            Zero => Some(*$self),
            QuadraticNonResidue => None,
            QuadraticResidue => {
                let mut z = $Self::qnr_to_t();
                let mut w = $self.pow($P::T_MINUS_ONE_DIV_TWO);
                let mut x = w * $self;
                let mut b = x * &w;

                let mut v = $P::TWO_ADICITY as usize;
                // t = self^t
                #[cfg(debug_assertions)]
                {
                    let mut check = b;
                    for _ in 0..(v - 1) {
                        check.square_in_place();
                    }
                    if !check.is_one() {
                        panic!("Input is not a square root, but it passed the QR test")
                    }
                }

                while !b.is_one() {
                    let mut k = 0usize;

                    let mut b2k = b;
                    while !b2k.is_one() {
                        // invariant: b2k = b^(2^k) after entering this loop
                        b2k.square_in_place();
                        k += 1;
                    }

                    let j = v - k - 1;
                    w = z;
                    for _ in 0..j {
                        w.square_in_place();
                    }

                    z = w.square();
                    b *= &z;
                    x *= &w;
                    v = k;
                }

                Some(x)
            },
        }
    }};
}

// Implements AddAssign on Self by deferring to an implementation on &Self
#[macro_export]
macro_rules! impl_additive_ops_from_ref {
    ($type: ident, $params: ident) => {
        #[allow(unused_qualifications)]
        impl<P: $params> core::ops::Add<Self> for $type<P> {
            type Output = Self;

            #[inline]
            fn add(self, other: Self) -> Self {
                let mut result = self;
                result.add_assign(&other);
                result
            }
        }

        #[allow(unused_qualifications)]
        impl<'a, P: $params> core::ops::Add<&'a mut Self> for $type<P> {
            type Output = Self;

            #[inline]
            fn add(self, other: &'a mut Self) -> Self {
                let mut result = self;
                result.add_assign(&*other);
                result
            }
        }

        #[allow(unused_qualifications)]
        impl<P: $params> core::ops::Sub<Self> for $type<P> {
            type Output = Self;

            #[inline]
            fn sub(self, other: Self) -> Self {
                let mut result = self;
                result.sub_assign(&other);
                result
            }
        }

        #[allow(unused_qualifications)]
        impl<'a, P: $params> core::ops::Sub<&'a mut Self> for $type<P> {
            type Output = Self;

            #[inline]
            fn sub(self, other: &'a mut Self) -> Self {
                let mut result = self;
                result.sub_assign(&*other);
                result
            }
        }

        #[allow(unused_qualifications)]
        impl<P: $params> core::iter::Sum<Self> for $type<P> {
            fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
                iter.fold(Self::zero(), core::ops::Add::add)
            }
        }

        #[allow(unused_qualifications)]
        impl<'a, P: $params> core::iter::Sum<&'a Self> for $type<P> {
            fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
                iter.fold(Self::zero(), core::ops::Add::add)
            }
        }

        #[allow(unused_qualifications)]
        impl<P: $params> core::ops::AddAssign<Self> for $type<P> {
            fn add_assign(&mut self, other: Self) {
                self.add_assign(&other)
            }
        }

        #[allow(unused_qualifications)]
        impl<P: $params> core::ops::SubAssign<Self> for $type<P> {
            fn sub_assign(&mut self, other: Self) {
                self.sub_assign(&other)
            }
        }

        #[allow(unused_qualifications)]
        impl<'a, P: $params> core::ops::AddAssign<&'a mut Self> for $type<P> {
            fn add_assign(&mut self, other: &'a mut Self) {
                self.add_assign(&*other)
            }
        }

        #[allow(unused_qualifications)]
        impl<'a, P: $params> core::ops::SubAssign<&'a mut Self> for $type<P> {
            fn sub_assign(&mut self, other: &'a mut Self) {
                self.sub_assign(&*other)
            }
        }
    };
}

// Implements AddAssign on Self by deferring to an implementation on &Self
#[macro_export]
macro_rules! impl_multiplicative_ops_from_ref {
    ($type: ident, $params: ident) => {
        #[allow(unused_qualifications)]
        impl<P: $params> core::ops::Mul<Self> for $type<P> {
            type Output = Self;

            #[inline]
            fn mul(self, other: Self) -> Self {
                let mut result = self;
                result.mul_assign(&other);
                result
            }
        }

        #[allow(unused_qualifications)]
        impl<P: $params> core::ops::Div<Self> for $type<P> {
            type Output = Self;

            #[inline]
            fn div(self, other: Self) -> Self {
                let mut result = self;
                result.div_assign(&other);
                result
            }
        }

        #[allow(unused_qualifications)]
        impl<'a, P: $params> core::ops::Mul<&'a mut Self> for $type<P> {
            type Output = Self;

            #[inline]
            fn mul(self, other: &'a mut Self) -> Self {
                let mut result = self;
                result.mul_assign(&*other);
                result
            }
        }

        #[allow(unused_qualifications)]
        impl<'a, P: $params> core::ops::Div<&'a mut Self> for $type<P> {
            type Output = Self;

            #[inline]
            fn div(self, other: &'a mut Self) -> Self {
                let mut result = self;
                result.div_assign(&*other);
                result
            }
        }

        #[allow(unused_qualifications)]
        impl<P: $params> core::iter::Product<Self> for $type<P> {
            fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
                iter.fold(Self::one(), core::ops::Mul::mul)
            }
        }

        #[allow(unused_qualifications)]
        impl<'a, P: $params> core::iter::Product<&'a Self> for $type<P> {
            fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
                iter.fold(Self::one(), Mul::mul)
            }
        }

        #[allow(unused_qualifications)]
        impl<P: $params> core::ops::MulAssign<Self> for $type<P> {
            fn mul_assign(&mut self, other: Self) {
                self.mul_assign(&other)
            }
        }

        #[allow(unused_qualifications)]
        impl<'a, P: $params> core::ops::DivAssign<&'a mut Self> for $type<P> {
            fn div_assign(&mut self, other: &'a mut Self) {
                self.div_assign(&*other)
            }
        }

        #[allow(unused_qualifications)]
        impl<'a, P: $params> core::ops::MulAssign<&'a mut Self> for $type<P> {
            fn mul_assign(&mut self, other: &'a mut Self) {
                self.mul_assign(&*other)
            }
        }

        #[allow(unused_qualifications)]
        impl<P: $params> core::ops::DivAssign<Self> for $type<P> {
            fn div_assign(&mut self, other: Self) {
                self.div_assign(&other)
            }
        }
    };
}
