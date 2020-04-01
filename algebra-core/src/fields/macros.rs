macro_rules! impl_mul_assign {
    ($n:expr, $nm:expr) => {
        // Declare register variables
        // unroll!(6, |k, kp, km| registers!(k, km));
        impl<'a, P: Fp384Parameters> MulAssign<&'a Self> for Fp384<P> {
            #[inline]
            #[unroll_for_loops]
            fn mul_assign(&mut self, other: &Self) {
                let mut c1 = 0u64;
                let mut c2 = 0u64;
                unroll!($n, |i, ip, im| mul_assign_outer_loop!(self, k, $nm, c1, c2));
                // Assign temp register values to self
                unroll!($n, |i, ip, im| (self.0).0[im] = format!("r{}", im));
                self.reduce();
            }
        }
    }
}

macro_rules! mul_assign_outer_loop {
    ($self:ident, $i:expr, $nm:expr, $c1:ident, $c2:ident) => {
        let zero_r = fa::mac(r[0], ($self.0).0[0], (other.0).0[$i], &mut carry1);
        let m = zero_r.wrapping_mul(P::INV);
        fa::mac_discard(r[0], m, P::MODULUS.0[0], &mut carry2);
        unroll!(6, |j, jp, jm| mul_assign_inner_loop!($self, $i, jp, $c1, $c2));
        let format!("r{}", $nm) = $c1 + $c2;
    }
}

macro_rules! mul_assign_inner_loop {
    ($self:ident, $k:expr, $jp:expr, $c1:ident, $c2:ident) => {
        let format!("r{}", $jp) = fa::mac_with_carry(format!("r{}", $jp), ($self.0).0[$jp], (other.0).0[$k], &mut $c2);
        let format!("r{}", $j) = fa::mac_with_carry(format!("r{}", $jp), k, P::MODULUS.0[$j], &mut $c1);
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

macro_rules! registers {
    ($i:expr, $im:expr) => {
        mashup! {
            m["r0"] = r
            m[format!("r{}", $i)] = format!("r{}", $im) r
        }
    }
}

macro_rules! unroll {
    (0, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {};
    (1, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ let $i: usize = 1; let $i_plus: usize = 2; let $i_minus: usize = 0; $s; }};
    (2, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(1, |$i, $i_plus, $i_minus| $s); let $i: usize = 2; let $i_plus: usize = 3; let $i_minus: usize = 1; $s; }};
    (3, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(2, |$i, $i_plus, $i_minus| $s); let $i: usize = 3; let $i_plus: usize = 4; let $i_minus: usize = 2; $s; }};
    (4, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(3, |$i, $i_plus, $i_minus| $s); let $i: usize = 4; let $i_plus: usize = 5; let $i_minus: usize = 3; $s; }};
    (5, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(4, |$i, $i_plus, $i_minus| $s); let $i: usize = 5; let $i_plus: usize = 6; let $i_minus: usize = 4; $s; }};
    (6, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(5, |$i, $i_plus, $i_minus| $s); let $i: usize = 6; let $i_plus: usize = 7; let $i_minus: usize = 5; $s; }};
    (7, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(6, |$i, $i_plus, $i_minus| $s); let $i: usize = 7; let $i_plus: usize = 8; let $i_minus: usize = 6; $s; }};
    (8, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(7, |$i, $i_plus, $i_minus| $s); let $i: usize = 8; let $i_plus: usize = 9; let $i_minus: usize = 7; $s; }};
    (9, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(8, |$i, $i_plus, $i_minus| $s); let $i: usize = 9; let $i_plus: usize = 10; let $i_minus: usize = 8; $s;}};
    (10, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(9, |$i, $i_plus, $i_minus| $s); let $i: usize = 10; let $i_plus: usize = 11; let $i_minus: usize = 9; $s;}};
    (11, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(10, |$i, $i_plus, $i_minus| $s); let $i: usize = 11; let $i_plus: usize = 12; let $i_minus: usize = 10; $s;}};
    (13, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(12, |$i, $i_plus, $i_minus| $s); let $i: usize = 13; let $i_plus: usize = 14; let $i_minus: usize = 12; $s;}};
    (12, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(11, |$i, $i_plus, $i_minus| $s); let $i: usize = 12; let $i_plus: usize = 13; let $i_minus: usize = 11; $s;}};
    (14, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(13, |$i, $i_plus, $i_minus| $s); let $i: usize = 14; let $i_plus: usize = 15; let $i_minus: usize = 13; $s;}};
    (15, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(14, |$i, $i_plus, $i_minus| $s); let $i: usize = 15; let $i_plus: usize = 16; let $i_minus: usize = 14; $s;}};
    (16, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(15, |$i, $i_plus, $i_minus| $s); let $i: usize = 16; let $i_plus: usize = 17; let $i_minus: usize = 15; $s;}};
    (17, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(16, |$i, $i_plus, $i_minus| $s); let $i: usize = 17; let $i_plus: usize = 18; let $i_minus: usize = 16; $s;}};
    (18, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(17, |$i, $i_plus, $i_minus| $s); let $i: usize = 18; let $i_plus: usize = 19; let $i_minus: usize = 17; $s;}};
    (19, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(18, |$i, $i_plus, $i_minus| $s); let $i: usize = 19; let $i_plus: usize = 20; let $i_minus: usize = 18; $s;}};
    (20, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(19, |$i, $i_plus, $i_minus| $s); let $i: usize = 20; let $i_plus: usize = 21; let $i_minus: usize = 19; $s;}};
    (21, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(20, |$i, $i_plus, $i_minus| $s); let $i: usize = 21; let $i_plus: usize = 22; let $i_minus: usize = 20; $s;}};
    (22, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(21, |$i, $i_plus, $i_minus| $s); let $i: usize = 22; let $i_plus: usize = 23; let $i_minus: usize = 21; $s;}};
    (23, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(22, |$i, $i_plus, $i_minus| $s); let $i: usize = 23; let $i_plus: usize = 24; let $i_minus: usize = 22; $s;}};
    (24, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(23, |$i, $i_plus, $i_minus| $s); let $i: usize = 24; let $i_plus: usize = 25; let $i_minus: usize = 23; $s;}};
    (25, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(24, |$i, $i_plus, $i_minus| $s); let $i: usize = 25; let $i_plus: usize = 26; let $i_minus: usize = 24; $s;}};
    (26, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(25, |$i, $i_plus, $i_minus| $s); let $i: usize = 26; let $i_plus: usize = 27; let $i_minus: usize = 25; $s;}};
    (27, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(26, |$i, $i_plus, $i_minus| $s); let $i: usize = 27; let $i_plus: usize = 28; let $i_minus: usize = 26; $s;}};
    (28, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(27, |$i, $i_plus, $i_minus| $s); let $i: usize = 28; let $i_plus: usize = 29; let $i_minus: usize = 27; $s;}};
    (29, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(28, |$i, $i_plus, $i_minus| $s); let $i: usize = 29; let $i_plus: usize = 30; let $i_minus: usize = 28; $s;}};
    (30, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(29, |$i, $i_plus, $i_minus| $s); let $i: usize = 30; let $i_plus: usize = 31; let $i_minus: usize = 29; $s;}};
    (31, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(30, |$i, $i_plus, $i_minus| $s); let $i: usize = 31; let $i_plus: usize = 32; let $i_minus: usize = 30; $s;}};
    (32, |$i:ident, $i_plus:ident, $i_minus:ident| $s:stmt) => {{ unroll!(31, |$i, $i_plus, $i_minus| $s); let $i: usize = 32; let $i_plus: usize = 33; let $i_minus: usize = 31; $s;}};
}
