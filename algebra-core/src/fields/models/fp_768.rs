use core::{
    cmp::{Ord, Ordering, PartialOrd},
    fmt::{Display, Formatter, Result as FmtResult},
    marker::PhantomData,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    str::FromStr,
};
use num_traits::{One, Zero};
use unroll::unroll_for_loops;

use crate::{
    biginteger::{arithmetic as fa, BigInteger as _BigInteger, BigInteger768 as BigInteger},
    bytes::{FromBytes, ToBytes},
    fields::{FftField, Field, FpParameters, LegendreSymbol, PrimeField, SquareRootField},
    io::{Read, Result as IoResult, Write},
};

impl_Fp!(Fp768, Fp768Parameters, 12);
