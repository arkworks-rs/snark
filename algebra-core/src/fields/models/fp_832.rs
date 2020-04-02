use crate::{
    biginteger::{arithmetic as fa, BigInteger as _BigInteger, BigInteger832 as BigInteger},
    bytes::{FromBytes, ToBytes},
    fields::{Field, FpParameters, LegendreSymbol, PrimeField, SquareRootField},
    io::{Read, Result as IoResult, Write},
};
use core::{
    cmp::{Ord, Ordering, PartialOrd},
    fmt::{Display, Formatter, Result as FmtResult},
    marker::PhantomData,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    str::FromStr,
};
use num_traits::{One, Zero};
use unroll::unroll_for_loops;

impl_Fp!(Fp832, Fp832Parameters, 13);
