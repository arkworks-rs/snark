use std::{
    cmp::{Ord, Ordering, PartialOrd},
    fmt::{Display, Formatter, Result as FmtResult},
    io::{Read, Result as IoResult, Write, Error as IoError, ErrorKind},
    marker::PhantomData,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    str::FromStr,
};

use crate::{
    biginteger::{
        arithmetic as fa, BigInteger as _BigInteger, BigInteger768 as BigInteger
    },
    bytes::{
        FromBytes, ToBytes
    },
    fields::{
        Field, FpParameters, LegendreSymbol, PrimeField, SquareRootField
    },
    MulShort, MulShortAssign,
};

use unroll::unroll_for_loops;

impl_Fp!(Fp768, Fp768Parameters, 12);