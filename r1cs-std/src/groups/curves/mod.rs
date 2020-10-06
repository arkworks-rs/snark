/// This module generically implements arithmetic for Short
/// Weierstrass elliptic curves by following the complete formulae of
/// [[Renes, Costello, Batina 2015]](https://eprint.iacr.org/2015/1060).
pub mod short_weierstrass;

/// This module generically implements arithmetic for Twisted
/// Edwards elliptic curves by following the complete formulae described in the
/// [EFD](https://www.hyperelliptic.org/EFD/g1p/auto-twisted.html).
pub mod twisted_edwards;
