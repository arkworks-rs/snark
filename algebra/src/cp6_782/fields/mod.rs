#[cfg(any(
    feature = "cp6_782",
    feature = "ed_on_cp6_782",
    feature = "ed_on_bw6_761"
))]
pub mod fr;
#[cfg(any(
    feature = "cp6_782",
    feature = "ed_on_cp6_782",
    feature = "ed_on_bw6_761"
))]
pub use self::fr::*;

#[cfg(feature = "cp6_782")]
pub mod fq;
#[cfg(feature = "cp6_782")]
pub use self::fq::*;

#[cfg(feature = "cp6_782")]
pub mod fq3;
#[cfg(feature = "cp6_782")]
pub use self::fq3::*;

#[cfg(feature = "cp6_782")]
pub mod fq6;
#[cfg(feature = "cp6_782")]
pub use self::fq6::*;

#[cfg(all(feature = "cp6_782", test))]
mod tests;
