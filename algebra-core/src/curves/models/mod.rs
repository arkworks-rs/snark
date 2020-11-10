use crate::fields::{Field, PrimeField, SquareRootField};

pub mod bls12;
pub mod bn;
pub mod bw6;
pub mod mnt4;
pub mod mnt6;

#[macro_use]
pub(crate) mod sw_batch_affine;
#[macro_use]
pub mod short_weierstrass_affine;
#[macro_use]
pub mod short_weierstrass_jacobian;
pub mod twisted_edwards_extended;

pub use short_weierstrass_jacobian::SWModelParameters;
pub use twisted_edwards_extended::TEModelParameters;

pub trait ModelParameters: Send + Sync + 'static {
    type BaseField: Field + SquareRootField;
    type ScalarField: PrimeField
        + SquareRootField
        + Into<<Self::ScalarField as PrimeField>::BigInt>
        + From<<Self::ScalarField as PrimeField>::BigInt>;
}

pub trait MontgomeryModelParameters: ModelParameters {
    const COEFF_A: Self::BaseField;
    const COEFF_B: Self::BaseField;

    type TEModelParameters: TEModelParameters<BaseField = Self::BaseField>;
}
