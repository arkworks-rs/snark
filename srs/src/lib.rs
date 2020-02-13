pub mod groupmap;
#[cfg(test)]
mod tests;

use algebra::curves::AffineCurve;
use core::fmt::Error;

pub trait CompressedSRS<G: AffineCurve> {
    fn create(length: usize, num_y_coords: usize) -> Self;
    fn decompress(c: Self) -> Vec<(G::BaseField, G::BaseField)>;
}
