pub mod groupmap;

use algebra::curves::AffineCurve;

pub trait CompressedSRS<G: AffineCurve> {
    fn create(length: usize, num_y_coords: usize) -> Self;
    fn decompress(c: Self) -> Vec<(G::BaseField, G::BaseField)>;
}
