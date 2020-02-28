mod groupmap;

#[derive(Copy, Clone, Debug)]
pub struct CompressedSRS<G: AffineCurve> {
    pub y_coords: Vec<G::BaseField>,
    pub x_bits: usize, // could be [bool; 2] ?
}

impl CompressedSRS {
    fn create(length: usize, num_y_coords: usize) -> Self {}
    fn decompress(c: Self) -> Vec<(G::BaseField, G::BaseField)> {}
}
