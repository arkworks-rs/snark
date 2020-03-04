mod algebra::AffineCurve;
mod groupmap::GroupMap;

#[derive(Copy, Clone, Debug)]
pub struct CompressedSRS<G: AffineCurve> {
    pub y_coords: Vec<G::BaseField>,
    pub x_bits: usize, // could be [bool; 2] ?
    pub y_bits: usize,
}

// first impl for num_y_coords = length
// because how do you ... should it be num_x_coords?
impl CompressedSRS {
    fn create(length: usize, num_y_coords: usize) -> Self {
        let map_params = GroupMap::setup();
        let ivec: Vec<_> = (0..length).filter_map(|f| f.ok())
            .map(|i| G::BaseField::From(i))
            .collect();
        let xyvec = GroupMap::batch_to_group_x(map_params, &ivec);


    }
    fn decompress(c: Self) -> Vec<(G::BaseField, G::BaseField)> {}
}
