trait BatchArithmetic<> {

    // Computes [-p, p, -3p, 3p, ..., -2^wp, 2^wp]
    pub fn batch_wnaf_tables<P>(w: usize, a: Vec<GroupAffine<P>>) -> Vec<Vec<G>>;

    // This function consumes the scalars
    // We can make this more generic in the future to use other than u16.
    pub fn batch_wnaf_opcode_recoding<BigInt: BigInteger>(
        mut scalars: Vec<BigInt>,
        w: usize
    ) -> Vec<Vec<Option<u16>>>;

    // This function consumes the second op as it mutates it in place
    // to prevent memory allocation
    pub fn batch_double_in_place_with_edge_cases<'a, F: Field, I, E>(op_iter: I) -> ()
    where I: Iterator<Item = &'a GroupAffine<E::P>>;

    pub fn batch_double_in_place<'a, F: Field, I, E>(op_iter: I) -> ()
    where I: Iterator<Item = &'a GroupAffine<E::P>>;

    pub fn batch_add_in_place_with_edge_cases<'a, F: Field, I, P>(op_iter: I)
    where I: Iterator<Item = (&'a GroupAffine<P>, GroupAffine<P>)>;

    pub fn batch_double_in_place<'a, F: Field, I, E>(op_iter: I) -> ()
    where I: Iterator<Item = &'a GroupAffine<E::P>>;

    pub fn batch_scalar_mul_in_place<G: AffineCurve, BigInt: BigInteger>(
        w: usize,
        mut points: Vec<G>,
        mut scalars: Vec<BigInt>,
    );

    pub fn batch_scalar_mul_in_place_glv<G: AffineCurve, BigInt: BigInteger>(
        w: usize,
        mut points: Vec<G>,
        mut scalars: Vec<BigInt>,
    );
}
