pub fn batch_scalar_mul_in_place<BigInt: BigInteger>(
    w: usize,
    mut points: Vec<GroupAffine<P>>,
    mut scalars: Vec<BigInt>,
) {
    let no_op = 1 << (w + 1); // noop is encoded as window_size
    let opcode_vectorised = batch_wnaf_opcode_recoding(w, scalars);
    let tables = batch_wnaf_tables(w, points);

    for opcode_row in opcode_vectorised.rev() {
        let double_iterator = opcode_row.zip(points)
            .filter(|op| op.is_some())
            .map(|op, p| p);

        batch_double_in_place_with_edge_cases(double_iterator);

        let add_iterator = opcode_row.zip(points, tables)
            .filter(|op| op.is_some() && op != no_op)
            .map(|op, p, t| (p, t[op]));

        batch_add_in_place_with_edge_cases(add_iterator);
    }
}

pub fn batch_scalar_mul_in_place_glv<BigInt: BigInteger>(
    w: usize,
    mut points: Vec<GroupAffine<P>>,
    mut scalars: Vec<BigInt>,
) {
    assert_eq!(points.len(), scalars.len());
    let batch_size = points.len();
    let mut k1 = Vec::with_capacity(batch_size);
    let mut k2 = Vec::with_capacity(batch_size);

    let p2 = points.map(|p| p.glv_endomorphism());
    batch_scalar_mul_in_place(w, points, k1);
    batch_scalar_mul_in_place(w, p2, k2);
    batch_add_in_place_with_edge_cases(points, p2);
}
