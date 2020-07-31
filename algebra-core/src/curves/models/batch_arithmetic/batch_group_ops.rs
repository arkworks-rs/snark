use zexe_algebra_core::fields::Field;
use zexe_algebra_core::curves::short_weierstrass_jacobian::GroupAffine;

// This implementation takes particular care to make most use of points fetched from memory
// And to reuse memory to prevent reallocations
// It is directly adapted from Aztec's code.

// https://github.com/AztecProtocol/barretenberg/blob/standardplonkjson/barretenberg/src/
// aztec/ecc/curves/bn254/scalar_multiplication/scalar_multiplication.cpp

pub fn batch_double_in_place_with_edge_cases<'a, F: Field, I, E>(op_iter: I) -> ()
where
    I: Iterator<Item = &'a GroupAffine<E::P>>,
{
    let mut inversion_tmp = F::one();
    let mut scratch_space = Vec::with_capacity(op_iter.size_hint().0);
    // We run two loops over the data separated by an inversion
    for a in op_iter {
        if !a.is_zero() {
            if a.y.is_zero() {
                a.infinity = true;
            } else {
                let x_sq = a.x.square();
                let x_sq_3 = x_sq.double_in_place() + x_sq; // numerator = 3x^2
                scratch_space.push(x_sq_3 * inversion_tmp); // 3x^2 * tmp
                inversion_tmp *= a.x.double(); // update tmp
            }
        }
    }

    inversion_tmp.invert().unwrap(); // this is always in Fp*

    for a in op_iter.rev() {
        if !a.is_zero() {
            let lambda = scratch_space.pop() * inversion_tmp;
            inversion_tmp *= a.x; // Remove the top layer of the denominator

            // x3 = l^2 + 2x
            let x3 = lambda.square_in_place() + a.x.double();
            // y3 = l*(x - x3) - y
            a.y = lambda * (a.x - x3) - a.y;
            a.x = x3;
        }
    }
}

// May not be secure...
pub fn batch_double_in_place<'a, F: Field, I, E>(op_iter: I) -> ()
where
    I: Iterator<Item = &'a GroupAffine<E::P>>,
{
    let mut inversion_tmp = F::one();
    let mut scratch_space = Vec::with_capacity(op_iter.size_hint().0);
    // We run two loops over the data separated by an inversion
    for a in op_iter {
        let x_sq = a.x.square();
        let x_sq_3 = x_sq.double_in_place() + x_sq; // numerator = 3x^2
        scratch_space.push(x_sq_3 * inversion_tmp); // 3x^2 * tmp
        inversion_tmp *= a.x.double(); // update tmp
    }

    inversion_tmp.invert().unwrap(); // this is always in Fp*

    for a in op_iter.rev() {
        let lambda = scratch_space.pop() * inversion_tmp;
        inversion_tmp *= a.x; // Remove the top layer of the denominator

        // x3 = l^2 + 2x
        let x3 = lambda.square_in_place() + a.x.double();
        // y3 = l*(x - x3) - y
        a.y = lambda * (a.x - x3) - a.y;
        a.x = x3;
    }
}

// This implementation takes particular care to make most use of points fetched from memory
// And to reuse memory to prevent reallocations
// It is directly adapted from Aztec's code.

// https://github.com/AztecProtocol/barretenberg/blob/standardplonkjson/barretenberg/src/
// aztec/ecc/curves/bn254/scalar_multiplication/scalar_multiplication.cpp

// This function consumes the second op as it mutates it in place
// to prevent memory allocation
pub fn batch_add_in_place_with_edge_cases<'a, F: Field, I, P>(op_iter: I)
where
    I: Iterator<Item = (&'a GroupAffine<P>, GroupAffine<P>)>,
{
    let mut inversion_tmp = F::one();
    let half = F::from_repr(P::MODULUS_MINUS_ONE_DIV_TWO) + F::one(); // (p + 1)/2 * 2 = 1
    // We run two loops over the data separated by an inversion
    for (a, b) in op_iter {
        if a.is_zero() || b.is_zero() {
            continue;
        } else if a.x == b.x {
            // double.
            // In our model, we consider self additions rare.
            // So we consider it inconsequential to make them more expensive
            // This costs 1 modular mul more than a standard squaring
            if a.y == b.y {
                let x_sq = b.x.square();
                b.x -= b.y; // x - y
                a.x = b.y.double(); // denominator = 2y
                a.y = x_sq.double_in_place() + x_sq; // numerator = 3x^2
                b.y -= half * a.y; // y - 3x^2/2
                a.y *= inversion_tmp; // 3x^2 * tmp
                inversion_tmp *= a.x; // update tmp
            } else {
                // No inversions take place if either operand is zero
                a.infinity = true;
                b.infinity = true;
            }
        } else {
            a.x -= b.x; // denominator = x1 - x2. We can recover x1 + x2 from this. Note this is never 0.
            a.y -= b.y; // numerator = y1 - y2
            a.y *= inversion_tmp; // (y1 - y2)*tmp
            inversion_tmp *= a.x // update tmp
        }
    }

    inversion_tmp.invert().unwrap(); // this is always in Fp*

    for (a, b) in op_iter.rev() {
        if a.is_zero() {
            a = b;
        } else if !b.is_zero() {
            let lambda = a.y * inversion_tmp;
            inversion_tmp *= a.x; // Remove the top layer of the denominator

            // x3 = l^2 + x1 + x2 or for squaring: 2y + l^2 + 2x - 2y = l^2 + 2x
            a.x += lambda.square_in_place() + b.x.double();
            // y3 = l*(x2 - x3) - y2 or for squaring: 3x^2/2y(x - y - x3) - (y - 3x^2/2) = l*(x - x3) - y
            a.y = lambda * (b.x - a.x) - b.y;
        }
    }
}

// This function consumes b_vec as it mutates it in place
// to prevent memory allocations
pub fn batch_add_in_place<'a, F: Field, I, P>(op_iter: I)
where
    I: Iterator<Item = (&'a GroupAffine<P>, GroupAffine<P>)>,
{
    let mut inversion_tmp = F::one();
    // We run two loops over the data separated by an inversion
    // let mut scratch_space = Vec::<AffineElement>::with_capacity(a_vec.len());
    for (a, b) in op_iter {
        a.x -= b.x; // denominator = x1 - x2. We can recover x1 + x2 from this. Note this is never 0.
        a.y -= b.y; // numerator = y1 - y2
        a.y *= inversion_tmp; // (y1 - y2)*tmp
        inversion_tmp *= a.x // update tmp
    }

    inversion_tmp.invert().unwrap(); // this is always in Fp*

    for (a, b) in op_iter.rev() {
        let lambda = a.y * inversion_tmp;
        inversion_tmp *= a.x; // Remove the top layer of the denominator
        a.x += lambda.square_in_place() + b.x.double();   // x3 = l^2 + x1 + x2
        a.y = lambda * (b.x - a.x) - b.y;     // y3 = l*(x2 - x3) - y2
    }
}
