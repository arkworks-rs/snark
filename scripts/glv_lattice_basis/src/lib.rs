extern crate algebra;
extern crate algebra_core;
extern crate num_traits;

mod arithmetic;

use algebra_core::{BigInteger, Field, PrimeField, ProjectiveCurve};
pub use arithmetic::*;
use num_traits::Zero;
use std::ops::Neg;

/// Takes data from two endomorphisms and sorts out which corresponds to which
fn which_endo<G: ProjectiveCurve>(
    base_roots: (G::BaseField, G::BaseField),
    scalar_roots: (G::ScalarField, G::ScalarField),
) -> (
    (G::BaseField, G::ScalarField),
    (G::BaseField, G::ScalarField),
) {
    // println!("{:?}, {:?}", base_roots, scalar_roots);
    let g = G::prime_subgroup_generator();

    let mut g_endo = g;
    *g_endo.get_x() *= &base_roots.0;

    let d1 = if g.mul(scalar_roots.0) == g_endo {
        (base_roots.0, scalar_roots.0)
    } else {
        let mut g_endo = g;
        *g_endo.get_x() *= &base_roots.1;
        assert!(g.mul(scalar_roots.0) == g_endo);

        (base_roots.1, scalar_roots.0)
    };

    let d2 = if g.mul(scalar_roots.1) == g_endo {
        (base_roots.0, scalar_roots.1)
    } else {
        let mut g_endo = g;
        *g_endo.get_x() *= &base_roots.1;
        assert!(g.mul(scalar_roots.1) == g_endo);

        (base_roots.1, scalar_roots.1)
    };

    (d1, d2)
}

fn cube_root_unity<F: Field, B: BigInteger>() -> (F, F) {
    let char = B::from_slice(F::characteristic());
    let deg = F::extension_degree();
    let mut modulus = char;
    for _ in 1..deg {
        modulus = B::mul_no_reduce_lo(&modulus.as_ref(), &char.as_ref());
    }

    modulus.sub_noborrow(&B::from(1));
    let (q, r) = div_with_remainder(modulus, B::from(3));
    assert!(r == B::from(0));

    let mut g = 2u32;
    let mut root1 = F::one();
    loop {
        if root1 != F::one() {
            break;
        }
        let x = F::from(g);
        root1 = x.pow(q);
        g += 1;
    }
    let root2 = root1 * root1;
    assert!(root1.pow(&[3]) == F::one());
    assert!(root2.pow(&[3]) == F::one());
    assert!(root1 != root2);

    (root1, root2)
}

fn get_endo_data<G: ProjectiveCurve, B: BigInteger>() -> (G::BaseField, G::ScalarField) {
    // println!("{:?}", which_endo::<G>(
    //     cube_root_unity::<G::BaseField>(),
    //     cube_root_unity::<G::ScalarField>(),
    // ));
    which_endo::<G>(
        cube_root_unity::<G::BaseField, B>(),
        cube_root_unity::<G::ScalarField, <G::ScalarField as PrimeField>::BigInt>(),
    )
    .1
}

pub fn print_glv_params<G: ProjectiveCurve, WideBigInt: BigInteger, BaseFieldBigInt: BigInteger>() {
    let (omega, lambda) = get_endo_data::<G, BaseFieldBigInt>();
    let g = G::prime_subgroup_generator();
    let mut g_endo = g;
    *g_endo.get_x() *= &omega;
    assert!(g.mul(lambda) == g_endo);

    println!("const OMEGA: Self::BaseField = {:?};", omega);
    let n = <G::ScalarField as PrimeField>::modulus();
    println!("const LAMBDA: Self::ScalarField = {:?};", lambda);

    let vecs = get_lattice_basis::<G::ScalarField>(n, lambda.into_repr());

    // We check that `(|B1| + 2) * (|B2| + 2) <  2n`
    // and `B_i^2 < 2n` e.g. `|B_i| < \sqrt{2n}$
    // We use this to prove some bounds later
    let wide_modulus = WideBigInt::from_slice(&n.as_ref()[..]);
    let two_modulus = WideBigInt::mul_no_reduce_lo(
        &wide_modulus.as_ref()[..],
        &WideBigInt::from(2).as_ref()[..],
    );

    let mut b1 = ((vecs.0).1).1;
    let mut b2 = ((vecs.1).1).1;
    let two = <G::ScalarField as PrimeField>::BigInt::from(2);
    let b1b1 = WideBigInt::mul_no_reduce(&b1.as_ref()[..], &b1.as_ref()[..]);
    let b2b2 = WideBigInt::mul_no_reduce(&b2.as_ref()[..], &b2.as_ref()[..]);

    b1.add_nocarry(&two);
    b2.add_nocarry(&two);
    let b1b2 = WideBigInt::mul_no_reduce(&b1.as_ref()[..], &b2.as_ref()[..]);

    assert!(b1b1 < two_modulus);
    assert!(b2b2 < two_modulus);
    assert!(b1b2 < two_modulus);

    for (i, vec) in [vecs.0, vecs.1].iter().enumerate() {
        let (s1, (flag, t1)) = vec;

        let mut t1_big = WideBigInt::from_slice(t1.as_ref());
        let n_big = WideBigInt::from_slice(n.as_ref());
        t1_big.muln(<G::ScalarField as PrimeField>::BigInt::NUM_LIMBS as u32 * 64);
        let (g1_big, _) = div_with_remainder::<WideBigInt>(t1_big, n_big);
        let g1 = <G::ScalarField as PrimeField>::BigInt::from_slice(g1_big.as_ref());

        println!("/// |round(B{} * R / n)|", i + 1);
        println!(
            "const Q{}: <Self::ScalarField as PrimeField>::BigInt = {:?};",
            ((i + 1) % 2) + 1,
            g1
        );
        println!(
            "const B{}: <Self::ScalarField as PrimeField>::BigInt = {:?};",
            i + 1,
            t1
        );
        println!("const B{}_IS_NEG: bool = {:?};", i + 1, flag);

        debug_assert_eq!(
            recompose_integer(
                G::ScalarField::from_repr(*s1).unwrap(),
                if !flag {
                    G::ScalarField::from_repr(*t1).unwrap()
                } else {
                    G::ScalarField::from_repr(*t1).unwrap().neg()
                },
                lambda
            ),
            G::ScalarField::zero()
        );
    }
    println!(
        "const R_BITS: u32 = {:?};",
        <G::ScalarField as PrimeField>::BigInt::NUM_LIMBS * 64
    );
}

// We work on arrays of size 3
// We assume that |E(F_q)| < R = 2^{ceil(limbs/2) * 64}
pub fn get_lattice_basis<F: PrimeField>(
    n: F::BigInt,
    lambda: F::BigInt,
) -> (
    (F::BigInt, (bool, F::BigInt)),
    (F::BigInt, (bool, F::BigInt)),
) {
    let mut r = [n, lambda, n];
    let one = F::one();
    let zero = F::zero();
    let mut t: [F; 3] = [zero, one, zero];
    let max_num_bits_lattice = (F::BigInt::from_slice(F::characteristic()).num_bits() - 1) / 2 + 1;

    // We can use an approximation as we are merely using a heuristic. We should
    // check that the parameters obtained from this heuristic satisfies the
    // required conditions separately.
    let sqrt_n = as_f64(n.as_ref()).sqrt();

    // println!("Log sqrtn: {}", sqrt_n.log2());

    let mut i = 0;
    // While r_i >= sqrt(n), we perform the extended euclidean algorithm so that
    // si*n + ti*lambda = ri then return the vectors (r_i, (sign(t_i), |t_i|)),
    // (r_i+1, (sign(t_i+1), |t_i+1|)) Notice this makes ri + (-ti)*lambda = 0
    // mod n, which is what we desire for our short lattice basis
    while as_f64(r[(i + 1) % 3].as_ref()) >= sqrt_n {
        // while i < 20 {
        let (q, rem): (F::BigInt, F::BigInt) =
            div_with_remainder::<F::BigInt>(r[i % 3], r[(i + 1) % 3]);
        r[(i + 2) % 3] = rem;
        let int_q = F::from_repr(q).unwrap();
        t[(i + 2) % 3] = t[i % 3] - int_q * (t[(i + 1) % 3]);

        i += 1;
    }
    let just_computed = (i + 1) % 3;
    // We reverse the signs due to s_i*n = r_i - t_i*LAMBDA
    let (neg_flag1, t1) = if t[just_computed].into_repr().num_bits() <= max_num_bits_lattice {
        (true, t[just_computed].into_repr())
    } else {
        (false, t[just_computed].neg().into_repr())
    };
    let vec_1 = (r[just_computed], (neg_flag1, t1));

    let prev = i % 3;
    let (neg_flag2, t2) = if t[prev].into_repr().num_bits() <= max_num_bits_lattice {
        (true, t[prev].into_repr())
    } else {
        (false, t[prev].neg().into_repr())
    };
    let vec_2 = (r[prev], (neg_flag2, t2));

    (vec_1, vec_2)
}

pub fn recompose_integer<F: PrimeField>(k1: F, k2: F, lambda: F) -> F {
    k1 + &(k2 * &lambda)
}

fn as_f64(bigint_ref: &[u64]) -> f64 {
    let mut n_float: f64 = 0.0;
    for (i, limb) in bigint_ref.iter().enumerate() {
        n_float += (*limb as f64) * 2f64.powf((i as f64) * 64f64)
    }
    n_float
}
