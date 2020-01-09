use super::Group;
use crate::{AffineCurve, PairingEngine};
use crate::fields::Field;
use crate::bytes::{ToCompressed, FromCompressed, FromBytes};
use crate::UniformRand;
use rand::SeedableRng;
use rand_xorshift::XorShiftRng;

pub fn group_test<G: Group>(a: G, mut b: G) {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    let zero = G::zero();
    let fr_zero = G::ScalarField::zero();
    let fr_one = G::ScalarField::one();
    let fr_two = fr_one + &fr_one;
    assert_eq!(zero, zero);
    assert_eq!(zero.is_zero(), true);
    assert_eq!(a.mul(&fr_one), a);
    assert_eq!(a.mul(&fr_two), a + &a);
    assert_eq!(a.mul(&fr_zero), zero);
    assert_eq!(a.mul(&fr_zero) - &a, -a);
    assert_eq!(a.mul(&fr_one) - &a, zero);
    assert_eq!(a.mul(&fr_two) - &a, a);

    // a == a
    assert_eq!(a, a);
    // a + 0 = a
    assert_eq!(a + &zero, a);
    // a - 0 = a
    assert_eq!(a - &zero, a);
    // a - a = 0
    assert_eq!(a - &a, zero);
    // 0 - a = -a
    assert_eq!(zero - &a, -a);
    // a.double() = a + a
    assert_eq!(a.double(), a + &a);
    // b.double() = b + b
    assert_eq!(b.double(), b + &b);
    // a + b = b + a
    assert_eq!(a + &b, b + &a);
    // a - b = -(b - a)
    assert_eq!(a - &b, -(b - &a));
    // (a + b) + a = a + (b + a)
    assert_eq!((a + &b) + &a, a + &(b + &a));
    // (a + b).double() = (a + b) + (b + a)
    assert_eq!((a + &b).double(), (a + &b) + &(b + &a));

    // Check that double_in_place and double give the same result
    let original_b = b;
    b.double_in_place();
    assert_eq!(original_b.double(), b);

    let fr_rand1 = G::ScalarField::rand(&mut rng);
    let fr_rand2 = G::ScalarField::rand(&mut rng);
    let a_rand1 = a.mul(&fr_rand1);
    let a_rand2 = a.mul(&fr_rand2);
    let fr_three = fr_two + &fr_rand1;
    let a_two = a.mul(&fr_two);
    assert_eq!(a_two, a.double(), "(a * 2)  != a.double()");
    let a_six = a.mul(&(fr_three * &fr_two));
    assert_eq!(a_two.mul(&fr_three), a_six, "(a * 2) * 3 != a * (2 * 3)");

    assert_eq!(
        a_rand1.mul(&fr_rand2),
        a_rand2.mul(&fr_rand1),
        "(a * r1) * r2 != (a * r2) * r1"
    );
    assert_eq!(
        a_rand2.mul(&fr_rand1),
        a.mul(&(fr_rand1 * &fr_rand2)),
        "(a * r2) * r1 != a * (r1 * r2)"
    );
    assert_eq!(
        a_rand1.mul(&fr_rand2),
        a.mul(&(fr_rand1 * &fr_rand2)),
        "(a * r1) * r2 != a * (r1 * r2)"
    );
}

pub fn compression_test<G: AffineCurve>(even: G, odd: G) {

    //Test correct compression/de-compression of a non-zero point with even y
    let even_compressed = even.compress();
    let even_len = even_compressed.len();

    let infinity_flag_set = bool::read([(even_compressed[even_len - 1] >> 7) & 1].as_ref()).unwrap();
    assert!(!infinity_flag_set);
    let parity_flag_set = bool::read([(even_compressed[even_len - 1] >> 6) & 1].as_ref()).unwrap();
    assert!(!parity_flag_set);

    let even_decompressed = G::decompress(even_compressed.clone()).unwrap();
    assert_eq!(even, even_decompressed);

    //Test correct compression/de-compression of a non-zero point with odd y
    let odd_compressed = odd.compress();
    let odd_len = odd_compressed.len();

    let infinity_flag_set = bool::read([(odd_compressed[odd_len - 1] >> 7) & 1].as_ref()).unwrap();
    assert!(!infinity_flag_set);
    let parity_flag_set = bool::read([(odd_compressed[odd_len - 1] >> 6) & 1].as_ref()).unwrap();
    assert!(parity_flag_set);

    let odd_decompressed = G::decompress(odd_compressed.clone()).unwrap();
    assert_eq!(odd, odd_decompressed);

    //Test correct compression/decompression of a zero point
    let z = G::zero();
    let z_compressed = z.compress();
    let z_len = z_compressed.len();

    let infinity_flag_set = bool::read([(z_compressed[z_len - 1] >> 7) & 1].as_ref()).unwrap();
    assert!(infinity_flag_set);
    // When the point is zero, parity flag is set to zero too.
    let parity_flag_set = bool::read([(z_compressed[z_len - 1] >> 6) & 1].as_ref()).unwrap();
    assert!(!parity_flag_set);

    let z_decompressed = G::decompress(z_compressed).unwrap();
    assert_eq!(z, z_decompressed);
}

pub fn gt_compression_test<E: PairingEngine>(even: E::Fqk, odd: E::Fqk)
{
    //Test correct compression/de-compression of a non-zero point with even c0
    let even_compressed = even.compress();
    let even_len = even_compressed.len();

    let parity_flag_set = bool::read([(even_compressed[even_len - 1] >> 7) & 1].as_ref()).unwrap();
    assert!(!parity_flag_set);

    let even_decompressed = E::Fqk::decompress(even_compressed.clone()).unwrap();
    assert_eq!(even, even_decompressed);

    //Test correct compression/de-compression of a non-zero point with odd c0
    let odd_compressed = odd.compress();
    let odd_len = odd_compressed.len();

    let parity_flag_set = bool::read([(odd_compressed[odd_len - 1] >> 7) & 1].as_ref()).unwrap();
    assert!(parity_flag_set);

    let odd_decompressed = E::Fqk::decompress(odd_compressed.clone()).unwrap();
    assert_eq!(odd, odd_decompressed);
}