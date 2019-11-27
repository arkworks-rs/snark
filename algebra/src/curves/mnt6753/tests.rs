use crate::{curves::{
    mnt6753::{G1Affine, G1Projective, G2Affine, G2Projective, MNT6},
    tests::curve_tests,
    AffineCurve, PairingEngine,
}, fields::mnt6753::fr753b::Fr, groups::tests::group_test, ToCompressed, FromCompressed, ToBytes, to_bytes, ProjectiveCurve, Field};
use rand;

#[test]
fn test_g1_projective_curve() {
    curve_tests::<G1Projective>();
}

#[test]
fn test_g1_projective_group() {
    let a: G1Projective = rand::random();
    let b: G1Projective = rand::random();
    group_test(a, b);
}

#[test]
fn test_g1_generator() {
    let generator = G1Affine::prime_subgroup_generator();
    assert!(generator.is_on_curve());
    assert!(generator.is_in_correct_subgroup_assuming_on_curve());
}

#[test]
fn test_g1_compression_decompression() {
    let a: G1Projective = rand::random();
    let a = a.into_affine();
    let a_len = to_bytes!(a).unwrap().len();

    //Test correct compression/de-compression of a random point
    let a_compressed = a.compress();
    let a_decompressed = G1Affine::decompress(a_compressed.clone()).unwrap();
    assert_eq!(a, a_decompressed);
    assert_eq!(a_compressed.len(), (a_len/2));

    //Test correct compression/decompression of a zero point
    let b = G1Affine::zero();
    let b_compressed = b.compress();
    let b_decompressed = G1Affine::decompress(b_compressed).unwrap();
    assert_eq!(b, b_decompressed);

    //Test wrong compression/de-compression of a random point by masking with 0s a random byte
    let mut a_compressed_modified = a_compressed.clone();
    loop {
        let index: usize = rand::random();
        let max_idx = (a_len/2) - 4;
        if a_compressed_modified[index % max_idx] != 0u8 {
            a_compressed_modified[index % max_idx] &= 0u8;
            break;
        }
    }
    let a_decompressed_modified = G1Affine::decompress(a_compressed_modified);
    if a_decompressed_modified.is_some() {
        //a_decompressed_modified should be or a point that doesn't belong to the group, or a valid
        //point but different from the original one.
        assert_ne!(a, a_decompressed_modified.unwrap());
    }
}

#[test]
fn test_g2_projective_curve() {
    curve_tests::<G2Projective>();
}

#[test]
fn test_g2_projective_group() {
    let a: G2Projective = rand::random();
    let b: G2Projective = rand::random();
    group_test(a, b);
}

#[test]
fn test_g2_generator() {
    let generator = G2Affine::prime_subgroup_generator();
    assert!(generator.is_on_curve());
    assert!(generator.is_in_correct_subgroup_assuming_on_curve());
}

#[test]
fn test_g2_compression_decompression() {
    let a: G2Projective = rand::random();
    let a = a.into_affine();
    let a_len = to_bytes!(a).unwrap().len();

    //Test correct compression/de-compression of a random point
    let a_compressed = a.compress();
    let a_decompressed = G2Affine::decompress(a_compressed.clone()).unwrap();
    assert_eq!(a, a_decompressed);
    assert_eq!(a_compressed.len(), (a_len/2));

    //Test correct compression/decompression of a zero point
    let b = G2Affine::zero();
    let b_compressed = b.compress();
    let b_decompressed = G2Affine::decompress(b_compressed).unwrap();
    assert_eq!(b, b_decompressed);

    //Test wrong compression/de-compression of a random point by masking with 0s a random byte
    let mut a_compressed_modified = a_compressed.clone();
    loop {
        let index: usize = rand::random();
        let max_idx = (a_len/2) - 4;
        if a_compressed_modified[index % max_idx] != 0u8 {
            a_compressed_modified[index % max_idx] &= 0u8;
            break;
        }
    }
    let a_decompressed_modified = G2Affine::decompress(a_compressed_modified);
    if a_decompressed_modified.is_some(){
        //a_decompressed_modified should be or a point that doesn't belong to the group, or a valid
        //point but different from the original one.
        assert_ne!(a, a_decompressed_modified.unwrap());
    }
}

#[test]
fn test_bilinearity() {
    use crate::fields::{mnt6753::fq6::Fq6, PrimeField};

    let a: G1Projective = rand::random();
    let b: G2Projective = rand::random();
    let s: Fr = rand::random();

    let sa = a * &s;
    let sb = b * &s;

    let ans1 = MNT6::pairing(sa, b);
    let ans2 = MNT6::pairing(a, sb);
    let ans3 = MNT6::pairing(a, b).pow(s.into_repr());

    assert_eq!(ans1, ans2);
    assert_eq!(ans2, ans3);

    assert_ne!(ans1, Fq6::one());
    assert_ne!(ans2, Fq6::one());
    assert_ne!(ans3, Fq6::one());

    assert_eq!(ans1.pow(Fr::characteristic()), Fq6::one());
    assert_eq!(ans2.pow(Fr::characteristic()), Fq6::one());
    assert_eq!(ans3.pow(Fr::characteristic()), Fq6::one());
}

#[test]
fn test_gt_compression(){
    use crate::fields::mnt6753::Fq6;

    let a: G1Projective = rand::random();
    let b: G2Projective = rand::random();
    let c = MNT6::pairing(a, b);

    //Compression/decompression works only if the element of Fq6 is the output of a pairing
    let c_len = to_bytes!(c).unwrap().len();
    let c_compressed = c.compress();
    let c_decompressed = Fq6::decompress(c_compressed.clone()).unwrap();
    assert_eq!(c_decompressed.pow(Fr::characteristic()), Fq6::one());
    assert_eq!(c, c_decompressed);
    assert_eq!(c_compressed.len(), c_len/2);

    //Negative case test
    let mut c_compressed_modified = c_compressed.clone();
    loop {
        let index: usize = rand::random();
        let max_idx = (c_len/2) - 4;
        if c_compressed_modified[index % max_idx] != 0u8 {
            c_compressed_modified[index % max_idx] &= 0u8;
            break;
        }
    }
    let c_decompressed_modified = Fq6::decompress(c_compressed_modified);
    if c_decompressed_modified.is_some() {
        assert_ne!(c, c_decompressed_modified.unwrap());
    }

    //Test that compression/decompression won't work for a Fq6 element which is not a pairing output
    let d: Fq6 = rand::random();
    let d_compressed = d.compress();
    let d_decompressed = Fq6::decompress(d_compressed);
    if d_decompressed.is_some(){
        assert_ne!(d, d_decompressed.unwrap());
    }
}