macro_rules! f_bench {
    // Use this for base fields
    ($f:ident, $f_type:ty, $f_repr:ident, $f_repr_type:ty, $field_ident:ident) => {
        field_common!($f, $f_type, $field_ident);
        sqrt!($f, $f_type, $field_ident);
        field_base!($f, $f_type, $f_repr, $f_repr_type, $field_ident);
    };
    // use this for intermediate fields
    (1, $f:ident, $f_type:ty, $field_ident:ident) => {
        field_common!($f, $f_type, $field_ident);
        sqrt!($f, $f_type, $field_ident);
    };
    // Use this for the full extension field Fqk
    (2, $f:ident, $f_type:ty, $field_ident:ident) => {
        field_common!($f, $f_type, $field_ident);
    };
}

macro_rules! field_common {
    ($f:ident, $f_type:ty, $field_ident:ident) => {
        paste::item! {
            #[bench]
            fn [<bench_ $field_ident _add_assign>](b: &mut ::test::Bencher) {
                const SAMPLES: usize = 1000;

                let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

                let v: Vec<_> = (0..SAMPLES)
                    .map(|_| ($f::rand(&mut rng), $f::rand(&mut rng)))
                    .collect();

                let mut count = 0;
                b.iter(|| {
                    let mut tmp = v[count].0;
                    n_fold!(tmp, v, add_assign, count);
                    count = (count + 1) % SAMPLES;
                    tmp
                });
            }

            #[bench]
            fn [<bench_ $field_ident _sub_assign>](b: &mut ::test::Bencher) {
                const SAMPLES: usize = 1000;

                let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

                let v: Vec<_> = (0..SAMPLES)
                    .map(|_| ($f::rand(&mut rng), $f::rand(&mut rng)))
                    .collect();

                let mut count = 0;
                b.iter(|| {
                    let mut tmp = v[count].0;
                    n_fold!(tmp, v, sub_assign, count);
                    count = (count + 1) % SAMPLES;
                    tmp
                });
            }

            #[bench]
            fn [<bench_ $field_ident _mul_assign>](b: &mut ::test::Bencher) {
                const SAMPLES: usize = 1000;

                let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

                let v: Vec<_> = (0..SAMPLES)
                    .map(|_| ($f::rand(&mut rng), $f::rand(&mut rng)))
                    .collect();

                let mut count = 0;
                b.iter(|| {
                    let mut tmp = v[count].0;
                    n_fold!(tmp, v, mul_assign, count);
                    count = (count + 1) % SAMPLES;
                    tmp
                });
            }

            #[bench]
            fn [<bench_ $field_ident _double>](b: &mut ::test::Bencher) {
                const SAMPLES: usize = 1000;

                let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

                let v: Vec<$f_type> = (0..SAMPLES).map(|_| $f::rand(&mut rng)).collect();

                let mut count = 0;
                b.iter(|| {
                    let mut tmp = v[count];
                    n_fold!(tmp, double_in_place);
                    count = (count + 1) % SAMPLES;
                    tmp
                });
            }

            #[bench]
            fn [<bench_ $field_ident _square>](b: &mut ::test::Bencher) {
                const SAMPLES: usize = 1000;

                let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

                let v: Vec<$f_type> = (0..SAMPLES).map(|_| $f::rand(&mut rng)).collect();

                let mut count = 0;
                b.iter(|| {
                    let mut tmp = v[count];
                    n_fold!(tmp, square_in_place);
                    count = (count + 1) % SAMPLES;
                    tmp
                });
            }

            #[bench]
            fn [<bench_ $field_ident _inverse>](b: &mut ::test::Bencher) {
                const SAMPLES: usize = 1000;

                let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

                let v: Vec<$f_type> = (0..SAMPLES).map(|_| $f::rand(&mut rng)).collect();

                let mut count = 0;
                b.iter(|| {
                    let tmp = v[count].inverse();
                    count = (count + 1) % SAMPLES;
                    tmp
                });
            }

            #[bench]
            fn [<bench_ $field_ident _deser>](b: &mut ::test::Bencher) {
                use algebra::{CanonicalSerialize, CanonicalDeserialize};
                const SAMPLES: usize = 1000;

                let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

                let mut num_bytes = 0;
                let v: Vec<_> = (0..SAMPLES).flat_map(|_| {
                    let mut bytes = Vec::with_capacity(1000);
                    let tmp = $f::rand(&mut rng);
                    tmp.serialize(&mut bytes).unwrap();
                    num_bytes = bytes.len();
                    bytes
                }).collect();

                let mut count = 0;
                b.iter(|| {
                    count = (count + 1) % SAMPLES;
                    let index = count * num_bytes;
                    $f_type::deserialize(&v[index..(index + num_bytes)]).unwrap()
                });
            }

            #[bench]
            fn [<bench_ $field_ident _ser>](b: &mut ::test::Bencher) {
                use algebra::CanonicalSerialize;
                const SAMPLES: usize = 1000;

                let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

                let v: Vec<$f_type> = (0..SAMPLES).map(|_| $f::rand(&mut rng)).collect();
                let mut bytes = Vec::with_capacity(1000);

                let mut count = 0;
                b.iter(|| {
                    let tmp = v[count];
                    count = (count + 1) % SAMPLES;
                    bytes.clear();
                    tmp.serialize(&mut bytes)

                });
            }

            #[bench]
            fn [<bench_ $field_ident _deser_unchecked>](b: &mut ::test::Bencher) {
                use algebra::{CanonicalSerialize, CanonicalDeserialize};
                const SAMPLES: usize = 1000;

                let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

                let mut num_bytes = 0;
                let v: Vec<_> = (0..SAMPLES).flat_map(|_| {
                    let mut bytes = Vec::with_capacity(1000);
                    let tmp = $f::rand(&mut rng);
                    tmp.serialize_unchecked(&mut bytes).unwrap();
                    num_bytes = bytes.len();
                    bytes
                }).collect();

                let mut count = 0;
                b.iter(|| {
                    count = (count + 1) % SAMPLES;
                    let index = count * num_bytes;
                    $f_type::deserialize_unchecked(&v[index..(index + num_bytes)]).unwrap()
                });
            }

            #[bench]
            fn [<bench_ $field_ident _ser_unchecked>](b: &mut ::test::Bencher) {
                use algebra::CanonicalSerialize;
                const SAMPLES: usize = 1000;

                let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

                let v: Vec<$f_type> = (0..SAMPLES).map(|_| $f::rand(&mut rng)).collect();
                let mut bytes = Vec::with_capacity(1000);

                let mut count = 0;
                b.iter(|| {
                    let tmp = v[count];
                    count = (count + 1) % SAMPLES;
                    bytes.clear();
                    tmp.serialize_unchecked(&mut bytes)

                });
            }
        }
    };
}

macro_rules! sqrt {
    ($f:ident, $f_type:ty, $field_ident:ident) => {
        paste::item! {
            #[bench]
            fn [<bench_ $field_ident _sqrt>](b: &mut ::test::Bencher) {
                const SAMPLES: usize = 1000;

                let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

                let v: Vec<$f_type> = (0..SAMPLES)
                    .map(|_| {
                        let mut tmp = $f::rand(&mut rng);
                        tmp.square_in_place();
                        tmp
                    })
                    .collect();

                let mut count = 0;
                b.iter(|| {
                    count = (count + 1) % SAMPLES;
                    v[count].sqrt()
                });
            }
        }
    };
}

macro_rules! field_base {
    ($f:ident, $f_type:ty, $f_repr:ident, $f_repr_type:ty, $field_ident:ident) => {
        paste::item! {
            #[bench]
            fn [<bench_ $field_ident _repr_add_nocarry>](b: &mut ::test::Bencher) {
                const SAMPLES: usize = 1000;

                let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

                let v: Vec<_> = (0..SAMPLES)
                    .map(|_| {
                        let mut tmp1 = $f_repr::rand(&mut rng);
                        let mut tmp2 = $f_repr::rand(&mut rng);
                        // Shave a few bits off to avoid overflow.
                        for _ in 0..3 {
                            tmp1.div2();
                            tmp2.div2();
                        }
                        (tmp1, tmp2)
                    })
                    .collect();

                let mut count = 0;
                b.iter(|| {
                    let mut tmp = v[count].0;
                    n_fold!(tmp, v, add_nocarry, count);
                    count = (count + 1) % SAMPLES;
                    tmp
                });
            }

            #[bench]
            fn [<bench_ $field_ident _repr_sub_noborrow>](b: &mut ::test::Bencher) {
                const SAMPLES: usize = 1000;

                let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

                let v: Vec<_> = (0..SAMPLES)
                    .map(|_| {
                        let tmp1 = $f_repr::rand(&mut rng);
                        let mut tmp2 = tmp1;
                        // Ensure tmp2 is smaller than tmp1.
                        for _ in 0..10 {
                            tmp2.div2();
                        }
                        (tmp1, tmp2)
                    })
                    .collect();

                let mut count = 0;
                b.iter(|| {
                    let mut tmp = v[count].0;
                    n_fold!(tmp, v, sub_noborrow, count);
                    count = (count + 1) % SAMPLES;
                    tmp
                });
            }

            #[bench]
            fn [<bench_ $field_ident _repr_num_bits>](b: &mut ::test::Bencher) {
                const SAMPLES: usize = 1000;

                let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

                let v: Vec<$f_repr_type> = (0..SAMPLES).map(|_| $f_repr::rand(&mut rng)).collect();

                let mut count = 0;
                b.iter(|| {
                    let tmp = v[count].num_bits();
                    count = (count + 1) % SAMPLES;
                    tmp
                });
            }

            #[bench]
            fn [<bench_ $field_ident _repr_mul2>](b: &mut ::test::Bencher) {
                const SAMPLES: usize = 1000;

                let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

                let v: Vec<$f_repr_type> = (0..SAMPLES).map(|_| $f_repr::rand(&mut rng)).collect();

                let mut count = 0;
                b.iter(|| {
                    let mut tmp = v[count];
                    n_fold!(tmp, mul2);
                    count = (count + 1) % SAMPLES;
                    tmp
                });
            }

            #[bench]
            fn [<bench_ $field_ident _repr_div2>](b: &mut ::test::Bencher) {
                const SAMPLES: usize = 1000;

                let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

                let v: Vec<$f_repr_type> = (0..SAMPLES).map(|_| $f_repr::rand(&mut rng)).collect();

                let mut count = 0;
                b.iter(|| {
                    let mut tmp = v[count];
                    n_fold!(tmp, div2);
                    count = (count + 1) % SAMPLES;
                    tmp
                });
            }

            #[bench]
            fn [<bench_ $field_ident _negate>](b: &mut ::test::Bencher) {
                const SAMPLES: usize = 1000;

                let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

                let v: Vec<$f_type> = (0..SAMPLES).map(|_| $f::rand(&mut rng)).collect();

                let mut count = 0;
                b.iter(|| {
                    let mut tmp = v[count];
                    tmp = -tmp;
                    count = (count + 1) % SAMPLES;
                    tmp
                });
            }

            #[bench]
            fn [<bench_ $field_ident _into_repr>](b: &mut ::test::Bencher) {
                const SAMPLES: usize = 1000;

                let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

                let v: Vec<$f_type> = (0..SAMPLES).map(|_| $f::rand(&mut rng)).collect();

                let mut count = 0;
                b.iter(|| {
                    count = (count + 1) % SAMPLES;
                    v[count].into_repr()
                });
            }

            #[bench]
            fn [<bench_ $field_ident _from_repr>](b: &mut ::test::Bencher) {
                const SAMPLES: usize = 1000;

                let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

                let v: Vec<$f_repr_type> = (0..SAMPLES)
                    .map(|_| $f::rand(&mut rng).into_repr())
                    .collect();

                let mut count = 0;
                b.iter(|| {
                    count = (count + 1) % SAMPLES;
                    $f::from(v[count])
                });
            }
        }
    };
}
