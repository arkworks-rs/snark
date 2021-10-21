use crate::{
    biginteger::BigInteger768 as BigInteger,
    field_new,
    fields::{
        fp3::{Fp3, Fp3Parameters},
        mnt6753::fq::Fq,
    },
};

pub type Fq3 = Fp3<Fq3Parameters>;

pub struct Fq3Parameters;

impl Fp3Parameters for Fq3Parameters {
    type Fp = Fq;

    // alpha = 11
    const NONRESIDUE: Fq = field_new!(
        Fq,
        BigInteger([
            0x4768931cfff9c7d4,
            0xc45e46d6ada96ca0,
            0x479b0bdb0b3c0107,
            0x362a089610f8d41b,
            0xdbafcec2c8a91aaf,
            0x78428b0ff9d96a06,
            0xf2e4472a9080c353,
            0xc9006ed33f0e971c,
            0x0794d9d10bdb7288,
            0x3c1e44cab5419e2c,
            0x49b5fc6c81f4560c,
            0x1c287777c30ba,
        ])
    );

    const TWO_ADICITY: u32 = 30;

    //t=(p^3-1)/2
    const T_MINUS_ONE_DIV_TWO: &'static [u64] = &[
        0xd6447f9d762cc94d,
        0xfc72f2d69c49b1dd,
        0x56524f8eca1d3e92,
        0x8f1633f602c3b2ae,
        0x45d5bebb37be973c,
        0x36b885fe0423c666,
        0x1b5aefa50853c03d,
        0x549ba23c3c70fa49,
        0xb323e0add7f13ec2,
        0x39c6bf6b757e6ec2,
        0x9017af105004645a,
        0x7d05c9b5544267a3,
        0xff83ee77adbe22f9,
        0xabe49e95ab5133f0,
        0xb98c227558b1b9e1,
        0xa54641bd1a4e20c8,
        0x52c5a4bad703a538,
        0xd4fd4c0c949ac98b,
        0x61c6203eb008385d,
        0xc65ed5664f9b95a9,
        0x55c4ecdf6ca7c4f5,
        0xc795504c013a1fb3,
        0xfc04ff3e3afea252,
        0xf2ae66577c689a10,
        0xaae48029a805f455,
        0xa827c78687948639,
        0x2f3433f22bf74542,
        0xc0f9bb9fe47134a2,
        0x98460e01b1baceca,
        0x54b654cc62afaea5,
        0x4116b8ae7f04bd20,
        0x43bcac41205e99c6,
        0x1abcd4f53d1d225e,
        0xbbcd53c3b60dd859,
        0xb10a9b0dc2128,
    ];

    // quadratic non-residue (c0+ 0* X + 0*X^2),
    // c0=1659781419976566415021745064391095587555604711122877233394175714744952\
    // 8858260736949611040045144870330052211081080633146517319686355350402298\
    // 7667807845209238605163122279088377413675555771794427110778974744994732\
    // 54624198506809678
    const QUADRATIC_NONRESIDUE_TO_T: (Fq, Fq, Fq) = (
        field_new!(
            Fq,
            BigInteger([
                0x2217cbfb0feb469c,
                0x68216255ea00e214,
                0xe1391d4fa199ab8,
                0x915ac7dfd6cbc927,
                0xdc90acba889c8eff,
                0x74d377e9be5fc824,
                0x7d23df9a20eabf7a,
                0x2891082620e9a3e6,
                0x820481f8ecaea6f8,
                0xf0b43af4e2ce8c2e,
                0x97cc7da5fef0c28a,
                0x6157c1dabadf,
            ])
        ),
        field_new!(Fq, BigInteger([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
        field_new!(Fq, BigInteger([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
    );

    const FROBENIUS_COEFF_FP3_C1: &'static [Fq] = &[
        //X^{q^0} = alpha^((q^0 - 1)/ 3)*X = 1*X
        field_new!(
            Fq,
            BigInteger([
                0xb99680147fff6f42,
                0x4eb16817b589cea8,
                0xa1ebd2d90c79e179,
                0xf725caec549c0da,
                0xab0c4ee6d3e6dad4,
                0x9fbca908de0ccb62,
                0x320c3bb713338498,
                0x598b4302d2f00a62,
                0x4074c9cbfd8ca621,
                0xfa47edb3865e88c,
                0x95455fb31ff9a195,
                0x7b479ec8e242,
            ])
        ),
        //X^{q^1} = alpha^((q^1 - 1)/ 3)*X
        field_new!(
            Fq,
            BigInteger([
                0x6b66f7b83f968680,
                0x1379b1ebf803e51e,
                0x9bb6f43b5282969c,
                0x3f64a98166c46a97,
                0x524a1cc56c78e977,
                0xf480725d1dc6e2f1,
                0xe660b05c89764d7d,
                0xe5b38512c92d9f5b,
                0xa75658e33e25f9f0,
                0xb4b96c948f0e9992,
                0xb8b087523d7db902,
                0x11d5033223a5d,
            ])
        ),
        //X^{q^2} = alpha^((q^2 - 1)/ 3)*X
        field_new!(
            Fq,
            BigInteger([
                0xb409ff15806a0a3f,
                0xec757f1362138688,
                0x9920baa7e003df81,
                0x6b08f346088b0f32,
                0x41955ee7e8c161eb,
                0x1e2f40c2cc85fb47,
                0x816438c604b12587,
                0xc8bef1104c8343cf,
                0x76ecc64a30e53860,
                0xf39babe0941b2dce,
                0xc22ca91d916b7315,
                0x2c2e5ba7a770,
            ])
        ),
    ];

    const FROBENIUS_COEFF_FP3_C2: &'static [Fq] = &[
        //(X^2)^{q^0} = alpha^(2(q^0 - 1)/ 3)*X^2 = 1*X^2
        field_new!(
            Fq,
            BigInteger([
                0xb99680147fff6f42,
                0x4eb16817b589cea8,
                0xa1ebd2d90c79e179,
                0xf725caec549c0da,
                0xab0c4ee6d3e6dad4,
                0x9fbca908de0ccb62,
                0x320c3bb713338498,
                0x598b4302d2f00a62,
                0x4074c9cbfd8ca621,
                0xfa47edb3865e88c,
                0x95455fb31ff9a195,
                0x7b479ec8e242,
            ])
        ),
        //(X^2)^{q^1} = alpha^(2(q^1 - 1)/ 3)*X^2
        field_new!(
            Fq,
            BigInteger([
                0xb409ff15806a0a3f,
                0xec757f1362138688,
                0x9920baa7e003df81,
                0x6b08f346088b0f32,
                0x41955ee7e8c161eb,
                0x1e2f40c2cc85fb47,
                0x816438c604b12587,
                0xc8bef1104c8343cf,
                0x76ecc64a30e53860,
                0xf39babe0941b2dce,
                0xc22ca91d916b7315,
                0x2c2e5ba7a770,
            ])
        ),
        //(X^2)^{q^2} = alpha^(2(q^2 - 1)/ 3)*X^2
        field_new!(
            Fq,
            BigInteger([
                0x6b66f7b83f968680,
                0x1379b1ebf803e51e,
                0x9bb6f43b5282969c,
                0x3f64a98166c46a97,
                0x524a1cc56c78e977,
                0xf480725d1dc6e2f1,
                0xe660b05c89764d7d,
                0xe5b38512c92d9f5b,
                0xa75658e33e25f9f0,
                0xb4b96c948f0e9992,
                0xb8b087523d7db902,
                0x11d5033223a5d,
            ])
        ),
    ];
}
