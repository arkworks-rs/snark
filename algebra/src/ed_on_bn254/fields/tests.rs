use crate::ed_on_bn254::{Fq, Fr};
use algebra_core::{
    biginteger::BigInteger256 as BigInteger,
    bytes::{FromBytes, ToBytes},
    fields::{Field, LegendreSymbol::*, SquareRootField},
    test_rng, One, Zero,
};

use crate::tests::fields::{field_test, primefield_test};

use core::str::FromStr;
use rand::Rng;

#[test]
fn test_fr() {
    let mut rng = test_rng();
    let a: Fr = rng.gen();
    let b: Fr = rng.gen();
    field_test(a, b);
    primefield_test::<Fr>();
}

#[test]
fn test_fq() {
    let mut rng = test_rng();
    let a: Fq = rng.gen();
    let b: Fq = rng.gen();
    field_test(a, b);
    primefield_test::<Fq>();
}

#[test]
fn test_fq_add() {
    let f1 = Fq::from_str(
        "18386742314266644595564329008376577163854043021652781768352795308532764650733",
    )
    .unwrap();
    let f2 = Fq::from_str(
        "39786307610986038981023499868190793548353538256264351797285876981647142458383",
    )
    .unwrap();
    let f3 = Fq::from_str(
        "14396564181574133132095017386052820535110852477085064878242263917028290117882",
    )
    .unwrap();
    assert!(!f1.is_zero());
    assert!(!f2.is_zero());
    assert!(!f3.is_zero());
    assert_eq!(f1 + &f2, f3);
}

#[test]
fn test_fq_add_one() {
    let f1 = Fq::from_str(
        "4946875394261337176810256604189376311946643975348516311606738923340201185904",
    )
    .unwrap();
    let f2 = Fq::from_str(
        "4946875394261337176810256604189376311946643975348516311606738923340201185905",
    )
    .unwrap();
    assert!(!f1.is_zero());
    assert!(!f2.is_zero());
    assert_eq!(f1 + &Fq::one(), f2);
}

#[test]
fn test_fq_mul() {
    let f1 = Fq::from_str(
        "24703123148064348394273033316595937198355721297494556079070134653139656190956",
    )
    .unwrap();
    let f2 = Fq::from_str(
        "38196797080882758914424853878212529985425118523754343117256179679117054302131",
    )
    .unwrap();
    let f3 = Fq::from_str(
        "1321267396236123309645330145349353750536542060403774171357889269349508194307",
    )
    .unwrap();
    assert!(!f1.is_zero());
    assert!(!f2.is_zero());
    assert!(!f3.is_zero());
    assert_eq!(f1 * &f2, f3);
}

#[test]
fn test_fq_triple_mul() {
    let f1 = Fq::from_str(
        "23834398828139479510988224171342199299644042568628082836691700490363123893905",
    )
    .unwrap();
    let f2 = Fq::from_str(
        "48343809612844640454129919255697536258606705076971130519928764925719046689317",
    )
    .unwrap();
    let f3 = Fq::from_str(
        "22704845471524346880579660022678666462201713488283356385810726260959369106033",
    )
    .unwrap();
    let f4 = Fq::from_str(
        "7747776931431194635550680695131420638163057297019399136408144301550822179875",
    )
    .unwrap();
    assert!(!f1.is_zero());
    assert!(!f2.is_zero());
    assert!(!f3.is_zero());
    assert_eq!(f1 * &f2 * &f3, f4);
}

#[test]
fn test_fq_div() {
    let f1 = Fq::from_str(
        "31892744363926593013886463524057935370302352424137349660481695792871889573091",
    )
    .unwrap();
    let f2 = Fq::from_str(
        "47695868328933459965610498875668250916462767196500056002116961816137113470902",
    )
    .unwrap();
    let f3 = Fq::from_str(
        "7301086967624450577859019086314322648061398679982346993011603220910508457334",
    )
    .unwrap();
    assert!(!f1.is_zero());
    assert!(!f2.is_zero());
    assert!(!f3.is_zero());
    assert_eq!(f1 / &f2, f3);
}

#[test]
fn test_fq_sub() {
    let f1 = Fq::from_str(
        "18695869713129401390241150743745601908470616448391638969502807001833388904079",
    )
    .unwrap();
    let f2 = Fq::from_str(
        "10105476028534616828778879109836101003805485072436929139123765141153277007373",
    )
    .unwrap();
    let f3 = Fq::from_str(
        "8590393684594784561462271633909500904665131375954709830379041860680111896706",
    )
    .unwrap();
    assert!(!f1.is_zero());
    assert!(!f2.is_zero());
    assert!(!f3.is_zero());
    assert_eq!(f1 - &f2, f3);
}

#[test]
fn test_fq_double_in_place() {
    let mut f1 = Fq::from_str(
        "29729289787452206300641229002276778748586801323231253291984198106063944136114",
    )
    .unwrap();
    let f3 = Fq::from_str(
        "15682093831225862156789646514039007320076873845630437896571987838976271280994",
    )
    .unwrap();
    assert!(!f1.is_zero());
    assert!(!f3.is_zero());
    f1.double_in_place();
    assert_eq!(f1, f3);
}

#[test]
fn test_fq_double_in_place_thrice() {
    let mut f1 = Fq::from_str(
        "32768907806651393940832831055386272949401004221411141755415956893066040832473",
    )
    .unwrap();
    let f3 = Fq::from_str(
        "21380590862979124081952185245260157621176025366712756262647409092194433207997",
    )
    .unwrap();
    assert!(!f1.is_zero());
    assert!(!f3.is_zero());
    f1.double_in_place();
    f1.double_in_place();
    f1.double_in_place();
    assert_eq!(f1, f3);
}

#[test]
fn test_fq_generate_random_ed_on_bn254_point() {
    let a = Fq::from_str("168700").unwrap();

    let d = Fq::from_str("168696").unwrap();
    let y = Fq::from_str(
        "19987327827845206670850937090314462639017692512983955920885166014935289314257",
    )
    .unwrap();
    let x2 = Fq::from_str(
        "2144239075372598103060889495211040948751593385312551803225522963913923559328",
    )
    .unwrap();

    let computed_y2 = y.square();
    let y2 = Fq::from_str(
        "11134206686211572308995578277928848431421308813024790181507137950838333998633",
    )
    .unwrap();
    assert_eq!(y2, computed_y2);

    let computed_dy2 = d * &computed_y2;
    let dy2 =
        Fq::from_str("345576003677591687256955722467813448317229128849323754147891993737799010947")
            .unwrap();
    assert_eq!(dy2, computed_dy2);

    let computed_divisor = computed_dy2 - a;
    let divisor =
        Fq::from_str("345576003677591687256955722467813448317229128849323754147891993737798842247")
            .unwrap();
    assert_eq!(divisor, computed_divisor);

    let computed_x2 = (computed_y2 - &Fq::one()) / &computed_divisor;
    assert_eq!(x2, computed_x2);

    let x = Fq::from_str(
        "4801447892755635304907919953550459075619191823587157449340656925102682829025",
    )
    .unwrap();
    let computed_x = computed_x2.sqrt().unwrap();
    assert_eq!(computed_x.square(), x2);
    assert_eq!(x, computed_x);

    fn add<'a>(curr: (Fq, Fq), other: &'a (Fq, Fq)) -> (Fq, Fq) {
        let y1y2 = curr.1 * &other.1;
        let x1x2 = curr.0 * &other.0;
        let a = Fq::from_str("168700").unwrap();
        let d = Fq::from_str("168696").unwrap();
        let dx1x2y1y2 = d * &y1y2 * &x1x2;

        let d1 = Fq::one() + &dx1x2y1y2;
        let d2 = Fq::one() - &dx1x2y1y2;

        let x1y2 = curr.0 * &other.1;
        let y1x2 = curr.1 * &other.0;

        let x = (x1y2 + &y1x2) / &d1;
        let y = (y1y2 - a * &x1x2) / &d2;

        (x, y)
    }

    let result = add((x, y), &(x, y));
    let result = add(result, &result);
    let result = add(result, &result);

    let point_x =
        Fq::from_str("380676173762867192861894055350059333852732198308367125138259398265363727587")
            .unwrap();
    let point_y = Fq::from_str(
        "8435074244857818446059206728316702149733931432112984450960434710303841866985",
    )
    .unwrap();
    assert_eq!((point_x, point_y), result);
}

#[test]
fn test_fq_square_in_place() {
    let mut f1 = Fq::from_str(
        "6060110850233386730847324622937480088943976359504617699731744947670229990461",
    )
    .unwrap();
    let f3 = Fq::from_str(
        "17018926051730832095053393285350575966874590491719897015583930476179087429684",
    )
    .unwrap();
    assert!(!f1.is_zero());
    assert!(!f3.is_zero());
    f1.square_in_place();
    assert_eq!(f1, f3);
}

#[test]
fn test_fq_sqrt() {
    let f1 = Fq::from_str(
        "5830207146824777307592559303161432403393380070279905260050870500920682305217",
    )
    .unwrap();
    let f3 = Fq::from_str(
        "2108183130040740552565127577293974960058698876185401671087892009247563211475",
    )
    .unwrap();
    assert_eq!(f1.sqrt().unwrap(), f3);
}

#[test]
fn test_fq_from_str() {
    let f1_from_repr = Fq::from(BigInteger([
        0xab8a2535947d1a77,
        0x9ba74cbfda0bbcda,
        0xe928b59724d60baf,
        0x1cccaaeb9bb1680a,
    ]));
    let f1 = Fq::from_str(
        "13026376210409056429264774981357153555336288129100724591327877625017068755575",
    )
    .unwrap();
    let f2_from_repr = Fq::from(BigInteger([
        0x97e9103775d2f35c,
        0xbe6756b6c587544b,
        0x6ee38c3afd88ef4b,
        0x2bacd150f540c677,
    ]));
    let f2 = Fq::from_str(
        "19754794831832707859764530223239420866832328728734160755396495950822165902172",
    )
    .unwrap();
    assert_eq!(f1_from_repr, f1);
    assert_eq!(f2_from_repr, f2);
}

#[test]
fn test_fq_legendre() {
    assert_eq!(QuadraticResidue, Fq::one().legendre());
    assert_eq!(Zero, Fq::zero().legendre());

    let e = BigInteger([
        0x2e8de1a676c03be8,
        0x73350d34fe25a560,
        0x7ea085919029688e,
        0x1d0868cb993cf28,
    ]);
    assert_eq!(QuadraticResidue, Fq::from(e).legendre());
    let e = BigInteger([
        0x891d8cc23c8d0706,
        0xe91800e007db2698,
        0xfff380321e9ac7a7,
        0x2659e28bd17eab6,
    ]);
    assert_eq!(QuadraticNonResidue, Fq::from(e).legendre());
}

#[test]
fn test_fq_bytes() {
    let f1_from_repr = Fq::from(BigInteger([
        0xab8a2535947d1a77,
        0x9ba74cbfda0bbcda,
        0xe928b59724d60baf,
        0x1cccaaeb9bb1680a,
    ]));

    let mut f1_bytes = [0u8; 32];
    f1_from_repr.write(f1_bytes.as_mut()).unwrap();

    let f1 = Fq::read(f1_bytes.as_ref()).unwrap();
    assert_eq!(f1_from_repr, f1);
}

#[test]
fn test_fr_add() {
    let f1 = Fr::from(BigInteger([
        0xccfc9a195e0f5c46,
        0xaed4874d13fb1285,
        0x27368f86ca2848eb,
        0x4f8adcfeb44fccc,
    ]));
    let f2 = Fr::from(BigInteger([
        0x661ff05bf8570851,
        0x1b171f4c59be97ef,
        0x5d2ce7f9b4d701f3,
        0x1e0e794623e0f68,
    ]));
    let f3 = Fr::from(BigInteger([
        0xcba9f2991d453da6,
        0x1eacb8e13498bc6a,
        0x4d596ec9aecf1fd3,
        0xcd0b95f15cd82f,
    ]));
    assert_eq!(f1 + &f2, f3);
}

#[test]
fn test_fr_mul() {
    let f1 = Fr::from(BigInteger([
        0xc2964d2dd5fb980f,
        0xbab64d599c57e496,
        0x39cae13e7d1d4f78,
        0x1aa995aa4de205c,
    ]));
    let f2 = Fr::from(BigInteger([
        0xc256e720cd43533b,
        0x3bfbadf6247e13bb,
        0x94c3d63a53714f63,
        0x10f8a7bf74efd57,
    ]));
    let f3 = Fr::from(BigInteger([
        0x5eac88be41e0e1fd,
        0x57aab36675b11e24,
        0x835582d896b4d13f,
        0x4808736e213036e,
    ]));
    assert_eq!(f1 * &f2, f3);
}
#[test]
fn test_fr_bytes() {
    let f1_from_repr = Fr::from(BigInteger([
        0xc81265fb4130fe0c,
        0xb308836c14e22279,
        0x699e887f96bff372,
        0x84ecc7e76c11ad,
    ]));

    let mut f1_bytes = [0u8; 32];
    f1_from_repr.write(f1_bytes.as_mut()).unwrap();

    let f1 = Fr::read(f1_bytes.as_ref()).unwrap();
    assert_eq!(f1_from_repr, f1);
}

#[test]
fn test_fr_from_str() {
    let f100_from_repr = Fr::from(BigInteger([0x64, 0, 0, 0]));
    let f100 = Fr::from_str("100").unwrap();
    assert_eq!(f100_from_repr, f100);
}
