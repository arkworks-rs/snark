use crate::{
    biginteger::BigInteger256 as BigInteger,
    bytes::{FromBytes, ToBytes},
    fields::{
        tests::{field_test, primefield_test},
        tweedle::{fq::Fq, fr::Fr},
        Field,
        LegendreSymbol::*,
        PrimeField, SquareRootField,
    },
};
use std::str::FromStr;

#[test]
fn test_tweedle_fr() {
    let a: Fr = rand::random();
    let b: Fr = rand::random();
    field_test(a, b);
    primefield_test::<Fr>();
}

#[test]
fn test_tweedle_fq() {
    let a: Fq = rand::random();
    let b: Fq = rand::random();
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
        "19786307610986038981023499868190793548353538256264351797285876981647142458383",
    )
    .unwrap();
    let f3 = Fq::from_str(
        "9225027615923634720695082624395393748885377621962700438691588327011328770299",
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
        "2946875394261337176810256604189376311946643975348516311606738923340201185904",
    )
    .unwrap();
    let f2 = Fq::from_str(
        "2946875394261337176810256604189376311946643975348516311606738923340201185905",
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
        "18196797080882758914424853878212529985425118523754343117256179679117054302131",
    )
    .unwrap();
    let f3 = Fq::from_str(
        "17200426197156579651574285400495776861464966086836221879823411702083205504671",
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
        "28343809612844640454129919255697536258606705076971130519928764925719046689317",
    )
    .unwrap();
    let f3 = Fq::from_str(
        "22704845471524346880579660022678666462201713488283356385810726260959369106033",
    )
    .unwrap();
    let f4 = Fq::from_str(
        "16214848981638217704332828370994953257254063208571075715629895139281735615384",
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
        "28343809612844640454129919255697536258606705076971130519928764925719046689317",
    )
    .unwrap();
    let f2 = Fq::from_str(
        "10565087619621820656254220984570060779538680183776893272687321117751927412212",
    )
    .unwrap();
    let f3 = Fq::from_str(
        "18524871784212090982875844178168704466262002646942949074484608991277082160887",
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
        "27729289787452206300641229002276778748586801323231253291984198106063944136114",
    )
    .unwrap();
    let f3 = Fq::from_str(
        "26510557265575363745389711752381580533851398990508073457021312248959309933411",
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
        "12768907806651393940832831055386272949401004221411141755415956893066040832473",
    )
    .unwrap();
    let f3 = Fq::from_str(
        "15307195525224004958984409686574252705241422803425834662486403255022591643333",
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
fn test_fq_square_in_place() {
    let mut f1 = Fq::from_str(
        "21864651240005695523200639428464570946052769938774601449735727714436878540682",
    )
    .unwrap();
    let f3 = Fq::from_str(
        "19581484219153942072047858331099766329275818131597627721493461404774907481642",
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
        "10875927553327821418567659853801220899541454800710193788767706167237535308236",
    )
    .unwrap();
    let f3 = Fq::from_str(
        "14305375816171855067009304014627046917501979889993032068503214923182596004917",
    )
    .unwrap();
    assert_eq!(f1.sqrt().unwrap(), f3);
}

#[test]
fn test_fq_from_str() {
    let f1_from_repr = Fq::from_repr(BigInteger([
        0x17a655e2b3cd9f8a,
        0xe98745acbc60cf8,
        0x8318964c0b265e48,
        0x3056f43c9ee293fd,
    ]));
    let f1 = Fq::from_str(
        "21864651240005695523200639428464570946052769938774601449735727714436878540682",
    )
    .unwrap();
    let f2_from_repr = Fq::from_repr(BigInteger([
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

    let e = Fq::from_str(
        "19754794831832707859764530223239420866832328728734160755396495950822165902172",
    )
    .unwrap();
    assert_eq!(QuadraticResidue, e.legendre());
    let e = Fq::from_str(
        "7894070009960485405056471228743059385328854667547937089962899125529157892247",
    )
    .unwrap();
    assert_eq!(QuadraticNonResidue, e.legendre());
}

#[test]
fn test_fq_bytes() {
    let f1_from_repr = Fq::from_repr(BigInteger([
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
    let f1 = Fr::from_str(
        "28343809612844640454129919255697536258606705076971130519928764925719046689317",
    )
    .unwrap();
    let f2 = Fr::from_str(
        "10565087619621820656254220984570060779538680183776893272687321117751927412212",
    )
    .unwrap();
    let f3 = Fr::from_str(
        "9960874923137412254491393988095620074823181604792704735842768974107331995672",
    )
    .unwrap();
    assert!(!f1.is_zero());
    assert!(!f2.is_zero());
    assert!(!f3.is_zero());
    assert_eq!(f1 + &f2, f3);
}

#[test]
fn test_fr_mul() {
    let f1 = Fr::from_str(
        "24703123148064348394273033316595937198355721297494556079070134653139656190956",
    )
    .unwrap();
    let f2 = Fr::from_str(
        "18196797080882758914424853878212529985425118523754343117256179679117054302131",
    )
    .unwrap();
    let f3 = Fr::from_str(
        "3313453718804050382618227313767372842847071020036631866735220624776250855218",
    )
    .unwrap();
    assert!(!f1.is_zero());
    assert!(!f2.is_zero());
    assert!(!f3.is_zero());
    assert_eq!(f1 * &f2, f3);
}

#[test]
fn test_fr_bytes() {
    let f1_from_repr = Fr::from_repr(BigInteger([
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
    let f100_from_repr = Fr::from_repr(BigInteger([0x64, 0, 0, 0]));
    let f100 = Fr::from_str("100").unwrap();
    assert_eq!(f100_from_repr, f100);
}

#[test]
#[ignore]
fn print_field() {
    println!("one: {:?}", Fq::one());
    println!("zero: {:?}", Fq::zero());
    println!(
        "256 in repr: {:?}",
        Fq::from_repr(BigInteger([0, 0, 1, 255]))
    );
    println!("256: {:?}", Fq::from_str("256").unwrap().into_repr());
    println!("1024: {:?}", Fq::from_str("1024").unwrap().into_repr());
    println!(
        "255 to bytes: {:?}",
        to_bytes![Fq::from_str("255").unwrap().into_repr()].unwrap()
    );
    println!(
        "256 to bytes: {:?}",
        to_bytes![Fq::from_str("256").unwrap().into_repr()].unwrap()
    );
    println!(
        "1023 to bytes: {:?}",
        to_bytes![Fq::from_str("1023").unwrap().into_repr()].unwrap()
    );
}
