use crate::{
    curves::{
        models::short_weierstrass_jacobian::GroupAffine,
        sw6::{
            g1::SW6G1Parameters, g2::SW6G2Parameters, G1Affine, G1Projective, G2Affine,
            G2Projective, SW6,
        },
        tests::{curve_tests, sw_curve_serialization_test},
        AffineCurve, PairingEngine,
    },
    groups::tests::group_test,
    CanonicalSerialize,
};
use num_traits::One;

#[test]
fn test_g1_projective_curve() {
    curve_tests::<G1Projective>();

    let byte_size = <GroupAffine<SW6G1Parameters> as CanonicalSerialize>::buffer_size();
    sw_curve_serialization_test::<SW6G1Parameters>(byte_size);
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
fn test_g2_projective_curve() {
    curve_tests::<G2Projective>();

    let byte_size = <GroupAffine<SW6G2Parameters> as CanonicalSerialize>::buffer_size();
    sw_curve_serialization_test::<SW6G2Parameters>(byte_size);
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
fn test_bilinearity() {
    use crate::fields::{
        sw6::{fq6::Fq6, fr::Fr},
        Field, PrimeField,
    };
    use rand;

    let a: G1Projective = rand::random();
    let b: G2Projective = rand::random();
    let s: Fr = rand::random();

    let sa = a * &s;
    let sb = b * &s;

    let ans1 = SW6::pairing(sa, b);
    let ans2 = SW6::pairing(a, sb);
    let ans3 = SW6::pairing(a, b).pow(s.into_repr());

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
#[ignore]
fn print_g1_generator() {
    use crate::fields::sw6::fq::Fq;

    let x: Fq = "5511163824921585887915590525772884263960974614921003940645351443740084257508990841338974915037175497689287870585840954231884082785026301437744745393958283053278991955159266640440849940136976927372133743626748847559939620888818486853646".parse().unwrap();
    let y: Fq = "7913123550914612057135582061699117755797758113868200992327595317370485234417808273674357776714522052694559358668442301647906991623400754234679697332299689255516547752391831738454121261248793568285885897998257357202903170202349380518443".parse().unwrap();

    println!("pub const G1_GENERATOR_X: Fq = Fq::new({});", x.0);
    println!("pub const G1_GENERATOR_Y: Fq = Fq::new({});", y.0);
}

#[test]
#[ignore]
fn print_g2_generator() {
    use crate::fields::sw6::fq::Fq;

    let x_c0: Fq = "13426761183630949215425595811885033211332897733228446437546263564078445562454176776915160094418980045665397361295624472103734543457352048745726512354895954850428989867542989474136256025045975283415690491751906307188562464175510373683338".parse().unwrap();
    let x_c1: Fq = "20471601555918880743198170952645906008198510944268658573129351735028343217532386920456705632337352161031960990613816401042894531220068552819818037605513359562118363589199569321421558696125646867661360498323171027455638052943806292028610".parse().unwrap();
    let x_c2: Fq = "3905053196875761830053608605277158152930144841844497593936739534395003062685449846381431331169369910535935138116320442345524758217411779027270883193856999691582831339845600938304719916501940381093815781408183227875600753651697934495980".parse().unwrap();

    let y_c0: Fq = "8567517639523571619872938228644013584947463594196306323477160496987712111576624702939472765993995586889532559039169098780892505598589581147768095093536988446010255611523736706017580686335404469207486594272103717837888228343074699140243".parse().unwrap();
    let y_c1: Fq = "3890537069205870914984502594450293167889863914413852788876350245583932846980126025043974070704295857226211547108005650399870458089721518559480870503159804530091559886149680718531004778697982910253701559194337987238111062202037698927752".parse().unwrap();
    let y_c2: Fq = "10936269922612615564271188303104593362724754284143779051599749016735041389483971486958818324356025479751246744831831158558101688599198721653921723013062333636402617118847009085485166284126970598561393411916461254016145116183331671450721".parse().unwrap();

    println!("pub const G2_GENERATOR_X_C0: Fq = Fq::new({});", x_c0.0);
    println!("pub const G2_GENERATOR_X_C1: Fq = Fq::new({});", x_c1.0);
    println!("pub const G2_GENERATOR_X_C2: Fq = Fq::new({});", x_c2.0);

    println!("pub const G2_GENERATOR_Y_C0: Fq = Fq::new({});", y_c0.0);
    println!("pub const G2_GENERATOR_Y_C1: Fq = Fq::new({});", y_c1.0);
    println!("pub const G2_GENERATOR_Y_C2: Fq = Fq::new({});", y_c2.0);
}
