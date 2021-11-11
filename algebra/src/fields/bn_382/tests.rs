use crate::{
    biginteger::{BigInteger, BigInteger384},
    fields::{
        bn_382::{
            Fq, Fq12, Fq12Parameters, Fq2, Fq2Parameters, Fq6, Fq6Parameters, FqParameters, Fr,
        },
        fp12_2over3over2::Fp12Parameters,
        fp6_3over2::Fp6Parameters,
        tests::{field_test, frobenius_test, primefield_test, sqrt_field_test},
        Field, Fp2Parameters, FpParameters, PrimeField, SquareRootField,
    },
    UniformRand,
};
use rand::SeedableRng;
use rand_xorshift::XorShiftRng;
use std::{
    cmp::Ordering,
    ops::{AddAssign, MulAssign},
};

pub(crate) const ITERATIONS: usize = 5;

#[test]
fn test_size_fr() {
    println!("{}", <Fr as PrimeField>::Params::MODULUS_BITS);
    println!("{}", std::mem::size_of::<<Fr as PrimeField>::BigInt>() * 8)
}

#[test]
fn test_bn_382_fr() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    for _ in 0..ITERATIONS {
        let a: Fr = UniformRand::rand(&mut rng);
        let b: Fr = UniformRand::rand(&mut rng);
        field_test(a, b);
        primefield_test::<Fr>();
        sqrt_field_test(b);
    }
}

#[test]
fn test_bn_382_fq() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    for _ in 0..ITERATIONS {
        let a: Fq = UniformRand::rand(&mut rng);
        let b: Fq = UniformRand::rand(&mut rng);
        field_test(a, b);
        primefield_test::<Fq>();
        sqrt_field_test(a);
    }
}

#[test]
fn test_bn_382_fq2() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    for _ in 0..ITERATIONS {
        let a: Fq2 = UniformRand::rand(&mut rng);
        let b: Fq2 = UniformRand::rand(&mut rng);
        field_test(a, b);
        sqrt_field_test(a);
    }
    frobenius_test::<Fq2, _>(Fq::characteristic(), 13);
}

#[test]
fn test_bn_382_fq6() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    for _ in 0..ITERATIONS {
        let g: Fq6 = UniformRand::rand(&mut rng);
        let h: Fq6 = UniformRand::rand(&mut rng);
        field_test(g, h);
    }
    frobenius_test::<Fq6, _>(Fq::characteristic(), 13);
}

#[test]
fn test_bn_382_fq12() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    for _ in 0..ITERATIONS {
        let g: Fq12 = UniformRand::rand(&mut rng);
        let h: Fq12 = UniformRand::rand(&mut rng);
        field_test(g, h);
    }
    frobenius_test::<Fq12, _>(Fq::characteristic(), 13);
}

#[test]
fn test_bn_382_negative_one() {
    let neg_one = Fq::new(BigInteger384([
        0x8,
        0xc0060c0c0,
        0xc1848c18180c00,
        0xa451b0a144d8480c,
        0x8a81e34d84edfc45,
        0x202449fed6c43c73,
    ]));
    assert_eq!(neg_one, -Fq::one());
}

#[test]
fn test_frob_coeffs() {
    let nqr = Fq::from_repr(BigInteger384::from(7));
    let cq = Fq::characteristic();
    let q = BigInteger384([cq[0], cq[1], cq[2], cq[3], cq[4], cq[5]]);

    let q_minus_1_over_2 = {
        let mut x = q.clone();
        let _ = x.sub_noborrow(&1.into());
        x.div2();
        x
    };
    assert_eq!(Fq2Parameters::FROBENIUS_COEFF_FP2_C1[0], Fq::one());
    assert_eq!(
        Fq2Parameters::FROBENIUS_COEFF_FP2_C1[1],
        // (q - 1) / 2
        nqr.pow(q_minus_1_over_2)
    );

    let nqr = Fq6Parameters::NONRESIDUE.clone();
    assert_eq!(Fq6Parameters::FROBENIUS_COEFF_FP6_C1[0], Fq2::one());
    assert_eq!(
        Fq6Parameters::FROBENIUS_COEFF_FP6_C1[1],
        // (q - 1) / 3
        nqr.pow([
            0,
            2147747848,
            9225641637177262208,
            10485338002510381824,
            2721697516351690029,
            865116730873272282,
        ])
    );
    assert_eq!(
        Fq6Parameters::FROBENIUS_COEFF_FP6_C1[2],
        // (q^2 - 1) / 3
        nqr.pow([
            0,
            4295495696,
            13843001656410866112,
            12619659656556060161,
            7863949179600421831,
            17824932913741859728,
            14106172767398838916,
            10392428039924848012,
            18240516690138460617,
            4172794130351852566,
            3052012245806222584,
            121716920077541807,
        ])
    );
    assert_eq!(
        Fq6Parameters::FROBENIUS_COEFF_FP6_C1[3],
        // (q^3 - 1) / 3
        nqr.pow([
            0,
            6443243544,
            13852080057700811712,
            13974760756070671107,
            8474791303092970207,
            15247111912366805117,
            12758405836849367900,
            12846794460241389429,
            13532760723228230634,
            8526073349132696454,
            16145725457143290934,
            3082578096989227526,
            6830165298950233162,
            17302433147824971781,
            12528394078520259374,
            13725539263945993487,
            13263483163850314236,
            17124866627198422,
        ])
    );
    assert_eq!(
        // (q^4 - 1) / 3
        Fq6Parameters::FROBENIUS_COEFF_FP6_C1[4],
        nqr.pow([
            0,
            8590991392,
            9252876841047099008,
            3675693021278299142,
            3349650391287462279,
            13911023695139878569,
            16275601568726772140,
            3014688273044934821,
            1472059661550192864,
            2035158628369053974,
            1713684101115859557,
            9284221196945340952,
            2683613555443221405,
            15729065835098504351,
            15984104543943309973,
            12723413988826450866,
            41170366029762151,
            16008086735292104253,
            15596425719292126299,
            15108390649563014334,
            52931762420676395,
            13732307747156784034,
            7818978216170078201,
            2409369681819975,
        ])
    );
    assert_eq!(
        // (q^5 - 1) / 3
        Fq6Parameters::FROBENIUS_COEFF_FP6_C1[5],
        nqr.pow([
            0,
            10738739240,
            45392006449728000,
            7740996319822131978,
            15478087213462929103,
            17472710613337942216,
            9044015039479477296,
            11586903911799256759,
            2813762230899659887,
            11104325984644937458,
            4302092015276168721,
            8234334058539386480,
            11962170479666254089,
            1461528360371107444,
            6162227182064005366,
            7232630698859231536,
            8412330371357627039,
            363274686941722124,
            4438052204811362316,
            15807890790833648603,
            17923236621871238366,
            10400515525131393666,
            14107173108564206483,
            5590046816246288581,
            16903738071061848386,
            5921390934571784652,
            6669210272393187091,
            1983999480480749653,
            857708666524832295,
            338984378100408,
        ])
    );

    assert_eq!(Fq6Parameters::FROBENIUS_COEFF_FP6_C2[0], Fq2::one());
    assert_eq!(
        Fq6Parameters::FROBENIUS_COEFF_FP6_C2[1],
        // (2 * q - 2) / 3
        nqr.pow([
            0,
            4295495696,
            4539200644972800,
            2523931931311212033,
            5443395032703380059,
            1730233461746544564,
        ])
    );
    assert_eq!(
        Fq6Parameters::FROBENIUS_COEFF_FP6_C2[2],
        // (2 * q^2 - 2) / 3
        nqr.pow([
            0,
            8590991392,
            9239259239112180608,
            6792575239402568707,
            15727898359200843663,
            17203121753774167840,
            9765601461088126217,
            2338112006140144409,
            18034289306567369619,
            8345588260703705133,
            6104024491612445168,
            243433840155083614,
        ])
    );
    assert_eq!(
        Fq6Parameters::FROBENIUS_COEFF_FP6_C2[3],
        // (2 * q^3 - 2) / 3
        nqr.pow([
            0,
            12886487088,
            9257416041692071808,
            9502777438431790599,
            16949582606185940415,
            12047479751024058618,
            7070067599989184185,
            7246844846773227243,
            8618777372746909653,
            17052146698265392909,
            13844706840577030252,
            6165156193978455053,
            13660330597900466324,
            16158122221940391946,
            6610044083330967133,
            9004334454182435359,
            8080222253991076857,
            34249733254396845,
        ])
    );
    assert_eq!(
        Fq6Parameters::FROBENIUS_COEFF_FP6_C2[4],
        // (2 * q^4 - 2) / 3
        nqr.pow([
            0,
            17181982784,
            59009608384646400,
            7351386042556598285,
            6699300782574924558,
            9375303316570205522,
            14104459063743992665,
            6029376546089869643,
            2944119323100385728,
            4070317256738107948,
            3427368202231719114,
            121698320181130288,
            5367227110886442811,
            13011387596487457086,
            13521465014177068331,
            7000083903943350117,
            82340732059524303,
            13569429396874656890,
            12746107364874700983,
            11770037225416477053,
            105863524841352791,
            9017871420604016452,
            15637956432340156403,
            4818739363639950,
        ])
    );
    assert_eq!(
        Fq6Parameters::FROBENIUS_COEFF_FP6_C2[5],
        // (2 * q^5 - 2) / 3
        nqr.pow([
            0,
            21477478480,
            90784012899456000,
            15481992639644263956,
            12509430353216306590,
            16498677152966332817,
            18088030078958954593,
            4727063749888961902,
            5627524461799319775,
            3761907895580323300,
            8604184030552337443,
            16468668117078772960,
            5477596885622956562,
            2923056720742214889,
            12324454364128010732,
            14465261397718463072,
            16824660742715254078,
            726549373883444248,
            8876104409622724632,
            13169037507957745590,
            17399729170032925117,
            2354286976553235717,
            9767602143418861351,
            11180093632492577163,
            15360732068414145156,
            11842781869143569305,
            13338420544786374182,
            3967998960961499306,
            1715417333049664590,
            677968756200816,
        ])
    );

    let nqr = Fq6Parameters::NONRESIDUE.clone();
    assert_eq!(Fq12Parameters::FROBENIUS_COEFF_FP12_C1[0], Fq2::one());
    assert_eq!(
        Fq12Parameters::FROBENIUS_COEFF_FP12_C1[1],
        nqr.pow([
            0,
            1073873924,
            4612820818588631104,
            14466041038109966720,
            1360848758175845014,
            432558365436636141,
        ])
    );
    assert_eq!(
        Fq12Parameters::FROBENIUS_COEFF_FP12_C1[2],
        nqr.pow([
            0,
            2147747848,
            16144872865060208864,
            15533201865132805888,
            3931974589800210915,
            8912466456870929864,
            7053086383699419458,
            14419586056817199814,
            9120258345069230308,
            2086397065175926283,
            10749378159757887100,
            60858460038770903,
        ])
    );
    assert_eq!(
        Fq12Parameters::FROBENIUS_COEFF_FP12_C1[3],
        nqr.pow(vec![
            0,
            3221621772,
            16149412065705181664,
            16210752414890111361,
            13460767688401260911,
            7623555956183402558,
            15602574955279459758,
            6423397230120694714,
            6766380361614115317,
            4263036674566348227,
            8072862728571645467,
            1541289048494613763,
            12638454686329892389,
            8651216573912485890,
            15487569076114905495,
            6862769631972996743,
            6631741581925157118,
            8562433313599211,
        ])
    );
    assert_eq!(
        Fq12Parameters::FROBENIUS_COEFF_FP12_C1[4],
        nqr.pow(vec![
            0,
            4295495696,
            4626438420523549504,
            11061218547493925379,
            10898197232498506947,
            6955511847569939284,
            17361172821218161878,
            1507344136522467410,
            736029830775096432,
            10240951351039302795,
            856842050557929778,
            13865482635327446284,
            10565178814576386510,
            17087904954404027983,
            7992052271971654986,
            15585079031268001241,
            9243957219869656883,
            17227415404500827934,
            7798212859646063149,
            16777567361636282975,
            26465881210338197,
            16089525910433167825,
            13132861144939814908,
            1204684840909987,
        ])
    );
    assert_eq!(
        Fq12Parameters::FROBENIUS_COEFF_FP12_C1[5],
        nqr.pow(vec![
            0,
            5369369620,
            22696003224864000,
            13093870196765841797,
            7739043606731464551,
            8736355306668971108,
            13745379556594514456,
            15016823992754404187,
            1406881115449829943,
            14775535029177244537,
            2151046007638084360,
            13340539066124469048,
            5981085239833127044,
            730764180185553722,
            3081113591032002683,
            12839687386284391576,
            4206165185678813519,
            181637343470861062,
            11442398139260456966,
            7903945395416824301,
            8961618310935619183,
            14423629799420472641,
            16276958591136879049,
            2795023408123144290,
            8451869035530924193,
            12184067504140668134,
            12557977173051369353,
            10215371777095150634,
            428854333262416147,
            169492189050204,
        ])
    );
    assert_eq!(
        Fq12Parameters::FROBENIUS_COEFF_FP12_C1[6],
        nqr.pow(vec![
            0,
            6443243544,
            2338184813809125152,
            16871233222817902855,
            9128410254440549604,
            6358414737107760671,
            7468729267714966206,
            13941217655799039844,
            13346431459341595049,
            7898852908234196435,
            15336792284145632846,
            11388163220812683193,
            3081093040835939111,
            3129346386773272968,
            9927228518296795807,
            7377526438469140102,
            10982033484834786477,
            4338217268182375900,
            9361642498942237193,
            95770243317576020,
            12381045952473046729,
            13583764370608885618,
            6914038083938247018,
            11758473285549676044,
            11115849729267060548,
            17409651029690199254,
            1556144995363291,
            201277993468444800,
            4475903768997462991,
            14292114947578464987,
            11232182361686567882,
            8023223058519243510,
            8460645903912770899,
            3728985165610438147,
            9557646917806612605,
            23846570632805,
        ])
    );
    assert_eq!(
        Fq12Parameters::FROBENIUS_COEFF_FP12_C1[7],
        nqr.pow(vec![
            0,
            7517117468,
            11572904852276332960,
            16955833485762150793,
            13861723677667078339,
            3351040880758239373,
            10187361698481445917,
            13133145996058314825,
            9880302656011032102,
            13564594818746138045,
            7715332987162020029,
            8507406103252913234,
            3841936539250661075,
            1805614599139944761,
            203070458119069734,
            6012523277532087348,
            17138495458943996700,
            554721240220003420,
            12440154572222299423,
            16155796765560152746,
            16906961404535980601,
            16988772101138510128,
            10387531528977690481,
            11239496299840526324,
            17694275296936601295,
            8678865958464098749,
            9320448299481179664,
            4686598267618731179,
            13742235150573415491,
            12374458970711564027,
            269339566263283418,
            10481813823058826132,
            717376942815119705,
            13666875719826635227,
            12280940418647842802,
            9519820163220694865,
            617514567176427309,
            12583439039711230168,
            17395426026137895032,
            8613762055003349927,
            10640527628973352346,
            3355074556131,
        ])
    );
    assert_eq!(
        Fq12Parameters::FROBENIUS_COEFF_FP12_C1[8],
        nqr.pow(vec![
            0,
            8590991392,
            9280112044916935808,
            7910196845710627852,
            14384733437153267341,
            13506577854283292941,
            10857690159372235189,
            272593318034492385,
            12923380669675470042,
            2100692327127267929,
            13230996614761816179,
            7168823705819996997,
            13956769771436729488,
            5681613690936661631,
            11503781881124304312,
            17948477454354511940,
            2585515360459221833,
            8259262986641375999,
            10186761771126765300,
            138677943505860314,
            5250511679047117534,
            14897331219308048187,
            569947668211672644,
            9595808377407591607,
            3673903986281427142,
            4870403510282402826,
            2096666344982093356,
            12051547529286607576,
            1844298388091022287,
            1306834614315497461,
            6474597614146232073,
            12878772894351383459,
            14872419472349185632,
            14433096506366908483,
            6489113408900444583,
            12661100796117668959,
            91512520663806584,
            1296368661047785873,
            10476259365810701958,
            2053707460771923724,
            6065289503083254206,
            6873892662994326422,
            12294887061544157066,
            1245036048899490368,
            878322480486140511,
            3487162457538806614,
            11345805986413709742,
            472039583826,
        ])
    );
    assert_eq!(
        Fq12Parameters::FROBENIUS_COEFF_FP12_C1[9],
        nqr.pow(vec![
            0,
            9664865316,
            13906550465440485312,
            2743593236484927887,
            15240256226051785162,
            5666101964554392148,
            15486630721293003604,
            11085304049122138569,
            3743538961979692360,
            10479969010368013395,
            8805888377599446479,
            1551988704758378836,
            14640760467567754482,
            17095114512289855362,
            17093001622846958786,
            7870427257463436135,
            5710915045644416575,
            4183002850442633052,
            2183746690979144293,
            10361570738316409762,
            1823519866094145509,
            1979140920109133656,
            3934637357680171129,
            483707244939023109,
            1565783218869661592,
            2070037913090854503,
            14241273410721161637,
            10692540910451051929,
            15363289899804154929,
            15065156723151968355,
            14866176115383228793,
            1754717951053774339,
            6219089168165776907,
            2843316444832961545,
            7769394471483342945,
            6209130940448592686,
            9895465624289118027,
            3768298914925623782,
            9749235554500066774,
            1141391204377577507,
            14056362858548820458,
            16607719785754580914,
            15673364900352983470,
            580200748086459451,
            10509892878537744119,
            15607998199247548635,
            12824174783078822258,
            6481915545385316101,
            1957565790411709196,
            349002659769060597,
            16415225026456770303,
            1404654246286235655,
            6356350610771785647,
            66413239101,
        ])
    );
    assert_eq!(
        Fq12Parameters::FROBENIUS_COEFF_FP12_C1[10],
        nqr.pow(vec![
            0,
            10738739240,
            7005476040137429856,
            14465292591906644755,
            14621431796216200706,
            14525857187455188078,
            10700440084621593392,
            12242072263566209620,
            8264012505466920242,
            14417507739553060693,
            16334656267535174813,
            5608731228249949498,
            13823809389348153905,
            7430045463521736622,
            2833765277650240876,
            1309907354854159809,
            5021927718196352462,
            15874794418293664475,
            15608826217061854163,
            1595251607803040137,
            8888388384778788,
            7831801103413923870,
            6005456231811099037,
            6096848355567875777,
            15850729980334142493,
            1624611499560916485,
            6184106485895126123,
            6387589785991100343,
            4441436714440535396,
            10674155509432621698,
            10042502358973557699,
            2576341273935999356,
            18156873608396522635,
            617515025121520734,
            9881222500003602426,
            4756903745958102601,
            2006577979583400419,
            8327703057715215960,
            15321419544975328394,
            12628506519071548335,
            13177222782269231712,
            2411912121646877137,
            10868769581993413381,
            16185388380426370530,
            5539608759663320000,
            10723559808687022239,
            14719067825964600572,
            460658681716792668,
            10091533243992712257,
            14530462048548757687,
            12200158948752460980,
            18189052240063574460,
            2913087417685331074,
            4736340601006000495,
            2372095410131630914,
            16516702885801480121,
            17997825935487833985,
            9257298735064255248,
            1536693001487705261,
            9343958598,
        ])
    );
    assert_eq!(
        Fq12Parameters::FROBENIUS_COEFF_FP12_C1[11],
        nqr.pow(vec![
            0,
            11812613164,
            7023632842717321056,
            744332624668717463,
            4371722958200983232,
            11636017688460936382,
            17152858786832518158,
            5682050630732907671,
            8727981457605556050,
            10057238357286604968,
            4003863022247557308,
            7604055838654582988,
            12261673555141771681,
            15883167835850577292,
            5714157687636561277,
            14618404407668008357,
            16402832752841580893,
            5972425114553822180,
            12155917574793055344,
            8808187807425270992,
            9211841961893389931,
            12725660675278824272,
            12414773954685858373,
            9792868552842682366,
            5719933945026251524,
            12537786363063843971,
            18187538106036260758,
            12740838136551074431,
            7969056988229557200,
            872948111533580865,
            17726711914655473767,
            770549698442050016,
            16298837429388465735,
            5178749236651908864,
            18138853862298625203,
            7793664090201224685,
            11193293738282492653,
            9873076932945272135,
            15807756096133577526,
            10083282791902955936,
            2195531347646939012,
            16238427073634898251,
            5493913095325849250,
            14170967139044366675,
            12805124499495589934,
            11391084290952387861,
            9364144152903456451,
            1899440368800537432,
            2322858144550417821,
            10240253012091294361,
            15065302116186703837,
            12448570265819329800,
            6343718079088714168,
            4380121812362314977,
            8816058002082775395,
            8311367630575408773,
            1256689681524789185,
            15712089605351878639,
            2370617092363650854,
            1292526687163685054,
            1968496925520971753,
            14509852769280912127,
            83170576256002463,
            15256590279050317540,
            18173816319152067096,
            1314640927,
        ])
    );
}

#[test]
fn test_fq_repr_from() {
    assert_eq!(
        BigInteger384::from(100),
        BigInteger384([100, 0, 0, 0, 0, 0])
    );
}

#[test]
fn test_fq_repr_is_odd() {
    assert!(!BigInteger384::from(0).is_odd());
    assert!(BigInteger384::from(0).is_even());
    assert!(BigInteger384::from(1).is_odd());
    assert!(!BigInteger384::from(1).is_even());
    assert!(!BigInteger384::from(324834872).is_odd());
    assert!(BigInteger384::from(324834872).is_even());
    assert!(BigInteger384::from(324834873).is_odd());
    assert!(!BigInteger384::from(324834873).is_even());
}

#[test]
fn test_fq_repr_is_zero() {
    assert!(BigInteger384::from(0).is_zero());
    assert!(!BigInteger384::from(1).is_zero());
    assert!(!BigInteger384([0, 0, 0, 0, 1, 0]).is_zero());
}

#[test]
fn test_fq_repr_div2() {
    let mut a = BigInteger384([
        0x8b0ad39f8dd7482a,
        0x147221c9a7178b69,
        0x54764cb08d8a6aa0,
        0x8519d708e1d83041,
        0x41f82777bd13fdb,
        0xf43944578f9b771b,
    ]);
    a.div2();
    assert_eq!(
        a,
        BigInteger384([
            0xc58569cfc6eba415,
            0xa3910e4d38bc5b4,
            0xaa3b265846c53550,
            0xc28ceb8470ec1820,
            0x820fc13bbde89fed,
            0x7a1ca22bc7cdbb8d,
        ])
    );
    for _ in 0..10 {
        a.div2();
    }
    assert_eq!(
        a,
        BigInteger384([
            0x6d31615a73f1bae9,
            0x54028e443934e2f1,
            0x82a8ec99611b14d,
            0xfb70a33ae11c3b06,
            0xe36083f04eef7a27,
            0x1e87288af1f36e,
        ])
    );
    for _ in 0..300 {
        a.div2();
    }
    assert_eq!(
        a,
        BigInteger384([0x7288af1f36ee3608, 0x1e8, 0x0, 0x0, 0x0, 0x0])
    );
    for _ in 0..50 {
        a.div2();
    }
    assert_eq!(a, BigInteger384([0x7a1ca2, 0x0, 0x0, 0x0, 0x0, 0x0]));
    for _ in 0..22 {
        a.div2();
    }
    assert_eq!(a, BigInteger384([0x1, 0x0, 0x0, 0x0, 0x0, 0x0]));
    a.div2();
    assert!(a.is_zero());
}

#[test]
fn test_fq_repr_divn() {
    let mut a = BigInteger384([
        0xaa5cdd6172847ffd,
        0x43242c06aed55287,
        0x9ddd5b312f3dd104,
        0xc5541fd48046b7e7,
        0x16080cf4071e0b05,
        0x1225f2901aea514e,
    ]);
    a.divn(0);
    assert_eq!(
        a,
        BigInteger384([
            0xaa5cdd6172847ffd,
            0x43242c06aed55287,
            0x9ddd5b312f3dd104,
            0xc5541fd48046b7e7,
            0x16080cf4071e0b05,
            0x1225f2901aea514e,
        ])
    );
    a.divn(1);
    assert_eq!(
        a,
        BigInteger384([
            0xd52e6eb0b9423ffe,
            0x21921603576aa943,
            0xceeead98979ee882,
            0xe2aa0fea40235bf3,
            0xb04067a038f0582,
            0x912f9480d7528a7,
        ])
    );
    a.divn(50);
    assert_eq!(
        a,
        BigInteger384([
            0x8580d5daaa50f54b,
            0xab6625e7ba208864,
            0x83fa9008d6fcf3bb,
            0x19e80e3c160b8aa,
            0xbe52035d4a29c2c1,
            0x244,
        ])
    );
    a.divn(130);
    assert_eq!(
        a,
        BigInteger384([
            0xa0fea40235bf3cee,
            0x4067a038f0582e2a,
            0x2f9480d7528a70b0,
            0x91,
            0x0,
            0x0,
        ])
    );
    a.divn(64);
    assert_eq!(
        a,
        BigInteger384([0x4067a038f0582e2a, 0x2f9480d7528a70b0, 0x91, 0x0, 0x0, 0x0])
    );
}

#[test]
fn test_fq_repr_mul2() {
    let mut a = BigInteger384::from(23712937547);
    a.mul2();
    assert_eq!(a, BigInteger384([0xb0acd6c96, 0x0, 0x0, 0x0, 0x0, 0x0]));
    for _ in 0..60 {
        a.mul2();
    }
    assert_eq!(
        a,
        BigInteger384([0x6000000000000000, 0xb0acd6c9, 0x0, 0x0, 0x0, 0x0])
    );
    for _ in 0..300 {
        a.mul2();
    }
    assert_eq!(
        a,
        BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0xcd6c960000000000])
    );
    for _ in 0..17 {
        a.mul2();
    }
    assert_eq!(
        a,
        BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x2c00000000000000])
    );
    for _ in 0..6 {
        a.mul2();
    }
    assert!(a.is_zero());
}

#[test]
fn test_fq_repr_num_bits() {
    let mut a = BigInteger384::from(0);
    assert_eq!(0, a.num_bits());
    a = BigInteger384::from(1);
    for i in 1..385 {
        assert_eq!(i, a.num_bits());
        a.mul2();
    }
    assert_eq!(0, a.num_bits());
}

#[test]
fn test_fq_repr_sub_noborrow() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut t = BigInteger384([
        0x827a4a08041ebd9,
        0x3c239f3dcc8f0d6b,
        0x9ab46a912d555364,
        0x196936b17b43910b,
        0xad0eb3948a5c34fd,
        0xd56f7b5ab8b5ce8,
    ]);
    t.sub_noborrow(&BigInteger384([
        0xc7867917187ca02b,
        0x5d75679d4911ffef,
        0x8c5b3e48b1a71c15,
        0x6a427ae846fd66aa,
        0x7a37e7265ee1eaf9,
        0x7c0577a26f59d5,
    ]));
    assert!(
        t == BigInteger384([
            0x40a12b8967c54bae,
            0xdeae37a0837d0d7b,
            0xe592c487bae374e,
            0xaf26bbc934462a61,
            0x32d6cc6e2b7a4a03,
            0xcdaf23e091c0313,
        ])
    );

    for _ in 0..1000 {
        let mut a = BigInteger384::rand(&mut rng);
        a.0[5] >>= 30;
        let mut b = a;
        for _ in 0..10 {
            b.mul2();
        }
        let mut c = b;
        for _ in 0..10 {
            c.mul2();
        }

        assert!(a < b);
        assert!(b < c);

        let mut csub_ba = c;
        csub_ba.sub_noborrow(&b);
        csub_ba.sub_noborrow(&a);

        let mut csub_ab = c;
        csub_ab.sub_noborrow(&a);
        csub_ab.sub_noborrow(&b);

        assert_eq!(csub_ab, csub_ba);
    }

    // Subtracting q+1 from q should produce -1 (mod 2**384)
    let mut qplusone = BigInteger384([
        0xb9feffffffffaaab,
        0x1eabfffeb153ffff,
        0x6730d2a0f6b0f624,
        0x64774b84f38512bf,
        0x4b1ba7b6434bacd7,
        0x1a0111ea397fe69a,
    ]);
    qplusone.sub_noborrow(&BigInteger384([
        0xb9feffffffffaaac,
        0x1eabfffeb153ffff,
        0x6730d2a0f6b0f624,
        0x64774b84f38512bf,
        0x4b1ba7b6434bacd7,
        0x1a0111ea397fe69a,
    ]));
    assert_eq!(
        qplusone,
        BigInteger384([
            0xffffffffffffffff,
            0xffffffffffffffff,
            0xffffffffffffffff,
            0xffffffffffffffff,
            0xffffffffffffffff,
            0xffffffffffffffff,
        ])
    );
}

#[test]
fn test_fq_repr_add_nocarry() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut t = BigInteger384([
        0x827a4a08041ebd9,
        0x3c239f3dcc8f0d6b,
        0x9ab46a912d555364,
        0x196936b17b43910b,
        0xad0eb3948a5c34fd,
        0xd56f7b5ab8b5ce8,
    ]);
    t.add_nocarry(&BigInteger384([
        0xc7867917187ca02b,
        0x5d75679d4911ffef,
        0x8c5b3e48b1a71c15,
        0x6a427ae846fd66aa,
        0x7a37e7265ee1eaf9,
        0x7c0577a26f59d5,
    ]));
    assert!(
        t == BigInteger384([
            0xcfae1db798be8c04,
            0x999906db15a10d5a,
            0x270fa8d9defc6f79,
            0x83abb199c240f7b6,
            0x27469abae93e1ff6,
            0xdd2fd2d4dfab6be,
        ])
    );

    // Test for the associativity of addition.
    for _ in 0..1000 {
        let mut a = BigInteger384::rand(&mut rng);
        let mut b = BigInteger384::rand(&mut rng);
        let mut c = BigInteger384::rand(&mut rng);

        // Unset the first few bits, so that overflow won't occur.
        a.0[5] >>= 3;
        b.0[5] >>= 3;
        c.0[5] >>= 3;

        let mut abc = a;
        abc.add_nocarry(&b);
        abc.add_nocarry(&c);

        let mut acb = a;
        acb.add_nocarry(&c);
        acb.add_nocarry(&b);

        let mut bac = b;
        bac.add_nocarry(&a);
        bac.add_nocarry(&c);

        let mut bca = b;
        bca.add_nocarry(&c);
        bca.add_nocarry(&a);

        let mut cab = c;
        cab.add_nocarry(&a);
        cab.add_nocarry(&b);

        let mut cba = c;
        cba.add_nocarry(&b);
        cba.add_nocarry(&a);

        assert_eq!(abc, acb);
        assert_eq!(abc, bac);
        assert_eq!(abc, bca);
        assert_eq!(abc, cab);
        assert_eq!(abc, cba);
    }

    // Adding 1 to (2^384 - 1) should produce zero
    let mut x = BigInteger384([
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
    ]);
    x.add_nocarry(&BigInteger384::from(1));
    assert!(x.is_zero());
}

#[test]
fn test_fq_add_assign() {
    {
        // Random number
        let mut tmp = Fq::new(BigInteger384([
            0x624434821df92b69,
            0x503260c04fd2e2ea,
            0xd9df726e0d16e8ce,
            0xfbcb39adfd5dfaeb,
            0x86b8a22b0c88b112,
            0x165a2ed809e4201b,
        ]));
        // Test that adding zero has no effect.
        tmp.add_assign(&Fq::new(BigInteger384::from(0)));
        assert_eq!(
            tmp,
            Fq::new(BigInteger384([
                0x624434821df92b69,
                0x503260c04fd2e2ea,
                0xd9df726e0d16e8ce,
                0xfbcb39adfd5dfaeb,
                0x86b8a22b0c88b112,
                0x165a2ed809e4201b,
            ]))
        );
    }

    // Test associativity

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..1000 {
        // Generate a, b, c and ensure (a + b) + c == a + (b + c).
        let a = Fq::rand(&mut rng);
        let b = Fq::rand(&mut rng);
        let c = Fq::rand(&mut rng);

        let mut tmp1 = a;
        tmp1.add_assign(&b);
        tmp1.add_assign(&c);

        let mut tmp2 = b;
        tmp2.add_assign(&c);
        tmp2.add_assign(&a);

        assert_eq!(tmp1, tmp2);
    }
}

#[test]
fn test_fq_mul_assign() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..1000000 {
        // Ensure that (a * b) * c = a * (b * c)
        let a = Fq::rand(&mut rng);
        let b = Fq::rand(&mut rng);
        let c = Fq::rand(&mut rng);

        let mut tmp1 = a;
        tmp1.mul_assign(&b);
        tmp1.mul_assign(&c);

        let mut tmp2 = b;
        tmp2.mul_assign(&c);
        tmp2.mul_assign(&a);

        assert_eq!(tmp1, tmp2);
    }

    for _ in 0..1000000 {
        // Ensure that r * (a + b + c) = r*a + r*b + r*c

        let r = Fq::rand(&mut rng);
        let mut a = Fq::rand(&mut rng);
        let mut b = Fq::rand(&mut rng);
        let mut c = Fq::rand(&mut rng);

        let mut tmp1 = a;
        tmp1.add_assign(&b);
        tmp1.add_assign(&c);
        tmp1.mul_assign(&r);

        a.mul_assign(&r);
        b.mul_assign(&r);
        c.mul_assign(&r);

        a.add_assign(&b);
        a.add_assign(&c);

        assert_eq!(tmp1, a);
    }
}

#[test]
fn test_fq_squaring() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..1000000 {
        // Ensure that (a * a) = a^2
        let a = Fq::rand(&mut rng);

        let mut tmp = a;
        tmp.square_in_place();

        let mut tmp2 = a;
        tmp2.mul_assign(&a);

        assert_eq!(tmp, tmp2);
    }
}

#[test]
fn test_fq_inverse() {
    assert!(Fq::zero().inverse().is_none());

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let one = Fq::one();

    for _ in 0..1000 {
        // Ensure that a * a^-1 = 1
        let mut a = Fq::rand(&mut rng);
        let ainv = a.inverse().unwrap();
        a.mul_assign(&ainv);
        assert_eq!(a, one);
    }
}

#[test]
fn test_fq_double_in_place() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..1000 {
        // Ensure doubling a is equivalent to adding a to itself.
        let mut a = Fq::rand(&mut rng);
        let mut b = a;
        b.add_assign(&a);
        a.double_in_place();
        assert_eq!(a, b);
    }
}

#[test]
fn test_fq_negate() {
    {
        let a = -Fq::zero();

        assert!(a.is_zero());
    }

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..1000 {
        // Ensure (a - (-a)) = 0.
        let mut a = Fq::rand(&mut rng);
        let b = -a;
        a.add_assign(&b);

        assert!(a.is_zero());
    }
}

#[test]
fn test_fq_pow() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for i in 0..1000 {
        // Exponentiate by various small numbers and ensure it consists with repeated
        // multiplication.
        let a = Fq::rand(&mut rng);
        let target = a.pow(&[i]);
        let mut c = Fq::one();
        for _ in 0..i {
            c.mul_assign(&a);
        }
        assert_eq!(c, target);
    }

    for _ in 0..1000 {
        // Exponentiating by the modulus should have no effect in a prime field.
        let a = Fq::rand(&mut rng);

        assert_eq!(a, a.pow(Fq::characteristic()));
    }
}

#[test]
fn test_fq_sqrt() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    assert_eq!(Fq::zero().sqrt().unwrap(), Fq::zero());

    for _ in 0..1000 {
        // Ensure sqrt(a^2) = a or -a
        let a = Fq::rand(&mut rng);
        let nega = -a;
        let mut b = a;
        b.square_in_place();

        let b = b.sqrt().unwrap();

        assert!(a == b || nega == b);
    }

    for _ in 0..1000 {
        // Ensure sqrt(a)^2 = a for random a
        let a = Fq::rand(&mut rng);

        if let Some(mut tmp) = a.sqrt() {
            tmp.square_in_place();

            assert_eq!(a, tmp);
        }
    }
}

#[test]
fn test_fq_num_bits() {
    assert_eq!(FqParameters::MODULUS_BITS, 382);
    assert_eq!(FqParameters::CAPACITY, 381);
}

#[test]
fn test_fq_ordering() {
    // BigInteger384's ordering is well-tested, but we still need to make sure the
    // Fq elements aren't being compared in Montgomery form.
    for i in 0..100 {
        assert!(Fq::from_repr(BigInteger384::from(i + 1)) > Fq::from_repr(BigInteger384::from(i)));
    }
}

// #[test]
// fn fq_repr_tests() {
//     ::tests::repr::random_repr_tests::<BigInteger384>();
// }

#[test]
fn test_fq_legendre() {
    use crate::fields::LegendreSymbol::*;

    assert_eq!(QuadraticResidue, Fq::one().legendre());
    assert_eq!(Zero, Fq::zero().legendre());

    assert_eq!(
        QuadraticResidue,
        Fq::from_repr(BigInteger384::from(2)).legendre()
    );
    assert_eq!(
        QuadraticResidue,
        Fq::from_repr(BigInteger384::from(4)).legendre()
    );
}

#[test]
fn test_fq2_ordering() {
    let mut a = Fq2::new(Fq::zero(), Fq::zero());

    let mut b = a.clone();

    assert!(a.cmp(&b) == Ordering::Equal);
    b.c0.add_assign(&Fq::one());
    assert!(a.cmp(&b) == Ordering::Less);
    a.c0.add_assign(&Fq::one());
    assert!(a.cmp(&b) == Ordering::Equal);
    b.c1.add_assign(&Fq::one());
    assert!(a.cmp(&b) == Ordering::Less);
    a.c0.add_assign(&Fq::one());
    assert!(a.cmp(&b) == Ordering::Less);
    a.c1.add_assign(&Fq::one());
    assert!(a.cmp(&b) == Ordering::Greater);
    b.c0.add_assign(&Fq::one());
    assert!(a.cmp(&b) == Ordering::Equal);
}

#[test]
fn test_fq2_basics() {
    assert_eq!(Fq2::new(Fq::zero(), Fq::zero(),), Fq2::zero());
    assert_eq!(Fq2::new(Fq::one(), Fq::zero(),), Fq2::one());
    assert!(Fq2::zero().is_zero());
    assert!(!Fq2::one().is_zero());
    assert!(!Fq2::new(Fq::zero(), Fq::one(),).is_zero());
}

#[test]
fn test_fq2_inverse() {
    assert!(Fq2::zero().inverse().is_none());

    let a = Fq2::new(
        Fq::from_repr(BigInteger384([
            0x57712f649f231bf3,
            0xdd7cddf2f4603dcb,
            0xf350693f560f254,
            0xd6bd3fdd0fc3b5b4,
            0xdb7c5e8258036911,
            0xe84c7a88c502271,
        ])),
        Fq::from_repr(BigInteger384([
            0xc04ec3bd1ede925e,
            0x54d5f95f26628255,
            0xf2ff26a3f286464b,
            0x72d8f68e3c11222b,
            0xfa5a1b93b0813f7a,
            0x168178c674fe5e1c,
        ])),
    );
    let a = a.inverse().unwrap();
    assert_eq!(
        a,
        Fq2::new(
            Fq::from_repr(BigInteger384([
                0x76e2a7728a451901,
                0x2717121041121420,
                0x2fdb1e07d9ac7a6e,
                0xb76e522717bbcdf5,
                0x5fec3d6b14481c4,
                0xc639eebd0f0f919
            ])),
            Fq::from_repr(BigInteger384([
                0xb219dc3c41a7fc1e,
                0xef737caadc1e,
                0xf8b3637027635f79,
                0x35f7aaba0c8c876a,
                0x51ab431de9fad2e5,
                0xa52a6c3c748fbec
            ])),
        )
    );
}

#[test]
fn test_fq2_addition() {
    let mut a = Fq2::new(
        Fq::from_repr(BigInteger384([
            0x2d0078036923ffc7,
            0x11e59ea221a3b6d2,
            0x8b1a52e0a90f59ed,
            0xb966ce3bc2108b13,
            0xccc649c4b9532bf3,
            0xf8d295b2ded9dc,
        ])),
        Fq::from_repr(BigInteger384([
            0x977df6efcdaee0db,
            0x946ae52d684fa7ed,
            0xbe203411c66fb3a5,
            0xb3f8afc0ee248cad,
            0x4e464dea5bcfd41e,
            0x12d1137b8a6a837,
        ])),
    );
    a.add_assign(&Fq2::new(
        Fq::from_repr(BigInteger384([
            0x619a02d78dc70ef2,
            0xb93adfc9119e33e8,
            0x4bf0b99a9f0dca12,
            0x3b88899a42a6318f,
            0x986a4a62fa82a49d,
            0x13ce433fa26027f5,
        ])),
        Fq::from_repr(BigInteger384([
            0x66323bf80b58b9b9,
            0xa1379b6facf6e596,
            0x402aef1fb797e32f,
            0x2236f55246d0d44d,
            0x4c8c1800eb104566,
            0x11d6e20e986c2085,
        ])),
    ));
    assert_eq!(
        a,
        Fq2::new(
            Fq::from_repr(BigInteger384([
                0x8e9a7adaf6eb0eb9,
                0xcb207e6b3341eaba,
                0xd70b0c7b481d23ff,
                0xf4ef57d604b6bca2,
                0x65309427b3d5d090,
                0x14c715d5553f01d2,
            ])),
            Fq::from_repr(BigInteger384([
                0xfdb032e7d9079a94,
                0x35a2809d15468d83,
                0xfe4b23317e0796d5,
                0xd62fa51334f560fa,
                0x9ad265eb46e01984,
                0x1303f3465112c8bc,
            ])),
        )
    );
}

#[test]
fn test_fq2_doubling() {
    let mut a = Fq2::new(
        Fq::from_repr(BigInteger384([
            0x2d0078036923ffc7,
            0x11e59ea221a3b6d2,
            0x8b1a52e0a90f59ed,
            0xb966ce3bc2108b13,
            0xccc649c4b9532bf3,
            0xf8d295b2ded9dc,
        ])),
        Fq::from_repr(BigInteger384([
            0x977df6efcdaee0db,
            0x946ae52d684fa7ed,
            0xbe203411c66fb3a5,
            0xb3f8afc0ee248cad,
            0x4e464dea5bcfd41e,
            0x12d1137b8a6a837,
        ])),
    );
    a.double_in_place();
    assert_eq!(
        a,
        Fq2::new(
            Fq::from_repr(BigInteger384([
                0x5a00f006d247ff8e,
                0x23cb3d4443476da4,
                0x1634a5c1521eb3da,
                0x72cd9c7784211627,
                0x998c938972a657e7,
                0x1f1a52b65bdb3b9,
            ])),
            Fq::from_repr(BigInteger384([
                0x2efbeddf9b5dc1b6,
                0x28d5ca5ad09f4fdb,
                0x7c4068238cdf674b,
                0x67f15f81dc49195b,
                0x9c8c9bd4b79fa83d,
                0x25a226f714d506e,
            ])),
        )
    );
}

#[test]
fn test_fq2_legendre() {
    use crate::fields::LegendreSymbol::*;

    assert_eq!(Zero, Fq2::zero().legendre());
    // i^2 = -1
    let mut m1 = -Fq2::one();
    assert_eq!(QuadraticResidue, m1.legendre());
    m1 = Fq6Parameters::mul_fp2_by_nonresidue(&m1);
    assert_eq!(QuadraticNonResidue, m1.legendre());
}

#[test]
fn test_fq6_mul_nonresidue() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let nqr = Fq6::new(Fq2::zero(), Fq2::one(), Fq2::zero());

    for _ in 0..1000 {
        let mut a = Fq6::rand(&mut rng);
        let mut b = a;
        a = Fq12Parameters::mul_fp6_by_nonresidue(&a);
        b.mul_assign(&nqr);

        assert_eq!(a, b);
    }
}

#[test]
fn test_fq6_mul_by_1() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..1000 {
        let c1 = Fq2::rand(&mut rng);
        let mut a = Fq6::rand(&mut rng);
        let mut b = a;

        a.mul_by_1(&c1);
        b.mul_assign(&Fq6::new(Fq2::zero(), c1, Fq2::zero()));

        assert_eq!(a, b);
    }
}

#[test]
fn test_fq6_mul_by_01() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..1000 {
        let c0 = Fq2::rand(&mut rng);
        let c1 = Fq2::rand(&mut rng);
        let mut a = Fq6::rand(&mut rng);
        let mut b = a;

        a.mul_by_01(&c0, &c1);
        b.mul_assign(&Fq6::new(c0, c1, Fq2::zero()));

        assert_eq!(a, b);
    }
}

#[test]
fn test_fq12_mul_by_014() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..1000 {
        let c0 = Fq2::rand(&mut rng);
        let c1 = Fq2::rand(&mut rng);
        let c5 = Fq2::rand(&mut rng);
        let mut a = Fq12::rand(&mut rng);
        let mut b = a;

        a.mul_by_014(&c0, &c1, &c5);
        b.mul_assign(&Fq12::new(
            Fq6::new(c0, c1, Fq2::zero()),
            Fq6::new(Fq2::zero(), c5, Fq2::zero()),
        ));

        assert_eq!(a, b);
    }
}
