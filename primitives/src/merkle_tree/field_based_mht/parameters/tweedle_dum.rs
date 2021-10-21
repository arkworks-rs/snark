use algebra::{biginteger::BigInteger256, field_new, fields::tweedle::Fq as TweedleFq};

use crate::{crh::poseidon::TweedleFqPoseidonHash, FieldBasedMerkleTreePrecomputedZeroConstants};

// PoseidonHash("This represents an empty Merkle Root for a TweedleDumPoseidonHash based Merkle Tree.")
pub const TWEEDLE_DUM_PHANTOM_MERKLE_ROOT: TweedleFq = field_new!(
    TweedleFq,
    BigInteger256([
        4904612841964010928,
        11732269297394565570,
        6035769393555604445,
        3097632584773363944
    ])
);

pub const TWEEDLE_DUM_MHT_POSEIDON_PARAMETERS: FieldBasedMerkleTreePrecomputedZeroConstants<
    'static,
    TweedleFqPoseidonHash,
> = FieldBasedMerkleTreePrecomputedZeroConstants {
    nodes: &[
        field_new!(TweedleFq, BigInteger256([0, 0, 0, 0])),
        field_new!(
            TweedleFq,
            BigInteger256([
                6139372262132429377,
                6616513606251009568,
                10180936183985522127,
                1871256090268734000
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                16581232002391324140,
                10128427957487402226,
                14615445748711189182,
                106768847887941616
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                7422282436897880714,
                2759687527703036765,
                18261471787718551385,
                4333102577557650426
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                8960034495897434649,
                5976620936990978514,
                11770448562828825262,
                420408715434692497
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                1639501920800041527,
                7047118941613450008,
                8439584256723729208,
                2340548282573108138
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                17474590584127356599,
                17005018111462329626,
                9038520322116564398,
                2842033007168063862
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                17907603573830123375,
                7485449917333291794,
                8497334770128174690,
                985616778997111667
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                8437238815112828116,
                7558012184601411834,
                4810203110390380299,
                3693304440212097843
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                15867381531984164088,
                14436926832167307336,
                7771858470711342912,
                4382274482182735339
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                6780535927106125369,
                893601027906750002,
                2168364268659532015,
                1398450762353999324
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                761179424970697872,
                18228584916778786433,
                3789686406673394224,
                315352039785223877
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                10506753105967550322,
                10479934490182553064,
                4711731144016154768,
                2619376882049526951
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                5491581052635507564,
                11596637832753061833,
                15074841229727970086,
                51841591242253653
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                10751695168434495081,
                11265029839378255208,
                2292162601217563466,
                2013346614746080729
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                1065409681300367041,
                14926577987640137998,
                13035113752352820427,
                777620396488025701
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                17671398566235914769,
                11145386535386684404,
                14965778970796212802,
                2548163089059524405
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                4147123003843975335,
                16801467547722941168,
                9567918343610890265,
                1232574073465013165
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                11272246842976262912,
                5714009969114434676,
                13201780159870357255,
                1761070557194703564
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                13412942200607847147,
                13612865033638762108,
                4989519727612199232,
                224999250619481366
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                16131593917109878013,
                2118354752416554125,
                2560734345425875426,
                265596394308507171
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                8350852403008315484,
                14236833013206637959,
                1704211568971383225,
                3223189361656905860
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                2280418127254463522,
                11760863143783095452,
                14517820451477144924,
                388576072386685144
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                14719543229856123753,
                7602837078929085760,
                13722828357997036314,
                3711455166874429393
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                14339280075298071695,
                9047408951566991594,
                2938258787476493148,
                154196581154821338
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                14855917590773565952,
                16320462695690950774,
                13608606206008960107,
                4021935482326785433
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                5664553553817478280,
                15018914933284007567,
                12910027197340055896,
                69535126663687267
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                12939490041852130828,
                13667911504789686466,
                6125435652755305576,
                1414802449081846718
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                8636287925162701366,
                8473826958806005476,
                11034325557957988141,
                2161057426287045498
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                3540680633243804891,
                11704818265718910048,
                10222301835080698341,
                1142989787151462434
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                4048987655390738064,
                8466268957972885142,
                600328630911781217,
                4059422043777361500
            ])
        ),
        field_new!(
            TweedleFq,
            BigInteger256([
                10396877418832343447,
                4626792818372145897,
                1925158989055659802,
                1880821489306116509
            ])
        ),
    ],
    merkle_arity: 2,
};

#[cfg(test)]
mod test {

    use super::{TWEEDLE_DUM_MHT_POSEIDON_PARAMETERS, TWEEDLE_DUM_PHANTOM_MERKLE_ROOT};
    use crate::{
        crh::TweedleFqPoseidonHash,
        merkle_tree::field_based_mht::parameters::{
            generate_mht_empty_nodes, generate_phantom_merkle_root_from_magic_string,
        },
        FieldBasedMerkleTreePrecomputedZeroConstants,
    };
    use algebra::{fields::tweedle::Fq, Field};

    #[ignore]
    #[test]
    fn test_generate_tweedle_fq_phantom_merkle_root() {
        let expected_root = generate_phantom_merkle_root_from_magic_string::<
            Fq,
            TweedleFqPoseidonHash,
        >(
            "This represents an empty Merkle Root for a TweedleDumPoseidonHash based Merkle Tree.",
        );
        assert_eq!(expected_root, TWEEDLE_DUM_PHANTOM_MERKLE_ROOT);
    }

    #[ignore]
    #[test]
    fn test_generate_binary_tweedle_fq_mht_empty_nodes() {
        let merkle_arity = 2;
        let max_height = 32;

        let empty_nodes = generate_mht_empty_nodes::<Fq, TweedleFqPoseidonHash>(
            merkle_arity,
            max_height,
            Fq::zero(),
        );
        assert_eq!(empty_nodes.len(), max_height);

        let params = FieldBasedMerkleTreePrecomputedZeroConstants::<TweedleFqPoseidonHash> {
            nodes: empty_nodes.as_slice(),
            merkle_arity,
        };
        assert_eq!(params, TWEEDLE_DUM_MHT_POSEIDON_PARAMETERS)
    }
}
