use algebra::{biginteger::BigInteger256, field_new, fields::tweedle::Fr as TweedleFr};

use crate::{crh::poseidon::TweedleFrPoseidonHash, FieldBasedMerkleTreePrecomputedZeroConstants};

// PoseidonHash("This represents an empty Merkle Root for a TweedleDeePoseidonHash based Merkle Tree.")
pub const TWEEDLE_DEE_PHANTOM_MERKLE_ROOT: TweedleFr = field_new!(
    TweedleFr,
    BigInteger256([
        3037328570344736808,
        7713672412249441083,
        9580110273020338221,
        816345156093001009
    ])
);

pub const TWEEDLE_DEE_MHT_POSEIDON_PARAMETERS: FieldBasedMerkleTreePrecomputedZeroConstants<
    'static,
    TweedleFrPoseidonHash,
> = FieldBasedMerkleTreePrecomputedZeroConstants {
    nodes: &[
        field_new!(TweedleFr, BigInteger256([0, 0, 0, 0])),
        field_new!(
            TweedleFr,
            BigInteger256([
                6846511105464766538,
                15768966942874777847,
                16388715769057780159,
                3605183713290623682
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                9222333104797974540,
                2988232145305907562,
                16209565825461578695,
                3126222989224963312
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                14417722119675398228,
                4278309788110750045,
                4043558729910385260,
                1385476922717649264
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                12405703624929027638,
                17686987702583161392,
                14818595643264832920,
                1298091960176016512
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                11220962421518165700,
                13583264328995303902,
                3004999268640918219,
                1836274747239137718
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                9143319041823283548,
                10625485209067256567,
                3101953621268315084,
                2784075795165174292
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                7856896111698860209,
                12274291086498461139,
                12254863429498589520,
                1157091047461829565
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                12781995487914321830,
                8909319778259688775,
                9744152270041314391,
                120486371160620658
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                16887679367162088832,
                3897372093715987906,
                1678885614393805654,
                3178520008167028611
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                14480369559671180193,
                7603632368518351521,
                15818547859043187272,
                2573528473272863098
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                18229837180404920553,
                15631813461238913726,
                1585236863667313179,
                1429740507100771895
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                7597046153975937193,
                7149755588864973802,
                12498868822806140200,
                3677394355095085380
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                9816128661650764975,
                9363556778543621599,
                16728536662023759852,
                4081025247655585206
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                395144311372928536,
                4050112514364960005,
                10671415354094218204,
                539401662144033117
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                7786395140808940918,
                9820478310629028028,
                637683664312732797,
                1223633973447551913
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                16337123754752253291,
                9321181451798296893,
                14228134238565922823,
                2055512899004156098
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                7184362807226572198,
                881197418286383995,
                12302845140236449565,
                3264750990323283693
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                16713308676519013057,
                6067655313659401642,
                16860473533228667971,
                1782908225150503344
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                16396585298471795251,
                13694615398856184131,
                3261625920828783531,
                3104500293296803657
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                3841871426667137165,
                4714827166841758205,
                10093050349190081639,
                3118264349018454210
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                6711494586375705172,
                1999689148281285955,
                10635255840256041308,
                1563702960800716865
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                6209055503806423570,
                14160757492870172822,
                2585173026931081813,
                45434336531305746
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                14967778432602946057,
                8401376301444944267,
                163291714506409423,
                1721295824263050447
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                15893452450907463361,
                15302997450958726322,
                2883783336883997986,
                192767118370846437
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                1635644627404959933,
                1389608663507338955,
                2055512431619577463,
                2018832189365185308
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                12430793728438799303,
                17460069178822475215,
                5065750370745188226,
                2955160172133821296
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                6809777662530062064,
                13028945633714685770,
                7708962218143035851,
                3805976379826616327
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                16227885774693608785,
                15582103355717546603,
                13023980185169825094,
                2292196543775082438
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                7124941499245031893,
                2481564590150069704,
                15682136442122527111,
                2927348525097878888
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                1642502257843298186,
                12569295554878034139,
                2314986831057740561,
                2218381406452172661
            ])
        ),
        field_new!(
            TweedleFr,
            BigInteger256([
                11027003307842316033,
                3331825336977729667,
                9686421639239063727,
                1895668150828362785
            ])
        ),
    ],
    merkle_arity: 2,
};

#[cfg(test)]
mod test {

    use super::{TWEEDLE_DEE_MHT_POSEIDON_PARAMETERS, TWEEDLE_DEE_PHANTOM_MERKLE_ROOT};
    use crate::{
        crh::TweedleFrPoseidonHash,
        merkle_tree::field_based_mht::parameters::{
            generate_mht_empty_nodes, generate_phantom_merkle_root_from_magic_string,
        },
        FieldBasedMerkleTreePrecomputedZeroConstants,
    };
    use algebra::{fields::tweedle::Fr, Field};

    #[ignore]
    #[test]
    fn test_generate_tweedle_fr_phantom_merkle_root() {
        let expected_root = generate_phantom_merkle_root_from_magic_string::<
            Fr,
            TweedleFrPoseidonHash,
        >(
            "This represents an empty Merkle Root for a TweedleDeePoseidonHash based Merkle Tree.",
        );
        assert_eq!(expected_root, TWEEDLE_DEE_PHANTOM_MERKLE_ROOT);
    }

    #[ignore]
    #[test]
    fn test_generate_binary_tweedle_fr_mht_empty_nodes() {
        let merkle_arity = 2;
        let max_height = 32;

        let empty_nodes = generate_mht_empty_nodes::<Fr, TweedleFrPoseidonHash>(
            merkle_arity,
            max_height,
            Fr::zero(),
        );
        assert_eq!(empty_nodes.len(), max_height);

        let params = FieldBasedMerkleTreePrecomputedZeroConstants::<TweedleFrPoseidonHash> {
            nodes: empty_nodes.as_slice(),
            merkle_arity,
        };
        assert_eq!(params, TWEEDLE_DEE_MHT_POSEIDON_PARAMETERS)
    }
}
