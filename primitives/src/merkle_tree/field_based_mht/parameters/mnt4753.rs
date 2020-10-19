use algebra::{
    fields::mnt4753::Fr as MNT4753Fr,
    biginteger::BigInteger768,
    field_new,
};

// PoseidonHash("This represents an empty Merkle Root for a MNT4753PoseidonHash based Merkle Tree.") padded with 0s
pub const MNT4753_PHANTOM_MERKLE_ROOT: MNT4753Fr =
    field_new!(MNT4753Fr, BigInteger768([
        8534937690304963668,
        5486292534803213323,
        1720870611961422927,
        11405914840660719672,
        7162329517212056783,
        11658292353137306079,
        17490588101047840223,
        12735752881395833110,
        11735157047413601083,
        6658060155531600932,
        1470933043432945054,
        312822709740712
    ]));

#[cfg(test)]
mod test {

    #[ignore]
    #[test]
    fn generate_mnt4753_phantom_merkle_root(){
        use algebra::{
            fields::mnt4753::Fr,
            FromBytes, PrimeField, FpParameters
        };
        use crate::crh::{
            MNT4PoseidonHash, FieldBasedHash,
        };
        use super::MNT4753_PHANTOM_MERKLE_ROOT;

        let field_size_in_bytes = (Fr::size_in_bits() + (<Fr as PrimeField>::Params::REPR_SHAVE_BITS as usize))/8;
        let magic_string = b"This represents an empty Merkle Root for a MNT4753PoseidonHash based Merkle Tree.";

        let mut hash_input = Vec::new();
        hash_input.extend_from_slice(magic_string);
        for _ in magic_string.len()..field_size_in_bytes { hash_input.push(0u8) }
        let hash_input_f = Fr::read(hash_input.as_slice()).unwrap();

        let hash = MNT4PoseidonHash::evaluate(&[hash_input_f]).unwrap();
        assert_eq!(hash, MNT4753_PHANTOM_MERKLE_ROOT);
    }
}


