pub trait ToBits {
    /// Serialize `self` into a bit vector.
    fn write_bits(&self) -> Vec<bool>;
}

pub trait FromBits: Sized {
    /// Reads `self` from `bits`
    fn read_bits(bits: Vec<bool>) -> Self;
}

pub trait ToCompressedBits {
    fn compress(&self) -> Vec<bool>;
}

pub trait FromCompressedBits: Sized  {
    fn decompress(compressed: Vec<bool>) -> Option<Self>;
}