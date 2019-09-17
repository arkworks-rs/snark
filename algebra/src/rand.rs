use rand::{Rng, distributions::{Distribution, Standard}};

pub trait Rand: Sized {
    fn rand<R>(rng: &mut R) -> Self
    where
        R: Rng,
        Standard: Distribution<Self>;
}
