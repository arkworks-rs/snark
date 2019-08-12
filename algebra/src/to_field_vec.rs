use crate::{
    curves::{
        models::{SWModelParameters, TEModelParameters},
        short_weierstrass_jacobian::{GroupAffine as SWAffine, GroupProjective as SWProjective},
        twisted_edwards_extended::{GroupAffine as TEAffine, GroupProjective as TEProjective},
        ProjectiveCurve,
    },
    Fp2, Fp2Parameters, FpParameters, Field, PrimeField,
};

type Error = Box<dyn std::error::Error>;

/// Types that can be converted to a vector of `F` elements. Useful for specifying
/// how public inputs to a constraint system should be represented inside
/// that constraint system.
pub trait ToConstraintField<F: Field> {
    fn to_field_elements(&self) -> Result<Vec<F>, Error>;
}

impl<F: PrimeField> ToConstraintField<F> for F {
    fn to_field_elements(&self) -> Result<Vec<F>, Error> {
        Ok(vec![*self])
    }
}

// Impl for base field
impl<F: Field> ToConstraintField<F> for [F] {
    #[inline]
    fn to_field_elements(&self) -> Result<Vec<F>, Error> {
        Ok(self.to_vec())
    }
}

impl<ConstraintF: Field> ToConstraintField<ConstraintF> for () {
    #[inline]
    fn to_field_elements(&self) -> Result<Vec<ConstraintF>, Error> {
        Ok(Vec::new())
    }
}

// Impl for Fp2<ConstraintF>
impl<P: Fp2Parameters> ToConstraintField<P::Fp> for Fp2<P>
{
    #[inline]
    fn to_field_elements(&self) -> Result<Vec<P::Fp>, Error> {
        let mut c0 = self.c0.to_field_elements()?;
        let c1 = self.c1.to_field_elements()?;
        c0.extend_from_slice(&c1);
        Ok(c0)
    }
}

impl<M: TEModelParameters, ConstraintF: Field> ToConstraintField<ConstraintF> for TEAffine<M>
where
    M::BaseField: ToConstraintField<ConstraintF>,
{
    #[inline]
    fn to_field_elements(&self) -> Result<Vec<ConstraintF>, Error> {
        let mut x_fe = self.x.to_field_elements()?;
        let y_fe = self.y.to_field_elements()?;
        x_fe.extend_from_slice(&y_fe);
        Ok(x_fe)
    }
}

impl<M: TEModelParameters, ConstraintF: Field> ToConstraintField<ConstraintF> for TEProjective<M>
where
    M::BaseField: ToConstraintField<ConstraintF>,
{
    #[inline]
    fn to_field_elements(&self) -> Result<Vec<ConstraintF>, Error> {
        let affine = self.into_affine();
        let mut x_fe = affine.x.to_field_elements()?;
        let y_fe = affine.y.to_field_elements()?;
        x_fe.extend_from_slice(&y_fe);
        Ok(x_fe)
    }
}

impl<M: SWModelParameters, ConstraintF: Field> ToConstraintField<ConstraintF> for SWAffine<M>
where
    M::BaseField: ToConstraintField<ConstraintF>,
{
    #[inline]
    fn to_field_elements(&self) -> Result<Vec<ConstraintF>, Error> {
        let mut x_fe = self.x.to_field_elements()?;
        let y_fe = self.y.to_field_elements()?;
        x_fe.extend_from_slice(&y_fe);
        Ok(x_fe)
    }
}

impl<M: SWModelParameters, ConstraintF: Field> ToConstraintField<ConstraintF> for SWProjective<M>
where
    M::BaseField: ToConstraintField<ConstraintF>,
{
    #[inline]
    fn to_field_elements(&self) -> Result<Vec<ConstraintF>, Error> {
        let affine = self.into_affine();
        let mut x_fe = affine.x.to_field_elements()?;
        let y_fe = affine.y.to_field_elements()?;
        x_fe.extend_from_slice(&y_fe);
        Ok(x_fe)
    }
}

impl<ConstraintF: PrimeField> ToConstraintField<ConstraintF> for [u8] {
    #[inline]
    fn to_field_elements(&self) -> Result<Vec<ConstraintF>, Error> {
        let max_size = <ConstraintF as PrimeField>::Params::CAPACITY / 8;
        let max_size = max_size as usize;
        let fes = self
            .chunks(max_size)
            .map(|chunk| {
                let mut chunk = chunk.to_vec();
                let len = chunk.len();
                for _ in len..(max_size + 1) {
                    chunk.push(0u8);
                }
                ConstraintF::read(chunk.as_slice())
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(fes)
    }
}

impl<ConstraintF: PrimeField> ToConstraintField<ConstraintF> for [u8; 32] {
    #[inline]
    fn to_field_elements(&self) -> Result<Vec<ConstraintF>, Error> {
        self.as_ref().to_field_elements()
    }
}
