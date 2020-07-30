use crate::{
    biginteger::BigInteger,
    curves::{
        models::{SWModelParameters, TEModelParameters},
        short_weierstrass_jacobian::{GroupAffine as SWAffine, GroupProjective as SWProjective},
        twisted_edwards_extended::{GroupAffine as TEAffine, GroupProjective as TEProjective},
    },
    Box, Field, FpParameters, PrimeField, Vec,
};

type Error = Box<dyn crate::Error>;

/// Types that can be converted to a vector of `F` elements. Useful for
/// specifying how public inputs to a constraint system should be represented
/// inside that constraint system.
pub trait ToConstraintField<F: Field> {
    fn to_field_elements(&self) -> Result<Vec<F>, Error>;
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

// Impl for boolean
impl<F: Field> ToConstraintField<F> for bool {
    #[inline]
    fn to_field_elements(&self) -> Result<Vec<F>, Error> {
        if *self {
            Ok(vec![F::one()])
        } else {
            Ok(vec![F::zero()])
        }
    }
}

impl<M: TEModelParameters> ToConstraintField<<M::BaseField as Field>::BaseRepresentationField>
    for TEAffine<M>
{
    #[inline]
    fn to_field_elements(
        &self,
    ) -> Result<Vec<<M::BaseField as Field>::BaseRepresentationField>, Error> {
        let mut x_fe = self.x.to_field_elements()?;
        let y_fe = self.y.to_field_elements()?;
        x_fe.extend_from_slice(&y_fe);
        Ok(x_fe)
    }
}

impl<M: TEModelParameters> ToConstraintField<<M::BaseField as Field>::BaseRepresentationField>
    for TEProjective<M>
{
    #[inline]
    fn to_field_elements(
        &self,
    ) -> Result<Vec<<M::BaseField as Field>::BaseRepresentationField>, Error> {
        TEAffine::from(*self).to_field_elements()
    }
}

impl<M: SWModelParameters> ToConstraintField<<M::BaseField as Field>::BaseRepresentationField>
    for SWAffine<M>
{
    #[inline]
    fn to_field_elements(
        &self,
    ) -> Result<Vec<<M::BaseField as Field>::BaseRepresentationField>, Error> {
        let mut x_fe = self.x.to_field_elements()?;
        let y_fe = self.y.to_field_elements()?;
        let infinity_fe = self.infinity.to_field_elements()?;
        x_fe.extend_from_slice(&y_fe);
        x_fe.extend_from_slice(&infinity_fe);
        Ok(x_fe)
    }
}

impl<M: SWModelParameters> ToConstraintField<<M::BaseField as Field>::BaseRepresentationField>
    for SWProjective<M>
{
    #[inline]
    fn to_field_elements(
        &self,
    ) -> Result<Vec<<M::BaseField as Field>::BaseRepresentationField>, Error> {
        SWAffine::from(*self).to_field_elements()
    }
}

impl<ConstraintF: PrimeField> ToConstraintField<ConstraintF> for [u8] {
    #[inline]
    fn to_field_elements(&self) -> Result<Vec<ConstraintF>, Error> {
        let max_size = <ConstraintF as PrimeField>::Params::CAPACITY / 8;
        let max_size = max_size as usize;
        let bigint_size = <ConstraintF as PrimeField>::BigInt::NUM_LIMBS * 8;
        let fes = self
            .chunks(max_size)
            .map(|chunk| {
                let mut chunk = chunk.to_vec();
                let len = chunk.len();
                for _ in len..bigint_size {
                    chunk.push(0u8);
                }
                ConstraintF::read(chunk.as_slice())
            })
            .collect::<Result<Vec<_>, _>>()
            .map_err(crate::SerializationError::from)
            .map_err(|e| Box::new(e))?;
        Ok(fes)
    }
}

impl<ConstraintF: PrimeField> ToConstraintField<ConstraintF> for [u8; 32] {
    #[inline]
    fn to_field_elements(&self) -> Result<Vec<ConstraintF>, Error> {
        self.as_ref().to_field_elements()
    }
}

impl<ConstraintF: PrimeField> ToConstraintField<ConstraintF> for Vec<u8> {
    #[inline]
    fn to_field_elements(&self) -> Result<Vec<ConstraintF>, Error> {
        let max_size = <ConstraintF as PrimeField>::Params::CAPACITY / 8;
        let max_size = max_size as usize;
        let bigint_size = <ConstraintF as PrimeField>::BigInt::NUM_LIMBS * 8;
        let fes = self
            .chunks(max_size)
            .map(|chunk| {
                let mut chunk = chunk.to_vec();
                let len = chunk.len();
                for _ in len..bigint_size {
                    chunk.push(0u8);
                }
                ConstraintF::read(chunk.as_slice())
            })
            .collect::<Result<Vec<_>, _>>()
            .map_err(crate::SerializationError::from)
            .map_err(|e| Box::new(e))?;
        Ok(fes)
    }
}
