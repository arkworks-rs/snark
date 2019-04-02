use crate::groups::{curves::short_weierstrass::AffineGadget, GroupGadget};
use algebra::{
    curves::bls12::{Bls12Parameters, G1Prepared, TwistType},
    fields::Field,
    BitIterator, PairingEngine, ProjectiveCurve,
};
use snark::{ConstraintSystem, SynthesisError};

use crate::{
    fields::{fp::FpGadget, fp2::Fp2Gadget, FieldGadget},
    uint8::UInt8,
    utils::{NEqGadget, ToBytesGadget},
};

use std::fmt::Debug;

pub mod bls12_377;

pub type G1Gadget<P, E> = AffineGadget<<P as Bls12Parameters>::G1Parameters, E, FpGadget<E>>;
pub type G2Gadget<P, E> = AffineGadget<<P as Bls12Parameters>::G2Parameters, E, Fp2G<P, E>>;

#[derive(Derivative)]
#[derivative(
    Clone(bound = "G1Gadget<P, E>: Clone"),
    Debug(bound = "G1Gadget<P, E>: Debug")
)]
pub struct G1PreparedGadget<P: Bls12Parameters<Fp = E::Fr>, E: PairingEngine>(pub G1Gadget<P, E>);

impl<P: Bls12Parameters<Fp = E::Fr>, E: PairingEngine> G1PreparedGadget<P, E> {
    pub fn get_value(&self) -> Option<G1Prepared<P>> {
        Some(G1Prepared::from_affine(
            self.0.get_value().unwrap().into_affine(),
        ))
    }

    pub fn from_affine<CS: ConstraintSystem<E>>(
        _cs: CS,
        q: &G1Gadget<P, E>,
    ) -> Result<Self, SynthesisError> {
        Ok(G1PreparedGadget(q.clone()))
    }
}

impl<P: Bls12Parameters<Fp = E::Fr>, E: PairingEngine> ToBytesGadget<E> for G1PreparedGadget<P, E> {
    #[inline]
    fn to_bytes<CS: ConstraintSystem<E>>(&self, mut cs: CS) -> Result<Vec<UInt8>, SynthesisError> {
        self.0.to_bytes(&mut cs.ns(|| "g_alpha to bytes"))
    }

    fn to_bytes_strict<CS: ConstraintSystem<E>>(
        &self,
        cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        self.to_bytes(cs)
    }
}

type Fp2G<P, E> = Fp2Gadget<<P as Bls12Parameters>::Fp2Params, E>;
type LCoeff<P, E> = (Fp2G<P, E>, Fp2G<P, E>);
#[derive(Derivative)]
#[derivative(
    Clone(bound = "Fp2Gadget<P::Fp2Params, E>: Clone"),
    Debug(bound = "Fp2Gadget<P::Fp2Params, E>: Debug")
)]
pub struct G2PreparedGadget<P: Bls12Parameters<Fp = E::Fr>, E: PairingEngine> {
    pub ell_coeffs: Vec<LCoeff<P, E>>,
}

impl<P: Bls12Parameters<Fp = E::Fr>, E: PairingEngine> ToBytesGadget<E> for G2PreparedGadget<P, E> {
    #[inline]
    fn to_bytes<CS: ConstraintSystem<E>>(&self, mut cs: CS) -> Result<Vec<UInt8>, SynthesisError> {
        let mut bytes = Vec::new();
        for (i, coeffs) in self.ell_coeffs.iter().enumerate() {
            let mut cs = cs.ns(|| format!("Iteration {}", i));
            bytes.extend_from_slice(&coeffs.0.to_bytes(&mut cs.ns(|| "c0"))?);
            bytes.extend_from_slice(&coeffs.1.to_bytes(&mut cs.ns(|| "c1"))?);
        }
        Ok(bytes)
    }

    fn to_bytes_strict<CS: ConstraintSystem<E>>(
        &self,
        cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        self.to_bytes(cs)
    }
}

impl<P: Bls12Parameters<Fp = E::Fr>, E: PairingEngine> G2PreparedGadget<P, E> {
    pub fn from_affine<CS: ConstraintSystem<E>>(
        mut cs: CS,
        q: &G2Gadget<P, E>,
    ) -> Result<Self, SynthesisError> {
        let two_inv = P::Fp::one().double().inverse().unwrap();
        let zero = G2Gadget::<P, E>::zero(cs.ns(|| "zero"))?;
        q.enforce_not_equal(cs.ns(|| "enforce not zero"), &zero)?;
        let mut ell_coeffs = vec![];
        let mut r = q.clone();

        for (j, i) in BitIterator::new(P::X).skip(1).enumerate() {
            let mut cs = cs.ns(|| format!("Iteration {}", j));
            ell_coeffs.push(Self::double(cs.ns(|| "double"), &mut r, &two_inv)?);

            if i {
                ell_coeffs.push(Self::add(cs.ns(|| "add"), &mut r, &q)?);
            }
        }

        Ok(Self { ell_coeffs })
    }

    fn double<CS: ConstraintSystem<E>>(
        mut cs: CS,
        r: &mut G2Gadget<P, E>,
        two_inv: &P::Fp,
    ) -> Result<LCoeff<P, E>, SynthesisError> {
        let a = r.y.inverse(cs.ns(|| "Inverse"))?;
        let mut b = r.x.square(cs.ns(|| "square x"))?;
        let b_tmp = b.clone();
        b.mul_by_fp_constant_in_place(cs.ns(|| "mul by two_inv"), two_inv)?;
        b.add_in_place(cs.ns(|| "compute b"), &b_tmp)?;

        let c = a.mul(cs.ns(|| "compute c"), &b)?;
        let d = r.x.double(cs.ns(|| "compute d"))?;
        let x3 = c.square(cs.ns(|| "c^2"))?.sub(cs.ns(|| "sub d"), &d)?;
        let e = c
            .mul(cs.ns(|| "c*r.x"), &r.x)?
            .sub(cs.ns(|| "sub r.y"), &r.y)?;
        let c_x3 = c.mul(cs.ns(|| "c*x_3"), &x3)?;
        let y3 = e.sub(cs.ns(|| "e = c * x3"), &c_x3)?;
        let mut f = c;
        f.negate_in_place(cs.ns(|| "c = -c"))?;
        r.x = x3;
        r.y = y3;
        match P::TWIST_TYPE {
            TwistType::M => Ok((e, f)),
            TwistType::D => Ok((f, e)),
        }
    }

    fn add<CS: ConstraintSystem<E>>(
        mut cs: CS,
        r: &mut G2Gadget<P, E>,
        q: &G2Gadget<P, E>,
    ) -> Result<LCoeff<P, E>, SynthesisError> {
        let a =
            q.x.sub(cs.ns(|| "q.x - r.x"), &r.x)?
                .inverse(cs.ns(|| "calc a"))?;
        let b = q.y.sub(cs.ns(|| "q.y - r.y"), &r.y)?;
        let c = a.mul(cs.ns(|| "compute c"), &b)?;
        let d = r.x.add(cs.ns(|| "r.x + q.x"), &q.x)?;
        let x3 = c.square(cs.ns(|| "c^2"))?.sub(cs.ns(|| "sub d"), &d)?;

        let e =
            r.x.sub(cs.ns(|| "r.x - x3"), &x3)?
                .mul(cs.ns(|| "c * (r.x - x3)"), &c)?;
        let y3 = e.sub(cs.ns(|| "calc y3"), &r.y)?;
        let g = c
            .mul(cs.ns(|| "c*r.x"), &r.x)?
            .sub(cs.ns(|| "calc g"), &r.y)?;
        let mut f = c;
        f.negate_in_place(cs.ns(|| "c = -c"))?;
        r.x = x3;
        r.y = y3;
        match P::TWIST_TYPE {
            TwistType::M => Ok((g, f)),
            TwistType::D => Ok((f, g)),
        }
    }
}
