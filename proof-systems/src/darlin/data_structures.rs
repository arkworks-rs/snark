use algebra::{AffineCurve, ToConstraintField, ToBits, ProjectiveCurve, UniformRand, serialize::*};
use marlin::Proof as MarlinProof;
use crate::darlin::accumulators::dlog::DLogItem;
use poly_commit::ipa_pc::{
    SuccinctCheckPolynomial, InnerProductArgPC,
    CommitterKey as DLogCommitterKey, Commitment,
};
use digest::Digest;
use rand::RngCore;

// Maybe later this will include more element as we will deferr algebraic checks over G1::BaseField
#[derive(Default, Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct FinalDarlinDeferredData<G1: AffineCurve, G2: AffineCurve> {
    pub(crate) previous_acc:       DLogItem<G2>,
    pub(crate) pre_previous_acc:   DLogItem<G1>,
}

impl<G1, G2> FinalDarlinDeferredData<G1, G2>
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    pub fn generate_random<R: RngCore, D: Digest>(
        rng: &mut R,
        committer_key_g1: &DLogCommitterKey<G1>,
        committer_key_g2: &DLogCommitterKey<G2>
    ) -> Self
    {
        // Generate valid accumulator over G1 starting from random xi_s
        let log_key_len_g1 = algebra::log2(committer_key_g1.comm_key.len());
        let random_xi_s_g1 = SuccinctCheckPolynomial::<G1::ScalarField>(vec![u128::rand(rng).into(); log_key_len_g1 as usize]);
        let g_final_g1 = InnerProductArgPC::<G1, D>::cm_commit(
            committer_key_g1.comm_key.as_slice(),
            random_xi_s_g1.compute_coeffs().as_slice(),
            None,
            None,
        );

        let acc_g1 = DLogItem::<G1> {
            g_final: Commitment::<G1> {comm: vec![g_final_g1.into_affine()], shifted_comm: None },
            xi_s: random_xi_s_g1
        };

        // Generate valid accumulator over G2 starting from random xi_s
        let log_key_len_g2 = algebra::log2(committer_key_g2.comm_key.len());
        let random_xi_s_g2 = SuccinctCheckPolynomial::<G2::ScalarField>(vec![u128::rand(rng).into(); log_key_len_g2 as usize]);

        let g_final_g2 = InnerProductArgPC::<G2, D>::cm_commit(
            committer_key_g2.comm_key.as_slice(),
            random_xi_s_g2.compute_coeffs().as_slice(),
            None,
            None,
        );

        let acc_g2 = DLogItem::<G2> {
            g_final: Commitment::<G2> {comm: vec![g_final_g2.into_affine()], shifted_comm: None },
            xi_s: random_xi_s_g2
        };

        // Return accumulators in deferred struct
        Self {
            previous_acc: acc_g2,
            pre_previous_acc: acc_g1
        }
    }
}

impl<G1, G2> ToConstraintField<G1::ScalarField> for FinalDarlinDeferredData<G1, G2>
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    fn to_field_elements(&self) -> Result<Vec<G1::ScalarField>, Box<dyn std::error::Error>> {
        let mut fes = Vec::new();

        // Convert previous_acc into G1::ScalarField field elements (matching field)

        let g_final_g2 = self.previous_acc.g_final.comm.clone();
        for c in g_final_g2.into_iter() {
            fes.append(&mut c.to_field_elements()?);
        }

        // Convert xi_s to G1::ScalarField elements, the field is not matching:
        // Since the xi_s will be only 128 bits long, we optimize the packing by
        // concatenating all the bits in a single array and (safely) read as
        // many native field elements as possible.
        let mut xi_s_bits = Vec::new();
        for fe in self.previous_acc.xi_s.0.clone().into_iter() {
            xi_s_bits.extend_from_slice(&fe.write_bits()[..128]);
        }
        fes.append(&mut xi_s_bits.to_field_elements()?);

        // Convert g_final_g1 to G1::ScalarField elements, the field is not matching:
        // serialize them all to bits and pack them safely into native field elements

        let g_final_g1 = self.pre_previous_acc.g_final.comm.clone();
        let mut g_final_g1_bits = Vec::new();
        for c in g_final_g1 {
            let c_fes = c.to_field_elements()?;
            for fe in c_fes {
                g_final_g1_bits.append(&mut fe.write_bits());
            }
        }
        fes.append(&mut g_final_g1_bits.to_field_elements()?);

        // Convert pre_previous_acc into G1::ScalarField field elements.
        // Even if the field is matching, we are still interested in
        // optimizing the final number of G1::ScalarField elements by
        // leveraging the fact that they are only 128 bits long.
        let mut xi_s_bits = Vec::new();
        for fe in self.pre_previous_acc.xi_s.0.clone().into_iter() {
            xi_s_bits.extend_from_slice(&fe.write_bits()[..128]);
        }
        fes.append(&mut xi_s_bits.to_field_elements()?);

        Ok(fes)
    }
}

#[derive(Derivative)]
#[derivative(Clone(bound = ""), Debug(bound = ""), Eq(bound = ""), PartialEq(bound = ""))]
#[derive(CanonicalSerialize, CanonicalDeserialize)]
/// FinalDarlinProof with two deferred DLOG accumulators.
pub struct FinalDarlinProof<G1: AffineCurve, G2: AffineCurve, D: Digest> {
    /// Full Marlin proof without deferred arithmetics in G1.
    pub proof:       MarlinProof<G1::ScalarField, InnerProductArgPC<G1, D>>,
    /// Deferred accumulators
    pub deferred:    FinalDarlinDeferredData<G1, G2>,
}