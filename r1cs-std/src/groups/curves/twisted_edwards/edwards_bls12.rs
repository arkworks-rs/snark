use crate::groups::curves::twisted_edwards::AffineGadget;
use algebra::curves::edwards_bls12::EdwardsParameters;
use algebra::fields::edwards_bls12::fq::Fq;

use crate::fields::edwards_bls12::FqGadget;

pub type EdwardsBlsGadget = AffineGadget<EdwardsParameters, Fq, FqGadget>;

#[cfg(test)]
mod test {
    use super::EdwardsBlsGadget as EdwardsG;
    use crate::{
        groups::curves::twisted_edwards::test::{edwards_constraint_costs, edwards_test},
        test_constraint_system::TestConstraintSystem,
    };
    use algebra::fields::edwards_bls12::fq::Fq;
    use algebra::curves::edwards_bls12::EdwardsParameters;

    #[test]
    fn edwards_constraint_costs_test() {
        let mut cs = TestConstraintSystem::<Fq>::new();
        edwards_constraint_costs::<_, EdwardsParameters, EdwardsG, _>(&mut cs);
        assert!(cs.is_satisfied());
    }

    #[test]
    fn edwards_bls12_gadget_test() {
        let mut cs = TestConstraintSystem::<Fq>::new();
        edwards_test::<_, EdwardsParameters, EdwardsG, _>(&mut cs);
        assert!(cs.is_satisfied());
    }
}
