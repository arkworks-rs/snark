use crate::groups::curves::twisted_edwards::AffineGadget;
use algebra::curves::jubjub::JubJubParameters;
use algebra::fields::jubjub::fq::Fq;

use crate::fields::jubjub::FqGadget;

pub type JubJubGadget = AffineGadget<JubJubParameters, Fq, FqGadget>;

#[cfg(test)]
mod test {
    use super::JubJubGadget as EdwardsG;
    use crate::{
        groups::curves::twisted_edwards::test::{edwards_constraint_costs, edwards_test},
        test_constraint_system::TestConstraintSystem,
    };
    use algebra::fields::jubjub::fq::Fq;
    use algebra::curves::jubjub::JubJubParameters as EdwardsParameters;

    #[test]
    fn edwards_constraint_costs_test() {
        let mut cs = TestConstraintSystem::<Fq>::new();
        edwards_constraint_costs::<_, EdwardsParameters, EdwardsG, _>(&mut cs);
        assert!(cs.is_satisfied());
    }

    #[test]
    fn jubjub_gadget_test() {
        let mut cs = TestConstraintSystem::<Fq>::new();
        edwards_test::<_, EdwardsParameters, EdwardsG, _>(&mut cs);
        assert!(cs.is_satisfied());
    }
}
