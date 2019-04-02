use crate::groups::curves::twisted_edwards::AffineGadget;
use algebra::curves::{bls12_381::Bls12_381, jubjub::JubJubParameters};

use crate::fields::jubjub::FqGadget;

pub type JubJubGadget = AffineGadget<JubJubParameters, Bls12_381, FqGadget>;

#[cfg(test)]
mod test {
    use super::JubJubGadget as EdwardsG;
    use crate::{
        groups::curves::twisted_edwards::test::{edwards_constraint_costs, edwards_test},
        test_constraint_system::TestConstraintSystem,
    };
    use algebra::curves::{bls12_381::Bls12_381, jubjub::JubJubParameters as EdwardsParameters};

    #[test]
    fn edwards_constraint_costs_test() {
        let mut cs = TestConstraintSystem::<Bls12_381>::new();
        edwards_constraint_costs::<_, EdwardsParameters, EdwardsG, _>(&mut cs);
        assert!(cs.is_satisfied());
    }

    #[test]
    fn edwards_sw6_gadget_test() {
        let mut cs = TestConstraintSystem::<Bls12_381>::new();
        edwards_test::<_, EdwardsParameters, EdwardsG, _>(&mut cs);
        assert!(cs.is_satisfied());
    }
}
