use crate::groups::curves::twisted_edwards::AffineGadget;
use algebra::curves::edwards_sw6::EdwardsParameters;
use algebra::fields::edwards_sw6::fq::Fq;

use crate::fields::edwards_sw6::FqGadget;

pub type EdwardsSWGadget = AffineGadget<EdwardsParameters, Fq, FqGadget>;

#[cfg(test)]
mod test {
    use super::EdwardsSWGadget as EdwardsG;
    use crate::{
        groups::curves::twisted_edwards::test::{edwards_constraint_costs, edwards_test},
        test_constraint_system::TestConstraintSystem,
    };
    use algebra::curves::edwards_sw6::EdwardsParameters;
    use algebra::fields::edwards_sw6::fq::Fq;

    #[test]
    fn edwards_constraint_costs_test() {
        let mut cs = TestConstraintSystem::<Fq>::new();
        edwards_constraint_costs::<_, EdwardsParameters, EdwardsG, _>(&mut cs);
        assert!(cs.is_satisfied());
    }

    #[test]
    fn edwards_sw6_gadget_test() {
        let mut cs = TestConstraintSystem::<Fq>::new();
        edwards_test::<_, EdwardsParameters, EdwardsG, _>(&mut cs);
        assert!(cs.is_satisfied());
    }
}
