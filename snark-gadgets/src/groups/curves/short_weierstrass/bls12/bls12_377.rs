use crate::groups::bls12::{
    G1Gadget as Bls12G1Gadget, G1PreparedGadget as Bls12G1PreparedGadget,
    G2Gadget as Bls12G2Gadget, G2PreparedGadget as Bls12G2PreparedGadget,
};
use algebra::curves::{bls12_377::Bls12_377Parameters, sw6::SW6};

pub type G1Gadget = Bls12G1Gadget<Bls12_377Parameters, SW6>;
pub type G2Gadget = Bls12G2Gadget<Bls12_377Parameters, SW6>;

pub type G1PreparedGadget = Bls12G1PreparedGadget<Bls12_377Parameters, SW6>;
pub type G2PreparedGadget = Bls12G2PreparedGadget<Bls12_377Parameters, SW6>;

#[cfg(test)]
mod test {
    use crate::fields::FieldGadget;
    use rand;

    use super::{G1Gadget, G2Gadget};
    use crate::{
        boolean::Boolean,
        groups::GroupGadget,
        test_constraint_system::TestConstraintSystem,
        utils::{AllocGadget, CondSelectGadget, EqGadget},
    };
    use algebra::{
        curves::{
            bls12_377::{G1Projective as G1, G2Projective as G2},
            sw6::SW6,
        },
        fields::bls12_377::Fr,
        AffineCurve, BitIterator, PrimeField, ProjectiveCurve,
    };
    use snark::ConstraintSystem;

    #[test]
    fn bls12_g1_constraint_costs() {
        use crate::boolean::AllocatedBit;

        let mut cs = TestConstraintSystem::<SW6>::new();

        let bit = AllocatedBit::alloc(&mut cs.ns(|| "bool"), || Ok(true))
            .unwrap()
            .into();

        let a: G1 = rand::random();
        let b: G1 = rand::random();
        let gadget_a = G1Gadget::alloc(&mut cs.ns(|| "a"), || Ok(a)).unwrap();
        let gadget_b = G1Gadget::alloc(&mut cs.ns(|| "b"), || Ok(b)).unwrap();
        let alloc_cost = cs.num_constraints();
        let _ = G1Gadget::conditionally_select(
            &mut cs.ns(|| "cond_select"),
            &bit,
            &gadget_a,
            &gadget_b,
        )
        .unwrap();
        let cond_select_cost = cs.num_constraints() - alloc_cost;

        let _ = gadget_a.add(&mut cs.ns(|| "ab"), &gadget_b).unwrap();
        let add_cost = cs.num_constraints() - cond_select_cost - alloc_cost;

        assert!(cs.is_satisfied());
        assert_eq!(
            cond_select_cost,
            <G1Gadget as CondSelectGadget<SW6>>::cost()
        );
        assert_eq!(add_cost, G1Gadget::cost_of_add());
    }

    #[test]
    fn bls12_g2_constraint_costs() {
        use crate::boolean::AllocatedBit;

        let mut cs = TestConstraintSystem::<SW6>::new();

        let bit = AllocatedBit::alloc(&mut cs.ns(|| "bool"), || Ok(true))
            .unwrap()
            .into();

        let a: G2 = rand::random();
        let b: G2 = rand::random();
        let gadget_a = G2Gadget::alloc(&mut cs.ns(|| "a"), || Ok(a)).unwrap();
        let gadget_b = G2Gadget::alloc(&mut cs.ns(|| "b"), || Ok(b)).unwrap();
        let alloc_cost = cs.num_constraints();
        let _ = G2Gadget::conditionally_select(
            &mut cs.ns(|| "cond_select"),
            &bit,
            &gadget_a,
            &gadget_b,
        )
        .unwrap();
        let cond_select_cost = cs.num_constraints() - alloc_cost;

        let _ = gadget_a.add(&mut cs.ns(|| "ab"), &gadget_b).unwrap();
        let add_cost = cs.num_constraints() - cond_select_cost - alloc_cost;

        assert!(cs.is_satisfied());
        assert_eq!(
            cond_select_cost,
            <G2Gadget as CondSelectGadget<SW6>>::cost()
        );
        assert_eq!(add_cost, G2Gadget::cost_of_add());
    }

    #[test]
    fn bls12_g1_gadget_test() {
        use rand::{Rand, SeedableRng, XorShiftRng};
        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        let mut cs = TestConstraintSystem::<SW6>::new();

        let a = G1::rand(&mut rng);
        let b = G1::rand(&mut rng);
        let a_affine = a.into_affine();
        let b_affine = b.into_affine();
        let mut gadget_a = G1Gadget::alloc(&mut cs.ns(|| "a"), || Ok(a)).unwrap();
        let gadget_b = G1Gadget::alloc(&mut cs.ns(|| "b"), || Ok(b)).unwrap();
        assert_eq!(gadget_a.x.value.unwrap(), a_affine.x);
        assert_eq!(gadget_a.y.value.unwrap(), a_affine.y);
        assert_eq!(gadget_b.x.value.unwrap(), b_affine.x);
        assert_eq!(gadget_b.y.value.unwrap(), b_affine.y);

        // Check addition
        let ab = a + &b;
        let ab_affine = ab.into_affine();
        let gadget_ab = gadget_a.add(&mut cs.ns(|| "ab"), &gadget_b).unwrap();
        let gadget_ba = gadget_b.add(&mut cs.ns(|| "ba"), &gadget_a).unwrap();
        gadget_ba
            .enforce_equal(&mut cs.ns(|| "b + a == a + b?"), &gadget_ab)
            .unwrap();

        let ab_val = gadget_ab
            .get_value()
            .expect("Doubling should be successful")
            .into_affine();
        assert_eq!(ab_val, ab_affine, "Result of addition is unequal");

        // Check doubling
        let aa = a.double();
        let aa_affine = aa.into_affine();
        gadget_a.double_in_place(&mut cs.ns(|| "2a")).unwrap();
        let aa_val = gadget_a
            .get_value()
            .expect("Doubling should be successful")
            .into_affine();
        assert_eq!(
            aa_val, aa_affine,
            "Gadget and native values are unequal after double."
        );

        // Check mul_bits
        let scalar = Fr::rand(&mut rng);
        let native_result = aa.into_affine().mul(scalar) + &b;
        let native_result = native_result.into_affine();

        let mut scalar: Vec<bool> = BitIterator::new(scalar.into_repr()).collect();
        // Get the scalar bits into little-endian form.
        scalar.reverse();
        let input = Vec::<Boolean>::alloc(cs.ns(|| "Input"), || Ok(scalar)).unwrap();
        let result = gadget_a
            .mul_bits(cs.ns(|| "mul_bits"), &gadget_b, input.iter())
            .unwrap();
        let result_val = result.get_value().unwrap().into_affine();
        assert_eq!(
            result_val, native_result,
            "gadget & native values are diff. after scalar mul"
        );

        if !cs.is_satisfied() {
            println!("{:?}", cs.which_is_unsatisfied().unwrap());
        }

        assert!(cs.is_satisfied());
    }

    #[test]
    fn bls12_g2_gadget_test() {
        let mut cs = TestConstraintSystem::<SW6>::new();

        let a: G2 = rand::random();
        let b: G2 = rand::random();
        let a_affine = a.into_affine();
        let b_affine = b.into_affine();

        let mut gadget_a = G2Gadget::alloc(&mut cs.ns(|| "a"), || Ok(a)).unwrap();
        let gadget_b = G2Gadget::alloc(&mut cs.ns(|| "b"), || Ok(b)).unwrap();
        assert_eq!(gadget_a.x.get_value().unwrap(), a_affine.x);
        assert_eq!(gadget_a.y.get_value().unwrap(), a_affine.y);
        assert_eq!(gadget_b.x.get_value().unwrap(), b_affine.x);
        assert_eq!(gadget_b.y.get_value().unwrap(), b_affine.y);

        let ab = a + &b;
        let ab_affine = ab.into_affine();
        let gadget_ab = gadget_a.add(&mut cs.ns(|| "ab"), &gadget_b).unwrap();
        let gadget_ba = gadget_b.add(&mut cs.ns(|| "ba"), &gadget_a).unwrap();
        gadget_ba
            .enforce_equal(&mut cs.ns(|| "b + a == a + b?"), &gadget_ab)
            .unwrap();
        assert_eq!(gadget_ab.x.get_value().unwrap(), ab_affine.x);
        assert_eq!(gadget_ab.y.get_value().unwrap(), ab_affine.y);

        let aa = a.double();
        let aa_affine = aa.into_affine();
        gadget_a.double_in_place(&mut cs.ns(|| "2a")).unwrap();

        assert_eq!(gadget_a.x.get_value().unwrap(), aa_affine.x);
        assert_eq!(gadget_a.y.get_value().unwrap(), aa_affine.y);

        if !cs.is_satisfied() {
            println!("{:?}", cs.which_is_unsatisfied().unwrap());
        }

        assert!(cs.is_satisfied());
    }
}
