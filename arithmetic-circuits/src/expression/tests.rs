use ark_bls12_377::{Fq, G1Affine};
use ark_bn254::Fr;
use ark_ec::short_weierstrass::Affine;
use ark_ff::{Field, UniformRand};
use ark_std::{collections::HashMap, ops::Deref, test_rng};
use itertools::Itertools;

use crate::{arithmetic_circuit::Node, expression::examples::*};

use super::{Expression, ExpressionInner};

#[test]
fn test_get_variables() {
    let circuit = generate_bls12_377_expression().to_arithmetic_circuit();

    assert_eq!(circuit.get_variable("x"), 4);
    assert_eq!(circuit.get_variable("y"), 0);

    let circuit = generate_lemniscate_expression().to_arithmetic_circuit();

    assert_eq!(circuit.get_variable("x"), 10);
    assert_eq!(circuit.get_variable("y"), 8);
}

#[test]
fn test_same_reference() {
    let mut f1 = Expression::<Fr>::variable("x");
    let mut f2 = Expression::<Fr>::variable("y");

    let original_f1 = f1.clone();

    for _ in 0..10 {
        let next = f1.clone() + f2.clone();
        f1 = f2;
        f2 = next;
    }

    let curious_result = f2 * original_f1.clone();

    match curious_result.0.deref() {
        ExpressionInner::Mul(_, right) => {
            assert_eq!(right.pointer(), original_f1.pointer())
        }
        _ => panic!("Expected a multiplication expression"),
    }
}

#[test]
fn test_addition() {
    let a = Expression::variable("x");
    let b = Expression::variable("y");

    let expression = a + b;
    let circuit = expression.to_arithmetic_circuit();

    assert_eq!(
        circuit.evaluate_with_labels(vec![("x", Fr::from(3)), ("y", Fr::from(5))]),
        Fr::from(8),
    );
}

#[test]
fn test_multiplication() {
    let a = Expression::variable("x");
    let b = Expression::variable("y");

    let expression = a * b;
    let circuit = expression.to_arithmetic_circuit();

    assert_eq!(
        circuit.evaluate_with_labels(vec![("x", Fr::from(3)), ("y", Fr::from(5))]),
        Fr::from(15),
    );
}

#[test]
fn test_subtraction() {
    let a = Expression::variable("x");
    let b = Expression::variable("y");

    let expression = a - b;
    let circuit = expression.to_arithmetic_circuit();

    assert_eq!(
        circuit.evaluate_with_labels(vec![("x", Fr::from(3)), ("y", Fr::from(5))]),
        Fr::from(-2),
    );
}

#[test]
fn test_some_operations() {
    let x_f = Fr::from(5);
    let y_f = Fr::from(3);

    let output = x_f.pow([3]) + (y_f - Fr::ONE).pow([11]) + Fr::from(13);

    let x_exp = Expression::constant(x_f);
    let y_exp = Expression::constant(y_f);

    let output_exp = 13 + x_exp.pow(3) + (y_exp - Fr::ONE).pow(11);
    let circ_output = output_exp.to_arithmetic_circuit().evaluate(vec![]);

    assert_eq!(output, circ_output);
}

/*
                                                Add
                    /                                                                   \
                  Add                                                                    Mul
        /                       \                                     /                                      \
    Constant(3)                 Mul                                  Add                                    Add
                        /              \                        /               \                   /                   \
                    Constant(2)        Mul                  Constant(3)         Mul             Constant(1)              Mul
                                    /       \                               /       \                             /             \
                                    |       |                         Constant(2)   |                         Constant(2)   Variable(1)
                                    |       |                                       |                                           |
                                    -------------------------------------------------                                           |
                                            |               |                                                                   |
                                            -------------------------------------------------------------------------------------
                                                            |                           |
                                                        Variable(0)                 Variable(1)
Original:

    16: Add(15, 5)
    15: Add(14, 13)
    14: Constant(3)
    13: Mul(12, 11)
    12: Constant(2)
    11: Mul(10, 9)
    10: Variable(0)
    9: Variable(1)
    8: Mul(4, 10)
    7: Add(6, 8)
    6: Constant(3)
    5: Mul(3, 7)
    4: Constant(2)
    3: Add(2, 1)
    2: Constant(1)
    1: Mul(0, 9)
    0: Constant(2)

After filtering:

    13: Add(12, 7) -> (3 + 2 * x * y) + (3 + 2 * x) * (1 + 2 * y) = 60
    12: Add(5, 11) -> 3 + 2 * x * y = 15
    11: Mul(0, 10) -> 2 * x * y = 12
    10: Mul(9, 8)  -> x * y = 6
    9: Variable(0) -> x = 3
    8: Variable(1) -> y = 2
    7: Mul(6, 3)   -> (3 + 2 * x) * (1 + 2 * y) = 45
    6: Add(5, 4)   -> 3 + 2 * x = 9
    5: Constant(3) -> 3
    4: Mul(0, 9)   -> 2 * x = 6
    3: Add(2, 1)   -> 1 + 2 * y = 5
    2: Constant(1) -> 1
    1: Mul(0, 8)   -> 2 * y = 4
    0: Constant(2) -> 2

*/
#[test]
fn test_to_arithmetic_circuit_1() {
    let x = Expression::variable("x");
    let y = Expression::variable("y");

    let expression = (3 + 2 * (x.clone() * y.clone())) + ((3 + 2 * x) * (1 + 2 * y));

    let circuit = expression.to_arithmetic_circuit();

    assert_eq!(
        circuit.nodes,
        vec![
            Node::Add(12, 7),
            Node::Add(5, 11),
            Node::Mul(0, 10),
            Node::Mul(9, 8),
            Node::Variable("x".to_string()),
            Node::Variable("y".to_string()),
            Node::Mul(6, 3),
            Node::Add(5, 4),
            Node::Constant(Fr::from(3)),
            Node::Mul(0, 9),
            Node::Add(2, 1),
            Node::Constant(Fr::ONE),
            Node::Mul(0, 8),
            Node::Constant(Fr::from(2)),
        ]
        .iter()
        .rev()
        .cloned()
        .collect::<Vec<_>>()
    );

    assert_eq!(
        circuit.constants,
        [(Fr::from(3), 5), (Fr::ONE, 2), (Fr::from(2), 0)]
            .iter()
            .cloned()
            .collect::<HashMap<_, _>>()
    );

    assert_eq!(
        circuit.evaluation_trace_with_labels(vec![("x", Fr::from(3)), ("y", Fr::from(2))], 13),
        vec![
            Some(Fr::from(60)),
            Some(Fr::from(15)),
            Some(Fr::from(12)),
            Some(Fr::from(6)),
            Some(Fr::from(3)),
            Some(Fr::from(2)),
            Some(Fr::from(45)),
            Some(Fr::from(9)),
            Some(Fr::from(3)),
            Some(Fr::from(6)),
            Some(Fr::from(5)),
            Some(Fr::from(1)),
            Some(Fr::from(4)),
            Some(Fr::from(2)),
        ]
        .iter()
        .rev()
        .cloned()
        .collect::<Vec<_>>()
    );
}

/*

                    Mul()
            /                   \
            Add()                 Add()
        /       \          /              \
        |        |      Variable(c)      Mul()
        |        |                     /       \
        -------------------------------        |
            |     |                            |
            |     ------------------------------
            |                   |
            Variable(a)         Variable(b)

    6: Mul(5, 2) -> 35
    5: Add(4, 3) -> 5
    4: Variable(a) -> 3
    3: Variable(b) -> 2
    2: Add(1, 0) -> 7
    1: Variable(c) -> 1
    0: Mul(4, 3) -> 6

*/

#[test]
fn test_to_arithmetic_circuit_2() {
    let a = Expression::variable("a");
    let b = Expression::variable("b");
    let c = Expression::variable("c");

    let expression = (a.clone() + b.clone()) * (c + a * b);

    let circuit = expression.to_arithmetic_circuit();

    assert_eq!(
        circuit.nodes,
        vec![
            Node::Mul(5, 2),
            Node::Add(4, 3),
            Node::Variable("a".to_string()),
            Node::Variable("b".to_string()),
            Node::Add(1, 0),
            Node::Variable("c".to_string()),
            Node::Mul(4, 3),
        ]
        .iter()
        .rev()
        .cloned()
        .collect::<Vec<_>>()
    );

    assert_eq!(circuit.constants, HashMap::new());

    assert_eq!(
        circuit.evaluation_trace_with_labels(
            vec![("a", Fr::from(3)), ("b", Fr::from(2)), ("c", Fr::from(1))],
            6
        ),
        vec![
            Some(Fr::from(6)),
            Some(Fr::from(1)),
            Some(Fr::from(7)),
            Some(Fr::from(2)),
            Some(Fr::from(3)),
            Some(Fr::from(5)),
            Some(Fr::from(35))
        ]
    );
}

#[test]
fn test_to_arithmetic_circuit_3() {
    let circuit = generate_3_by_3_determinant_expression().to_arithmetic_circuit();

    let values = (0..3)
        .cartesian_product(0..3)
        .map(|(i, j)| {
            (
                format!("x_{}_{}", i, j),
                Fr::from((3 * i + j) * (3 * i + j)),
            )
        })
        .collect::<Vec<_>>();

    let vars = [values, vec![("det".to_string(), Fr::from(-216))]].concat();

    assert_eq!(
        circuit.evaluate_with_labels(vars.iter().map(|(k, v)| (k.as_str(), *v)).collect()),
        Fr::ONE,
    );
}

#[test]
fn test_to_arithmetic_circuit_4() {
    let circuit = generate_bls12_377_expression().to_arithmetic_circuit();
    let Affine { x, y, .. } = G1Affine::rand(&mut test_rng());

    assert_eq!(
        circuit.evaluate_with_labels(vec![("x", x), ("y", y)]),
        Fq::ONE
    );
}

#[test]
fn test_to_arithmetic_circuit_5() {
    let circuit = generate_lemniscate_expression().to_arithmetic_circuit();

    assert_eq!(
        circuit.evaluate_with_labels(vec![("x", Fr::from(8)), ("y", Fr::from(4))]),
        Fr::ONE
    );
}
