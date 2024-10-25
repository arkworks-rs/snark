use super::ArithmeticCircuit;
use ark_bls12_377::Fq as FqBLS;
use ark_bn254::Fr as FrBN;
use ark_ff::Field;

// Defining equation of BLS12-377: y^2 = x^3 + 1 (over Fq)
pub fn generate_bls12_377_circuit() -> ArithmeticCircuit<FqBLS> {
    let mut circuit = ArithmeticCircuit::new();

    // Ligero circuits must start with a constant 1
    let one = circuit.constant(FqBLS::ONE);

    let x = circuit.new_variable_with_label("x");
    let y = circuit.new_variable_with_label("y");

    let y_squared = circuit.pow(y, 2);
    let minus_y_squared = circuit.minus(y_squared);
    let x_cubed = circuit.pow(x, 3);

    // Ligero will prove x^3 + 1 - y^2 + 1 = 1 Note that one could compute the
    // left-hand side as x^3 + 2 - y^2 in order to save one addition gate
    circuit.add_nodes([x_cubed, one, minus_y_squared, one]);
    circuit

    // n_i = 2, s = 8

    // Original circuit
    //     0: Constant(1)
    //     1: Variable
    //     2: Variable
    //     3: node(2) * node(2)
    //     4: Constant(21888242871839275222246405745257275088696311157297823662689037894645226208582)
    //     5: node(4) * node(3)
    //     6: node(1) * node(1)
    //     7: node(6) * node(1)
    //     8: node(7) + node(0)
    //     9: node(8) + node(5)
    //     10: node(9) + node(0)
}

/// (x^2 + y^2)^2 - 120x^2 + 80y^2 + 1 = 1
pub fn generate_lemniscate_circuit() -> ArithmeticCircuit<FrBN> {
    let mut circuit = ArithmeticCircuit::new();

    // Ligero circuits must start with a constant 1
    let one = circuit.constant(FrBN::ONE);

    let x = circuit.new_variable();
    let y = circuit.new_variable();

    let a = circuit.constant(FrBN::from(120));
    let b = circuit.constant(FrBN::from(80));

    let x_2 = circuit.mul(x, x);
    let y_2 = circuit.mul(y, y);

    let a_x_2 = circuit.mul(a, x_2);
    let b_y_2 = circuit.mul(b, y_2);
    let minus_a_x_2 = circuit.minus(a_x_2);

    let x_2_plus_y_2 = circuit.add(x_2, y_2);
    let b_y_2_minus_a_x_2 = circuit.add(b_y_2, minus_a_x_2);

    let x_2_plus_y_2_2 = circuit.mul(x_2_plus_y_2, x_2_plus_y_2);

    circuit.add_nodes([x_2_plus_y_2_2, b_y_2_minus_a_x_2, one]);
    circuit
}

pub fn generate_3_by_3_determinant_circuit() -> ArithmeticCircuit<FrBN> {
    let mut circuit = ArithmeticCircuit::new();

    // Ligero circuits must start with a constant 1
    let one = circuit.constant(FrBN::ONE);

    let vars = circuit.new_variables(9);
    let det = circuit.new_variable();

    let aei = circuit.mul_nodes([vars[0], vars[4], vars[8]]);
    let bfg = circuit.mul_nodes([vars[1], vars[5], vars[6]]);
    let cdh = circuit.mul_nodes([vars[2], vars[3], vars[7]]);

    let ceg = circuit.mul_nodes([vars[2], vars[4], vars[6]]);
    let bdi = circuit.mul_nodes([vars[1], vars[3], vars[8]]);
    let afh = circuit.mul_nodes([vars[0], vars[5], vars[7]]);

    let sum1 = circuit.add_nodes([aei, bfg, cdh]);
    let sum2 = circuit.add_nodes([ceg, bdi, afh]);

    let minus_sum2 = circuit.minus(sum2);
    let minus_det = circuit.minus(det);

    circuit.add_nodes([sum1, minus_sum2, minus_det, one]);
    circuit
}
