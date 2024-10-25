use super::Expression;
use ark_bls12_377::Fq;
use ark_bn254::Fr;

pub fn generate_bls12_377_expression() -> Expression<Fq> {
    let x = Expression::variable("x");
    let y = Expression::variable("y");

    1 + (1 + x.pow(3) - y.pow(2))
}

/// (x^2 + y^2)^2 - 120x^2 + 80y^2 + 1 = 1
pub fn generate_lemniscate_expression() -> Expression<Fr> {
    let x = Expression::variable("x");
    let y = Expression::variable("y");

    1 + (x.clone().pow(2) + y.clone().pow(2)).pow(2) - 120 * x.pow(2) + 80 * y.pow(2)
}

pub fn generate_3_by_3_determinant_expression() -> Expression<Fr> {
    let matrix = (0..3)
        .map(|i| {
            (0..3)
                .map(|j| Expression::variable(&format!("x_{i}_{j}")))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let possitive_diagonal = (0..3)
        .map(|k| {
            vec![0, 4, 8]
                .into_iter()
                .zip(0..3)
                .map(|(j, i)| matrix[i][(j + k) % 3].clone())
                .product()
        })
        .sum::<Expression<Fr>>();

    let negative_diagonal = (0..3)
        .map(|k| {
            vec![2, 4, 6]
                .into_iter()
                .zip(0..3)
                .map(|(j, i)| matrix[i][(j + k) % 3].clone())
                .product()
        })
        .sum::<Expression<Fr>>();

    let det = Expression::variable("det");

    1 + (possitive_diagonal - negative_diagonal - det)
}
