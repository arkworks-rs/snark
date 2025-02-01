use ark_ff::Field;
use ark_std::vec::Vec;
/// A sparse representation of constraint matrices.
pub type Matrix<F> = Vec<Vec<(F, usize)>>;

/// Transpose a matrix of field elements.
pub fn transpose<F: Field>(matrix: &Matrix<F>, num_col: usize) -> Matrix<F> {
    // Initialize the transposed matrix with empty vectors
    let mut transposed: Matrix<F> = vec![Vec::new(); num_col];

    // Iterate through each row and each element in the row
    for (row_index, row) in matrix.iter().enumerate() {
        for &(value, col_index) in row {
            // Add the element to the new row (which is originally a column) in the
            // transposed matrix
            transposed[col_index].push((value, row_index));
        }
    }

    // Return the transposed matrix
    transposed
}

/// Multiply a matrix by a vector.
pub fn mat_vec_mul<F: Field>(matrix: &Matrix<F>, vector: &[F]) -> Vec<F> {
    let mut output: Vec<F> = Vec::new();
    for row in matrix {
        let mut sum: F = F::zero();
        for (value, col) in row {
            sum += vector[*col] * value;
        }
        output.push(sum);
    }
    output
}
