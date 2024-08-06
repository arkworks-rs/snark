use ark_ff::Field;
use ark_std::vec::Vec;
/// A sparse representation of constraint matrices.
pub type Matrix<F> = Vec<Vec<(F, usize)>>;

/// Transpose a matrix of field elements.
pub fn transpose<F: Field>(matrix: &Matrix<F>) -> Matrix<F> {
    // First, find the maximum column index to know the size of the transposed
    // matrix
    let max_cols = matrix
        .iter()
        .flat_map(|row| row.iter().map(|&(_, col)| col + 1))
        .max()
        .unwrap_or(0);

    // Initialize the transposed matrix with empty vectors
    let mut transposed: Matrix<F> = vec![Vec::new(); max_cols];

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
