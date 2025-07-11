/// Assertion that a matrix is symmetric.
#[macro_export]
macro_rules! assert_symmetric_mat {
    ($matrix:expr) => {
        assert_eq!(
            $matrix.nrows(),
            $matrix.ncols(),
            "Matrix is not symmetric: {} rows != {} cols",
            $matrix.nrows(),
            $matrix.ncols()
        );
    };
}

/// Assertion that two matrices have the same number of rows.
#[macro_export]
macro_rules! assert_nrows {
    ($matrix1:expr, $matrix2:expr) => {
        assert_eq!(
            $matrix1.nrows(),
            $matrix2.nrows(),
            "Matrices have different number of rows: {} != {}",
            $matrix1.nrows(),
            $matrix2.nrows()
        );
    };
}
