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
