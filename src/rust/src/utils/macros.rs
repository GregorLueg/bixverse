///////////////////
// Matrix macros //
///////////////////

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

/// Assertion that two matrices have the same dimensions (rows and columns).
#[macro_export]
macro_rules! assert_same_dims {
    ($matrix1:expr, $matrix2:expr) => {
        assert_eq!(
            ($matrix1.nrows(), $matrix1.ncols()),
            ($matrix2.nrows(), $matrix2.ncols()),
            "Matrices have different dimensions: {}x{} != {}x{}",
            $matrix1.nrows(),
            $matrix1.ncols(),
            $matrix2.nrows(),
            $matrix2.ncols()
        );
    };
}

///////////////////
// Vector macros //
///////////////////

/// Assertion that two vectors have the same length.
#[macro_export]
macro_rules! assert_same_len_vec2 {
    ($vec1:expr, $vec2:expr) => {
        let len1 = $vec1.len();
        let len2 = $vec2.len();

        assert!(
            len1 == len2,
            "Vectors have different lengths: {} != {}",
            len1,
            len2,
        );
    };
}

/// Assertion that three vectors have the same length.
#[macro_export]
macro_rules! assert_same_len_vec3 {
    ($vec1:expr, $vec2:expr, $vec3:expr) => {
        let len1 = $vec1.len();
        let len2 = $vec2.len();
        let len3 = $vec3.len();

        assert!(
            len1 == len2 && len2 == len3,
            "Vectors have different lengths: {} != {} != {}",
            len1,
            len2,
            len3
        );
    };
}
