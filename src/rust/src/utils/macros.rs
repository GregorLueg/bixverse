///////////////////
// Matrix macros //
///////////////////

/// Assertion that a matrix is symmetric.
#[macro_export]
macro_rules! assert_symmetric_mat {
    ($matrix:expr) => {{
        let nrows = $matrix.nrows();
        let ncols = $matrix.ncols();

        assert_eq!(
            nrows, ncols,
            "Matrix is not square: {} rows != {} cols",
            nrows, ncols
        );

        let check_size = nrows.min(10);

        for i in 0..check_size {
            for j in (i + 1)..check_size {
                let val_ij = $matrix.get(i, j);
                let val_ji = $matrix.get(j, i);
                assert!(
                    (val_ij - val_ji).abs() < 1e-10,
                    "Matrix not symmetric at ({}, {}): {} != {}",
                    i,
                    j,
                    val_ij,
                    val_ji
                );
            }
        }
    }};
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

/// Assertion that all vectors have the same length.
#[macro_export]
macro_rules! assert_same_len {
    ($($vec:expr),+ $(,)?) => {
        {
            let lengths: Vec<usize> = vec![$($vec.len()),+];
            let first_len = lengths[0];

            if !lengths.iter().all(|&len| len == first_len) {
                panic!(
                    "Vectors have different lengths: {:?}",
                    lengths
                );
            }
        }
    };
}
