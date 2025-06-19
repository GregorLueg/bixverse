use ndarray::{array, Array1, Array2, Axis, s};
use ndarray_linalg::Solve;

fn blockpivot(m1: Array2<f64>, m2: Array2<f64>) {
    let m1_t_m1 = m1.t().dot(&m1);
    let m1_t_m2 = m1.t().dot(&m2);

    let shape = m1_t_m2.shape();
    let (n, k) = (shape[0], shape[1]);
    let max_iter = n * 5;

    
    let x = Array2::<f64>::zeros((n, k));
    // TODO - do we need to clone here?
    let y = -m1_t_m2.clone();

    let p_bar = 3;
    let mut pass_set = Array2::<bool>::from_elem((n, k), false);
    let mut p_vec = Array1::<u32>::from_elem( k, p_bar);
    let mut ninf_vec = Array1::<u32>::from_elem( k, n as u32 + 1);

    let not_opt_set = y.map(|x| *x < 0.0);
    let infea_set = Array2::<bool>::from_elem((n, k), false);

    let not_good = not_opt_set.mapv(|x|x as u32).sum_axis(Axis(0)); 

    let not_opt_colset  = not_good.map(|x| *x > 0);
    let not_opt_cols: Vec<usize> = non_zero_indexes(&not_opt_colset);

    let mut big_iter = 0;
    let mut success = true;

    while !not_opt_cols.is_empty() {
        // TODO - shift this to end of loop?
        big_iter += 1;

        if big_iter > max_iter {
            success = false;
            break
        }

        let cols_set_1 = &not_opt_colset & (
            &not_good.iter()
                .zip(ninf_vec.iter())
                .map(|(&x, &y)| x < y)
                .collect::<Array1<bool>>()
            );
        

        let temp_1 = &not_opt_colset & (
            &not_good.iter()
                .zip(ninf_vec.iter())
                .map(|(&x, &y)| x >= y)
                .collect::<Array1<bool>>()
            );

        let temp_2 = p_vec.map(|x| *x >= 1);

        let cols_set_2 = &temp_1 & &temp_2;
        let cols_set_3 = &temp_1 & &temp_2.map(|x| !x);

        let cols_1 = non_zero_indexes(&cols_set_1);
        let cols_2 = non_zero_indexes(&cols_set_2);
        let cols_3 = non_zero_indexes(&cols_set_3);

        if !cols_1.is_empty() {
            for &i in &cols_1 {
                p_vec[i] = p_bar;
                ninf_vec[i] = not_good[i];
            }

            let cols_set_1 = cols_set_1.broadcast(not_opt_set.shape()).unwrap().to_owned();

            let true_set = &not_opt_set & &cols_set_1;
            let false_set = &infea_set & &cols_set_1;

            true_set
                .iter()
                .zip(pass_set.iter_mut())
                .for_each(|(&b, p)| if b { *p = true; });

            false_set
                .iter()
                .zip(pass_set.iter_mut())
                .for_each(|(&b, p)| if b { *p = false; });

        }
        
        if !cols_2.is_empty() {
            for &i in &cols_2 {
                p_vec[i] = p_vec[i] - 1;
            }

            let cols_set_2 = cols_set_2.broadcast(not_opt_set.shape()).unwrap().to_owned();

            let true_set = &not_opt_set & &cols_set_2;
            let false_set = &infea_set & &cols_set_2;

            true_set
                .iter()
                .zip(pass_set.iter_mut())
                .for_each(|(&b, p)| if b { *p = true; });

            false_set
                .iter()
                .zip(pass_set.iter_mut())
                .for_each(|(&b, p)| if b { *p = false; });

        } 
        
        if !cols_3.is_empty() {
            for &col in &cols_3 {
                let candi_set= &not_opt_set.column(col) | &infea_set.column(col);
                let to_change = *non_zero_indexes(&candi_set).iter().max().unwrap();

                for i in 0..pass_set.shape()[0] {
                    pass_set[[i, to_change]] = !pass_set[[i, to_change]];
                }
            }

        }

        //let m1_t_m2_slice = m1_t_m2.select(Axis(1), &not_opt_cols).to_owned();
        // let Z = m1_t_m1.solve(&m1_t_m2_slice).unwrap();

    } 
}

fn non_zero_indexes(array: &Array1<bool>) -> Vec<usize> {
    array
        .iter()
        .enumerate()
        .filter_map(|(i, &not_opt)| if not_opt { Some(i) } else { None })
        .collect::<Vec<usize>>()
}


#[cfg(test)]
mod tests {
    use super::*;
    use faer::{mat, Mat, Col};
    use faer::linalg::solvers::Solve;
    use faer::linalg::solvers::PartialPivLu;

    #[test]
    fn test_group_cols() {
        let a = array![
            [0.0, 1.0, 2.0],
            [9.0, 7.0, 6.0],
            [3.0, 4.0, 5.0],
        ];

        let w = array![
            [3.0, 1.0],
            [5.0, 3.0],
            [2.0, 4.0]
        ];
        
        let a = mat![
            [38.0, 26.0],
            [26.0, 26.0],
        ];

        // Define a right-hand side vector b
        let b = mat![
            [51.0, 46.0, 46.0],
            [39.0, 38.0, 40.0]
        ];

        // Compute the LU decomposition
        let lu = PartialPivLu::new(a.as_ref());

        // Solve the system
        let x = lu.solve(&b);
        println!("{:?}", x);

        // blockpivot(w, a);
        

    }
}