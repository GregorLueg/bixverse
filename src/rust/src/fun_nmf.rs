
use faer::{Mat, mat};

fn blockpivot(m1: Mat<f64>, m2: Mat<f64>) {
    let m1_t_m1 = m1.transpose() * &m1;
    let m1_t_m2 = m1.transpose() * &m2;

    let (n, k) = m1_t_m2.shape();
    let max_iter = n * 5;

    let x = Mat::<f64>::zeros(n, k);
    let y = -m1_t_m2;

    // Use vector as faer only supports f64/f32 matrices.
    let mut pass_set = vec![vec![false; k]; n];

    let p_bar = 3;
    let mut p_vec = vec![p_bar; k];
    let mut ninf_vec = vec![(n + 1) as i32; k];

    let mut not_opt_set =  vec![vec![false; k]; n];
    for i in 0..n {
        for j in 0..k {
            not_opt_set[i][j] = y[(i, j)] < 0.0
        }
    }
    let infea_set = vec![vec![false; k]; n];

    let mut not_good = vec![0 as i32; k];
    for row in not_opt_set {
        for (j, &val) in row.iter().enumerate() {
            if val {
                not_good[j] += 1;
            }
        }
    }
    let not_opt_colset: Vec<bool> = not_good.iter().map(|&x| x > 0).collect();
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

        let not_good_le_ninf= &not_good
            .iter()
            .zip(ninf_vec.iter())
            .map(|(&x, &y)| x < y)
            .collect::<Vec<bool>>();

        let not_good_gr_ninf = &not_good
            .iter()
            .zip(ninf_vec.iter())
            .map(|(&x, &y)| x >= y)
            .collect::<Vec<bool>>();

        let cols_set_1 = &not_opt_colset
            .iter()
            .zip(not_good_le_ninf)
            .map(|(&x, &y)| x & y)
            .collect::<Vec<bool>>();

        let temp_1 = &not_opt_colset
            .iter()
            .zip(not_good_gr_ninf)
            .map(|(&x, &y)| x & y)
            .collect::<Vec<bool>>();

        let temp_2 = p_vec.iter().map(|&x| x >= 1).collect::<Vec<bool>>();

        let cols_set_2 = &temp_1
            .iter()
            .zip(&temp_2)
            .map(|(&x, y)| x & y)
            .collect::<Vec<bool>>();

        let cols_set_3 = &temp_1
            .iter()
            .zip(&temp_2)
            .map(|(&x, y)| x & !y)
            .collect::<Vec<bool>>();
       
        let cols_1 = non_zero_indexes(cols_set_1);
        let cols_2 = non_zero_indexes(cols_set_2);
        let cols_3 = non_zero_indexes(cols_set_3);

        println!("{:?}", cols_1);
        println!("{:?}", cols_2);
        println!("{:?}", cols_3);
 
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
        
        /* 
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

        let m1_t_m2_slice = m1_t_m2.select(Axis(1), &not_opt_cols);
        println!("m1_t_m1 shape: {:?}", m1_t_m1.shape());
        println!("m1_t_m2_slice shape: {:?}", m1_t_m2_slice.shape());


        */
    } 
}

fn non_zero_indexes(bool_vec: &Vec<bool>) -> Vec<usize> {
    bool_vec
        .iter()
        .enumerate()
        .filter_map(|(i, &not_opt)| if not_opt { Some(i) } else { None })
        .collect::<Vec<usize>>()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_group_cols() {
        let a = mat![
            [0.0, 1.0, 2.0],
            [9.0, 7.0, 6.0],
            [3.0, 4.0, 5.0],
        ];

        let w = mat![
            [3.0, 1.0],
            [5.0, 3.0],
            [2.0, 4.0]
        ];

        blockpivot(w, a);

    }
}