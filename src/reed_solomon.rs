use crate::field::Fp;
use crate::polynomial::Polynomial;

pub fn decode_poly(wD: Polynomial, borrows: &[u8], t: Fp) -> Polynomial {
    // Compute syndromes poly
    let sD = compute_syndroms_poly(&wD, t);

    println!("Syndromes Poly: {:?}", sD);
    unimplemented!()
}

fn compute_syndroms_poly(wd: &Polynomial, t: Fp) -> Polynomial {
    let mut result: Vec<Fp> = vec![];
    for i in 1..=(2 * t.0) {
        result.push(wd.evaluate(i.into()))
    }
    Polynomial::from_fp_coefficients_vec(result)
}

fn compute_borrows_localizator_poly(borrows: &[u8]) -> Polynomial {
    let mut borrow_locators = vec![];

    for borrow_position in borrows {
        borrow_locators.push(Polynomial::from_coefficients_slice(&[1, *borrow_position]));
    }

    borrow_locators.iter().product()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn syndroms_comp() {
        let wD = Polynomial::from_coefficients_slice(&[
            27, 1, 24, 25, 1, 27, 21, 12, 3, 4, 5, 8, 9, 30, 3, 0, 3, 0, 22, 24, 12, 22, 30, 30,
            28, 23, 1, 1, 30, 0, 4,
        ]);
        let t: Fp = 3.into();
        assert!(wD.degree() == 30);
        let sD = compute_syndroms_poly(&wD, t);
        println!("Syndromes Poly: {:?}", sD);
    }

    #[test]
    fn compute_borrows_localizator_poly_test() {
        let wD = Polynomial::from_coefficients_slice(&[
            27, 1, 24, 25, 1, 27, 21, 12, 3, 4, 5, 8, 9, 30, 3, 0, 3, 0, 22, 24, 12, 22, 30, 30,
            28, 23, 1, 1, 30, 0, 4,
        ]);
        let t: Fp = 3.into();
        assert!(wD.degree() == 30);
        let sD = compute_syndroms_poly(&wD, t);
        let lD = compute_borrows_localizator_poly(&[13, 16]);
        println!("Borrow Localizators Poly: {:?}", lD);
    }
}
