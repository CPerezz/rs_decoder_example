use crate::field::Fp;
use crate::polynomial::Polynomial;
use std::convert::TryInto;

pub fn decode_poly(wD: Polynomial, borrows: &[u8], t: Fp) -> Polynomial {
    // Compute syndromes poly
    let sD = compute_syndroms_poly(&wD, t);
    println!("Syndromes Poly: {:?}", sD);
    // Compute borrows poly
    let lD = compute_borrows_localizator_poly(borrows);
    println!("Borrows Poly: {:?}", lD);
    // Compute modified syndrome poly
    let tD = compute_modified_syndrome_poly(sD, lD, t);
    println!("Modif Syndromes Poly: {:?}", tD);
    // Compute evaluator poly
    let evD = compute_evaluator_poly(tD, t, ((borrows.len() - 1) as u8).try_into().unwrap());
    println!("Evaluator Poly: {:?}", evD);
    Polynomial::zero()
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

fn compute_modified_syndrome_poly(sD: Polynomial, lD: Polynomial, t: Fp) -> Polynomial {
    // Compute x^2t poly
    let mut d_2t_poly = vec![Fp::zero(); 2 * t.0 as usize + 2];
    d_2t_poly[2 * t.0 as usize + 1] = Fp::one();

    // Modified syndrome poly = L(D) * S(D) % D^2t
    // Take the residue of (L(D) * S(D)) / D^2t
    (Polynomial::from_fp_coefficients_vec(d_2t_poly) / (lD * sD)).1
}

/// Computes the evaluator poly by applying the gdc to (D^2t, t(D)) until
/// if finds a residue with degree < (t + #borrows/2).
fn compute_evaluator_poly(tD: Polynomial, t: Fp, num_borrows: u8) -> Polynomial {
    // Compute x^2t poly
    let mut d_2t_poly = vec![Fp::zero(); 2 * t.0 as usize + 2];
    d_2t_poly[2 * t.0 as usize + 1] = Fp::one();
    let d_2t_poly = Polynomial::from_fp_coefficients_vec(d_2t_poly);

    // Execute custom GCD until residue with degree < (t + #borrows/2).
    d_2t_poly.gcd_for_eavluator(tD, t, num_borrows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decodif_test() {
        let wD = Polynomial::from_coefficients_slice(&[
            27, 1, 24, 25, 1, 27, 21, 12, 3, 4, 5, 8, 9, 30, 3, 0, 3, 0, 22, 24, 12, 22, 30, 30,
            28, 23, 1, 1, 30, 0, 4,
        ]);
        let borrows_positions = &[13, 16];
        let t = Fp(3);
        decode_poly(wD, borrows_positions, t);
    }

    #[test]
    fn syndroms_comp() {
        let wD = Polynomial::from_coefficients_slice(&[
            27, 1, 24, 25, 1, 27, 21, 12, 3, 4, 5, 8, 9, 30, 3, 0, 3, 0, 22, 24, 12, 22, 30, 30,
            28, 23, 1, 1, 30, 0, 4,
        ]);
        let t: Fp = 3.into();
        assert!(wD.degree() == 30);
        let sD = compute_syndroms_poly(&wD, t);
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
    }
}
