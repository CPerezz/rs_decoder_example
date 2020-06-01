use crate::field::Fp;
use crate::polynomial::Polynomial;
use std::convert::TryInto;

pub fn decode_poly(wD: Polynomial, borrows: &[u8], t: Fp) -> Polynomial {
    // Compute syndromes poly
    let sD = compute_syndroms_poly(&wD, t);
    println!("Syndromes Poly: {:?}", sD);
    // Compute borrows localizator poly
    let lD = compute_borrows_localizator_poly(borrows);
    println!("Borrows Poly: {:?}", lD);
    // Compute modified syndrome poly
    let tD = compute_modified_syndrome_poly(sD, lD.clone(), t);
    println!("Modif Syndromes Poly: {:?}", tD);
    // Compute evaluator poly
    let evD = compute_evaluator_poly(
        tD.clone(),
        t,
        ((borrows.len() - 1) as u8).try_into().unwrap(),
    );
    println!("Evaluator Poly: {:?}", evD);
    // Compute localizator poly
    let l_hat_t = compute_localizator_poly(tD, t);
    println!("Localizator Poly: {:?}", l_hat_t.clone());
    // Chien search over l_hat_d (Localizator poly)
    // These are stored starting from the biggest degree term
    let error_positions = chien_search(l_hat_t.clone());

    println!("Error positions: {:?}", error_positions);
    // Calculate L'(D) poly
    let l_prime_d = localizator_prime_poly(lD, l_hat_t);
    println!("L'(D) Poly: {:?}", l_prime_d.clone());

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

/// Computes the localizator poly by applying the Extended GCD
/// algorigthm.
fn compute_localizator_poly(tD: Polynomial, t: Fp) -> Polynomial {
    // Compute x^2t poly
    let mut d_2t_poly = vec![Fp::zero(); 2 * t.0 as usize + 2];
    d_2t_poly[2 * t.0 as usize + 1] = Fp::one();
    let d_2t_poly = Polynomial::from_fp_coefficients_vec(d_2t_poly);

    // Extended GCD between tD & x^2t
    d_2t_poly.extended_gcd(tD).0
}

/// Evaluates the localizator poly on points 0... until it finds
/// as many 0's as the degree of the localizator error.
///
/// This is indeed because the localizator error degree is already
/// telling us how many errors does the recieved message has. So we know
/// that we need to be evaluating the poly until two evaluations give 0.
/// The error_positions are returned starting by the most significant
/// monomial of the polynomial.
fn chien_search(base_poly: Polynomial) -> Vec<u8> {
    let mut error_positions = vec![];
    let mut i = 1u8;

    while error_positions.len() < base_poly.degree() {
        let eval = base_poly.evaluate(Fp(i));
        if eval == Fp::zero() {
            error_positions.push(i);
        };
        i += 1;
    }

    error_positions
}

/// Computes the L(D) poly by multiplying the borrow localizator times
/// the localizator_poly and derivating the result.
fn localizator_prime_poly(lD: Polynomial, l_hat_d: Polynomial) -> Polynomial {
    let lD_final = lD * l_hat_d;

    // Return ld_final's derivate
    lD_final.derivate()
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
