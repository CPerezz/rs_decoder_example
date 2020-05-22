//! Implementation of Polynomial working in GF(32).
//!
//! `MODULUS` is set as a constant so you can work in
//! whatever GF you want.
//!

use crate::field::Fp;
use std::ops::{Add, Div, Mul, Neg, Shl, Shr, Sub};
/// Stores a polynomial in coefficient form.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Polynomial {
    /// The coefficient of `x^i` is stored at location `i` in `self.coeffs`.
    pub coeffs: Vec<Fp>,
}

pub const MODULUS: Fp = Fp(31u8);

impl Polynomial {
    /// Returns the zero polynomial.
    pub fn zero() -> Self {
        Self { coeffs: Vec::new() }
    }

    /// Checks if the given polynomial is zero.
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty() || self.coeffs.iter().all(|coeff| coeff == &Fp::zero())
    }

    /// Constructs a new polynomial from a list of coefficients.
    pub fn from_coefficients_slice(coeffs: &[u8]) -> Self {
        Self::from_coefficients_vec(coeffs.to_vec())
    }

    /// Constructs a new polynomial from a list of coefficients.
    pub fn from_coefficients_vec(coeffs: Vec<u8>) -> Self {
        let mut result = Self {
            coeffs: coeffs.iter().map(|elem| Fp(*elem)).collect(),
        };
        // While there are zeros at the end of the coefficient vector, pop them off.
        result.truncate_leading_zeros();

        result
    }

    /// Constructs a new polynomial from a list of coefficients.
    pub fn from_fp_coefficients_vec(coeffs: Vec<Fp>) -> Self {
        let mut result = Self {
            coeffs: coeffs.iter().map(|elem| *elem).collect(),
        };
        // While there are zeros at the end of the coefficient vector, pop them off.
        result.truncate_leading_zeros();

        result
    }

    /// Returns the degree of the polynomial.
    pub fn degree(&self) -> usize {
        if self.is_zero() {
            0
        } else {
            assert!(self
                .coeffs
                .last()
                .map_or(false, |coeff| *coeff != Fp::zero()));
            self.coeffs.len() - 1
        }
    }

    fn truncate_leading_zeros(&mut self) {
        while self
            .coeffs
            .last()
            .map_or(false, |coeff| coeff == &Fp::zero())
        {
            self.coeffs.pop();
        }
    }

    /// Evaluates `self` at the given `point` in the field.
    pub fn evaluate(&self, point: Fp) -> Fp {
        if self.is_zero() {
            0u64;
        }
        let mut powers_of_point = vec![Fp::one()];
        let mut cur = point;
        for _ in 0..self.degree() {
            powers_of_point.push(cur);
            cur = cur * point;
        }
        assert_eq!(powers_of_point.len(), self.coeffs.len());

        powers_of_point
            .iter()
            .zip(self.coeffs.iter())
            .map(|(power, coeff)| *power * *coeff)
            .sum()
    }
}

// Resulting poly will have the same degree or less than the bigger poly that took
// part in the addition.
impl Add<&Polynomial> for &Polynomial {
    type Output = Polynomial;
    fn add(self, _rhs: &Polynomial) -> Polynomial {
        match (self.is_zero(), _rhs.is_zero()) {
            (true, false) => return _rhs.clone(),
            (false, true) => return self.clone(),
            _ => (),
        };
        let mut result = if self.degree() >= _rhs.degree() {
            let mut result = self.clone();
            for (a, b) in result.coeffs.iter_mut().zip(&_rhs.coeffs) {
                *a = *a + *b
            }
            result
        } else {
            let mut result = _rhs.clone();
            for (a, b) in result.coeffs.iter_mut().zip(&self.coeffs) {
                *a = *a + *b
            }
            result
        };
        result.truncate_leading_zeros();
        result
    }
}

impl Add<Polynomial> for Polynomial {
    type Output = Polynomial;
    fn add(self, _rhs: Polynomial) -> Polynomial {
        &self + &_rhs
    }
}

impl Sub<Polynomial> for Polynomial {
    type Output = Polynomial;
    fn sub(self, _rhs: Polynomial) -> Polynomial {
        let mut result = if self.is_zero() {
            let mut result = _rhs.clone();
            for coeff in &mut result.coeffs {
                *coeff = -(*coeff);
            }
            result
        } else if _rhs.is_zero() {
            self.clone()
        } else if self.degree() >= _rhs.degree() {
            let mut result = self.clone();
            for (a, b) in result.coeffs.iter_mut().zip(&_rhs.coeffs) {
                *a = *a - *b
            }
            result
        } else {
            let mut result = self.clone();
            result.coeffs.resize(_rhs.coeffs.len(), Fp::zero());
            for (a, b) in result.coeffs.iter_mut().zip(&_rhs.coeffs) {
                *a = *a - *b;
            }
            result
        };
        result.truncate_leading_zeros();
        result
    }
}

/// Perform a naive n^2 multiplication of `self` by `_rhs`.
impl Mul<&Polynomial> for &Polynomial {
    type Output = Polynomial;
    fn mul(self, _rhs: &Polynomial) -> Polynomial {
        if self.is_zero() || _rhs.is_zero() {
            Polynomial::zero()
        } else {
            let mut result = vec![Fp::zero(); self.degree() + _rhs.degree() + 1];
            for (i, self_coeff) in self.coeffs.iter().enumerate() {
                for (j, _rhs_coeff) in _rhs.coeffs.iter().enumerate() {
                    result[i + j] = result[i + j] + (*self_coeff * *_rhs_coeff);
                }
            }
            Polynomial::from_fp_coefficients_vec(result)
        }
    }
}

impl Mul<Polynomial> for Polynomial {
    type Output = Polynomial;
    fn mul(self, _rhs: Polynomial) -> Polynomial {
        &self * &_rhs
    }
}

// For shifting to the right, we just shift the coeffs position of the poly.
impl Shr<usize> for &Polynomial {
    type Output = Polynomial;
    fn shr(self, shifter: usize) -> Polynomial {
        let mut result = vec![Fp::zero(); self.degree() - shifter];
        result
            .iter_mut()
            .zip(self.coeffs.iter().skip(shifter))
            .for_each(|(res, shift_coeff)| *res = *shift_coeff);
        Polynomial::from_fp_coefficients_vec(result)
    }
}

// Remember that Shifting-left by 2 is like multiplying by x^2
impl Shl<usize> for Polynomial {
    type Output = Polynomial;
    fn shl(self, shifter: usize) -> Polynomial {
        // Generate x^shifter poly.
        let mut shifter_term = vec![Fp::zero(); shifter + 1];
        shifter_term[shifter] = Fp::one();

        self * Polynomial::from_fp_coefficients_vec(shifter_term)
    }
}

/// Implements long polynomial division algorithm
impl Div<&Polynomial> for &Polynomial {
    type Output = (Polynomial, Polynomial);
    fn div(self, divisor: &Polynomial) -> (Polynomial, Polynomial) {
        if self.is_zero() {
            (Polynomial::zero(), Polynomial::zero())
        } else if divisor.is_zero() {
            panic!("Dividing by zero polynomial")
        } else if self.degree() < divisor.degree() {
            (Polynomial::zero(), self.clone().into())
        } else {
            // Now we know that self.degree() >= divisor.degree();
            let mut quotient = vec![Fp::zero(); self.degree() - divisor.degree() + 1];
            let mut remainder: Polynomial = self.clone();
            // Can unwrap here because we know self is not zero.
            let divisor_leading_inv = divisor.coeffs.last().unwrap().inverse(MODULUS);
            while !remainder.is_zero() && remainder.degree() >= divisor.degree() {
                let cur_q_coeff = *remainder.coeffs.last().unwrap() * divisor_leading_inv;
                let cur_q_degree = remainder.degree() - divisor.degree();
                quotient[cur_q_degree] = cur_q_coeff;

                for (i, div_coeff) in divisor
                    .coeffs
                    .iter()
                    .cloned()
                    .enumerate()
                    .collect::<Vec<(usize, Fp)>>()
                {
                    remainder.coeffs[cur_q_degree + i] =
                        remainder.coeffs[cur_q_degree + i] - (cur_q_coeff * div_coeff);
                }
                while let Some(true) = remainder.coeffs.last().map(|c| c == &Fp::zero()) {
                    remainder.coeffs.pop();
                }
                /*panic!(
                    "quot: {:?}, rem: {:?}, leading_div: {:?}, invers: {:?}",
                    quotient,
                    remainder,
                    divisor.coeffs.last().unwrap(),
                    divisor_leading_inv
                );*/
            }
            (Polynomial::from_fp_coefficients_vec(quotient), remainder)
        }
    }
}

impl Div<Polynomial> for Polynomial {
    type Output = (Polynomial, Polynomial);
    fn div(self, divisor: Polynomial) -> (Polynomial, Polynomial) {
        &self / &divisor
    }
}

use core::borrow::Borrow;
use core::iter::Product;

impl<T> Product<T> for Polynomial
where
    T: Borrow<Polynomial>,
{
    fn product<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Polynomial::from_coefficients_slice(&[1]), |acc, item| {
            acc * item.borrow().clone()
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_poly_addition() {
        let p1 = Polynomial::from_coefficients_slice(&[1, 15, 18, 25, 30]);

        let p2 = Polynomial::from_coefficients_slice(&[30, 15, 2, 0, 25]);

        let expected_res = Polynomial::from_coefficients_slice(&[0, 30, 20, 25, 24]);
        assert!((p1 + p2) == expected_res);
    }

    #[test]
    fn test_poly_subtraction() {
        let p1 = Polynomial::from_coefficients_slice(&[1, 15, 18, 25, 30]);

        let p2 = Polynomial::from_coefficients_slice(&[30, 15, 2, 0, 25]);

        let expected_res = Polynomial::from_coefficients_slice(&[2, 0, 16, 25, 5]);
        assert!((p1 - p2) == expected_res);
    }

    #[test]
    fn test_poly_mul() {
        let p1 = Polynomial::from_coefficients_slice(&[1, 15, 18, 25, 30]);

        let p2 = Polynomial::from_coefficients_slice(&[30, 15, 0, 0, 25]);

        let expected_res = Polynomial::from_coefficients_slice(&[30, 0, 21, 28, 29, 19, 16, 5, 6]);

        assert!(expected_res.degree() == (p1.degree() + p2.degree()));
        assert_eq!(p1 * p2, expected_res);
    }

    #[test]
    fn test_poly_div_by_one() {
        let p1 = Polynomial::from_coefficients_slice(&[1, 15, 18, 25, 30]);
        let p2 = Polynomial::from_coefficients_slice(&[1]);

        assert_eq!((p1.clone() / p2.clone()).0, p1.clone());
        assert_eq!((p1.clone() / p1.clone()).0, p2.clone());
        assert_eq!((p1.clone() / p2.clone()).1, Polynomial::zero());
    }

    #[test]
    fn test_big_poly_div() {
        let p1 = Polynomial::from_coefficients_slice(&[1, 15, 18, 25, 3]);
        let p2 = Polynomial::from_coefficients_slice(&[29, 27, 18, 30]);
        let quot = Polynomial::from_coefficients_slice(&[14, 28]);
        let residue = Polynomial::from_coefficients_slice(&[29, 3, 2]);

        assert_eq!((quot.clone() * p2.clone()) + residue.clone(), p1.clone());

        assert_eq!((p1.clone() / p2.clone()).0, quot);
        assert_eq!((p1.clone() / p2.clone()).1, residue);
    }
}
