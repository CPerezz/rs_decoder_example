use crate::ops::mod_inv;
use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct Fp(pub u8);

pub const MODULUS: u16 = 31;

impl Add<Fp> for Fp {
    type Output = Fp;
    fn add(self, _rhs: Fp) -> Self {
        Fp(((self.0 as u16 + _rhs.0 as u16) % MODULUS) as u8)
    }
}

impl Sub<Fp> for Fp {
    type Output = Fp;
    fn sub(self, _rhs: Fp) -> Self {
        self + Fp(MODULUS as u8 - _rhs.0)
    }
}

impl Mul<Fp> for Fp {
    type Output = Fp;
    fn mul(self, _rhs: Fp) -> Self {
        Fp(((self.0 as i16 * _rhs.0 as i16) % MODULUS as i16) as u8)
    }
}

impl Div<Fp> for Fp {
    type Output = Fp;
    fn div(self, _rhs: Fp) -> Self {
        self * Fp(mod_inv(_rhs.0 as i16, MODULUS as i16))
    }
}

impl Neg for Fp {
    type Output = Fp;
    fn neg(self) -> Fp {
        Fp(MODULUS as u8 - self.0)
    }
}

impl From<u8> for Fp {
    fn from(num: u8) -> Self {
        Fp(num)
    }
}

use core::borrow::Borrow;
use core::iter::Sum;

impl<T> Sum<T> for Fp
where
    T: Borrow<Fp>,
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Self::zero(), |acc, item| acc + *item.borrow())
    }
}

impl Fp {
    pub fn zero() -> Fp {
        Fp(0u8)
    }

    pub fn one() -> Fp {
        Fp(1u8)
    }

    pub fn inverse(&self, modulus: Fp) -> Fp {
        let modulus = modulus.0 as i16;
        let mut mn = (modulus as i16, self.0 as i16);
        let mut xy = (0, 1);
        while mn.1 != 0 {
            xy = (xy.1, xy.0 - (mn.0 / mn.1) * xy.1);
            mn = (mn.1, mn.0 % mn.1);
        }
        while xy.0 < 0 {
            xy.0 += modulus;
        }
        // This is guaranteed to be correct since
        // the result is always positive and since
        // we work with u8, a i16 can be collapsed into
        // it.
        Fp(xy.0 as u8)
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn field_addition() {
        let zero = Fp::zero();
        let a = Fp::from(30);
        let b = Fp::from(1);
        let c = Fp::from(27);

        assert_eq!(zero + b, b);
        assert_eq!(a + b, zero);
        assert_eq!(a + c, 26.into())
    }

    #[test]
    fn field_subctraction() {
        let zero = Fp::zero();
        let a = Fp::from(30);
        let b = Fp::from(1);
        let c = Fp::from(27);

        assert_eq!(zero - b, a);
        assert_eq!(a - b, 29.into());
        assert_eq!(a - c, 3.into());
        assert_eq!(c - a, 28.into());
    }

    #[test]
    fn field_mul() {
        let zero = Fp::zero();
        let a = Fp::from(30);
        let b = Fp::from(1);
        let c = Fp::from(27);

        assert_eq!(zero * b, zero);
        assert_eq!(a * b, a);
        assert_eq!(a * c, 4.into());
        assert_eq!(c * 25.into(), 24.into());
    }

    #[test]
    fn field_div() {
        let zero = Fp::zero();
        let a = Fp::from(30);
        let b = Fp::from(1);
        let c = Fp::from(27);

        assert_eq!(zero / b, zero);
        assert_eq!(a / b, a);
        assert_eq!(a / c, 8.into());
        assert_eq!(c / 25.into(), 11.into());
    }
}
