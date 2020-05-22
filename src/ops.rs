//! Modular operation helpers

pub(crate) fn mod_inv(a: i16, module: i16) -> u8 {
    let mut mn = (module, a);
    let mut xy = (0, 1);
    while mn.1 != 0 {
        xy = (xy.1, xy.0 - (mn.0 / mn.1) * xy.1);
        mn = (mn.1, mn.0 % mn.1);
    }
    while xy.0 < 0 {
        xy.0 += module;
    }
    // This is guaranteed to be correct since
    // the result is always positive and since
    // we work with u8, a i16 can be collapsed into
    // it.
    xy.0 as u8
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mod_inv_test() {
        assert_eq!(mod_inv(30i16, 31i16), 30u8);
    }
}
