//! Operators for finding greatest common divisor.

use crate::{
    arch::word::Word,
    buffer::Buffer,
    div, gcd,
    ibig::IBig,
    memory::MemoryAllocation,
    ubig::{Repr::*, UBig},
};

impl UBig {
    /// Compute the greatest common divisor between self and the other operand
    ///
    /// # Example
    /// ```
    /// # use ibig::ubig;
    /// assert_eq!(ubig!(12).gcd(&ubig!(18)), ubig!(6));
    /// ```
    ///
    /// Panics if two oprands are both zero.
    #[inline]
    pub fn gcd(&self, rhs: &UBig) -> UBig {
        match (self.repr(), rhs.repr()) {
            (Small(word0), Small(word1)) => UBig::gcd_word(*word0, *word1),
            (Small(word0), Large(buffer1)) => UBig::gcd_large_word(buffer1, *word0),
            (Large(buffer0), Small(word1)) => UBig::gcd_large_word(buffer0, *word1),
            (Large(buffer0), Large(buffer1)) => UBig::gcd_large(buffer0.clone(), buffer1.clone()),
        }
    }

    /// Compute the greatest common divisor between self and the other operand, and return
    /// both the common divisor `g` and the Bézout coefficients.
    ///
    /// # Example
    /// ```
    /// # use ibig::{ibig, ubig};
    /// assert_eq!(ubig!(12).extended_gcd(&ubig!(18)), (ubig!(6), ibig!(-1), ibig!(1)));
    /// ```
    ///
    /// # Panics
    /// 
    /// Panics if two operands are both zero.
    #[inline]
    pub fn extended_gcd(&self, rhs: &UBig) -> (UBig, IBig, IBig) {
        match (self.clone().into_repr(), rhs.clone().into_repr()) {
            (Small(word0), Small(word1)) => UBig::extended_gcd_word(word0, word1),
            (Large(buffer0), Small(word1)) => UBig::extended_gcd_large_word(buffer0, word1),
            (Small(word0), Large(buffer1)) => {
                let (g, s, t) = UBig::extended_gcd_large_word(buffer1, word0);
                (g, t, s)
            }
            (Large(buffer0), Large(buffer1)) => UBig::extended_gcd_large(buffer0, buffer1),
        }
    }

    /// Perform gcd on two `Word`s.
    #[inline]
    fn gcd_word(a: Word, b: Word) -> UBig {
        if a == 0 || b == 0 {
            if a == 0 && b == 0 {
                UBig::panic_gcd_with_0();
            }
            return UBig::from_word(a | b);
        }

        UBig::from_word(gcd::gcd_word_by_word(a, b))
    }

    /// Perform extended gcd on two `Word`s.
    #[inline]
    fn extended_gcd_word(a: Word, b: Word) -> (UBig, IBig, IBig) {
        if a == 0 && b == 0 {
            UBig::panic_gcd_with_0();
        }

        // cases where a = 0 or b = 0 can be correctly handled by extended euclidean algorithm
        let (r, s, t) = gcd::xgcd_word_by_word(a, b, false);
        (UBig::from_word(r), IBig::from(s), IBig::from(t))
    }

    /// Perform gcd on a large number with a `Word`.
    #[inline]
    fn gcd_large_word(buffer: &Buffer, rhs: Word) -> UBig {
        if rhs == 0 {
            return buffer.clone().into();
        }

        // reduce the large number
        let small = div::rem_by_word(buffer, rhs);
        if small == 0 {
            return UBig::from_word(rhs);
        }

        UBig::from_word(gcd::gcd_word_by_word(small, rhs))
    }

    /// Perform extended gcd on a large number with a `Word`.
    #[inline]
    fn extended_gcd_large_word(mut buffer: Buffer, rhs: Word) -> (UBig, IBig, IBig) {
        if rhs == 0 {
            return (buffer.into(), IBig::from(1u8), IBig::from(0u8));
        }

        // reduce the large number
        let rem = div::div_by_word_in_place(&mut buffer, rhs);
        if rem == 0 {
            return (UBig::from_word(rhs), IBig::from(0u8), IBig::from(1u8));
        }

        let (r, s, t) = gcd::xgcd_word_by_word(rhs, rem, false);
        let new_t = -t * IBig::from(UBig::from(buffer)) + s;
        (UBig::from_word(r), IBig::from(t), new_t)
    }

    /// Perform gcd on two large numbers.
    #[inline]
    fn gcd_large(mut lhs: Buffer, mut rhs: Buffer) -> UBig {
        let len = gcd::gcd_in_place(&mut lhs, &mut rhs);
        lhs.truncate(len);
        lhs.into()
    }

    /// Perform extended gcd on two large numbers.
    #[inline]
    fn extended_gcd_large(mut lhs: Buffer, mut rhs: Buffer) -> (UBig, IBig, IBig) {
        let res_len = lhs.len().min(rhs.len());
        let mut buffer = Buffer::allocate(res_len);
        buffer.push_zeros(res_len);

        let mut allocation =
            MemoryAllocation::new(gcd::memory_requirement_exact(lhs.len(), rhs.len()));
        let mut memory = allocation.memory();

        let (lhs_sign, rhs_sign) =
            gcd::xgcd_in_place(&mut lhs, &mut rhs, &mut buffer, false, &mut memory);
        (
            buffer.into(),
            IBig::from_sign_magnitude(lhs_sign, rhs.into()),
            IBig::from_sign_magnitude(rhs_sign, lhs.into()),
        )
    }

    #[inline]
    fn panic_gcd_with_0() -> ! {
        panic!("the greatest common divisor is not defined between zeros!")
    }
}
