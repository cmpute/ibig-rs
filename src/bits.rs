use crate::{
    buffer::Buffer,
    ibig::IBig,
    primitive::{double_word, Word, WORD_BITS},
    ubig::{Repr::*, UBig},
};
use core::{
    mem,
    ops::{BitAnd, BitAndAssign, BitOr, BitOrAssign, BitXor, BitXorAssign},
};

impl UBig {
    /// Returns true if the `n`-th bit is set.
    ///
    /// # Examples
    ///
    /// ```
    /// # use ibig::prelude::*;
    /// assert_eq!(ubig!(0b10010).bit(1), true);
    /// assert_eq!(ubig!(0b10010).bit(3), false);
    /// assert_eq!(ubig!(0b10010).bit(100), false);
    /// ```
    #[inline]
    pub fn bit(&self, n: usize) -> bool {
        match self.repr() {
            Small(word) => n < WORD_BITS as usize && word & (1 as Word) << n != 0,
            Large(buffer) => {
                let idx = n / WORD_BITS as usize;
                idx < buffer.len() && buffer[idx] & (1 as Word) << (n % WORD_BITS as usize) != 0
            }
        }
    }

    /// Set the `n`-th bit.
    ///
    /// # Examples
    ///
    /// ```
    /// # use ibig::prelude::*;
    /// let mut a = ubig!(0b100);
    /// a.set_bit(0);
    /// assert_eq!(a, ubig!(0b101));
    /// a.set_bit(10);
    /// assert_eq!(a, ubig!(0b10000000101));
    /// ```
    #[inline]
    pub fn set_bit(&mut self, n: usize) {
        match mem::take(self).into_repr() {
            Small(word) => {
                if n < WORD_BITS as usize {
                    *self = UBig::from_word(word | 1 << n)
                } else {
                    *self = UBig::with_bit_word_slow(word, n)
                }
            }
            Large(buffer) => *self = UBig::with_bit_large(buffer, n),
        }
    }

    fn with_bit_word_slow(word: Word, n: usize) -> UBig {
        debug_assert!(n >= WORD_BITS as usize);
        let idx = n / WORD_BITS as usize;
        let mut buffer = Buffer::allocate(idx + 1);
        buffer.push(word);
        buffer.extend((1..idx).map(|_| 0));
        buffer.push(1 << n % WORD_BITS as usize);
        buffer.into()
    }

    fn with_bit_large(mut buffer: Buffer, n: usize) -> UBig {
        let idx = n / WORD_BITS as usize;
        if idx < buffer.len() {
            buffer[idx] |= 1 << n % WORD_BITS as usize;
        } else {
            buffer.ensure_capacity(idx + 1);
            buffer.push_zeros(idx - buffer.len());
            buffer.push(1 << n % WORD_BITS as usize);
        }
        buffer.into()
    }

    /// Clear the `n`-th bit.
    ///
    /// # Examples
    ///
    /// ```
    /// # use ibig::prelude::*;
    /// let mut a = ubig!(0b101);
    /// a.clear_bit(0);
    /// assert_eq!(a, ubig!(0b100));
    /// ```
    #[inline]
    pub fn clear_bit(&mut self, n: usize) {
        match mem::take(self).into_repr() {
            Small(word) => {
                if n < WORD_BITS as usize {
                    *self = UBig::from_word(word & !(1 << n))
                }
            }
            Large(buffer) => *self = UBig::without_bit_large(buffer, n),
        }
    }

    fn without_bit_large(mut buffer: Buffer, n: usize) -> UBig {
        let idx = n / WORD_BITS as usize;
        if idx < buffer.len() {
            buffer[idx] &= !(1 << n % WORD_BITS as usize);
        }
        buffer.into()
    }

    /// Returns the number of trailing zeros in the binary representation.
    ///
    /// In other words, it is the smallest `n` such that 2 to the power of `n` divides the number.
    ///
    /// For 0, it returns `None`.
    ///
    /// # Examples
    ///
    /// ```
    /// # use ibig::prelude::*;
    /// assert_eq!(ubig!(17).trailing_zeros(), Some(0));
    /// assert_eq!(ubig!(48).trailing_zeros(), Some(4));
    /// assert_eq!(ubig!(0b101000000).trailing_zeros(), Some(6));
    /// assert_eq!(ubig!(0).trailing_zeros(), None);
    /// ```
    #[inline]
    pub fn trailing_zeros(&self) -> Option<usize> {
        match self.repr() {
            Small(0) => None,
            Small(word) => Some(word.trailing_zeros() as usize),
            Large(buffer) => Some(UBig::trailing_zeros_large(buffer)),
        }
    }

    fn trailing_zeros_large(words: &[Word]) -> usize {
        debug_assert!(*words.last().unwrap() != 0);

        for (idx, word) in words.iter().enumerate() {
            if *word != 0 {
                return idx * WORD_BITS as usize + word.trailing_zeros() as usize;
            }
        }
        panic!("trailing_zeros_large(0)")
    }

    /// Integer logarithm base 2.
    ///
    /// Returns the floor of the logarithm base 2 of the number.
    /// In other words, it is the position of the highest 1 bit in the binary representation.
    ///
    /// For 0, it returns `None`.
    ///
    /// # Examples
    ///
    /// ```
    /// # use ibig::prelude::*;
    /// assert_eq!(ubig!(17).ilog2(), Some(4));
    /// assert_eq!(ubig!(0b101000000).ilog2(), Some(8));
    /// assert_eq!(ubig!(0).ilog2(), None);
    /// ```
    #[inline]
    pub fn ilog2(&self) -> Option<usize> {
        match self.repr() {
            Small(0) => None,
            Small(word) => Some((WORD_BITS - 1 - word.leading_zeros()) as usize),
            Large(buffer) => Some(
                buffer.len() * WORD_BITS as usize
                    - 1
                    - buffer.last().unwrap().leading_zeros() as usize,
            ),
        }
    }

    /// True if the number is a power of 2.
    ///
    /// # Examples
    ///
    /// ```
    /// # use ibig::prelude::*;
    /// assert_eq!(ubig!(0).is_power_of_two(), false);
    /// assert_eq!(ubig!(8).is_power_of_two(), true);
    /// assert_eq!(ubig!(9).is_power_of_two(), false);
    /// ```
    #[inline]
    pub fn is_power_of_two(&self) -> bool {
        match self.repr() {
            Small(word) => word.is_power_of_two(),
            Large(buffer) => UBig::is_power_of_two_large(buffer),
        }
    }

    fn is_power_of_two_large(words: &[Word]) -> bool {
        debug_assert!(*words.last().unwrap() != 0);

        words[..words.len() - 1].iter().all(|x| *x == 0) && words.last().unwrap().is_power_of_two()
    }
}

impl IBig {
    /// Returns the number of trailing zeros in the two's complement binary representation.
    ///
    /// In other words, it is the smallest `n` such that 2 to the power of `n` divides the number.
    ///
    /// For 0, it returns `None`.
    ///
    /// # Examples
    ///
    /// ```
    /// # use ibig::prelude::*;
    /// assert_eq!(ibig!(17).trailing_zeros(), Some(0));
    /// assert_eq!(ibig!(-48).trailing_zeros(), Some(4));
    /// assert_eq!(ibig!(-0b101000000).trailing_zeros(), Some(6));
    /// assert_eq!(ibig!(0).trailing_zeros(), None);
    /// ```
    #[inline]
    pub fn trailing_zeros(&self) -> Option<usize> {
        self.magnitude().trailing_zeros()
    }
}

/// Next power of two.
///
/// # Examples
/// ```
/// # use ibig::prelude::*;
/// assert_eq!(ubig!(5).next_power_of_two(), ubig!(8));
/// ```
pub trait NextPowerOfTwo {
    type Output;

    fn next_power_of_two(self) -> Self::Output;
}

impl NextPowerOfTwo for UBig {
    type Output = UBig;

    #[inline]
    fn next_power_of_two(self) -> UBig {
        match self.into_repr() {
            Small(word) => match word.checked_next_power_of_two() {
                Some(p) => UBig::from_word(p),
                None => UBig::from(double_word(0, 1)),
            },
            Large(buffer) => UBig::next_power_of_two_large(buffer),
        }
    }
}

impl NextPowerOfTwo for &UBig {
    type Output = UBig;

    #[inline]
    fn next_power_of_two(self) -> UBig {
        self.clone().next_power_of_two()
    }
}

impl UBig {
    fn next_power_of_two_large(mut buffer: Buffer) -> UBig {
        debug_assert!(*buffer.last().unwrap() != 0);

        let n = buffer.len();
        let mut iter = buffer[..n - 1].iter_mut().skip_while(|x| **x == 0);

        let carry = match iter.next() {
            None => 0,
            Some(x) => {
                *x = 0;
                for x in iter {
                    *x = 0;
                }
                1
            }
        };

        let last = buffer.last_mut().unwrap();
        match last
            .checked_add(carry)
            .and_then(|x| x.checked_next_power_of_two())
        {
            Some(p) => *last = p,
            None => {
                *last = 0;
                buffer.ensure_capacity(n + 1);
                buffer.push(1);
            }
        }

        buffer.into()
    }
}

impl BitAnd<UBig> for UBig {
    type Output = UBig;

    #[inline]
    fn bitand(self, rhs: UBig) -> UBig {
        match (self.into_repr(), rhs.into_repr()) {
            (Small(word0), Small(word1)) => UBig::from_word(word0 & word1),
            (Small(word0), Large(buffer1)) => UBig::from_word(word0 & buffer1.first().unwrap()),
            (Large(buffer0), Small(word1)) => UBig::from_word(buffer0.first().unwrap() & word1),
            (Large(buffer0), Large(buffer1)) => {
                if buffer0.len() <= buffer1.len() {
                    UBig::bitand_large(buffer0, &buffer1)
                } else {
                    UBig::bitand_large(buffer1, &buffer0)
                }
            }
        }
    }
}

impl BitAnd<&UBig> for UBig {
    type Output = UBig;

    #[inline]
    fn bitand(self, rhs: &UBig) -> UBig {
        match self.into_repr() {
            Small(word0) => match rhs.repr() {
                Small(word1) => UBig::from_word(word0 & word1),
                Large(buffer1) => UBig::from_word(word0 & buffer1.first().unwrap()),
            },
            Large(buffer0) => match rhs.repr() {
                Small(word1) => UBig::from_word(buffer0.first().unwrap() & word1),
                Large(buffer1) => UBig::bitand_large(buffer0, &buffer1),
            },
        }
    }
}

impl BitAnd<UBig> for &UBig {
    type Output = UBig;

    #[inline]
    fn bitand(self, rhs: UBig) -> UBig {
        rhs.bitand(self)
    }
}

impl BitAnd<&UBig> for &UBig {
    type Output = UBig;

    #[inline]
    fn bitand(self, rhs: &UBig) -> UBig {
        match (self.repr(), rhs.repr()) {
            (Small(word0), Small(word1)) => UBig::from_word(word0 & word1),
            (Small(word0), Large(buffer1)) => UBig::from_word(word0 & buffer1.first().unwrap()),
            (Large(buffer0), Small(word1)) => UBig::from_word(buffer0.first().unwrap() & word1),
            (Large(buffer0), Large(buffer1)) => {
                if buffer0.len() <= buffer1.len() {
                    UBig::bitand_large(buffer0.clone(), buffer1)
                } else {
                    UBig::bitand_large(buffer1.clone(), buffer0)
                }
            }
        }
    }
}

impl BitAndAssign<UBig> for UBig {
    #[inline]
    fn bitand_assign(&mut self, rhs: UBig) {
        *self = mem::take(self) & rhs;
    }
}

impl BitAndAssign<&UBig> for UBig {
    #[inline]
    fn bitand_assign(&mut self, rhs: &UBig) {
        *self = mem::take(self) & rhs;
    }
}

impl UBig {
    fn bitand_large(mut buffer: Buffer, rhs: &[Word]) -> UBig {
        if buffer.len() > rhs.len() {
            buffer.truncate(rhs.len());
        }
        for (x, y) in buffer.iter_mut().zip(rhs.iter()) {
            *x &= *y;
        }
        buffer.into()
    }
}

impl BitOr<UBig> for UBig {
    type Output = UBig;

    #[inline]
    fn bitor(self, rhs: UBig) -> UBig {
        match (self.into_repr(), rhs.into_repr()) {
            (Small(word0), Small(word1)) => UBig::from_word(word0 | word1),
            (Small(word0), Large(buffer1)) => UBig::bitor_large_word(buffer1, word0),
            (Large(buffer0), Small(word1)) => UBig::bitor_large_word(buffer0, word1),
            (Large(buffer0), Large(buffer1)) => {
                if buffer0.len() >= buffer1.len() {
                    UBig::bitor_large(buffer0, &buffer1)
                } else {
                    UBig::bitor_large(buffer1, &buffer0)
                }
            }
        }
    }
}

impl BitOr<&UBig> for UBig {
    type Output = UBig;

    #[inline]
    fn bitor(self, rhs: &UBig) -> UBig {
        match self.into_repr() {
            Small(word0) => match rhs.repr() {
                Small(word1) => UBig::from_word(word0 | word1),
                Large(buffer1) => UBig::bitor_large_word(buffer1.clone(), word0),
            },
            Large(buffer0) => match rhs.repr() {
                Small(word1) => UBig::bitor_large_word(buffer0, *word1),
                Large(buffer1) => UBig::bitor_large(buffer0, buffer1),
            },
        }
    }
}

impl BitOr<UBig> for &UBig {
    type Output = UBig;

    #[inline]
    fn bitor(self, rhs: UBig) -> UBig {
        rhs.bitor(self)
    }
}

impl BitOr<&UBig> for &UBig {
    type Output = UBig;

    #[inline]
    fn bitor(self, rhs: &UBig) -> UBig {
        match (self.repr(), rhs.repr()) {
            (Small(word0), Small(word1)) => UBig::from_word(word0 | word1),
            (Small(word0), Large(buffer1)) => UBig::bitor_large_word(buffer1.clone(), *word0),
            (Large(buffer0), Small(word1)) => UBig::bitor_large_word(buffer0.clone(), *word1),
            (Large(buffer0), Large(buffer1)) => {
                if buffer0.len() >= buffer1.len() {
                    UBig::bitor_large(buffer0.clone(), buffer1)
                } else {
                    UBig::bitor_large(buffer1.clone(), buffer0)
                }
            }
        }
    }
}

impl BitOrAssign<UBig> for UBig {
    #[inline]
    fn bitor_assign(&mut self, rhs: UBig) {
        *self = mem::take(self) | rhs;
    }
}

impl BitOrAssign<&UBig> for UBig {
    #[inline]
    fn bitor_assign(&mut self, rhs: &UBig) {
        *self = mem::take(self) | rhs;
    }
}

impl UBig {
    fn bitor_large_word(mut buffer: Buffer, rhs: Word) -> UBig {
        debug_assert!(buffer.len() >= 2);

        *buffer.first_mut().unwrap() |= rhs;
        buffer.into()
    }

    fn bitor_large(mut buffer: Buffer, rhs: &[Word]) -> UBig {
        for (x, y) in buffer.iter_mut().zip(rhs.iter()) {
            *x |= *y;
        }
        if rhs.len() > buffer.len() {
            buffer.ensure_capacity(rhs.len());
            buffer.extend(&rhs[buffer.len()..]);
        }
        buffer.into()
    }
}

impl BitXor<UBig> for UBig {
    type Output = UBig;

    #[inline]
    fn bitxor(self, rhs: UBig) -> UBig {
        match (self.into_repr(), rhs.into_repr()) {
            (Small(word0), Small(word1)) => UBig::from_word(word0 ^ word1),
            (Small(word0), Large(buffer1)) => UBig::bitxor_large_word(buffer1, word0),
            (Large(buffer0), Small(word1)) => UBig::bitxor_large_word(buffer0, word1),
            (Large(buffer0), Large(buffer1)) => {
                if buffer0.len() >= buffer1.len() {
                    UBig::bitxor_large(buffer0, &buffer1)
                } else {
                    UBig::bitxor_large(buffer1, &buffer0)
                }
            }
        }
    }
}

impl BitXor<&UBig> for UBig {
    type Output = UBig;

    #[inline]
    fn bitxor(self, rhs: &UBig) -> UBig {
        match self.into_repr() {
            Small(word0) => match rhs.repr() {
                Small(word1) => UBig::from_word(word0 ^ word1),
                Large(buffer1) => UBig::bitxor_large_word(buffer1.clone(), word0),
            },
            Large(buffer0) => match rhs.repr() {
                Small(word1) => UBig::bitxor_large_word(buffer0, *word1),
                Large(buffer1) => UBig::bitxor_large(buffer0, buffer1),
            },
        }
    }
}

impl BitXor<UBig> for &UBig {
    type Output = UBig;

    #[inline]
    fn bitxor(self, rhs: UBig) -> UBig {
        rhs.bitxor(self)
    }
}

impl BitXor<&UBig> for &UBig {
    type Output = UBig;

    #[inline]
    fn bitxor(self, rhs: &UBig) -> UBig {
        match (self.repr(), rhs.repr()) {
            (Small(word0), Small(word1)) => UBig::from_word(word0 ^ word1),
            (Small(word0), Large(buffer1)) => UBig::bitxor_large_word(buffer1.clone(), *word0),
            (Large(buffer0), Small(word1)) => UBig::bitxor_large_word(buffer0.clone(), *word1),
            (Large(buffer0), Large(buffer1)) => {
                if buffer0.len() >= buffer1.len() {
                    UBig::bitxor_large(buffer0.clone(), buffer1)
                } else {
                    UBig::bitxor_large(buffer1.clone(), buffer0)
                }
            }
        }
    }
}

impl BitXorAssign<UBig> for UBig {
    #[inline]
    fn bitxor_assign(&mut self, rhs: UBig) {
        *self = mem::take(self) ^ rhs;
    }
}

impl BitXorAssign<&UBig> for UBig {
    #[inline]
    fn bitxor_assign(&mut self, rhs: &UBig) {
        *self = mem::take(self) ^ rhs;
    }
}

impl UBig {
    fn bitxor_large_word(mut buffer: Buffer, rhs: Word) -> UBig {
        debug_assert!(buffer.len() >= 2);

        *buffer.first_mut().unwrap() ^= rhs;
        buffer.into()
    }

    fn bitxor_large(mut buffer: Buffer, rhs: &[Word]) -> UBig {
        for (x, y) in buffer.iter_mut().zip(rhs.iter()) {
            *x ^= *y;
        }
        if rhs.len() > buffer.len() {
            buffer.ensure_capacity(rhs.len());
            buffer.extend(&rhs[buffer.len()..]);
        }
        buffer.into()
    }
}

/// Bitwise AND NOT operation.
///
/// `x.and_not(y)` is equivalent to `x & !y`. For `UBig` the latter is not a valid expression
/// because the `!` operator is not defined.
///
/// # Examples
///
/// ```
/// # use ibig::prelude::*;
/// assert_eq!(ubig!(0xff).and_not(ubig!(0x1111)), ubig!(0xee));
/// ```
pub trait AndNot<Rhs = Self> {
    type Output;

    fn and_not(self, rhs: Rhs) -> Self::Output;
}

impl AndNot<UBig> for UBig {
    type Output = UBig;

    #[inline]
    fn and_not(self, rhs: UBig) -> UBig {
        match (self.into_repr(), rhs.into_repr()) {
            (Small(word0), Small(word1)) => UBig::from_word(word0 & !word1),
            (Small(word0), Large(buffer1)) => UBig::from_word(word0 & !buffer1.first().unwrap()),
            (Large(buffer0), Small(word1)) => UBig::and_not_large_word(buffer0, word1),
            (Large(buffer0), Large(buffer1)) => UBig::and_not_large(buffer0, &buffer1),
        }
    }
}

impl AndNot<&UBig> for UBig {
    type Output = UBig;

    #[inline]
    fn and_not(self, rhs: &UBig) -> UBig {
        match self.into_repr() {
            Small(word0) => match rhs.repr() {
                Small(word1) => UBig::from_word(word0 & !word1),
                Large(buffer1) => UBig::from_word(word0 & !buffer1.first().unwrap()),
            },
            Large(buffer0) => match rhs.repr() {
                Small(word1) => UBig::and_not_large_word(buffer0, *word1),
                Large(buffer1) => UBig::and_not_large(buffer0, buffer1),
            },
        }
    }
}

impl AndNot<UBig> for &UBig {
    type Output = UBig;

    #[inline]
    fn and_not(self, rhs: UBig) -> UBig {
        match self.repr() {
            Small(word0) => match rhs.into_repr() {
                Small(word1) => UBig::from_word(word0 & !word1),
                Large(buffer1) => UBig::from_word(word0 & !buffer1.first().unwrap()),
            },
            Large(buffer0) => match rhs.into_repr() {
                Small(word1) => UBig::and_not_large_word(buffer0.clone(), word1),
                Large(buffer1) => UBig::and_not_large(buffer0.clone(), &buffer1),
            },
        }
    }
}

impl AndNot<&UBig> for &UBig {
    type Output = UBig;

    #[inline]
    fn and_not(self, rhs: &UBig) -> UBig {
        match (self.repr(), rhs.repr()) {
            (Small(word0), Small(word1)) => UBig::from_word(word0 & !word1),
            (Small(word0), Large(buffer1)) => UBig::from_word(word0 & !buffer1.first().unwrap()),
            (Large(buffer0), Small(word1)) => UBig::and_not_large_word(buffer0.clone(), *word1),
            (Large(buffer0), Large(buffer1)) => UBig::and_not_large(buffer0.clone(), buffer1),
        }
    }
}

impl UBig {
    fn and_not_large_word(mut buffer: Buffer, rhs: Word) -> UBig {
        debug_assert!(buffer.len() >= 2);

        *buffer.first_mut().unwrap() &= !rhs;
        buffer.into()
    }

    fn and_not_large(mut buffer: Buffer, rhs: &[Word]) -> UBig {
        for (x, y) in buffer.iter_mut().zip(rhs.iter()) {
            *x &= !*y;
        }
        buffer.into()
    }
}
