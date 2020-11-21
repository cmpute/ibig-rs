//! Conversions of UBig to and from primitive integer types.

use super::{
    UBig,
    allocate::allocate_words,
    word::{bit_size, Word, WORD_BITS},
};
use core::convert::TryFrom;

macro_rules! impl_from_unsigned {
    ($t:ty) => {
        impl From<$t> for UBig {
            fn from(x: $t) -> UBig {
                if bit_size::<$t>() <= WORD_BITS || x <= Word::MAX as $t {
                    UBig::from_word(x as Word)
                } else {
                    let n = (bit_size::<$t>() - x.leading_zeros() as usize + (WORD_BITS - 1))
                        / WORD_BITS;
                    let mut words = allocate_words(n);
                    let mut remaining_bits = x;
                    // Makes the shift non-constant to silence error for smaller bit sizes where
                    // we never reach this loop.
                    let shift = WORD_BITS;
                    for _ in 0..n {
                        words.push(remaining_bits as Word);
                        remaining_bits >>= shift;
                    }
                    debug_assert!(remaining_bits == 0);
                    UBig::from_words_normalized_correct_capacity(words)
                }
            }
        }
    };
}

impl_from_unsigned!(u8);
impl_from_unsigned!(u16);
impl_from_unsigned!(u32);
impl_from_unsigned!(u64);
impl_from_unsigned!(u128);
impl_from_unsigned!(usize);

impl From<bool> for UBig {
    fn from(b: bool) -> UBig {
        u8::from(b).into()
    }
}

impl From<char> for UBig {
    fn from(c: char) -> UBig {
        u32::from(c).into()
    }
}

macro_rules! impl_from_signed {
    ($t:ty as $u:ty) => {
        impl TryFrom<$t> for UBig {
            type Error = <$u as TryFrom<$t>>::Error;

            fn try_from(x: $t) -> Result<UBig, Self::Error> {
                let y = <$u as TryFrom<$t>>::try_from(x)?;
                Ok(y.into())
            }
        }
    };
}

impl_from_signed!(i8 as u8);
impl_from_signed!(i16 as u16);
impl_from_signed!(i32 as u32);
impl_from_signed!(i64 as u64);
impl_from_signed!(i128 as u128);
impl_from_signed!(isize as usize);

#[cfg(test)]
mod tests {
    use super::*;

    // TODO: test for all word sizes by not using from_words
    #[cfg(target_pointer_width = "64")]
    #[test]
    fn test_from_unsigned() {
        assert_eq!(UBig::from(0xf1u8), UBig::from_word(0xf1));
        assert_eq!(UBig::from(0xf123u16), UBig::from_word(0xf123));
        assert_eq!(UBig::from(0xf1234567u32), UBig::from_word(0xf1234567));
        assert_eq!(
            UBig::from(0xf123456701234567u64),
            UBig::from_word(0xf123456701234567)
        );
        assert_eq!(
            UBig::from(0xf1234567012345670123456701234567u128),
            UBig::from_words(vec![0x0123456701234567, 0xf123456701234567])
        );
        assert_eq!(UBig::from(5u128), UBig::from_word(5));
        assert_eq!(UBig::from(5usize), UBig::from_word(5));
    }

    #[test]
    fn test_from_bool() {
        assert_eq!(UBig::from(false), UBig::from(0u8));
        assert_eq!(UBig::from(true), UBig::from(1u8));
    }

    #[test]
    fn test_from_char() {
        assert_eq!(UBig::from('a'), UBig::from(0x61u8));
        assert_eq!(UBig::from('Ł'), UBig::from(0x141u16));
    }

    #[test]
    fn test_from_signed() {
        assert!(UBig::try_from(-5i32).is_err());
        assert_eq!(UBig::try_from(5i32), Ok(UBig::from(5u32)));
        assert_eq!(UBig::try_from(5i128 << 120), Ok(UBig::from(5u128 << 120)));
    }
}
