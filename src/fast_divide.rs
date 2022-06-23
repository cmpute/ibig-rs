//! Divide by a prearranged Word quickly using multiplication by the reciprocal.

use crate::{
    arch::word::{DoubleWord, Word},
    assert::{assert_in_const_fn, debug_assert_in_const_fn},
    math,
    primitive::{double_word, extend_word, split_double_word, WORD_BITS},
};

/// Divide a Word by a prearranged divisor.
///
/// Granlund, Montgomerry "Division by Invariant Integers using Multiplication"
/// Algorithm 4.1.
#[derive(Clone, Copy)]
pub(crate) struct FastDivideSmall {
    // 2 <= divisor < 2^N, N = WORD_BITS
    divisor: Word,

    // Let n = ceil(log_2(divisor))
    // 2^(n-1) < divisor <= 2^n
    // shift = n - 1
    shift: u32,

    // m = floor(B * 2^n / divisor) + 1 - B, where B = 2^N
    m: Word,
}

impl FastDivideSmall {
    #[inline]
    pub(crate) const fn new(divisor: Word) -> Self {
        assert_in_const_fn(divisor > 1);
        let n = math::ceil_log_2_word(divisor);

        // Calculate:
        // m = floor(B * 2^n / divisor) + 1 - B
        // m >= B + 1 - B >= 1
        // m <= B * 2^n / (2^(n-1) + 1) + 1 - B
        //    = (B * 2^n + 2^(n-1) + 1) / (2^(n-1) + 1) - B
        //    = B * (2^n + 2^(n-1-N) + 2^-N) / (2^(n-1)+1) - B
        //    < B * (2^n + 2^1) / (2^(n-1)+1) - B
        //    = B
        // So m fits in a Word.
        //
        // Note:
        // divisor * (B + m) = divisor * floor(B * 2^n / divisor + 1)
        // = B * 2^n + k, 1 <= k <= divisor

        // m = floor(B * (2^n-1 - (divisor-1)) / divisor) + 1
        let (lo, _hi) = split_double_word(
            double_word(0, math::ones_word(n) - (divisor - 1)) / extend_word(divisor),
        );
        // assert!(_hi == 0);
        FastDivideSmall {
            divisor,
            shift: n - 1,
            m: lo + 1,
        }
    }

    /// ( a / divisor, a % divisor)
    #[inline]
    pub(crate) fn div_rem(&self, a: Word) -> (Word, Word) {
        // q = floor( (B + m) * a / (B * 2^n) )
        //
        // Remember that divisor * (B + m) = B * 2^n + k, 1 <= k <= 2^n
        //
        // (B + m) * a / (B * 2^n)
        // = a / divisor * (B * 2^n + k) / (B * 2^n)
        // = a / divisor + k * a / (divisor * B * 2^n)
        // On one hand, this is >= a / divisor
        // On the other hand, this is:
        // <= a / divisor + 2^n * (B-1) / (2^n * B) / divisor
        // < (a + 1) / divisor
        //
        // Therefore the floor is always the exact quotient.

        // t = m * n / B
        let (_, t) = split_double_word(extend_word(self.m) * extend_word(a));
        // q = (t + a) / 2^n = (t + (a - t)/2) / 2^(n-1)
        let q = (t + ((a - t) >> 1)) >> self.shift;
        let r = a - q * self.divisor;
        (q, r)
    }

    #[inline]
    pub(crate) const fn dummy() -> Self {
        FastDivideSmall {
            divisor: 0,
            shift: 0,
            m: 0,
        }
    }
}

/// Divide a DoubleWord by a prearranged divisor.
///
/// Assumes quotient fits in a Word.
///
/// Möller, Granlund, "Improved division by invariant integers", Algorithm 4.
#[derive(Clone, Copy)]
pub(crate) struct FastDivideNormalized {
    // Top bit must be 1.
    divisor: Word,

    // floor ((B^2 - 1) / divisor) - B, where B = 2^WORD_BITS
    m: Word,
}

impl FastDivideNormalized {
    /// Calculate the inverse m > 0 of a normalized divisor (fit in a word), such that
    /// 
    /// (m + B) * divisor = B^2 - k for some 1 <= k <= divisor
    /// 
    #[inline]
    pub(crate) const fn invert_word(divisor: Word) -> Word {
        let (m, _hi) = split_double_word(DoubleWord::MAX / extend_word(divisor));
        assert_in_const_fn(_hi == 1);
        m
    }

    /// Initialize from a given normalized divisor.
    ///
    /// divisor must have top bit of 1
    #[inline]
    pub(crate) const fn new(divisor: Word) -> Self {
        assert_in_const_fn(divisor.leading_zeros() == 0);
        Self { divisor, m: Self::invert_word(divisor) }
    }

    #[inline]
    pub(crate) const fn div_rem_word(&self, a: Word) -> (Word, Word) {
        if a < self.divisor {
            (0, a)
        } else {
            (1, a - self.divisor)
        }
    }

    /// (a / divisor, a % divisor)
    /// The result must fit in a single word.
    #[inline]
    pub(crate) const fn div_rem(&self, a: DoubleWord) -> (Word, Word) {
        let (a_lo, a_hi) = split_double_word(a);
        debug_assert_in_const_fn!(a_hi < self.divisor);

        // Approximate quotient is (m + B) * a / B^2 ~= (m * a/B + a)/B.
        // This is q1 below.
        // This doesn't overflow because a_hi < Word::MAX.
        let (q0, q1) = split_double_word(extend_word(self.m) * extend_word(a_hi) + a);

        // q = q1 + 1 is our first approximation, but calculate mod B.
        // r = a - q * d
        let q = q1.wrapping_add(1);
        let r = a_lo.wrapping_sub(q.wrapping_mul(self.divisor));

        // Theorem: max(-d, q0+1-B) <= r < max(B-d, q0)
        // Proof:
        // r = a - q * d = a - q1 * d - d
        // = a - (q1 * B + q0 - q0) * d/B - d
        // = a - (m * a_hi + a - q0) * d/B - d
        // = a - ((m+B) * a_hi + a_lo - q0) * d/B - d
        // = a - ((B^2-k)/d * a_hi + a_lo - q0) * d/B - d
        // = a - B * a_hi + (a_hi * k - a_lo * d + q0 * d) / B - d
        // = (a_hi * k + a_lo * (B - d) + q0 * d) / B - d
        //
        // r >= q0 * d / B - d
        // r >= -d
        // r >= d/B (q0 - B) > q0-B
        // r >= max(-d, q0+1-B)
        //
        // r < (d * d + B * (B-d) + q0 * d) / B - d
        // = (B-d)^2 / B + q0 * d / B
        // = (1 - d/B) * (B-d) + (d/B) * q0
        // <= max(B-d, q0)
        // QED

        // if r mod B > q0 { q -= 1; r += d; }
        //
        // Consider two cases:
        // a) r >= 0:
        // Then r = r mod B > q0, hence r < B-d. Adding d will not overflow r.
        // b) r < 0:
        // Then r mod B = r-B > q0, and r >= -d, so adding d will make r non-negative.
        // In either case, this will result in 0 <= r < B.

        // In a branch-free way:
        // decrease = 0xffff.fff = -1 if r mod B > q0, 0 otherwise.
        let (_, decrease) = split_double_word(extend_word(q0).wrapping_sub(extend_word(r)));
        let mut q = q.wrapping_add(decrease);
        let mut r = r.wrapping_add(decrease & self.divisor);

        // At this point 0 <= r < B, i.e. 0 <= r < 2d.
        // the following fix step is unlikely to happen
        if r > self.divisor {
            q += 1;
            r -= self.divisor;
        }

        (q, r)
    }

    #[inline]
    pub(crate) const fn dummy() -> Self {
        FastDivideNormalized { divisor: 0, m: 0 }
    }
}

/// Divide a 3-Word by a prearranged DoubleWord divisor.
///
/// Assumes quotient fits in a Word.
///
/// Möller, Granlund, "Improved division by invariant integers"
/// Algorithm 5.
#[derive(Clone, Copy)]
pub(crate) struct FastDivideNormalized2 {
    // Top bit must be 1.
    divisor: DoubleWord,

    // floor ((B^3 - 1) / divisor) - B, where B = 2^WORD_BITS
    m: Word,
}


impl FastDivideNormalized2 {
    /// Calculate the inverse m > 0 of a normalized divisor (fit in a DoubleWord), such that
    /// 
    /// (m + B) * divisor = B^3 - k for some 1 <= k <= divisor
    /// 
    /// Möller, Granlund, "Improved division by invariant integers", Algorithm 6.
    #[inline]
    pub(crate) const fn invert_double_word(divisor: DoubleWord) -> Word {
        let (d0, d1) = split_double_word(divisor);
        let mut v = FastDivideNormalized::invert_word(d1);
        // then B^2 - d1 <= (B + v)d1 < B^2

        let (mut p, c) = d1.wrapping_mul(v).overflowing_add(d0);
        if c {
            v -= 1;
            if p >= d1 {
                v -= 1;
                p -= d1;
            }
            p = p.wrapping_sub(d1);
        }
        // then B^2 - d1 <= (B + v)d1 + d0 < B^2

        let (t0, t1) = split_double_word(extend_word(v) * extend_word(d0));
        let (p, c) = p.overflowing_add(t1);
        if c {
            v -= 1;
            if double_word(t0, p) >= divisor {
                v -= 1;
            }
        }

        v
    }

    /// Initialize from a given normalized divisor.
    ///
    /// divisor must have top bit of 1
    #[inline]
    pub(crate) const fn new(divisor: DoubleWord) -> Self {
        assert_in_const_fn(divisor.leading_zeros() == 0);
        Self { divisor, m: Self::invert_double_word(divisor) }
    }

    /// The input a is arranged as (lo, mi & hi)
    /// The output is the single word quotient and double word remainder
    pub(crate) const fn div_rem(&self, a: (Word, DoubleWord)) -> (Word, DoubleWord) {
        let (a0, a12) = a;
        let (a1, a2) = split_double_word(a12);
        let (d0, d1) = split_double_word(self.divisor);

        let (q0, q1) = split_double_word(extend_word(self.m) * extend_word(a2) + a12);
        let r1 = a1.wrapping_sub(q1.wrapping_mul(d1));
        let t = extend_word(d0) * extend_word(q1);
        let r = double_word(a0, r1).wrapping_sub(t).wrapping_sub(self.divisor);

        // The first guess of quotient is q1 + 1
        // if r1 >= q0 { r += d; } else { q1 += 1; }
        // In a branch-free way:
        // decrease = 0xffff.fff = -1 if r1 >= q0, 0 otherwise
        let (_, r1) = split_double_word(r);
        let (_, decrease) =
            split_double_word(extend_word(q0).wrapping_sub(extend_word(r1)));
        let mut q1 = q1.wrapping_sub(!decrease);
        let mut r = r.wrapping_add(double_word(decrease, decrease) & self.divisor);

        // the following fix step is unlikely to happen
        if r >= self.divisor {
            q1 += 1;
            r -= self.divisor;
        }

        (q1, r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitive::WORD_BITS;
    use rand::prelude::*;

    #[test]
    fn test_fast_divide_small() {
        let mut rng = StdRng::seed_from_u64(1);
        for _ in 0..1000000 {
            let d_bits = rng.gen_range(2..=WORD_BITS);
            let max_d = math::ones(d_bits);
            let d = rng.gen_range(max_d / 2 + 1..=max_d);
            let fast_div = FastDivideSmall::new(d);
            let n = rng.gen();
            let (q, r) = fast_div.div_rem(n);
            assert_eq!(q, n / d);
            assert_eq!(r, n % d);
        }
    }

    #[test]
    fn test_fast_divide_normalized() {
        let mut rng = StdRng::seed_from_u64(1);
        for _ in 0..1000000 {
            let d = rng.gen_range(Word::MAX / 2 + 1..=Word::MAX);
            let q = rng.gen();
            let r = rng.gen_range(0..d);
            let a = extend_word(q) * extend_word(d) + extend_word(r);
            let fast_div = FastDivideNormalized::new(d);
            assert_eq!(fast_div.div_rem(a), (q, r));
        }
    }

    #[test]
    fn test_fast_divide_normalized2() {
        let mut rng = StdRng::seed_from_u64(1);
        for _ in 0..100000 {
            let d = rng.gen_range(DoubleWord::MAX / 2 + 1..=DoubleWord::MAX);
            let r = rng.gen_range(0..d);
            let q = rng.gen();

            let (d0, d1) = split_double_word(d);
            let (r0, r1) = split_double_word(r);
            let (a0, c) = split_double_word(extend_word(q) * extend_word(d0) + extend_word(r0));
            let a12 = extend_word(q) * extend_word(d1) + extend_word(r1) + extend_word(c);

            let fast_div = FastDivideNormalized2::new(d);
            assert_eq!(fast_div.div_rem((a0, a12)), (q, r), "failed at {:?} / {}", (a0, a12), d);
        }
    }
}
