use ibig::{
    ops::Gcd,
    ubig, UBig,
};
use rand::prelude::*;

fn random_ubig<R>(bits: usize, rng: &mut R) -> UBig
where
    R: Rng + ?Sized,
{
    rng.gen_range(ubig!(1) << (bits - 1)..ubig!(1) << bits)
}

pub(crate) fn egcd(mut a: UBig, mut b: UBig) -> UBig {
    // find common factors of 2
    let shift = (&a | &b).trailing_zeros().unwrap();
    a >>= a.trailing_zeros().unwrap();
    b >>= b.trailing_zeros().unwrap();

    // the binary GCD algorithm
    while a != b {
        if &a > &b {
            a -= &b;
            a >>= a.trailing_zeros().unwrap();
        } else {
            b -= &a;
            b >>= b.trailing_zeros().unwrap();
        }
    }
    a << shift
}

fn bench() {
    let mut s = ubig!(0);
    let mut rng = StdRng::seed_from_u64(3);
    let bits = 10usize.pow(5);
    for _ in 0..4 {
        let a = random_ubig(bits, &mut rng);
        let b = random_ubig(bits, &mut rng);
        // let g = a.gcd(b);
        let g = egcd(a, b);
        s += g;
    }

    println!("sum: {}", s);
}

fn main() {
    bench()
}