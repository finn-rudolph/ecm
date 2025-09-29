use crate::coords::{MontgomeryPoint, Point, ProjectivePoint};
use crate::sieve;

use rug::{Integer, rand::RandState};

pub fn ecm(n: &Integer, b1: usize, b2: usize, rng: &mut RandState) -> Option<Integer> {
    assert!(b1 < b2);

    let p = ProjectivePoint::new_curve(n, rng);

    log::info!("using curve {}", p.curve());
    log::info!("using point {}", p);

    let q = match stage_1(n, b1, p) {
        Ok(factor) => return Some(factor),
        Err(point) => point,
    };

    None
}

fn stage_1<T: Point>(n: &Integer, b1: usize, mut p: T) -> Result<Integer, T> {
    let mut g = Integer::from(1);

    for prime in sieve::primes(b1) {
        let mut prime_power = prime;
        while prime_power <= b1 as u64 {
            p = p.mul(prime);
            g = (g * p.z()) % n;
            prime_power *= prime;
        }
    }

    log::info!("stage 1 finished");

    g = g.gcd(n);
    if 1 < g && g < *n {
        Ok(g)
    } else {
        if g == *n {
            log::warn!("gcd after stage 1 equals n");
        }
        Err(p)
    }
}

fn montgomery_stage_2(n: &Integer, b2: usize, mut p: MontgomeryPoint) -> Option<Integer> {}
