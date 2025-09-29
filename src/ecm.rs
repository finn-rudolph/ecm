use crate::coords::{MontgomeryPoint, Point, ProjectivePoint};
use crate::sieve;

use rug::Complete;
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
            log::warn!("gcd equals n at the end of stage 1");
        }
        Err(p)
    }
}

fn montgomery_stage_2(
    n: &Integer,
    b1: usize,
    b2: usize,
    d: usize,
    p: MontgomeryPoint,
) -> Option<Integer> {
    let mut multiples_of_2p: Vec<MontgomeryPoint> = Vec::with_capacity(d + 1);
    let mut xz: Vec<Integer> = Vec::with_capacity(d + 1);

    // primes are always odd, so we only need even differences
    multiples_of_2p.push(MontgomeryPoint::origin(p.curve_rc().clone()));
    multiples_of_2p.push(p.mul(2));
    multiples_of_2p.push(multiples_of_2p[1].mul(2));
    xz.push(Integer::from(0));
    xz.push((multiples_of_2p[1].x() * multiples_of_2p[1].z()).complete() % n);
    xz.push((multiples_of_2p[2].x() * multiples_of_2p[2].z()).complete() % n);

    for i in 3..=d {
        multiples_of_2p.push(
            multiples_of_2p[i - 1]
                .add_with_known_difference(&multiples_of_2p[1], &multiples_of_2p[i - 2]),
        );
        xz[i] = (multiples_of_2p[i].x() * multiples_of_2p[i].z()).complete();
    }

    let mut g = Integer::from(1);
    let b = b1 as u64 - 1;
    let t = p.mul(b - 2 * d as u64);
    let r = p.mul(b);

    for r in (b as usize..b2).step_by(2 * d) {}

    if 1 < g && g < *n {
        Some(g)
    } else {
        if g == *n {
            log::warn!("gcd equals n at the end of stage 2")
        }
        None
    }
}
