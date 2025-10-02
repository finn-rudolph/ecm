use crate::coords::{MontgomeryPoint, Point};
use crate::sieve::Sieve;

use rug::Complete;
use rug::Integer;

pub fn ecm(
    n: &Integer,
    b1: usize,
    b2: usize,
    d: usize,
    p: MontgomeryPoint,
    sieve: &Sieve,
) -> Option<Integer> {
    log::info!("using curve {}", p.curve());
    log::info!("using initial point {}", p);

    let q = match stage_1(n, b1, p, sieve) {
        Ok(factor) => return Some(factor),
        Err(point) => point,
    };

    log::debug!("beginning stage 2 from {}", q);

    montgomery_stage_2(n, b1, b2, d, q, sieve)
}

fn stage_1<T: Point>(n: &Integer, b1: usize, mut p: T, sieve: &Sieve) -> Result<Integer, T> {
    let mut g = Integer::from(1);

    for prime in sieve.primes(1, b1) {
        let mut prime_power = *prime;
        while prime_power <= b1 {
            p = p.mul(*prime as u64);
            g = (g * p.z()) % n;
            prime_power *= prime;
        }
    }

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
    sieve: &Sieve,
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
        xz.push((multiples_of_2p[i].x() * multiples_of_2p[i].z()).complete());
    }

    let mut g = Integer::from(1);
    let b = b1 - 1;
    assert!(b & 1 == 1); // so we hit primes
    let mut previous_base = p.mul((b - 2 * d) as u64);
    let mut base = p.mul(b as u64); // r * p, the point we move 2d steps forward

    let primes = sieve.primes(b, b2);
    let mut prime_iter = primes.iter();

    for r in (b..b2).step_by(2 * d) {
        let base_xz = (base.x() * base.z()).complete() % n;

        let mut q = match prime_iter.next() {
            Some(q) => *q,
            None => break,
        };

        while q <= (r + 2 * d) {
            let diff_to_base = (q - r) >> 1;
            g = (g
                * (((base.x() - multiples_of_2p[diff_to_base].x()).complete()
                    * (base.z() + multiples_of_2p[diff_to_base].z()).complete())
                    % n
                    - &base_xz
                    + &xz[diff_to_base]))
                % n;

            q = match prime_iter.next() {
                Some(q) => *q,
                None => break,
            };
        }

        let curr_base = base;
        base = curr_base.add_with_known_difference(&multiples_of_2p[d], &previous_base);
        previous_base = curr_base;
    }

    g = g.gcd(n);

    if 1 < g && g < *n {
        Some(g)
    } else {
        if g == *n {
            log::warn!("gcd equals n at the end of stage 2")
        }
        None
    }
}
