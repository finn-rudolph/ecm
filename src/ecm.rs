use crate::point::Point;
use crate::sieve;

use crate::{curve::WeierstrassCurve, curve_selection};
use rug::Complete;
use rug::{Integer, rand::RandState};

fn ecm(n: &Integer, b1: usize, rng: &mut RandState) -> Option<Integer> {
    let mut curve = WeierstrassCurve::default();
    let mut point = curve_selection::weierstrass(n, rng, &mut curve);

    for u in sieve::primes(b1) {
        let mut v = u;
        while v <= b1 as u64 {
            point = point.mul(u);
            let g = n.gcd_ref(&point.z).complete();
            if 1 < g && &g < n {
                return Some(g);
            }
            v *= u;
        }
    }

    None
}
