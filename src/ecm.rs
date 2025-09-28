use crate::{curve::WeierstrassCurve, curve_selection};
use crate::point::ProjPoint;
use rug::{Complete, Integer, rand::RandState};

fn random_curve<'c>(n: &Integer, rng: &mut RandState, curve: &mut WeierstrassCurve) -> Point<'c> {}

fn ecm(n: &Integer, b1: u64, rng: &mut RandState) {
    let mut curve = WeierstrassCurve::new();
    curve_selection::random()

    for u in primes(b1) {
        let mut v = u;
        while v <= b1 {
            p *= u;
            v *= u;
        }
    }
}
