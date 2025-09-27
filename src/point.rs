use crate::curve::WeierstrassCurve;
use std::ops;

use rug::{Complete, Integer};

trait Point<'c>: ops::Mul<u64, Output = Self> {
    fn add(rhs: Self) -> Self;
    fn double() -> Self;
}

struct ProjPoint<'c> {
    x: Integer,
    y: Integer,
    z: Integer,
    curve: &'c WeierstrassCurve,
}

// impl ops::Mul<u64> for Point {
//     type Output = ProjPoint;
//
//     fn mul(self, rhs: u64) -> Self::Output {}
// }

impl<'c> Point<'c> for ProjPoint<'c> {
    fn add(&self, rhs: &ProjPoint<'c>) -> Self {
        let n = &self.curve.n;

        let x2z1 = (&rhs.x * &self.z).complete();
        let x1z2 = (&self.x * &rhs.z).complete();
        let y2z1 = (&rhs.y * &self.z).complete();
        let y1z2 = (&self.y * &rhs.z).complete();

        let alpha = (&x2z1 - &x1z2).complete() % n;
        let beta = (x2z1 + &x1z2) % n;
        let gamma = (&y2z1 - &y1z2).complete() % n;
        let delta = (y2z1 + &y1z2) % n;

        let zeta = (&self.z * &rhs.z).complete() % n;

        let alpha_sq = alpha.clone().square() % n;
        let gamma_sq = gamma.clone().square() % n;

        let x =
            (&alpha * (((&gamma_sq * &zeta).complete() - (&alpha_sq * &beta).complete()) % n)) % n;

        let alpha_cb = (&alpha_sq * alpha) % n;

        let alpha_cb_delta = (&alpha_cb * delta) % n;
        let z = (alpha_cb * &zeta) % n;
        let y = halve_mod(
            (gamma * ((3 * alpha_sq * beta - gamma_sq * zeta) % n) - alpha_cb_delta) % n,
            &n,
        );

        ProjPoint {
            x,
            y,
            z,
            curve: self.curve,
        }
    }

    fn double() -> Self {
        todo!()
    }
}

fn halve_mod(x: Integer, n: &Integer) -> Integer {
    assert!(n.is_odd());

    if x.is_even() { x >> 1 } else { (x + n) >> 1 }
}
