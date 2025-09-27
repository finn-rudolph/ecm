use crate::curve::WeierstrassCurve;
use std::{ops, process::Output};

use rug::{Complete, Integer};

trait Point<'a, 'c>
where
    Self: 'a
        + Sized
        + ops::Mul<u64>
        + ops::MulAssign<u64>
        + ops::Add
        + ops::Add<&'a Self>
        + ops::Sub
        + ops::Sub<&'a Self>
        + ops::AddAssign,
    &'a Self: ops::Mul<u64> + ops::Add + ops::Add<&'a Self> + ops::Sub + ops::Sub<&'a Self>,
{
    fn double(&self) -> Self;
}

#[derive(Clone)]
struct ProjPoint<'c> {
    x: Integer,
    y: Integer,
    z: Integer,
    curve: &'c WeierstrassCurve,
}

impl<'c> ProjPoint<'c> {
    fn origin(curve: &'c WeierstrassCurve) -> ProjPoint<'c> {
        ProjPoint {
            x: Integer::from(0),
            y: Integer::from(1),
            z: Integer::from(0),
            curve,
        }
    }
}

impl<'c> ops::Mul<u64> for &ProjPoint<'c> {
    type Output = ProjPoint<'c>;

    fn mul(self, n: u64) -> Self::Output {
        if n == 0 {
            return ProjPoint::origin(self.curve);
        }

        let mut q = self.clone();
        let mut m = 3 * n;
        let b = 64 - m.leading_zeros();
        for i in (1..=b - 2).rev() {
            q = q.double();
            let (m_i, n_i) = ((m >> i) & 1, (n >> i) & 1);
            if (m_i, n_i) == (1, 0) {
                q += self;
            } else if (m_i, n_i) == (0, 1) {
                q -= self;
            }
        }
        q
    }
}

impl<'c> ops::Mul<u64> for ProjPoint<'c> {
    type Output = ProjPoint<'c>;

    fn mul(self, rhs: u64) -> Self::Output {
        (&self) * rhs
    }
}

impl<'a, 'c> ops::Add<&'a ProjPoint<'c>> for &ProjPoint<'c> {
    type Output = ProjPoint<'c>;

    fn add(self, rhs: &ProjPoint<'c>) -> Self::Output {
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
}

fn halve_mod(x: Integer, n: &Integer) -> Integer {
    assert!(n.is_odd());

    if x.is_even() { x >> 1 } else { (x + n) >> 1 }
}
