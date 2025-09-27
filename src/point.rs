use crate::curve::WeierstrassCurve;
use std::ops::{Add, Mul, Neg, Sub};

use rug::{Complete, Integer};

trait Point<'a, 'c>
where
    Self: 'a + Sized + Mul<u64> + Add + Add<&'a Self> + Sub + Sub<&'a Self> + Neg,
    &'a Self: Mul<u64> + Add + Add<&'a Self> + Sub + Sub<&'a Self> + Neg,
{
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

impl<'c> Mul<u64> for &ProjPoint<'c> {
    type Output = ProjPoint<'c>;

    fn mul(self, n: u64) -> Self::Output {
        if n == 0 {
            return ProjPoint::origin(self.curve);
        } else if n == 2 {
        }

        let mut q = self.clone();
        let m = 3 * n;
        let b = 64 - m.leading_zeros();
        for i in (1..=b - 2).rev() {
            q = q * 2;
            let (m_i, n_i) = ((m >> i) & 1, (n >> i) & 1);
            if (m_i, n_i) == (1, 0) {
                q = q + self;
            } else if (m_i, n_i) == (0, 1) {
                q = q - self;
            }
        }
        q
    }
}

impl<'c> Mul<u64> for ProjPoint<'c> {
    type Output = ProjPoint<'c>;

    fn mul(self, rhs: u64) -> Self::Output {
        (&self) * rhs
    }
}

impl<'c> Add<&ProjPoint<'c>> for &ProjPoint<'c> {
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
            n,
        );

        ProjPoint {
            x,
            y,
            z,
            curve: self.curve,
        }
    }
}

impl<'c> Neg for &ProjPoint<'c> {
    type Output = ProjPoint<'c>;

    fn neg(self) -> Self::Output {
        ProjPoint {
            x: self.x.clone(),
            y: (&self.y).neg().complete(),
            z: self.z.clone(),
            curve: self.curve,
        }
    }
}

fn halve_mod(x: Integer, n: &Integer) -> Integer {
    assert!(n.is_odd());

    if x.is_even() { x >> 1 } else { (x + n) >> 1 }
}

macro_rules! forward_ref_binop {
    (impl $imp:ident, $method:ident for $t:ty, $u:ty) => {
        impl<'a, 'c> $imp<$u> for &'a $t {
            type Output = <$t as $imp<$u>>::Output;

            #[inline]
            fn $method(self, rhs: $u) -> <$t as $imp<$u>>::Output {
                $imp::$method(self, &rhs)
            }
        }

        impl<'a, 'c> $imp<&'a $u> for $t {
            type Output = <$t as $imp<$u>>::Output;

            #[inline]
            fn $method(self, rhs: &'a $u) -> <$t as $imp<$u>>::Output {
                $imp::$method(&self, rhs)
            }
        }

        impl<'c> $imp<$u> for $t {
            type Output = $t;

            #[inline]
            fn $method(self, rhs: $u) -> <$t as $imp<$u>>::Output {
                $imp::$method(&self, &rhs)
            }
        }
    };
}

macro_rules! sub_from_add_neg {
    ($t:ty) => {
        impl<'a, 'c> Sub<&'a $t> for &$t {
            type Output = $t;

            #[inline]
            fn sub(self, rhs: &'a $t) -> Self::Output {
                self + (-rhs)
            }
        }
    };
}

forward_ref_binop! {impl Add, add for ProjPoint<'c>, ProjPoint<'c> }
sub_from_add_neg! {ProjPoint<'c>}
forward_ref_binop! {impl Sub, sub for ProjPoint<'c>, ProjPoint<'c> }
