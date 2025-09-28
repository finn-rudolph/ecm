use crate::curve::{Curve, WeierstrassCurve};
use std::ops::{Add, Mul, Neg, Sub};

use rug::{Complete, Integer};

pub trait Point<'a, 'c>
where
    for<'b> Self: 'a + Sized + Mul<u64> + Add + Add<&'b Self> + Sub + Sub<&'b Self> + Neg,
    for<'b> &'a Self:
        Mul<u64> + Add<&'b Self> + Add<&'b Self> + Sub<&'b Self> + Sub<&'b Self> + Neg,
{
    type CurveType: Curve;

    fn curve(&self) -> &'c Self::CurveType;
}

#[derive(Clone)]
pub struct ProjPoint<'c> {
    pub x: Integer,
    pub y: Integer,
    pub z: Integer,
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

impl<'a, 'c> Point<'a, 'c> for ProjPoint<'c>
where
    Self: 'a,
{
    type CurveType = WeierstrassCurve;

    fn curve(&self) -> &'c WeierstrassCurve {
        self.curve
    }
}

impl<'c> Mul<u64> for &ProjPoint<'c> {
    type Output = ProjPoint<'c>;

    fn mul(self, n: u64) -> Self::Output {
        if n == 0 {
            return ProjPoint::origin(self.curve);
        } else if n == 2 {
            let n = &self.curve.n;

            let lambda: Integer = ((&self.x * &self.y).complete() << 1) % n;
            let nu: Integer = ((&self.y * &self.z).complete() << 1) % n;
            let mu: Integer =
                (3 * self.x.clone().square() + &self.curve.a * self.z.clone().square()) % n;

            let mu_sq = mu.clone().square();
            let nu_sq = nu.clone().square() % n;
            let lambda_nu: Integer = lambda * &nu;
            let two_lambda_nu: Integer = Complete::complete(&lambda_nu << 1);

            let y_lhs = mu * ((&two_lambda_nu + lambda_nu - &mu_sq) % n);
            let x: Integer = (&nu * (&mu_sq - two_lambda_nu)) % n;
            let z = &nu_sq * nu;
            let y_rhs = (nu_sq * &self.y) % n;
            let y = (y_lhs - ((y_rhs * &self.y) << 1)) % n;

            return ProjPoint {
                x,
                y,
                z,
                curve: self.curve,
            };
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

macro_rules! forward_lhs_ref_binop {
    (impl $imp:ident, $method:ident for $t:ty, $u:ty) => {
        impl<'c> $imp<$u> for $t {
            type Output = $t;

            #[inline]
            fn $method(self, rhs: $u) -> <$t as $imp<$u>>::Output {
                $imp::$method(&self, rhs)
            }
        }
    };
}

macro_rules! forward_ref_binop {
    (impl $imp:ident, $method:ident for $t:ty, $u:ty) => {
        impl<'a, 'c> $imp<&'a $u> for $t {
            type Output = <$t as $imp<$u>>::Output;

            #[inline]
            fn $method(self, rhs: &'a $u) -> <$t as $imp<$u>>::Output {
                $imp::$method(&self, rhs)
            }
        }

        impl<'a, 'c> $imp<$u> for &'a $t {
            type Output = <$t as $imp<$u>>::Output;

            #[inline]
            fn $method(self, rhs: $u) -> <$t as $imp<$u>>::Output {
                $imp::$method(self, &rhs)
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

macro_rules! forward_ref_unop {
    (impl $imp: ident, $method:ident for $t:ty) => {
        impl<'c> $imp for $t {
            type Output = $t;

            #[inline]
            fn $method(self) -> Self::Output {
                $imp::$method(&self)
            }
        }
    };
}

macro_rules! gen_sub_from_add_neg {
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

forward_lhs_ref_binop! {impl Mul, mul for ProjPoint<'c>, u64}
forward_ref_binop! {impl Add, add for ProjPoint<'c>, ProjPoint<'c> }
gen_sub_from_add_neg! {ProjPoint<'c>}
forward_ref_binop! {impl Sub, sub for ProjPoint<'c>, ProjPoint<'c> }
forward_ref_unop! {impl Neg, neg for ProjPoint<'c>}
