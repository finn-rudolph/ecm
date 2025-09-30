use crate::coords::{Curve, Point};
use std::{fmt::Display, ops::Neg, rc::Rc};

use rug::{Complete, Integer};

// y^2z = x^3 + axz^2 + bz^3
#[derive(Clone, Debug)]
pub struct WeierstrassCurve {
    pub n: Integer,
    pub a: Integer,
    pub b: Integer,
}

impl Curve for WeierstrassCurve {
    fn n(&self) -> &Integer {
        &self.n
    }
}

impl Display for WeierstrassCurve {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "y^2 = x^3 + {}x + {} (Weierstrass form)",
            &self.a, &self.b
        )
    }
}

#[derive(Clone, Debug)]
pub struct ProjectivePoint {
    x: Integer,
    y: Integer,
    z: Integer,
    curve: Rc<WeierstrassCurve>,
}

impl ProjectivePoint {
    fn add(&self, rhs: &Self) -> Self {
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

        ProjectivePoint {
            x,
            y,
            z,
            curve: self.curve.clone(),
        }
    }

    fn neg(&self) -> Self {
        ProjectivePoint {
            x: self.x.clone(),
            y: (&self.y).neg().complete(),
            z: self.z.clone(),
            curve: self.curve.clone(),
        }
    }

    fn sub(&self, rhs: &Self) -> Self {
        self.add(&rhs.neg())
    }
}

impl Display for ProjectivePoint {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[{} : {} : {}]", &self.x, &self.y, &self.z)
    }
}

impl PartialEq for ProjectivePoint {
    fn eq(&self, rhs: &Self) -> bool {
        // TODO: this currently only works for fields.
        if !Rc::ptr_eq(self.curve_rc(), rhs.curve_rc()) {
            false
        } else {
            let n = self.n();
            (&self.x * &rhs.z).complete() % n == (&self.z * &rhs.x).complete() % n
                && (&self.y * &rhs.z).complete() % n == (&self.z * &rhs.y).complete() % n
        }
    }
}

impl Eq for ProjectivePoint {}

impl Point for ProjectivePoint {
    type CurveType = WeierstrassCurve;

    fn origin(curve: Rc<WeierstrassCurve>) -> ProjectivePoint {
        ProjectivePoint {
            x: Integer::from(0),
            y: Integer::from(1),
            z: Integer::from(0),
            curve,
        }
    }

    fn new_curve(n: &Integer, rng: &mut rug::rand::RandState) -> Self {
        loop {
            let x = Integer::random_below_ref(n, rng).complete();
            let y = Integer::random_below_ref(n, rng).complete();
            let a = Integer::random_below_ref(n, rng).complete();
            let b = (y.clone().square() - &x * (x.clone().square() + &a)) % n;

            let discriminant: Integer = ((a.clone().square() * &a) << 2) + 27 * b.clone().square();

            if discriminant.gcd(n) == 1 {
                let curve = Rc::new(WeierstrassCurve { n: n.clone(), a, b });
                return ProjectivePoint {
                    x,
                    y,
                    z: Integer::from(1),
                    curve,
                };
            }
        }
    }

    fn mul(&self, k: u64) -> Self {
        if k == 0 {
            return ProjectivePoint::origin(self.curve.clone());
        } else if k == 2 {
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

            return ProjectivePoint {
                x,
                y,
                z,
                curve: self.curve.clone(),
            };
        }

        let mut q = self.clone();
        let m = 3 * k;
        let b = 64 - m.leading_zeros();
        for i in (1..=b - 2).rev() {
            q = q.mul(2);
            let (m_i, k_i) = ((m >> i) & 1, (k >> i) & 1);
            if (m_i, k_i) == (1, 0) {
                q = q.add(self);
            } else if (m_i, k_i) == (0, 1) {
                q = q.sub(self);
            }
        }
        q
    }

    fn x(&self) -> &Integer {
        &self.x
    }

    fn z(&self) -> &Integer {
        &self.z
    }

    fn curve_rc(&self) -> &Rc<WeierstrassCurve> {
        &self.curve
    }
}

fn halve_mod(x: Integer, n: &Integer) -> Integer {
    assert!(n.is_odd());

    if x.is_even() { x >> 1 } else { (x + n) >> 1 }
}
