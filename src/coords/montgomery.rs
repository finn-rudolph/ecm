use std::{fmt::Display, rc::Rc};

use rug::{Complete, Integer};

use crate::coords::{Curve, Point};

// a = 1, b = 0 implicitly
struct MontgomeryCurve {
    n: Integer,
    c: Integer,
}

impl Curve for MontgomeryCurve {}

impl Display for MontgomeryCurve {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "y^2 = x^3 + {}x^2 + x (Montgomery form)", &self.c)
    }
}

#[derive(Clone)]
struct MontgomeryPoint {
    x: Integer,
    z: Integer,
    curve: Rc<MontgomeryCurve>,
}

impl MontgomeryPoint {
    fn add_with_diff(&self, rhs: &MontgomeryPoint, diff: &MontgomeryPoint) -> MontgomeryPoint {
        let n = &self.curve.n;

        let x1x2_minus_z1z2 = ((&self.x * &rhs.x).complete() - (&self.z * &rhs.z)) % n;
        let x = (diff.z() * (x1x2_minus_z1z2.square() % n)) % n;
        let x1z2_minus_x2z1 = ((&self.x * &rhs.z).complete() - (&self.z * &rhs.x)) % n;
        let z = (diff.x() * (x1z2_minus_x2z1.square() % n)) % n;

        MontgomeryPoint {
            x,
            z,
            curve: self.curve.clone(),
        }
    }
}

impl Point for MontgomeryPoint {
    type CurveType = MontgomeryCurve;

    fn new_curve(n: &Integer, rng: &mut rug::rand::RandState) -> Self {
        todo!()
    }

    fn mul(&self, n: u64) -> Self {
        if n == 0 {
            return MontgomeryPoint {
                x: Integer::from(0),
                z: Integer::from(0),
                curve: self.curve.clone(),
            };
        } else if n == 1 {
            return self.clone();
        } else if n == 2 {
            let x1_sq = self.x.square_ref().complete();
            let z1_sq = self.z.square_ref().complete();
            let x1_z1 = (&self.x * &self.z).complete() % n;

            let z_intermed = ((&x1_sq + (&self.curve.c * &x1_z1)).complete() + &z1_sq) % n;
            let z = ((x1_z1 * z_intermed) << 2) % n;
            let x = ((x1_sq - z1_sq) % n).square() % n;

            return MontgomeryPoint {
                x,
                z,
                curve: self.curve.clone(),
            };
        }

        let mut u = self.clone();
        let mut t = self.mul(2);

        let b = 64 - n.leading_zeros();
        for i in (0..=b - 2).rev() {
            if (n >> i) & 1 == 1 {
                u = t.add_with_diff(&u, self);
                t = t.mul(2);
            } else {
                t = u.add_with_diff(&t, self);
                u = u.mul(2);
            }
        }

        if n & 1 == 1 {
            u.add_with_diff(&t, self)
        } else {
            u.mul(2)
        }
    }

    fn curve(&self) -> &Self::CurveType {
        &self.curve
    }

    fn x(&self) -> &Integer {
        &self.x
    }

    fn z(&self) -> &Integer {
        &self.z
    }
}

impl Display for MontgomeryPoint {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[{} : ? : {}]", &self.x, &self.z)
    }
}
