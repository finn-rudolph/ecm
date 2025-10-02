use std::{fmt::Display, rc::Rc};

use rug::{Complete, Integer};

use crate::coords::{Curve, Point};

// y^2z = x^3 + cx^2z + xz^2
#[derive(Clone, Debug)]
pub struct MontgomeryCurve {
    n: Integer,
    c: Integer,
}

impl MontgomeryCurve {}

impl Curve for MontgomeryCurve {
    fn n(&self) -> &Integer {
        &self.n
    }
}

impl Display for MontgomeryCurve {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "y^2 = x^3 + {}x^2 + x (Montgomery form)", &self.c)
    }
}

#[derive(Clone, Debug)]
pub struct MontgomeryPoint {
    x: Integer,
    z: Integer,
    curve: Rc<MontgomeryCurve>,
}

impl MontgomeryPoint {
    pub fn add_with_known_difference(
        &self,
        rhs: &MontgomeryPoint,
        difference: &MontgomeryPoint,
    ) -> MontgomeryPoint {
        let n = &self.curve.n;

        let x1x2_minus_z1z2 = ((&self.x * &rhs.x).complete() - (&self.z * &rhs.z)) % n;
        let x = (difference.z() * (x1x2_minus_z1z2.square() % n)) % n;
        let x1z2_minus_x2z1 = ((&self.x * &rhs.z).complete() - (&self.z * &rhs.x)) % n;
        let z = (difference.x() * (x1z2_minus_x2z1.square() % n)) % n;

        MontgomeryPoint {
            x,
            z,
            curve: self.curve.clone(),
        }
    }

    pub fn from_sigma(n: &Integer, sigma: &Integer) -> Option<Self> {
        if *sigma == 0 || *sigma == 1 || *sigma == 2 {
            return None;
        }

        let u: Integer = (sigma.square_ref().complete() - 5) % n;
        let v: Integer = (sigma << 4u32).complete();
        let u_cb: Integer = u.pow_mod_ref(&Integer::from(3), n).unwrap().complete();
        let u_cb_v = (&u_cb * &v).complete();

        let g = n.gcd_ref(&u_cb_v).complete();
        if g > 1 {
            log::warn!(
                "discovered nontrivial divisor of n during curve construction: {}",
                g
            );
            return None;
        }

        let u_cb_v_inv = (u_cb_v << 2u32).invert(n).unwrap();
        let c = ((((&v - &u).complete().pow_mod(&Integer::from(3), n).unwrap() * (3 * u + &v))
            % n)
            * u_cb_v_inv
            - 2)
            % n;

        let curve = Rc::new(MontgomeryCurve { n: n.clone(), c });

        Some(MontgomeryPoint {
            x: u_cb,
            z: v.pow_mod_ref(&Integer::from(3), n).unwrap().complete(),
            curve,
        })
    }
}

impl Point for MontgomeryPoint {
    type CurveType = MontgomeryCurve;

    fn origin(curve: Rc<Self::CurveType>) -> Self {
        MontgomeryPoint {
            x: Integer::from(0),
            z: Integer::from(0),
            curve,
        }
    }

    // The returned point may not lie on y^2 = x^3 + cx^2 + x, but instead on a twist.
    // We do not explicitly keep track of the twist.
    fn random(n: &Integer, rng: &mut rug::rand::RandState) -> Self {
        loop {
            let sigma = n.random_below_ref(rng).complete();
            if let Some(point) = MontgomeryPoint::from_sigma(n, &sigma) {
                log::debug!("σ = {}", sigma);
                return point;
            }
        }
    }

    fn mul(&self, k: u64) -> Self {
        if k == 0 {
            return MontgomeryPoint {
                x: Integer::from(0),
                z: Integer::from(0),
                curve: self.curve.clone(),
            };
        } else if k == 1 {
            return self.clone();
        } else if k == 2 {
            let n = self.n();

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

        let b = 64 - k.leading_zeros();
        for i in (0..=b - 2).rev() {
            if (k >> i) & 1 == 1 {
                u = t.add_with_known_difference(&u, self);
                t = t.mul(2);
            } else {
                t = u.add_with_known_difference(&t, self);
                u = u.mul(2);
            }
        }

        if k & 1 == 1 {
            u.add_with_known_difference(&t, self)
        } else {
            u.mul(2)
        }
    }

    fn curve_rc(&self) -> &Rc<Self::CurveType> {
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

#[cfg(test)]
mod test {
    use rug::{integer::IntegerExt64, rand::RandState};

    use crate::coords::{ProjectivePoint, projective::WeierstrassCurve};

    use super::*;

    fn points_equal(p: &MontgomeryPoint, q: &ProjectivePoint) -> bool {
        todo!()
    }

    fn sqrt_mod_p(a: &Integer, p: &Integer) -> Integer {
        assert!(p.mod_u64(4) == 3);
        todo!()
    }

    // fn eq(&self, rhs: &Self) -> bool {
    //     // TODO: this currently only works for fields.
    //     if !Rc::ptr_eq(self.curve_rc(), rhs.curve_rc()) {
    //         false
    //     } else {
    //         let n = self.n();
    //         (&self.x * &rhs.z).complete() % n == (&self.z * &rhs.x).complete() % n
    //             && (&self.y * &rhs.z).complete() % n == (&self.z * &rhs.y).complete() % n
    //     }
    // }

    #[test]
    fn test_mul_2() {
        // NOTE: over composite moduli, need the implementations be equivalent?
        // (different addition chains)
        let mut n = Integer::from(u128::MAX).next_prime();
        let mut rng = RandState::new();

        for _ in 0..42 {
            // we have to compute a finite field sqrt to obtain the y coordinate
            while n.mod_u64(4) == 1 {
                n = (n + 1u32).next_prime();
            }

            for _ in 0..17 {
                let (p, y_sq) = loop {
                    let p = MontgomeryPoint::random(&n, &mut rng);
                    if p.z().is_zero() {
                        continue;
                    }
                    let y_sq =
                        (p.x() * (1u32 + (p.x() * (&p.curve().c + p.x()).complete()) % &n)) % &n;

                    if y_sq.legendre(&n) != 1i32 {
                        break (p, y_sq);
                    }
                };

                let y = y_sq.pow_mod(&((&n + 1u32).complete() >> 2), &n).unwrap();

                // x -> x - c/3
                // y -> y

                let c_sq = p.curve().c.square_ref().complete();
                let one_third = Integer::from(3).invert(&n).unwrap();
                let minus_c_over_3 = (&c_sq * -one_third) % &n;
                let a = Integer::from(1) + &minus_c_over_3;
                let mut b = minus_c_over_3.clone();
                let c_sq_over_9 = minus_c_over_3.square_ref().complete();
                b += &c_sq_over_9;
                b += c_sq_over_9 * &minus_c_over_3;

                let proj_p =
                    ProjectivePoint::new(n.clone(), a, b, p.x() - minus_c_over_3, y, 1.into());
                // assert!(points_equal(&p.mul(2), &p_proj.mul(2)));
            }
            n = (n + 1u32).next_prime();
        }
    }
}
