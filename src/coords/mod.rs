mod montgomery;
mod projective;

use std::{fmt::Display, rc::Rc};

use rug::{Integer, rand::RandState};

pub trait Curve: Display {
    fn n(&self) -> &Integer;
}

pub trait Point: Display + Clone {
    type CurveType: Curve;

    fn origin(curve: Rc<Self::CurveType>) -> Self;

    fn random_curve(n: &Integer, rng: &mut RandState) -> Self;

    fn mul(&self, k: u64) -> Self;

    fn curve(&self) -> &Self::CurveType {
        self.curve_rc()
    }

    fn curve_rc(&self) -> &Rc<Self::CurveType>;

    fn n(&self) -> &Integer {
        self.curve().n()
    }

    fn x(&self) -> &Integer;

    fn z(&self) -> &Integer;
}

pub use montgomery::MontgomeryPoint;
pub use projective::ProjectivePoint;
