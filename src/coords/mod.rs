mod montgomery;
mod projective;

use std::{fmt::Display, rc::Rc};

use rug::{Integer, rand::RandState};

pub trait Curve: Display {}

pub trait Point: Display + Clone {
    type CurveType: Curve;

    fn origin(curve: Rc<Self::CurveType>) -> Self;

    fn new_curve(n: &Integer, rng: &mut RandState) -> Self;

    fn mul(&self, n: u64) -> Self;

    fn curve(&self) -> &Self::CurveType {
        self.curve_rc()
    }

    fn curve_rc(&self) -> &Rc<Self::CurveType>;

    fn x(&self) -> &Integer;

    fn z(&self) -> &Integer;
}

pub use montgomery::MontgomeryPoint;
pub use projective::ProjectivePoint;
