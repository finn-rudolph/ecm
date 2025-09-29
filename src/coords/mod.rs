mod montgomery;
mod projective;

use std::fmt::Display;

use rug::{Integer, rand::RandState};

pub trait Curve: Display {}

pub trait Point: Display + Clone {
    type CurveType: Curve;

    fn new_curve(n: &Integer, rng: &mut RandState) -> Self;

    fn mul(&self, n: u64) -> Self;

    fn curve(&self) -> &Self::CurveType;

    fn x(&self) -> &Integer;

    fn z(&self) -> &Integer;
}

pub use montgomery::MontgomeryPoint;
pub use projective::ProjectivePoint;
