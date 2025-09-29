mod projective;

use std::fmt::Display;

use rug::{Integer, rand::RandState};

pub use projective::ProjectivePt;

pub trait Curve: Display {}

pub trait Point: Display {
    type CurveType: Curve;

    fn new_curve(n: &Integer, rng: &mut RandState) -> Self;

    fn add(&self, rhs: &Self) -> Self;

    fn sub(&self, rhs: &Self) -> Self;

    fn neg(&self) -> Self;

    fn mul(&self, n: u64) -> Self;

    fn curve(&self) -> &Self::CurveType;

    fn x(&self) -> &Integer;

    fn z(&self) -> &Integer;
}
