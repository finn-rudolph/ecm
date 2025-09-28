use rug::{Complete, Integer, rand::RandState};

pub trait Curve {}

pub struct WeierstrassCurve {
    pub n: Integer,
    pub a: Integer,
    pub b: Integer,
}

impl Curve for WeierstrassCurve {}
