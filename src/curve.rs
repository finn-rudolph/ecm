use rug::Integer;

pub trait Curve {}

#[derive(Default, Clone)]
pub struct WeierstrassCurve {
    pub n: Integer,
    pub a: Integer,
    pub b: Integer,
}

impl Curve for WeierstrassCurve {}
