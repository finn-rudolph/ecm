use std::fmt::Display;

use rug::Integer;

pub trait Curve {}

#[derive(Default, Clone)]
pub struct WeierstrassCurve {
    pub n: Integer,
    pub a: Integer,
    pub b: Integer,
}

impl Curve for WeierstrassCurve {}

impl Display for WeierstrassCurve {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "curve y^2 = x^3 + {} * x + {} mod {} (Weierstrass format)",
            &self.a, &self.b, &self.n
        )
    }
}
