use crate::curve::WeierstrassCurve;
use crate::point::ProjPoint;
use rug::{Complete, Integer, rand::RandState};

fn random_curve<'c>(n: &Integer) -> Point<'c> {}

fn ecm(n: &Integer, b1: u64, rng: &mut RandState) {
    let curve = WeierstrassCurve::new();
    loop {
        let x = Integer::random_below_ref(n, rng).complete();
        let y = Integer::random_below_ref(n, rng).complete();
        let a = Integer::random_below_ref(n, rng).complete();
        let b = (y.clone().square() - &x * (x.clone().square() + &a)) % n;

        let discriminant: Integer = ((a.clone().square() * &a) << 2) + 27 * b.clone().square();

        if discriminant.gcd(n) == 1 {
            let curve = WeierstrassCurve { n: n.clone(), a, b };
            let p = ProjPoint {
                x,
                y,
                z: Integer::from(1),
                curve: &curve,
            };
            break;
        }
    }

    for u in primes(b1) {
        let mut v = u;
        while v <= b1 {
            p *= u;
            v *= u;
        }
    }
}
