mod curve;
mod curve_selection;
mod ecm;
mod point;
mod sieve;

use clap::Parser;
use rug::{Integer, rand::RandState};

#[derive(Parser)]
struct Args {
    #[arg(short, long, help = "integer to be factored")]
    n: Integer,

    #[arg(long, help = "stage 1 smoothness bound")]
    b1: usize,

    #[arg(long, help = "stage 2 smoothness bound")]
    b2: usize,
}

fn main() {
    let args = Args::parse();

    let mut rng = RandState::new();
    match ecm::ecm(&args.n, args.b1, args.b2, &mut rng) {
        Some(factor) => println!("factor found: {factor}"),
        None => println!("no factor found"),
    }
}
