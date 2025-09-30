mod coords;
mod ecm;
mod sieve;

use std::process;

use clap::Parser;
use rug::{Integer, rand::RandState};

#[derive(Parser)]
struct Args {
    #[arg(short, long, help = "integer to be factored")]
    n: Integer,

    #[arg(long, help = "stage 1 smoothness bound")]
    b1: usize,

    #[arg(long, help = "stage 2 smoothness bound")]
    b2: Option<usize>,

    #[arg(short, long, help = "stage 2 step size")]
    d: Option<usize>,
}

impl Args {
    fn validate(&self) {
        if self.b2.is_some() && self.b1 > self.b2.unwrap() {
            eprint!(
                "error: argument for `--b2` must be greater or equal than the argument for `--b2`"
            );
            process::exit(1);
        }
    }

    fn choose_defaults(&mut self) {
        const B2B1_RATIO: usize = 127;

        if self.b2.is_none() {
            self.b2 = Some(self.b1 * B2B1_RATIO);
        }
        if self.d.is_none() {
            self.d = Some(B2B1_RATIO);
        }
    }
}

fn main() {
    let mut args = Args::parse();
    args.validate();
    args.choose_defaults();

    let mut rng = RandState::new();

    match ecm::ecm(
        &args.n,
        args.b1,
        args.b2.unwrap(),
        args.d.unwrap(),
        &mut rng,
    ) {
        Some(factor) => println!("factor found: {factor}"),
        None => println!("no factor found"),
    }
}
