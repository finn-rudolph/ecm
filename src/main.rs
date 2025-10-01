mod coords;
mod ecm;
mod sieve;

use std::env;
use std::process;
use std::time::SystemTime;
use std::time::UNIX_EPOCH;

use clap::Parser;
use rug::{Integer, rand::RandState};

use crate::coords::MontgomeryPoint;
use crate::coords::Point;

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

    #[arg(long, help = "curve selection paramter")]
    sigma: Option<Integer>,
}

impl Args {
    fn validate(&self) {
        if self.b2.is_some() && self.b1 > self.b2.unwrap() {
            eprintln!(
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
    if env::var("RUST_LOG").is_err() {
        unsafe { env::set_var("RUST_LOG", "info") }
    }
    env_logger::init();

    let mut args = Args::parse();
    args.validate();
    args.choose_defaults();

    if args.sigma.is_some() {
        match ecm::ecm(
            &args.n,
            args.b1,
            args.b2.unwrap(),
            args.d.unwrap(),
            match MontgomeryPoint::new_curve(&args.n, &args.sigma.unwrap()) {
                Some(point) => point,
                None => {
                    eprintln!("error: invalid curve parameter σ");
                    process::exit(1);
                }
            },
        ) {
            Some(factor) => println!("factor found: {factor}"),
            None => println!("no factor found"),
        }
        return;
    }

    let mut rng = RandState::new();
    rng.seed(
        &SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
            .into(),
    );

    match ecm::ecm(
        &args.n,
        args.b1,
        args.b2.unwrap(),
        args.d.unwrap(),
        MontgomeryPoint::random_curve(&args.n, &mut rng),
    ) {
        Some(factor) => println!("factor found: {factor}"),
        None => println!("no factor found"),
    }
}
