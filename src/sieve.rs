pub struct Sieve {
    primes: Vec<usize>,
}

impl Sieve {
    pub fn new(upper_bound: usize) -> Self {
        let mut is_prime = vec![true; upper_bound + 1];
        let mut primes: Vec<usize> = Vec::new();

        for i in 2..=upper_bound {
            if is_prime[i] {
                primes.push(i);
                for j in (i * i..=upper_bound).step_by(i) {
                    is_prime[j] = false;
                }
            }
        }

        Sieve { primes }
    }

    pub fn primes(&self, lower_bound: usize, upper_bound: usize) -> &[usize] {
        assert!(lower_bound < upper_bound);

        let bs_begin = self.primes.binary_search(&lower_bound);
        let begin = if let Ok(i) = bs_begin {
            i
        } else {
            bs_begin.unwrap_err()
        };

        let bs_end = self.primes.binary_search(&(upper_bound + 1));
        let end = if let Ok(i) = bs_end {
            i
        } else {
            bs_end.unwrap_err()
        };

        &self.primes[begin..end]
    }
}
