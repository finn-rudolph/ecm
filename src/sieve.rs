pub fn primes(upper_bound: usize) -> Vec<u64> {
    let mut is_prime = vec![true; upper_bound + 1];
    let mut primes: Vec<u64> = Vec::new();

    for i in 2..=upper_bound {
        if is_prime[i] {
            primes.push(i as u64);
            for j in (i * i..=upper_bound).step_by(i) {
                is_prime[j] = false;
            }
        }
    }

    primes
}
