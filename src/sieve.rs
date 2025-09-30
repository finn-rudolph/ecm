pub fn primes(lower_bound: usize, upper_bound: usize) -> Vec<usize> {
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

    let bs_result = primes.binary_search(&lower_bound);
    let begin = bs_result.unwrap_or(bs_result.unwrap_err());

    primes.into_iter().skip(begin).collect()
}
