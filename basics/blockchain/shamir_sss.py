import random

# ──────────────────────────────────────────────────────────────────────────────
# We operate over a finite field GF(p), where p is prime. All arithmetic
# (addition, multiplication, inversion) is done modulo p. In practice, you pick
# a prime larger than any possible secret. Here we choose the 256-bit secp256k1 prime
# (commonly used in Bitcoin/Ethereum elliptic curves), just as an example.
PRIME = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
# ──────────────────────────────────────────────────────────────────────────────


def _mod_inverse(a: int, p: int = PRIME) -> int:
    """
    Compute the modular inverse of 'a' under prime 'p', i.e., find x such that:
        (a * x) % p == 1.
    We use the Extended Euclidean Algorithm.
    
    Parameters:
    - a: the integer whose inverse we want
    - p: the prime modulus
    
    Returns:
    - The modular inverse of a modulo p.
    """
    # Initialize variables for the extended Euclidean algorithm
    lm, hm = 1, 0           # 'lm' and 'hm' will hold intermediate coefficients
    low, high = a % p, p    # 'low' = a mod p, 'high' = p

    # While we haven’t reduced to gcd = 1
    while low > 1:
        ratio = high // low           # integer division
        nm = hm - lm * ratio          # update for 'hm' vs 'lm'
        new = high - low * ratio      # update gcd pair
        # Shift variables down one step
        hm, lm = lm, nm
        high, low = low, new

    # At this point, 'lm' is the inverse of 'a' mod p (might be negative)
    return lm % p  # ensure positive result


def _eval_polynomial(coeffs: list[int], x: int, p: int = PRIME) -> int:
    """
    Evaluate a polynomial at point x, modulo p.
    The polynomial is specified by 'coeffs', where:
        f(x) = coeffs[0] + coeffs[1]*x + coeffs[2]*x^2 + ... + coeffs[d]*x^d (mod p)
    
    Parameters:
    - coeffs: list of coefficients [a0, a1, a2, ...]
    - x: the x-value at which to evaluate
    - p: the prime modulus
    
    Returns:
    - f(x) mod p
    """
    total = 0
    # Horner’s method is also possible, but this loop is fine for small degree k−1
    for exponent, coeff in enumerate(coeffs):
        # Compute coeff * (x^exponent) mod p, then add to total
        total = (total + coeff * pow(x, exponent, p)) % p
    return total


def generate_shares(secret: int, k: int, n: int, p: int = PRIME) -> list[tuple[int, int]]:
    """
    Split 'secret' into 'n' shares with threshold 'k' using a random polynomial
    of degree (k−1). Return a list of n points (x_i, y_i).
    
    Steps:
      1. Pick random coefficients a1, a2, ..., a_{k-1} uniformly from [0, p-1].
         Let a0 = secret.
         So the polynomial is:
            f(x) = a0 + a1*x + a2*x^2 + ... + a_{k-1}*x^(k-1)   (mod p)
      2. For each i = 1 to n, compute:
            x_i = i
            y_i = f(i) mod p
      3. Return the list: [(1, y1), (2, y2), …, (n, y_n)].
    
    Parameters:
    - secret: the integer you want to split (must be < p)
    - k: threshold number of shares required to reconstruct
    - n: total number of shares to create (n ≥ k)
    - p: prime modulus
    
    Returns:
    - A list of n tuples: (x_i, y_i), with x_i = i and y_i = f(i).
    """
    if k > n:
        raise ValueError("Threshold k cannot be greater than number of shares n")
    
    # Build coefficient list: [a0, a1, a2, ..., a_{k-1}]
    # a0 is the secret; a1..a_{k-1} are random in [0, p−1]
    coeffs = [secret] + [random.randrange(0, p) for _ in range(k - 1)]
    
    shares: list[tuple[int, int]] = []
    # Generate n shares by evaluating at x = 1, 2, ..., n
    for i in range(1, n + 1):
        x_i = i
        y_i = _eval_polynomial(coeffs, x_i, p)  # f(i) mod p
        shares.append((x_i, y_i))
    
    return shares


def reconstruct_secret(shares: list[tuple[int, int]], p: int = PRIME) -> int:
    """
    Given a subset of at least k shares [(x1, y1), (x2, y2), ..., (xk, yk)],
    reconstruct the original secret = f(0) using Lagrange interpolation in GF(p).
    
    The formula for Lagrange interpolation at x = 0 is:
       secret = f(0) = sum_{j=1..k} [ y_j * L_j(0) ]  (mod p),
    where:
       L_j(0) = ∏_{m=1..k, m≠j} (0 − x_m) * inv(x_j − x_m)  (mod p).
    
    Parameters:
    - shares: list of (x_i, y_i) pairs; length must be ≥ threshold k.
    - p: prime modulus
    
    Returns:
    - The reconstructed secret (an integer < p).
    """
    if len(shares) < 1:
        raise ValueError("Need at least one share to attempt reconstruction")
    
    secret = 0
    # Iterate over each share (xj, yj)
    for j, (xj, yj) in enumerate(shares):
        # Compute the Lagrange basis L_j(0) step by step:
        numerator = 1
        denominator = 1
        for m, (xm, _) in enumerate(shares):
            if m == j:
                continue
            # Multiply (0 - x_m) = -xm  into numerator
            numerator = (numerator * (-xm)) % p
            # Multiply (x_j - x_m) into denominator
            denominator = (denominator * (xj - xm)) % p
        
        # Now L_j(0) = numerator * inv(denominator) mod p
        lagrange_coeff = (numerator * _mod_inverse(denominator, p)) % p
        
        # Accumulate yj * L_j(0) into secret
        secret = (secret + (yj * lagrange_coeff)) % p
    
    return secret


# ──────────────────────────────────────────────────────────────────────────────
# Demo / Example Usage
# If you run this file directly (python shamir_sss.py), it will:
# 1) Choose a sample secret
# 2) Generate n shares with a threshold k
# 3) Reconstruct the secret from the first k shares
# 4) Print results and verify correctness
# ──────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    # 1) Define a secret (must be less than PRIME). Here we pick a large integer.
    secret_int = 12345678901234567890

    # 2) Choose threshold k and total number of shares n
    k = 3  # need any 3 shares to reconstruct
    n = 5  # create a total of 5 shares

    # 3) Generate the 5 shares
    shares_list = generate_shares(secret_int, k, n)
    print("Generated Shares (x, y):")
    for s in shares_list:
        print("  ", s)

    # 4) Choose any k = 3 shares (we’ll take the first 3 for demonstration)
    subset = shares_list[:k]

    # 5) Reconstruct the secret using those 3 shares
    reconstructed = reconstruct_secret(subset)
    print(f"\nReconstructed Secret from first {k} shares: {reconstructed}")

    # 6) Verify that it matches original
    assert reconstructed == secret_int, "Reconstruction failed!"
    print("✅ Secret successfully reconstructed!")
