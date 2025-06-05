# Blockchain Task: Shamir’s Secret Sharing (Polynomic Vault)

This folder contains a Python implementation of Shamir’s Secret Sharing Scheme,
which splits a secret into multiple shares so that only a threshold number of shares
can reconstruct the original secret.

## Files

- **shamir_sss.py**  
  - Implements all core logic for SSS, including:
    1. `generate_shares(secret: int, k: int, n: int)`  
       • Creates a random polynomial of degree (k−1) where the constant term is `secret`.  
       • Evaluates that polynomial at x = 1, 2, …, n to produce n shares `(x, y)`.  
       • Returns a list of n `(x_i, y_i)` pairs.
    2. `reconstruct_secret(shares: List[(x, y)])`  
       • Uses Lagrange interpolation in GF(p) to recover `f(0)`, which is the secret.  
    3. Demo code under `if __name__ == "__main__":`  
       • Shows generating 5 shares (3-of-5 scheme) and reconstructing the secret from any 3 of them.

## How to Run Locally

1. **Clone the repository** (if not already done):
   ```bash
   git clone https://github.com/<your-username>/chokkajosephabhishek_spider_task1.git
