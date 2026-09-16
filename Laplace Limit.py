"""Module for mathematical computation and analysis."""

import argparse

import mpmath


def compute_laplace_limit(target_digits: int) -> mpmath.mpf:
    """
    Computes the Laplace Limit using the Newton-Raphson method.
    The Laplace Limit L is the root of the equation:
        x * exp(sqrt(1 + x^2)) = 1 + sqrt(1 + x^2)
        
    We define f(x) = x * exp(sqrt(1 + x^2)) - sqrt(1 + x^2) - 1.
    And its derivative f'(x):
        Let y = sqrt(1 + x^2)
        f'(x) = exp(y) + (x^2 * exp(y)) / y - x / y
    
    To maintain an O(1) memory footprint and optimize calculation speed, we
    employ dynamic precision-stepping. We start at a low precision and double
    it at each Newton step because Newton's method converges quadratically.
    This limits CPU time spent on highly inaccurate early guesses.
    """
    # Safety margin for precision
    max_dps = target_digits + 50
    current_dps = 50
    
    mpmath.mp.dps = current_dps
    
    # Initial guess for the Laplace limit (~0.6627434193)
    x = mpmath.mpf("0.66274341934918158097")
    
    # Quadratic convergence: double the working precision each step
    while current_dps < max_dps:
        current_dps = min(current_dps * 2, max_dps)
        mpmath.mp.dps = current_dps
        
        y = mpmath.sqrt(1 + x*x)
        ey = mpmath.exp(y)
        
        f = x * ey - y - 1
        df = ey + (x * x * ey) / y - x / y
        
        x = x - f / df
        
    # Perform one final iteration at the maximum precision to guarantee accuracy
    y = mpmath.sqrt(1 + x*x)
    ey = mpmath.exp(y)
    
    f = x * ey - y - 1
    df = ey + (x * x * ey) / y - x / y
    
    x = x - f / df
    
    return x


def main():
    """Entry point — parse arguments and run the main computation.
    
    """
    parser = argparse.ArgumentParser(description="Calculate the Laplace Limit to N digits.")
    parser.add_argument("-n", "--digits", type=int, default=1000000,
                        help="Number of digits to compute (default: 1000)")
    args = parser.parse_args()
    
    target_digits = args.digits
    print(f"[*] Calculating Laplace Limit to {target_digits} digits...")
    
    x = compute_laplace_limit(target_digits)
    
    # Format the string with extra digits to allow for safe string truncation
    mpmath.mp.dps = target_digits + 50
    
    # Convert to string, ensuring fixed-point notation
    raw_str = mpmath.nstr(x, target_digits + 20, min_fixed=-mpmath.inf, max_fixed=mpmath.inf)
    
    if "." in raw_str:
        integer_part, fractional_part = raw_str.split(".")
    else:
        integer_part, fractional_part = raw_str, ""
        
    # Strict truncation at target_digits (compliant with OEIS)
    truncated_fraction = fractional_part[:target_digits]
    final_value = f"{integer_part}.{truncated_fraction}"
    
    # 1) Save continuous raw digits
    out_txt = f"Laplace_Limit_{target_digits}_digits.txt"
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(final_value)
    print(f"[*] Saved exact truncated value to {out_txt}")
    
    # 2) Save b-file for OEIS (A033259)
    b_file = "b033259.txt"
    with open(b_file, "w", encoding="utf-8") as f:
        for i, digit in enumerate(truncated_fraction, start=1):
            f.write(f"{i} {digit}\n")
    print(f"[*] Saved OEIS b-file to {b_file}")

if __name__ == "__main__":
    main()
