#!/usr/bin/env python3
"""
Black-Scholes Reference Tests using scipy.

This script generates reference values for testing the Move implementation
against the industry-standard Black-Scholes formulas in scipy.

Usage:
    python test_black_scholes.py

All values use WAD scaling (10^18) to match the Move implementation.
"""

import numpy as np
from scipy.stats import norm
from dataclasses import dataclass
from typing import Tuple

# WAD scale factor: 10^18
WAD = 10**18


@dataclass
class BlackScholesResult:
    """Black-Scholes calculation results."""
    d1: float
    d2: float
    call_price: float
    put_price: float
    delta_call: float
    delta_put: float
    gamma: float
    vega: float
    theta_call: float
    theta_put: float
    rho_call: float
    rho_put: float


def black_scholes(
    S: float,
    K: float,
    T: float,
    r: float,
    sigma: float
) -> BlackScholesResult:
    """
    Calculate Black-Scholes option pricing and Greeks.
    
    Args:
        S: Spot price
        K: Strike price
        T: Time to expiry (years)
        r: Risk-free rate
        sigma: Volatility
        
    Returns:
        BlackScholesResult with all values
    """
    # d1 and d2
    d1 = (np.log(S / K) + (r + sigma**2 / 2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    
    # Option prices
    call = S * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)
    put = K * np.exp(-r * T) * norm.cdf(-d2) - S * norm.cdf(-d1)
    
    # Greeks
    delta_call = norm.cdf(d1)
    delta_put = norm.cdf(d1) - 1
    
    gamma = norm.pdf(d1) / (S * sigma * np.sqrt(T))
    vega = S * norm.pdf(d1) * np.sqrt(T)
    
    theta_call = (
        -S * norm.pdf(d1) * sigma / (2 * np.sqrt(T))
        - r * K * np.exp(-r * T) * norm.cdf(d2)
    )
    theta_put = (
        -S * norm.pdf(d1) * sigma / (2 * np.sqrt(T))
        + r * K * np.exp(-r * T) * norm.cdf(-d2)
    )
    
    rho_call = K * T * np.exp(-r * T) * norm.cdf(d2)
    rho_put = -K * T * np.exp(-r * T) * norm.cdf(-d2)
    
    return BlackScholesResult(
        d1=d1,
        d2=d2,
        call_price=call,
        put_price=put,
        delta_call=delta_call,
        delta_put=delta_put,
        gamma=gamma,
        vega=vega,
        theta_call=theta_call,
        theta_put=theta_put,
        rho_call=rho_call,
        rho_put=rho_put,
    )


def to_wad(value: float) -> int:
    """Convert float to WAD-scaled integer."""
    return int(value * WAD)


def from_wad(value: int) -> float:
    """Convert WAD-scaled integer to float."""
    return value / WAD


def test_atm_option():
    """Test ATM option: S = K = 100, T = 1, r = 5%, σ = 20%."""
    print("=" * 60)
    print("TEST: ATM Option (S=100, K=100, T=1, r=5%, σ=20%)")
    print("=" * 60)
    
    result = black_scholes(S=100, K=100, T=1, r=0.05, sigma=0.2)
    
    print(f"\nd1 = {result.d1:.10f}")
    print(f"d2 = {result.d2:.10f}")
    print(f"\nCall Price = {result.call_price:.10f}")
    print(f"Put Price  = {result.put_price:.10f}")
    print(f"\nPut-Call Parity Check:")
    print(f"  C - P = {result.call_price - result.put_price:.10f}")
    print(f"  S - K*e^(-rT) = {100 - 100 * np.exp(-0.05):.10f}")
    
    print(f"\nGreeks:")
    print(f"  Delta (call) = {result.delta_call:.10f}")
    print(f"  Delta (put)  = {result.delta_put:.10f}")
    print(f"  Gamma        = {result.gamma:.10f}")
    print(f"  Vega         = {result.vega:.10f}")
    print(f"  Theta (call) = {result.theta_call:.10f}")
    print(f"  Theta (put)  = {result.theta_put:.10f}")
    print(f"  Rho (call)   = {result.rho_call:.10f}")
    print(f"  Rho (put)    = {result.rho_put:.10f}")
    
    print(f"\nWAD-scaled values for Move tests:")
    print(f"  d1 = {to_wad(result.d1)}")
    print(f"  d2 = {to_wad(result.d2)}")
    print(f"  Call = {to_wad(result.call_price)}")
    print(f"  Put  = {to_wad(result.put_price)}")
    print(f"  Delta (call) = {to_wad(result.delta_call)}")
    print(f"  Gamma = {to_wad(result.gamma)}")
    print(f"  Vega = {to_wad(result.vega)}")
    
    return result


def test_itm_call():
    """Test ITM call: S = 120, K = 100."""
    print("\n" + "=" * 60)
    print("TEST: ITM Call (S=120, K=100, T=1, r=5%, σ=20%)")
    print("=" * 60)
    
    result = black_scholes(S=120, K=100, T=1, r=0.05, sigma=0.2)
    
    print(f"\nd1 = {result.d1:.10f}")
    print(f"d2 = {result.d2:.10f}")
    print(f"\nCall Price = {result.call_price:.10f}")
    print(f"Put Price  = {result.put_price:.10f}")
    print(f"Intrinsic Value = {max(120 - 100, 0)}")
    print(f"Time Value = {result.call_price - max(120 - 100, 0):.10f}")
    
    return result


def test_otm_call():
    """Test OTM call: S = 80, K = 100."""
    print("\n" + "=" * 60)
    print("TEST: OTM Call (S=80, K=100, T=1, r=5%, σ=20%)")
    print("=" * 60)
    
    result = black_scholes(S=80, K=100, T=1, r=0.05, sigma=0.2)
    
    print(f"\nd1 = {result.d1:.10f}")
    print(f"d2 = {result.d2:.10f}")
    print(f"\nCall Price = {result.call_price:.10f}")
    print(f"Put Price  = {result.put_price:.10f}")
    
    return result


def test_deep_itm():
    """Test deep ITM call: S = 200, K = 100."""
    print("\n" + "=" * 60)
    print("TEST: Deep ITM Call (S=200, K=100, T=1, r=5%, σ=20%)")
    print("=" * 60)
    
    result = black_scholes(S=200, K=100, T=1, r=0.05, sigma=0.2)
    
    print(f"\nCall Price = {result.call_price:.10f}")
    print(f"Forward Diff = {200 - 100 * np.exp(-0.05):.10f}")
    print(f"Delta (call) = {result.delta_call:.10f} (should be ~1)")
    
    return result


def test_deep_otm():
    """Test deep OTM call: S = 50, K = 100."""
    print("\n" + "=" * 60)
    print("TEST: Deep OTM Call (S=50, K=100, T=1, r=5%, σ=20%)")
    print("=" * 60)
    
    result = black_scholes(S=50, K=100, T=1, r=0.05, sigma=0.2)
    
    print(f"\nCall Price = {result.call_price:.10f}")
    print(f"Delta (call) = {result.delta_call:.10f} (should be ~0)")
    
    return result


def test_short_expiry():
    """Test short expiry: T = 0.25 (3 months)."""
    print("\n" + "=" * 60)
    print("TEST: Short Expiry (S=100, K=100, T=0.25, r=5%, σ=20%)")
    print("=" * 60)
    
    result = black_scholes(S=100, K=100, T=0.25, r=0.05, sigma=0.2)
    
    print(f"\nCall Price = {result.call_price:.10f}")
    print(f"Put Price  = {result.put_price:.10f}")
    print(f"Theta (call) = {result.theta_call:.10f} (should be more negative)")
    
    return result


def test_high_volatility():
    """Test high volatility: σ = 50%."""
    print("\n" + "=" * 60)
    print("TEST: High Volatility (S=100, K=100, T=1, r=5%, σ=50%)")
    print("=" * 60)
    
    result = black_scholes(S=100, K=100, T=1, r=0.05, sigma=0.5)
    
    print(f"\nCall Price = {result.call_price:.10f}")
    print(f"Put Price  = {result.put_price:.10f}")
    print(f"Vega = {result.vega:.10f} (should be higher)")
    
    return result


def verify_put_call_parity(S: float, K: float, T: float, r: float, sigma: float) -> bool:
    """Verify put-call parity: C - P = S - K*e^(-rT)."""
    result = black_scholes(S, K, T, r, sigma)
    
    lhs = result.call_price - result.put_price
    rhs = S - K * np.exp(-r * T)
    
    diff = abs(lhs - rhs)
    tolerance = 1e-10
    
    return diff < tolerance


def generate_move_test_vectors():
    """Generate test vectors for Move unit tests."""
    print("\n" + "=" * 60)
    print("MOVE TEST VECTORS")
    print("=" * 60)
    
    test_cases = [
        # (S, K, T, r, sigma, description)
        (100, 100, 1.0, 0.05, 0.2, "ATM standard"),
        (120, 100, 1.0, 0.05, 0.2, "ITM call"),
        (80, 100, 1.0, 0.05, 0.2, "OTM call"),
        (100, 100, 0.25, 0.05, 0.2, "Short expiry"),
        (100, 100, 1.0, 0.05, 0.4, "Higher vol"),
        (100, 100, 1.0, 0.10, 0.2, "Higher rate"),
    ]
    
    print("\nconst SCALE: u256 = 1_000_000_000_000_000_000;\n")
    
    for S, K, T, r, sigma, desc in test_cases:
        result = black_scholes(S, K, T, r, sigma)
        
        print(f"// {desc}: S={S}, K={K}, T={T}, r={r*100}%, σ={sigma*100}%")
        print(f"// Call = {result.call_price:.6f}, Put = {result.put_price:.6f}")
        print(f"let spot = {int(S)} * SCALE;")
        print(f"let strike = {int(K)} * SCALE;")
        print(f"let time = {to_wad(T)}; // {T}")
        print(f"let rate = {to_wad(r)}; // {r*100}%")
        print(f"let vol = {to_wad(sigma)}; // {sigma*100}%")
        print(f"// Expected call: {to_wad(result.call_price)}")
        print(f"// Expected put:  {to_wad(result.put_price)}")
        print()


def main():
    """Run all tests."""
    print("Black-Scholes Reference Implementation")
    print("Using scipy.stats.norm for CDF/PDF calculations")
    print()
    
    # Run test cases
    test_atm_option()
    test_itm_call()
    test_otm_call()
    test_deep_itm()
    test_deep_otm()
    test_short_expiry()
    test_high_volatility()
    
    # Verify put-call parity for various cases
    print("\n" + "=" * 60)
    print("PUT-CALL PARITY VERIFICATION")
    print("=" * 60)
    
    test_cases = [
        (100, 100, 1, 0.05, 0.2),
        (120, 100, 1, 0.05, 0.2),
        (80, 100, 1, 0.05, 0.2),
        (100, 100, 0.25, 0.05, 0.2),
        (100, 100, 1, 0.10, 0.3),
    ]
    
    for S, K, T, r, sigma in test_cases:
        parity_ok = verify_put_call_parity(S, K, T, r, sigma)
        status = "✓" if parity_ok else "✗"
        print(f"{status} S={S}, K={K}, T={T}, r={r}, σ={sigma}")
    
    # Generate Move test vectors
    generate_move_test_vectors()
    
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print("All reference calculations completed.")
    print("Use the WAD-scaled values above for Move unit tests.")
    print("Tolerance for Move tests: < 0.1% for prices, < 0.5% for Greeks")


if __name__ == "__main__":
    main()
