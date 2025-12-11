# Black-Scholes Option Pricing Library

A Sui Move implementation of the Black-Scholes option pricing model for European options.

> ⚠️ **Experimental**: This is a learning/research project. It has not been audited and should not be used in production without thorough review.

## Core Dependency

This package is built on top of the [**move-gaussian**](https://github.com/Evan-Kim2028/move-gaussian) library, which provides:
- Normal distribution CDF/PDF for option pricing
- Transcendental functions (ln, exp, sqrt) for Black-Scholes calculations
- SignedWad fixed-point arithmetic for negative values (d₁, d₂, Greeks)

| Gaussian Library | |
|------------------|---|
| **Repository** | [github.com/Evan-Kim2028/move-gaussian](https://github.com/Evan-Kim2028/move-gaussian) |
| **Package (Testnet)** | [`0x66f9087a3d9ae3fe07a5f3c1475d503f1b0ea508d3b83b73b0b8637b57629f7f`](https://suiscan.xyz/testnet/object/0x66f9087a3d9ae3fe07a5f3c1475d503f1b0ea508d3b83b73b0b8637b57629f7f) |
| **Version Used** | v0.9.0 |

## Deployment

### Testnet (v0.2.0 - Latest)

| Field | Value |
|-------|-------|
| **Package ID** | [`0x1637ddc0495a8833ebd580224dad7154dfb33477f73d2c7fb41e2b350efa55b3`](https://suiscan.xyz/testnet/object/0x1637ddc0495a8833ebd580224dad7154dfb33477f73d2c7fb41e2b350efa55b3) |
| **Version** | v0.2.0 |
| **Modules** | `bs_events`, `d_values`, `entry`, `greeks`, `option_pricing` |

### Live Example Transactions

**1. Option Pricing (ATM Call/Put)**

[`CdAxPyw1T7tF4xMPpfVqVhJMDL4Xy6zeyC24YeQxpjJt`](https://suiscan.xyz/testnet/tx/CdAxPyw1T7tF4xMPpfVqVhJMDL4Xy6zeyC24YeQxpjJt)

```
Function: entry::price_option
Inputs:   S=$100, K=$100, T=1yr, r=5%, σ=20%

Emitted OptionPriceEvent:
├── call_price:      $10.45
├── put_price:       $5.57
└── parity_verified: true ✓
```

**2. Full Analysis (Prices + Greeks)**

[`48TFYV87TXRJMUuCzoMZ4T5CLVsFgQoT1fptR2w7NXPv`](https://suiscan.xyz/testnet/tx/48TFYV87TXRJMUuCzoMZ4T5CLVsFgQoT1fptR2w7NXPv)

```
Function: entry::full_analysis_call
Inputs:   S=$100, K=$100, T=1yr, r=5%, σ=20%

Emitted GreeksEvent:
├── delta (Δ):  0.637  (sensitivity to spot)
├── gamma (Γ):  0.019  (rate of change of delta)
├── vega (ν):   37.52  (sensitivity to vol)
├── theta (θ): -6.41   (time decay/year)
└── rho (ρ):    53.23  (sensitivity to rate)
```

## Overview

This package implements the Black-Scholes formula for European option pricing on Sui. It depends on the [gaussian](https://github.com/Evan-Kim2028/move-gaussian) library for normal distribution functions (CDF, PDF) and transcendental math (ln, exp, sqrt).

### What's Implemented

- **Option Pricing**: European call and put prices
- **Greeks**: Delta, Gamma, Vega, Theta, Rho
- **Put-Call Parity**: Verification that C - P = S - Ke^(-rT)
- **Entry Functions**: On-chain callable functions that emit events

### What's Been Tested

- 83 Move unit tests covering core functionality
- Compared against scipy's `scipy.stats.norm` for reference values
- Basic edge cases (ATM, ITM, OTM options)
- Put-call parity verification

### Known Limitations

- **European options only** - no early exercise
- **No dividends** - assumes no dividend payments  
- **Constant volatility** - single implied vol, no smile/skew
- **Domain bounds** - underlying gaussian CDF is clamped to ±6σ
- **Fixed-point precision** - WAD (10^18) arithmetic has inherent rounding

## Installation

Add to your `Move.toml`:

```toml
[dependencies]
black_scholes = { git = "https://github.com/Evan-Kim2028/move-black-scholes.git", rev = "v0.2.0" }
```

> **Note**: This will automatically pull in the [gaussian](https://github.com/Evan-Kim2028/move-gaussian) dependency (v0.9.0).

## Usage

All values use WAD scaling (multiply by 10^18):

```move
use black_scholes::option_pricing;

// S = 100, K = 100, T = 1 year, r = 5%, σ = 20%
let spot = 100_000_000_000_000_000_000;      // 100 * 1e18
let strike = 100_000_000_000_000_000_000;    // 100 * 1e18  
let time = 1_000_000_000_000_000_000;        // 1.0 * 1e18
let rate = 50_000_000_000_000_000;           // 0.05 * 1e18
let vol = 200_000_000_000_000_000;           // 0.20 * 1e18

let (call, put) = option_pricing::option_prices(spot, strike, time, rate, vol);
```

### Greeks

```move
use black_scholes::greeks;

let delta = greeks::delta_call(spot, strike, time, rate, vol);
let gamma = greeks::gamma(spot, strike, time, rate, vol);
let vega = greeks::vega(spot, strike, time, rate, vol);

// All at once
let (delta, gamma, vega, theta, rho) = greeks::all_greeks_call(spot, strike, time, rate, vol);
```

## Black-Scholes Formula

For reference, the standard Black-Scholes formulas:

```
d₁ = [ln(S/K) + (r + σ²/2)T] / (σ√T)
d₂ = d₁ - σ√T

Call = S·Φ(d₁) - K·e^(-rT)·Φ(d₂)
Put  = K·e^(-rT)·Φ(-d₂) - S·Φ(-d₁)
```

Where Φ(x) is the standard normal CDF.

## Modules

| Module | Description |
|--------|-------------|
| `d_values` | Computes d₁ and d₂ |
| `option_pricing` | Call/put prices, intrinsic values, parity check |
| `greeks` | Delta, Gamma, Vega, Theta, Rho |
| `entry` | Entry functions that emit events |
| `bs_events` | Event struct definitions |

## Entry Functions

These can be called directly via transactions:

| Function | Emits |
|----------|-------|
| `price_option` | `OptionPriceEvent` |
| `compute_greeks_call` | `GreeksEvent` |
| `compute_greeks_put` | `GreeksEvent` |
| `full_analysis_call` | Both events |
| `full_analysis_put` | Both events |

See [`DEPLOYMENTS.toml`](./DEPLOYMENTS.toml) for full API reference.

## Building & Testing

```bash
sui move build
sui move test
```

## License

MIT

## Disclaimer

This software is provided "as is" without warranty. It is experimental code for learning purposes. Do your own research and testing before using in any context where correctness matters.
