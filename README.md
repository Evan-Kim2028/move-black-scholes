# Black-Scholes Option Pricing Library

A Sui Move implementation of the Black-Scholes option pricing model for European options.

> ⚠️ **Experimental**: This is a learning/research project. It has not been audited and should not be used in production without thorough review.

## Deployment

### Testnet

| Field | Value |
|-------|-------|
| **Package ID** | `0xd33496acb857637f0498079aacd56cedae09c0a0fbed287239840c54f8e9719c` |
| **Version** | v0.1.0 |
| **Modules** | `bs_events`, `d_values`, `entry`, `greeks`, `option_pricing` |
| **Depends On** | [gaussian](https://github.com/Evan-Kim2028/move-gaussian) v0.7.0 |
| **Explorer** | [View on SuiVision](https://testnet.suivision.xyz/package/0xd33496acb857637f0498079aacd56cedae09c0a0fbed287239840c54f8e9719c) |

### Verified On-Chain

The package has been deployed and called on Sui testnet.

**Example Transaction**: [View on SuiVision](https://testnet.suivision.xyz/txblock/7BKaNWq6s5yfxsy2tYH5DHqspWvYpXK5Bh3ZpBkYKH8p)

```
Function: entry::price_option
Inputs:   S=100, K=100, T=1yr, r=5%, σ=20%

Emitted OptionPriceEvent:
├── call_price:      ~10.45
├── put_price:       ~5.57
└── parity_verified: true
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
black_scholes = { git = "https://github.com/Evan-Kim2028/move-black-scholes.git", rev = "v0.1.0" }
```

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

## Dependencies

This package uses the [gaussian](https://github.com/Evan-Kim2028/move-gaussian) library for:

- `transcendental::ln_wad`, `exp_wad`, `sqrt_wad`
- `normal_forward::cdf_standard`, `pdf_standard`
- `signed_wad` for signed fixed-point math

## License

MIT

## Disclaimer

This software is provided "as is" without warranty. It is experimental code for learning purposes. Do your own research and testing before using in any context where correctness matters.
