/// Entry functions for on-chain Black-Scholes calculations.
/// 
/// # Overview
/// 
/// This module provides `entry` functions that can be called directly from
/// transactions. Each function emits events for off-chain indexing.
/// 
/// # Functions
/// 
/// - `price_option` - Compute call and put prices
/// - `compute_greeks_call` - Compute all Greeks for a call option
/// - `compute_greeks_put` - Compute all Greeks for a put option
/// - `compute_d_values` - Compute d₁ and d₂ values
/// - `compute_intrinsic` - Compute intrinsic values
/// 
/// # Example Transaction
/// 
/// ```bash
/// sui client call --package $PACKAGE_ID --module entry --function price_option \
///     --args 100000000000000000000 100000000000000000000 1000000000000000000 \
///            50000000000000000 200000000000000000
/// ```
/// 
/// Parameters (all WAD-scaled):
/// - spot = 100 * 10^18 (S = $100)
/// - strike = 100 * 10^18 (K = $100)
/// - time = 1 * 10^18 (T = 1 year)
/// - rate = 0.05 * 10^18 (r = 5%)
/// - volatility = 0.2 * 10^18 (σ = 20%)
module black_scholes::entry {
    use black_scholes::option_pricing;
    use black_scholes::greeks;
    use black_scholes::d_values;
    use black_scholes::bs_events;
    use gaussian::signed_wad;

    // === Entry Functions ===

    /// Compute European call and put option prices.
    /// 
    /// Emits an `OptionPriceEvent` with the results.
    /// 
    /// # Arguments
    /// * `spot` - Current asset price S (WAD-scaled)
    /// * `strike` - Strike price K (WAD-scaled)
    /// * `time_to_expiry` - Time T in years (WAD-scaled)
    /// * `rate` - Risk-free rate r (WAD-scaled)
    /// * `volatility` - Volatility σ (WAD-scaled)
    /// 
    /// # Events
    /// Emits `OptionPriceEvent` with call/put prices and parity verification.
    entry fun price_option(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
        ctx: &TxContext,
    ) {
        let (call_price, put_price) = option_pricing::option_prices(
            spot, strike, time_to_expiry, rate, volatility
        );
        
        let parity_verified = option_pricing::verify_put_call_parity(
            spot, strike, time_to_expiry, rate, volatility
        );

        bs_events::emit_option_price(
            spot,
            strike,
            time_to_expiry,
            rate,
            volatility,
            call_price,
            put_price,
            parity_verified,
            ctx.sender(),
        );
    }

    /// Compute all Greeks for a call option.
    /// 
    /// Emits a `GreeksEvent` with delta, gamma, vega, theta, and rho.
    /// 
    /// # Arguments
    /// Same as `price_option`.
    /// 
    /// # Events
    /// Emits `GreeksEvent` with all five Greeks.
    entry fun compute_greeks_call(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
        ctx: &TxContext,
    ) {
        let (delta, gamma, vega, theta, rho) = greeks::all_greeks_call(
            spot, strike, time_to_expiry, rate, volatility
        );

        bs_events::emit_greeks(
            spot,
            strike,
            time_to_expiry,
            rate,
            volatility,
            true, // is_call
            delta,
            false, // call delta is always positive
            gamma,
            vega,
            signed_wad::abs(&theta),
            signed_wad::is_negative(&theta),
            rho,
            false, // call rho is always positive
            ctx.sender(),
        );
    }

    /// Compute all Greeks for a put option.
    /// 
    /// Emits a `GreeksEvent` with delta, gamma, vega, theta, and rho.
    /// 
    /// # Arguments
    /// Same as `price_option`.
    /// 
    /// # Events
    /// Emits `GreeksEvent` with all five Greeks.
    entry fun compute_greeks_put(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
        ctx: &TxContext,
    ) {
        let (delta, gamma, vega, theta, rho) = greeks::all_greeks_put(
            spot, strike, time_to_expiry, rate, volatility
        );

        bs_events::emit_greeks(
            spot,
            strike,
            time_to_expiry,
            rate,
            volatility,
            false, // is_call = false (put)
            signed_wad::abs(&delta),
            signed_wad::is_negative(&delta),
            gamma,
            vega,
            signed_wad::abs(&theta),
            signed_wad::is_negative(&theta),
            signed_wad::abs(&rho),
            signed_wad::is_negative(&rho),
            ctx.sender(),
        );
    }

    /// Compute d₁ and d₂ values used in Black-Scholes.
    /// 
    /// Useful for debugging and verification.
    /// 
    /// # Arguments
    /// Same as `price_option`.
    /// 
    /// # Events
    /// Emits `DValuesEvent` with d₁, d₂, and σ√T.
    entry fun compute_d_values_entry(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
        ctx: &TxContext,
    ) {
        let (d1, d2) = d_values::compute_d_values(
            spot, strike, time_to_expiry, rate, volatility
        );
        let vol_sqrt_t = d_values::compute_vol_sqrt_t(time_to_expiry, volatility);

        bs_events::emit_d_values(
            spot,
            strike,
            time_to_expiry,
            rate,
            volatility,
            signed_wad::abs(&d1),
            signed_wad::is_negative(&d1),
            signed_wad::abs(&d2),
            signed_wad::is_negative(&d2),
            vol_sqrt_t,
            ctx.sender(),
        );
    }

    /// Compute intrinsic values for call and put options.
    /// 
    /// # Arguments
    /// * `spot` - Current asset price S (WAD-scaled)
    /// * `strike` - Strike price K (WAD-scaled)
    /// 
    /// # Events
    /// Emits `IntrinsicValueEvent` with call and put intrinsic values.
    entry fun compute_intrinsic(
        spot: u256,
        strike: u256,
        ctx: &TxContext,
    ) {
        let call_intrinsic = option_pricing::call_intrinsic(spot, strike);
        let put_intrinsic = option_pricing::put_intrinsic(spot, strike);

        bs_events::emit_intrinsic_value(
            spot,
            strike,
            call_intrinsic,
            put_intrinsic,
            ctx.sender(),
        );
    }

    /// Batch compute: prices + Greeks for a call option.
    /// 
    /// More gas-efficient than calling separately.
    /// Emits both `OptionPriceEvent` and `GreeksEvent`.
    entry fun full_analysis_call(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
        ctx: &TxContext,
    ) {
        // Price the option
        let (call_price, put_price) = option_pricing::option_prices(
            spot, strike, time_to_expiry, rate, volatility
        );
        let parity_verified = option_pricing::verify_put_call_parity(
            spot, strike, time_to_expiry, rate, volatility
        );

        bs_events::emit_option_price(
            spot,
            strike,
            time_to_expiry,
            rate,
            volatility,
            call_price,
            put_price,
            parity_verified,
            ctx.sender(),
        );

        // Compute Greeks
        let (delta, gamma, vega, theta, rho) = greeks::all_greeks_call(
            spot, strike, time_to_expiry, rate, volatility
        );

        bs_events::emit_greeks(
            spot,
            strike,
            time_to_expiry,
            rate,
            volatility,
            true,
            delta,
            false,
            gamma,
            vega,
            signed_wad::abs(&theta),
            signed_wad::is_negative(&theta),
            rho,
            false,
            ctx.sender(),
        );
    }

    /// Batch compute: prices + Greeks for a put option.
    /// 
    /// More gas-efficient than calling separately.
    /// Emits both `OptionPriceEvent` and `GreeksEvent`.
    entry fun full_analysis_put(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
        ctx: &TxContext,
    ) {
        // Price the option
        let (call_price, put_price) = option_pricing::option_prices(
            spot, strike, time_to_expiry, rate, volatility
        );
        let parity_verified = option_pricing::verify_put_call_parity(
            spot, strike, time_to_expiry, rate, volatility
        );

        bs_events::emit_option_price(
            spot,
            strike,
            time_to_expiry,
            rate,
            volatility,
            call_price,
            put_price,
            parity_verified,
            ctx.sender(),
        );

        // Compute Greeks
        let (delta, gamma, vega, theta, rho) = greeks::all_greeks_put(
            spot, strike, time_to_expiry, rate, volatility
        );

        bs_events::emit_greeks(
            spot,
            strike,
            time_to_expiry,
            rate,
            volatility,
            false,
            signed_wad::abs(&delta),
            signed_wad::is_negative(&delta),
            gamma,
            vega,
            signed_wad::abs(&theta),
            signed_wad::is_negative(&theta),
            signed_wad::abs(&rho),
            signed_wad::is_negative(&rho),
            ctx.sender(),
        );
    }
}
