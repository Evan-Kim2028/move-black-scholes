/// Edge case tests for Black-Scholes implementation.
/// 
/// These tests cover boundary conditions and stress scenarios that
/// could reveal numerical issues in the implementation.
#[test_only]
module black_scholes::edge_case_tests {
    use black_scholes::option_pricing;
    use black_scholes::d_values;
    use black_scholes::greeks;
    use gaussian::signed_wad;

    const SCALE: u256 = 1_000_000_000_000_000_000;

    // ============================================================
    // EXTREME MONEYNESS TESTS
    // ============================================================

    #[test]
    /// Very deep ITM call: S = 1000, K = 100 (10x ITM)
    /// Call should be approximately S - K*e^(-rT) (pure intrinsic + small time value)
    fun test_very_deep_itm_call() {
        let spot = 1000 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let call = option_pricing::call_price(spot, strike, time, rate, vol);
        let discount = option_pricing::discount_factor(rate, time);
        let forward_diff = spot - (strike * discount) / SCALE;

        // Call should be very close to forward difference for deep ITM
        // Allow 1% tolerance
        assert!(call >= forward_diff * 99 / 100, 0);
        assert!(call <= forward_diff * 102 / 100, 1); // Slightly above due to time value
    }

    #[test]
    /// Very deep OTM call: S = 10, K = 100 (0.1x, 90% OTM)
    /// Call should be nearly zero
    fun test_very_deep_otm_call() {
        let spot = 10 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let call = option_pricing::call_price(spot, strike, time, rate, vol);
        
        // Should be less than 0.01% of spot
        assert!(call < spot / 10000, 0);
    }

    #[test]
    /// Symmetric deep ITM put: S = 10, K = 100
    fun test_very_deep_itm_put() {
        let spot = 10 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let put = option_pricing::put_price(spot, strike, time, rate, vol);
        let discount = option_pricing::discount_factor(rate, time);
        let discounted_strike = (strike * discount) / SCALE;

        // Deep ITM put ≈ K*e^(-rT) - S
        let expected_min = discounted_strike - spot;
        assert!(put >= expected_min * 95 / 100, 0);
    }

    // ============================================================
    // EXTREME TIME TESTS
    // ============================================================

    #[test]
    /// Very short expiry: T = 1 day = 1/365 years
    fun test_one_day_expiry() {
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE / 365;  // 1 day
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let (call, put) = option_pricing::option_prices(spot, strike, time, rate, vol);
        
        // Very short expiry ATM options should have small time value
        // Both should be positive but small
        assert!(call > 0, 0);
        assert!(put > 0, 1);
        assert!(call < 5 * SCALE, 2); // Less than $5
        assert!(put < 5 * SCALE, 3);
    }

    #[test]
    /// Very long expiry: T = 10 years
    fun test_ten_year_expiry() {
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = 10 * SCALE;  // 10 years
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let call = option_pricing::call_price(spot, strike, time, rate, vol);
        
        // Long expiry = more time value
        // 10-year ATM call with 20% vol should be substantial
        assert!(call > 30 * SCALE, 0); // Should be > $30
    }

    // ============================================================
    // EXTREME VOLATILITY TESTS
    // ============================================================

    #[test]
    /// Very low volatility: σ = 1%
    fun test_very_low_vol() {
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 100;  // 1%

        let (call, put) = option_pricing::option_prices(spot, strike, time, rate, vol);
        
        // Low vol = low time value
        // Call should be close to forward value
        assert!(call > 0, 0);
        assert!(put > 0, 1);
        
        // Put-call parity should still hold
        assert!(option_pricing::verify_put_call_parity(spot, strike, time, rate, vol), 2);
    }

    #[test]
    /// Very high volatility: σ = 100%
    fun test_very_high_vol() {
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE;  // 100%

        let (call, put) = option_pricing::option_prices(spot, strike, time, rate, vol);
        
        // High vol = high time value
        assert!(call > 30 * SCALE, 0); // Should be substantial
        assert!(put > 25 * SCALE, 1);
        
        // Put-call parity should still hold
        assert!(option_pricing::verify_put_call_parity(spot, strike, time, rate, vol), 2);
    }

    // ============================================================
    // EXTREME RATE TESTS
    // ============================================================

    #[test]
    /// Zero interest rate
    fun test_zero_rate() {
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = 0;  // 0%
        let vol = SCALE / 5;

        let (call, put) = option_pricing::option_prices(spot, strike, time, rate, vol);
        
        // With r=0, discount = 1, so C - P = S - K
        // For ATM (S = K), C should equal P
        let diff = if (call >= put) { call - put } else { put - call };
        
        // Should be very close (< 0.1%)
        assert!(diff < SCALE / 1000, 0);
    }

    #[test]
    /// High interest rate: r = 20%
    fun test_high_rate() {
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 5;  // 20%
        let vol = SCALE / 5;

        let (call, put) = option_pricing::option_prices(spot, strike, time, rate, vol);
        
        // High rate benefits calls (forward > spot)
        // Call should be worth more than put for ATM
        assert!(call > put, 0);
        
        // Put-call parity should still hold
        assert!(option_pricing::verify_put_call_parity(spot, strike, time, rate, vol), 1);
    }

    // ============================================================
    // EXTREME PRICE TESTS
    // ============================================================

    #[test]
    /// Very small prices: S = K = 0.01
    fun test_small_prices() {
        let spot = SCALE / 100;  // 0.01
        let strike = SCALE / 100;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let (call, put) = option_pricing::option_prices(spot, strike, time, rate, vol);
        
        // Should still work and be proportional
        assert!(call > 0, 0);
        assert!(put > 0, 1);
        
        // Put-call parity
        assert!(option_pricing::verify_put_call_parity(spot, strike, time, rate, vol), 2);
    }

    #[test]
    /// Large prices: S = K = 100,000
    fun test_large_prices() {
        let spot = 100000 * SCALE;
        let strike = 100000 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let (call, _put) = option_pricing::option_prices(spot, strike, time, rate, vol);
        
        // Should scale linearly with price
        // ATM 100k option should be ~1000x the 100 option
        let (call_100, _) = option_pricing::option_prices(
            100 * SCALE, 100 * SCALE, time, rate, vol
        );
        
        let ratio = (call * 100) / call_100;
        // Should be approximately 100000 (within 1%)
        assert!(ratio > 99000 && ratio < 101000, 0);
    }

    // ============================================================
    // D-VALUE EDGE CASES
    // ============================================================

    #[test]
    /// d-values at extreme moneyness should be clamped reasonably
    fun test_d_values_extreme_itm() {
        let spot = 1000 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let (d1, d2) = d_values::compute_d_values(spot, strike, time, rate, vol);
        
        // Both should be very positive for deep ITM
        assert!(!signed_wad::is_negative(&d1), 0);
        assert!(!signed_wad::is_negative(&d2), 1);
        
        // d1 should be > 6 (at CDF limit)
        let d1_mag = signed_wad::abs(&d1);
        assert!(d1_mag > 6 * SCALE, 2);
    }

    #[test]
    /// d-values at extreme OTM should be very negative
    fun test_d_values_extreme_otm() {
        let spot = 10 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let (d1, d2) = d_values::compute_d_values(spot, strike, time, rate, vol);
        
        // Both should be very negative for deep OTM
        assert!(signed_wad::is_negative(&d1), 0);
        assert!(signed_wad::is_negative(&d2), 1);
    }

    // ============================================================
    // GREEKS EDGE CASES
    // ============================================================

    #[test]
    /// Delta should approach 1 for very deep ITM
    fun test_delta_approaches_one() {
        let spot = 500 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let delta = greeks::delta_call(spot, strike, time, rate, vol);
        
        // Should be > 0.999
        assert!(delta > 999 * SCALE / 1000, 0);
    }

    #[test]
    /// Delta should approach 0 for very deep OTM
    fun test_delta_approaches_zero() {
        let spot = 20 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let delta = greeks::delta_call(spot, strike, time, rate, vol);
        
        // Should be < 0.001
        assert!(delta < SCALE / 1000, 0);
    }

    #[test]
    /// Gamma should be very small for deep ITM/OTM
    fun test_gamma_small_at_extremes() {
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let gamma_atm = greeks::gamma(100 * SCALE, strike, time, rate, vol);
        let gamma_itm = greeks::gamma(200 * SCALE, strike, time, rate, vol);
        let gamma_otm = greeks::gamma(50 * SCALE, strike, time, rate, vol);
        
        // ATM gamma should be much larger
        assert!(gamma_atm > gamma_itm * 2, 0);
        assert!(gamma_atm > gamma_otm * 2, 1);
    }

    // ============================================================
    // NUMERICAL STABILITY TESTS
    // ============================================================

    #[test]
    /// Test that very small time doesn't cause division issues
    fun test_small_time_stability() {
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE / 1000;  // 0.001 years (~8 hours)
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        // Should not panic
        let (call, put) = option_pricing::option_prices(spot, strike, time, rate, vol);
        
        assert!(call > 0 || call == 0, 0); // Just verify it returns
        assert!(put > 0 || put == 0, 1);
    }

    #[test]
    /// Verify symmetry: swapping S and K should swap call and put values
    /// (adjusted for discounting)
    fun test_price_symmetry() {
        let price_a = 80 * SCALE;
        let price_b = 120 * SCALE;
        let time = SCALE;
        let rate = 0;  // Use 0 rate for cleaner symmetry
        let vol = SCALE / 5;

        let call_a = option_pricing::call_price(price_a, price_b, time, rate, vol);
        let put_b = option_pricing::put_price(price_b, price_a, time, rate, vol);

        // With r=0: Call(S=80, K=120) should equal Put(S=120, K=80)
        let diff = if (call_a >= put_b) { call_a - put_b } else { put_b - call_a };
        
        // Allow 1% tolerance
        assert!(diff < call_a / 100 + SCALE / 100, 0);
    }
}
