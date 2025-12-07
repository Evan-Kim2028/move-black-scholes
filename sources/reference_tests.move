/// Reference tests comparing Move implementation to scipy Black-Scholes.
/// 
/// These tests verify accuracy against scipy.stats.norm calculations.
/// All expected values come from test_black_scholes.py.
#[test_only]
module black_scholes::reference_tests {
    use black_scholes::option_pricing;
    use black_scholes::d_values;
    use black_scholes::greeks;
    use gaussian::signed_wad;

    // === Constants ===
    
    const SCALE: u256 = 1_000_000_000_000_000_000;
    
    // Tolerances
    const PRICE_TOLERANCE_BPS: u256 = 100;  // 1% tolerance
    const GREEK_TOLERANCE_BPS: u256 = 500;  // 5% tolerance for Greeks

    // === Helper Functions ===

    fun within_tolerance(actual: u256, expected: u256, tolerance_bps: u256): bool {
        let mut tolerance = (expected * tolerance_bps) / 10000;
        if (tolerance == 0) { tolerance = 1 };
        
        if (actual >= expected) {
            actual - expected <= tolerance
        } else {
            expected - actual <= tolerance
        }
    }

    // === ATM Reference Tests ===

    #[test]
    /// scipy reference: S=100, K=100, T=1, r=5%, σ=20%
    /// Expected: call=10.4506, put=5.5735
    fun test_atm_call_scipy_reference() {
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;  // 5%
        let vol = SCALE / 5;    // 20%

        let call = option_pricing::call_price(spot, strike, time, rate, vol);
        
        // Expected: 10.4506 WAD = 10_450_600_000_000_000_000
        let expected = 10_450_600_000_000_000_000u256;
        assert!(within_tolerance(call, expected, PRICE_TOLERANCE_BPS), 0);
    }

    #[test]
    fun test_atm_put_scipy_reference() {
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let put = option_pricing::put_price(spot, strike, time, rate, vol);
        
        // Expected: 5.5735 WAD = 5_573_500_000_000_000_000
        let expected = 5_573_500_000_000_000_000u256;
        assert!(within_tolerance(put, expected, PRICE_TOLERANCE_BPS), 0);
    }

    #[test]
    /// scipy reference: d1=0.35, d2=0.15
    fun test_atm_d_values_scipy_reference() {
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let (d1, d2) = d_values::compute_d_values(spot, strike, time, rate, vol);
        
        // d1 ≈ 0.35
        let d1_expected = 350_000_000_000_000_000u256;
        let d1_actual = signed_wad::abs(&d1);
        assert!(within_tolerance(d1_actual, d1_expected, 200), 0); // 2% tolerance
        assert!(!signed_wad::is_negative(&d1), 1);
        
        // d2 ≈ 0.15
        let d2_expected = 150_000_000_000_000_000u256;
        let d2_actual = signed_wad::abs(&d2);
        assert!(within_tolerance(d2_actual, d2_expected, 200), 2);
        assert!(!signed_wad::is_negative(&d2), 3);
    }

    // === ITM Reference Tests ===

    #[test]
    /// scipy reference: S=120, K=100, T=1, r=5%, σ=20%
    /// Expected: call=26.17, put=1.29
    fun test_itm_call_scipy_reference() {
        let spot = 120 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let call = option_pricing::call_price(spot, strike, time, rate, vol);
        
        // Expected: 26.17 WAD
        let expected = 26_170_000_000_000_000_000u256;
        assert!(within_tolerance(call, expected, PRICE_TOLERANCE_BPS), 0);
    }

    #[test]
    fun test_itm_put_scipy_reference() {
        let spot = 120 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let put = option_pricing::put_price(spot, strike, time, rate, vol);
        
        // Expected: 1.29 WAD
        let expected = 1_290_000_000_000_000_000u256;
        assert!(within_tolerance(put, expected, PRICE_TOLERANCE_BPS), 0);
    }

    // === OTM Reference Tests ===

    #[test]
    /// scipy reference: S=80, K=100, T=1, r=5%, σ=20%
    /// Expected: call=1.86, put=16.98
    fun test_otm_call_scipy_reference() {
        let spot = 80 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let call = option_pricing::call_price(spot, strike, time, rate, vol);
        
        // Expected: 1.86 WAD
        let expected = 1_860_000_000_000_000_000u256;
        assert!(within_tolerance(call, expected, PRICE_TOLERANCE_BPS), 0);
    }

    #[test]
    fun test_otm_put_scipy_reference() {
        let spot = 80 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let put = option_pricing::put_price(spot, strike, time, rate, vol);
        
        // Expected: 16.98 WAD
        let expected = 16_980_000_000_000_000_000u256;
        assert!(within_tolerance(put, expected, PRICE_TOLERANCE_BPS), 0);
    }

    // === Short Expiry Reference Tests ===

    #[test]
    /// scipy reference: S=100, K=100, T=0.25, r=5%, σ=20%
    /// Expected: call=4.61, put=3.37
    fun test_short_expiry_call_scipy_reference() {
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE / 4;  // 0.25 years
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let call = option_pricing::call_price(spot, strike, time, rate, vol);
        
        // Expected: 4.61 WAD
        let expected = 4_610_000_000_000_000_000u256;
        assert!(within_tolerance(call, expected, PRICE_TOLERANCE_BPS), 0);
    }

    #[test]
    fun test_short_expiry_put_scipy_reference() {
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE / 4;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let put = option_pricing::put_price(spot, strike, time, rate, vol);
        
        // Expected: 3.37 WAD
        let expected = 3_370_000_000_000_000_000u256;
        assert!(within_tolerance(put, expected, PRICE_TOLERANCE_BPS), 0);
    }

    // === High Volatility Reference Tests ===

    #[test]
    /// scipy reference: S=100, K=100, T=1, r=5%, σ=40%
    /// Expected: call=18.02, put=13.15
    fun test_high_vol_call_scipy_reference() {
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = 2 * SCALE / 5;  // 40%

        let call = option_pricing::call_price(spot, strike, time, rate, vol);
        
        // Expected: 18.02 WAD
        let expected = 18_020_000_000_000_000_000u256;
        assert!(within_tolerance(call, expected, PRICE_TOLERANCE_BPS), 0);
    }

    #[test]
    fun test_high_vol_put_scipy_reference() {
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = 2 * SCALE / 5;

        let put = option_pricing::put_price(spot, strike, time, rate, vol);
        
        // Expected: 13.15 WAD
        let expected = 13_150_000_000_000_000_000u256;
        assert!(within_tolerance(put, expected, PRICE_TOLERANCE_BPS), 0);
    }

    // === High Rate Reference Tests ===

    #[test]
    /// scipy reference: S=100, K=100, T=1, r=10%, σ=20%
    /// Expected: call=13.27, put=3.75
    fun test_high_rate_call_scipy_reference() {
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 10;  // 10%
        let vol = SCALE / 5;

        let call = option_pricing::call_price(spot, strike, time, rate, vol);
        
        // Expected: 13.27 WAD
        let expected = 13_270_000_000_000_000_000u256;
        assert!(within_tolerance(call, expected, PRICE_TOLERANCE_BPS), 0);
    }

    #[test]
    fun test_high_rate_put_scipy_reference() {
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 10;
        let vol = SCALE / 5;

        let put = option_pricing::put_price(spot, strike, time, rate, vol);
        
        // Expected: 3.75 WAD
        let expected = 3_750_000_000_000_000_000u256;
        assert!(within_tolerance(put, expected, PRICE_TOLERANCE_BPS), 0);
    }

    // === Greeks Reference Tests ===

    #[test]
    /// scipy reference: delta_call=0.6368
    fun test_atm_delta_call_scipy_reference() {
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let delta = greeks::delta_call(spot, strike, time, rate, vol);
        
        // Expected: 0.6368 = 636_800_000_000_000_000
        let expected = 636_800_000_000_000_000u256;
        assert!(within_tolerance(delta, expected, GREEK_TOLERANCE_BPS), 0);
    }

    #[test]
    /// scipy reference: gamma=0.0188
    fun test_atm_gamma_scipy_reference() {
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let gamma = greeks::gamma(spot, strike, time, rate, vol);
        
        // Expected: 0.0188 = 18_800_000_000_000_000
        let expected = 18_800_000_000_000_000u256;
        assert!(within_tolerance(gamma, expected, GREEK_TOLERANCE_BPS), 0);
    }

    #[test]
    /// scipy reference: vega=37.52
    fun test_atm_vega_scipy_reference() {
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let vega = greeks::vega(spot, strike, time, rate, vol);
        
        // Expected: 37.52 = 37_520_000_000_000_000_000
        let expected = 37_520_000_000_000_000_000u256;
        assert!(within_tolerance(vega, expected, GREEK_TOLERANCE_BPS), 0);
    }

    // === Put-Call Parity Reference Tests ===

    #[test]
    fun test_put_call_parity_atm() {
        assert!(option_pricing::verify_put_call_parity(
            100 * SCALE, 100 * SCALE, SCALE, SCALE / 20, SCALE / 5
        ), 0);
    }

    #[test]
    fun test_put_call_parity_itm() {
        assert!(option_pricing::verify_put_call_parity(
            120 * SCALE, 100 * SCALE, SCALE, SCALE / 20, SCALE / 5
        ), 0);
    }

    #[test]
    fun test_put_call_parity_otm() {
        assert!(option_pricing::verify_put_call_parity(
            80 * SCALE, 100 * SCALE, SCALE, SCALE / 20, SCALE / 5
        ), 0);
    }

    #[test]
    fun test_put_call_parity_high_vol() {
        assert!(option_pricing::verify_put_call_parity(
            100 * SCALE, 100 * SCALE, SCALE, SCALE / 20, SCALE / 2
        ), 0);
    }

    #[test]
    fun test_put_call_parity_short_expiry() {
        assert!(option_pricing::verify_put_call_parity(
            100 * SCALE, 100 * SCALE, SCALE / 4, SCALE / 20, SCALE / 5
        ), 0);
    }

    // === Boundary Tests ===

    #[test]
    /// Deep ITM call delta should be ~1
    fun test_deep_itm_call_delta() {
        let delta = greeks::delta_call(
            200 * SCALE, 100 * SCALE, SCALE, SCALE / 20, SCALE / 5
        );
        // Should be > 0.99
        assert!(delta > 990_000_000_000_000_000, 0);
    }

    #[test]
    /// Deep OTM call delta should be ~0
    fun test_deep_otm_call_delta() {
        let delta = greeks::delta_call(
            50 * SCALE, 100 * SCALE, SCALE, SCALE / 20, SCALE / 5
        );
        // Should be < 0.01
        assert!(delta < 10_000_000_000_000_000, 0);
    }

    #[test]
    /// Call price should increase with spot
    fun test_call_monotonic_spot() {
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let call_80 = option_pricing::call_price(80 * SCALE, strike, time, rate, vol);
        let call_100 = option_pricing::call_price(100 * SCALE, strike, time, rate, vol);
        let call_120 = option_pricing::call_price(120 * SCALE, strike, time, rate, vol);

        assert!(call_80 < call_100, 0);
        assert!(call_100 < call_120, 1);
    }

    #[test]
    /// Call price should decrease with strike
    fun test_call_monotonic_strike() {
        let spot = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let call_k80 = option_pricing::call_price(spot, 80 * SCALE, time, rate, vol);
        let call_k100 = option_pricing::call_price(spot, 100 * SCALE, time, rate, vol);
        let call_k120 = option_pricing::call_price(spot, 120 * SCALE, time, rate, vol);

        assert!(call_k80 > call_k100, 0);
        assert!(call_k100 > call_k120, 1);
    }

    #[test]
    /// Call price should increase with volatility
    fun test_call_monotonic_vol() {
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;

        let call_vol10 = option_pricing::call_price(spot, strike, time, rate, SCALE / 10);
        let call_vol20 = option_pricing::call_price(spot, strike, time, rate, SCALE / 5);
        let call_vol40 = option_pricing::call_price(spot, strike, time, rate, 2 * SCALE / 5);

        assert!(call_vol10 < call_vol20, 0);
        assert!(call_vol20 < call_vol40, 1);
    }

    #[test]
    /// Call price should increase with time (for ATM)
    fun test_call_monotonic_time() {
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let call_3m = option_pricing::call_price(spot, strike, SCALE / 4, rate, vol);
        let call_6m = option_pricing::call_price(spot, strike, SCALE / 2, rate, vol);
        let call_1y = option_pricing::call_price(spot, strike, SCALE, rate, vol);

        assert!(call_3m < call_6m, 0);
        assert!(call_6m < call_1y, 1);
    }
}
