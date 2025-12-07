/// Computes d₁ and d₂ values for Black-Scholes formula.
/// 
/// # Overview
/// 
/// The d-values are intermediate calculations used in both option pricing
/// and Greeks. This module computes them once to avoid redundant work.
/// 
/// # Formulas
/// 
/// ```
/// d₁ = [ln(S/K) + (r + σ²/2)T] / (σ√T)
/// d₂ = d₁ - σ√T
/// ```
/// 
/// Where:
/// - S = Spot price
/// - K = Strike price  
/// - T = Time to expiry (years)
/// - r = Risk-free rate
/// - σ = Volatility
module black_scholes::d_values {
    use gaussian::transcendental;
    use gaussian::signed_wad::{Self, SignedWad};

    // === Constants ===

    /// WAD scale factor: 10^18
    const SCALE: u256 = 1_000_000_000_000_000_000;

    // === Errors ===

    /// Spot price must be positive
    const ESpotZero: u64 = 100;
    
    /// Strike price must be positive
    const EStrikeZero: u64 = 101;
    
    /// Time to expiry must be positive
    const ETimeZero: u64 = 102;
    
    /// Volatility must be positive
    const EVolatilityZero: u64 = 103;

    // === Public Functions ===

    /// Compute d₁ and d₂ values for Black-Scholes.
    /// 
    /// # Arguments
    /// * `spot` - Current asset price S (WAD-scaled)
    /// * `strike` - Strike price K (WAD-scaled)
    /// * `time_to_expiry` - Time T in years (WAD-scaled)
    /// * `rate` - Risk-free rate r (WAD-scaled, e.g., 0.05 = 5%)
    /// * `volatility` - Annualized volatility σ (WAD-scaled, e.g., 0.2 = 20%)
    /// 
    /// # Returns
    /// * `(SignedWad, SignedWad)` - (d₁, d₂) values
    /// 
    /// # Example
    /// ```move
    /// // S = 100, K = 100, T = 1 year, r = 5%, σ = 20%
    /// let spot = 100_000_000_000_000_000_000;
    /// let strike = 100_000_000_000_000_000_000;
    /// let time = 1_000_000_000_000_000_000;
    /// let rate = 50_000_000_000_000_000;
    /// let vol = 200_000_000_000_000_000;
    /// 
    /// let (d1, d2) = d_values::compute_d_values(spot, strike, time, rate, vol);
    /// // d1 ≈ 0.35, d2 ≈ 0.15
    /// ```
    public fun compute_d_values(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
    ): (SignedWad, SignedWad) {
        // Validate inputs
        assert!(spot > 0, ESpotZero);
        assert!(strike > 0, EStrikeZero);
        assert!(time_to_expiry > 0, ETimeZero);
        assert!(volatility > 0, EVolatilityZero);

        // Step 1: Compute ln(S/K)
        // S/K = (spot * SCALE) / strike
        let s_over_k = (spot * SCALE) / strike;
        let ln_s_over_k = transcendental::ln_wad(s_over_k);

        // Step 2: Compute σ²/2
        // vol^2 = (volatility * volatility) / SCALE
        let vol_squared = (volatility * volatility) / SCALE;
        // vol^2 / 2
        let half_vol_squared = vol_squared / 2;

        // Step 3: Compute (r + σ²/2) * T
        // (rate + half_vol_squared) * time_to_expiry / SCALE
        let rate_plus_half_vol_sq = rate + half_vol_squared;
        let drift_term = (rate_plus_half_vol_sq * time_to_expiry) / SCALE;
        let drift_signed = signed_wad::from_wad(drift_term);

        // Step 4: Compute numerator = ln(S/K) + (r + σ²/2)T
        let numerator = signed_wad::add(&ln_s_over_k, &drift_signed);

        // Step 5: Compute σ√T (denominator)
        let sqrt_t = transcendental::sqrt_wad(time_to_expiry);
        let vol_sqrt_t = (volatility * sqrt_t) / SCALE;

        // Step 6: d₁ = numerator / (σ√T)
        let vol_sqrt_t_signed = signed_wad::from_wad(vol_sqrt_t);
        let d1 = signed_wad::div_wad(&numerator, &vol_sqrt_t_signed);

        // Step 7: d₂ = d₁ - σ√T
        let d2 = signed_wad::sub(&d1, &vol_sqrt_t_signed);

        (d1, d2)
    }

    /// Compute only d₁ (useful for some Greeks calculations).
    /// 
    /// # Arguments
    /// Same as `compute_d_values`.
    /// 
    /// # Returns
    /// * `SignedWad` - d₁ value only
    public fun compute_d1(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
    ): SignedWad {
        let (d1, _) = compute_d_values(spot, strike, time_to_expiry, rate, volatility);
        d1
    }

    /// Compute σ√T (volatility times square root of time).
    /// 
    /// This is a common term used in d₂ = d₁ - σ√T and Vega calculations.
    /// 
    /// # Arguments
    /// * `time_to_expiry` - Time T in years (WAD-scaled)
    /// * `volatility` - Annualized volatility σ (WAD-scaled)
    /// 
    /// # Returns
    /// * `u256` - σ√T (WAD-scaled)
    public fun compute_vol_sqrt_t(
        time_to_expiry: u256,
        volatility: u256,
    ): u256 {
        let sqrt_t = transcendental::sqrt_wad(time_to_expiry);
        (volatility * sqrt_t) / SCALE
    }

    // === Test Functions ===

    #[test]
    fun test_d_values_atm() {
        // ATM option: S = K = 100, T = 1, r = 5%, σ = 20%
        // Expected: d1 ≈ 0.35, d2 ≈ 0.15
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE; // 1 year
        let rate = SCALE / 20; // 5% = 0.05
        let vol = SCALE / 5;   // 20% = 0.2

        let (d1, d2) = compute_d_values(spot, strike, time, rate, vol);

        // d1 should be positive and around 0.35
        assert!(!signed_wad::is_negative(&d1), 0);
        let d1_mag = signed_wad::abs(&d1);
        // d1 ≈ 0.35 = 350_000_000_000_000_000
        // Allow 5% tolerance: [0.33, 0.37]
        assert!(d1_mag > 330_000_000_000_000_000, 1);
        assert!(d1_mag < 370_000_000_000_000_000, 2);

        // d2 should be positive and around 0.15
        assert!(!signed_wad::is_negative(&d2), 3);
        let d2_mag = signed_wad::abs(&d2);
        // d2 ≈ 0.15 = 150_000_000_000_000_000
        // Allow 5% tolerance: [0.13, 0.17]
        assert!(d2_mag > 130_000_000_000_000_000, 4);
        assert!(d2_mag < 170_000_000_000_000_000, 5);

        // d1 - d2 should equal σ√T = 0.2
        let diff = signed_wad::sub(&d1, &d2);
        let diff_mag = signed_wad::abs(&diff);
        // σ√T = 0.2 * 1 = 0.2 = 200_000_000_000_000_000
        assert!(diff_mag > 195_000_000_000_000_000, 6);
        assert!(diff_mag < 205_000_000_000_000_000, 7);
    }

    #[test]
    fun test_d_values_itm_call() {
        // ITM call: S = 120, K = 100
        // d1 and d2 should be larger (more positive)
        let spot = 120 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let (d1, _d2) = compute_d_values(spot, strike, time, rate, vol);

        // d1 should be significantly positive (> 1.0)
        assert!(!signed_wad::is_negative(&d1), 0);
        let d1_mag = signed_wad::abs(&d1);
        assert!(d1_mag > SCALE, 1); // d1 > 1.0
    }

    #[test]
    fun test_d_values_otm_call() {
        // OTM call: S = 80, K = 100
        // d1 and d2 should be negative
        let spot = 80 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let (_d1, d2) = compute_d_values(spot, strike, time, rate, vol);

        // d2 should be negative for OTM
        assert!(signed_wad::is_negative(&d2), 0);
    }

    #[test]
    fun test_vol_sqrt_t() {
        // σ = 20%, T = 1 year
        // σ√T = 0.2 * 1 = 0.2
        let time = SCALE;
        let vol = SCALE / 5;

        let result = compute_vol_sqrt_t(time, vol);
        
        // Should be approximately 0.2
        assert!(result > 195_000_000_000_000_000, 0);
        assert!(result < 205_000_000_000_000_000, 1);
    }

    #[test]
    fun test_vol_sqrt_t_quarter_year() {
        // σ = 40%, T = 0.25 years (3 months)
        // σ√T = 0.4 * 0.5 = 0.2
        let time = SCALE / 4; // 0.25
        let vol = 2 * SCALE / 5; // 0.4

        let result = compute_vol_sqrt_t(time, vol);
        
        // Should be approximately 0.2
        assert!(result > 195_000_000_000_000_000, 0);
        assert!(result < 205_000_000_000_000_000, 1);
    }

    #[test]
    #[expected_failure(abort_code = ESpotZero)]
    fun test_zero_spot_fails() {
        compute_d_values(0, SCALE, SCALE, SCALE / 20, SCALE / 5);
    }

    #[test]
    #[expected_failure(abort_code = EStrikeZero)]
    fun test_zero_strike_fails() {
        compute_d_values(SCALE, 0, SCALE, SCALE / 20, SCALE / 5);
    }

    #[test]
    #[expected_failure(abort_code = ETimeZero)]
    fun test_zero_time_fails() {
        compute_d_values(SCALE, SCALE, 0, SCALE / 20, SCALE / 5);
    }

    #[test]
    #[expected_failure(abort_code = EVolatilityZero)]
    fun test_zero_vol_fails() {
        compute_d_values(SCALE, SCALE, SCALE, SCALE / 20, 0);
    }
}
