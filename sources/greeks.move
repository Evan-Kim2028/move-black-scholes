/// Option Greeks - sensitivities of option prices to various parameters.
/// 
/// # Overview
/// 
/// The "Greeks" measure how option prices change with respect to:
/// - Delta (Δ): Change in underlying price
/// - Gamma (Γ): Rate of change of delta
/// - Vega (ν): Change in volatility
/// - Theta (Θ): Time decay
/// - Rho (ρ): Change in interest rate
/// 
/// # Formulas
/// 
/// For a call option:
/// - Delta = Φ(d₁)
/// - Gamma = φ(d₁) / (S·σ·√T)
/// - Vega = S·φ(d₁)·√T
/// - Theta = -(S·φ(d₁)·σ)/(2√T) - r·K·e^(-rT)·Φ(d₂)
/// - Rho = K·T·e^(-rT)·Φ(d₂)
/// 
/// For a put option:
/// - Delta = Φ(d₁) - 1
/// - Gamma = same as call
/// - Vega = same as call
/// - Theta = -(S·φ(d₁)·σ)/(2√T) + r·K·e^(-rT)·Φ(-d₂)
/// - Rho = -K·T·e^(-rT)·Φ(-d₂)
module black_scholes::greeks {
    use gaussian::transcendental;
    use gaussian::normal_forward;
    use gaussian::signed_wad::{Self, SignedWad};
    use black_scholes::d_values;

    // === Constants ===

    /// WAD scale factor: 10^18
    const SCALE: u256 = 1_000_000_000_000_000_000;

    // === Errors ===

    /// Spot price must be positive
    const ESpotZero: u64 = 300;
    
    /// Strike price must be positive
    const EStrikeZero: u64 = 301;
    
    /// Time to expiry must be positive
    const ETimeZero: u64 = 302;
    
    /// Volatility must be positive
    const EVolatilityZero: u64 = 303;

    // === Delta ===

    /// Compute delta for a call option.
    /// 
    /// # Formula
    /// ```
    /// Delta_call = Φ(d₁)
    /// ```
    /// 
    /// Delta represents the change in option price for a $1 change in 
    /// the underlying. Call delta is always in [0, 1].
    /// 
    /// # Arguments
    /// * `spot` - Current asset price S (WAD-scaled)
    /// * `strike` - Strike price K (WAD-scaled)
    /// * `time_to_expiry` - Time T in years (WAD-scaled)
    /// * `rate` - Risk-free rate r (WAD-scaled)
    /// * `volatility` - Annualized volatility σ (WAD-scaled)
    /// 
    /// # Returns
    /// * `u256` - Delta in [0, SCALE] representing [0, 1]
    public fun delta_call(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
    ): u256 {
        validate_inputs(spot, strike, time_to_expiry, volatility);
        
        let d1 = d_values::compute_d1(spot, strike, time_to_expiry, rate, volatility);
        normal_forward::cdf_standard(&d1)
    }

    /// Compute delta for a put option.
    /// 
    /// # Formula
    /// ```
    /// Delta_put = Φ(d₁) - 1
    /// ```
    /// 
    /// Put delta is always in [-1, 0], returned as SignedWad.
    /// 
    /// # Returns
    /// * `SignedWad` - Delta (negative for puts)
    public fun delta_put(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
    ): SignedWad {
        let call_delta = delta_call(spot, strike, time_to_expiry, rate, volatility);
        // Put delta = Call delta - 1
        // Since call_delta is in [0, 1], result is in [-1, 0]
        if (call_delta >= SCALE) {
            signed_wad::zero()
        } else {
            signed_wad::new(SCALE - call_delta, true) // negative
        }
    }

    // === Gamma ===

    /// Compute gamma (same for both call and put).
    /// 
    /// # Formula
    /// ```
    /// Gamma = φ(d₁) / (S·σ·√T)
    /// ```
    /// 
    /// Gamma measures how fast delta changes. It's highest for ATM options.
    /// 
    /// # Returns
    /// * `u256` - Gamma (WAD-scaled)
    public fun gamma(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
    ): u256 {
        validate_inputs(spot, strike, time_to_expiry, volatility);
        
        let d1 = d_values::compute_d1(spot, strike, time_to_expiry, rate, volatility);
        let phi_d1 = normal_forward::pdf_standard(&d1); // φ(d₁)
        
        // Denominator: S·σ·√T
        let sqrt_t = transcendental::sqrt_wad(time_to_expiry);
        let s_vol_sqrt_t = (spot * volatility / SCALE) * sqrt_t / SCALE;
        
        // Gamma = φ(d₁) / (S·σ·√T)
        if (s_vol_sqrt_t == 0) {
            0
        } else {
            (phi_d1 * SCALE) / s_vol_sqrt_t
        }
    }

    // === Vega ===

    /// Compute vega (same for both call and put).
    /// 
    /// # Formula
    /// ```
    /// Vega = S·φ(d₁)·√T
    /// ```
    /// 
    /// Vega measures sensitivity to volatility changes.
    /// Note: This is vega per 1 unit change in vol, not per 1% change.
    /// 
    /// # Returns
    /// * `u256` - Vega (WAD-scaled)
    public fun vega(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
    ): u256 {
        validate_inputs(spot, strike, time_to_expiry, volatility);
        
        let d1 = d_values::compute_d1(spot, strike, time_to_expiry, rate, volatility);
        let phi_d1 = normal_forward::pdf_standard(&d1);
        
        let sqrt_t = transcendental::sqrt_wad(time_to_expiry);
        
        // Vega = S·φ(d₁)·√T
        (spot * phi_d1 / SCALE) * sqrt_t / SCALE
    }

    // === Theta ===

    /// Compute theta for a call option.
    /// 
    /// # Formula
    /// ```
    /// Theta_call = -(S·φ(d₁)·σ)/(2√T) - r·K·e^(-rT)·Φ(d₂)
    /// ```
    /// 
    /// Theta measures time decay. It's typically negative (options lose value over time).
    /// This returns theta per year; divide by 365 for daily theta.
    /// 
    /// # Returns
    /// * `SignedWad` - Theta (usually negative)
    public fun theta_call(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
    ): SignedWad {
        validate_inputs(spot, strike, time_to_expiry, volatility);
        
        let (d1, d2) = d_values::compute_d_values(spot, strike, time_to_expiry, rate, volatility);
        let phi_d1 = normal_forward::pdf_standard(&d1);
        let cdf_d2 = normal_forward::cdf_standard(&d2);
        
        let sqrt_t = transcendental::sqrt_wad(time_to_expiry);
        
        // Term 1: -(S·φ(d₁)·σ)/(2√T)
        // Numerator: S·φ(d₁)·σ
        let term1_num = (spot * phi_d1 / SCALE) * volatility / SCALE;
        // Denominator: 2√T
        let term1_denom = 2 * sqrt_t;
        let term1 = if (term1_denom > 0) {
            (term1_num * SCALE) / term1_denom
        } else {
            0
        };
        
        // Term 2: r·K·e^(-rT)·Φ(d₂)
        let r_times_t = (rate * time_to_expiry) / SCALE;
        let neg_r_t = signed_wad::new(r_times_t, true);
        let discount = transcendental::exp_wad(&neg_r_t);
        let term2 = (rate * strike / SCALE) * discount / SCALE * cdf_d2 / SCALE;
        
        // Theta = -term1 - term2
        let total = term1 + term2;
        signed_wad::new(total, true) // negative
    }

    /// Compute theta for a put option.
    /// 
    /// # Formula
    /// ```
    /// Theta_put = -(S·φ(d₁)·σ)/(2√T) + r·K·e^(-rT)·Φ(-d₂)
    /// ```
    /// 
    /// # Returns
    /// * `SignedWad` - Theta (can be positive for deep ITM puts)
    public fun theta_put(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
    ): SignedWad {
        validate_inputs(spot, strike, time_to_expiry, volatility);
        
        let (d1, d2) = d_values::compute_d_values(spot, strike, time_to_expiry, rate, volatility);
        let phi_d1 = normal_forward::pdf_standard(&d1);
        let neg_d2 = signed_wad::negate(&d2);
        let cdf_neg_d2 = normal_forward::cdf_standard(&neg_d2);
        
        let sqrt_t = transcendental::sqrt_wad(time_to_expiry);
        
        // Term 1: -(S·φ(d₁)·σ)/(2√T) - always negative
        let term1_num = (spot * phi_d1 / SCALE) * volatility / SCALE;
        let term1_denom = 2 * sqrt_t;
        let term1 = if (term1_denom > 0) {
            (term1_num * SCALE) / term1_denom
        } else {
            0
        };
        
        // Term 2: r·K·e^(-rT)·Φ(-d₂) - positive contribution
        let r_times_t = (rate * time_to_expiry) / SCALE;
        let neg_r_t = signed_wad::new(r_times_t, true);
        let discount = transcendental::exp_wad(&neg_r_t);
        let term2 = (rate * strike / SCALE) * discount / SCALE * cdf_neg_d2 / SCALE;
        
        // Theta = -term1 + term2
        if (term2 >= term1) {
            signed_wad::from_wad(term2 - term1)
        } else {
            signed_wad::new(term1 - term2, true)
        }
    }

    // === Rho ===

    /// Compute rho for a call option.
    /// 
    /// # Formula
    /// ```
    /// Rho_call = K·T·e^(-rT)·Φ(d₂)
    /// ```
    /// 
    /// Rho measures sensitivity to interest rate changes.
    /// This returns rho per 1 unit change in r, not per 1%.
    /// 
    /// # Returns
    /// * `u256` - Rho (WAD-scaled, positive for calls)
    public fun rho_call(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
    ): u256 {
        validate_inputs(spot, strike, time_to_expiry, volatility);
        
        let (_, d2) = d_values::compute_d_values(spot, strike, time_to_expiry, rate, volatility);
        let cdf_d2 = normal_forward::cdf_standard(&d2);
        
        let r_times_t = (rate * time_to_expiry) / SCALE;
        let neg_r_t = signed_wad::new(r_times_t, true);
        let discount = transcendental::exp_wad(&neg_r_t);
        
        // Rho = K·T·e^(-rT)·Φ(d₂)
        (strike * time_to_expiry / SCALE) * discount / SCALE * cdf_d2 / SCALE
    }

    /// Compute rho for a put option.
    /// 
    /// # Formula
    /// ```
    /// Rho_put = -K·T·e^(-rT)·Φ(-d₂)
    /// ```
    /// 
    /// # Returns
    /// * `SignedWad` - Rho (negative for puts)
    public fun rho_put(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
    ): SignedWad {
        validate_inputs(spot, strike, time_to_expiry, volatility);
        
        let (_, d2) = d_values::compute_d_values(spot, strike, time_to_expiry, rate, volatility);
        let neg_d2 = signed_wad::negate(&d2);
        let cdf_neg_d2 = normal_forward::cdf_standard(&neg_d2);
        
        let r_times_t = (rate * time_to_expiry) / SCALE;
        let neg_r_t = signed_wad::new(r_times_t, true);
        let discount = transcendental::exp_wad(&neg_r_t);
        
        // Rho = -K·T·e^(-rT)·Φ(-d₂)
        let magnitude = (strike * time_to_expiry / SCALE) * discount / SCALE * cdf_neg_d2 / SCALE;
        signed_wad::new(magnitude, true) // negative
    }

    // === All Greeks ===

    /// Compute all Greeks for a call option in one call.
    /// 
    /// More gas efficient than computing each Greek separately.
    /// 
    /// # Returns
    /// * `(u256, u256, u256, SignedWad, u256)` - (delta, gamma, vega, theta, rho)
    public fun all_greeks_call(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
    ): (u256, u256, u256, SignedWad, u256) {
        validate_inputs(spot, strike, time_to_expiry, volatility);
        
        // Compute d values once
        let (d1, d2) = d_values::compute_d_values(spot, strike, time_to_expiry, rate, volatility);
        
        // CDF and PDF values
        let cdf_d1 = normal_forward::cdf_standard(&d1);
        let cdf_d2 = normal_forward::cdf_standard(&d2);
        let pdf_d1 = normal_forward::pdf_standard(&d1);
        
        // Common terms
        let sqrt_t = transcendental::sqrt_wad(time_to_expiry);
        let r_times_t = (rate * time_to_expiry) / SCALE;
        let neg_r_t = signed_wad::new(r_times_t, true);
        let discount = transcendental::exp_wad(&neg_r_t);
        
        // Delta = Φ(d₁)
        let delta = cdf_d1;
        
        // Gamma = φ(d₁) / (S·σ·√T)
        let s_vol_sqrt_t = (spot * volatility / SCALE) * sqrt_t / SCALE;
        let gamma_val = if (s_vol_sqrt_t > 0) {
            (pdf_d1 * SCALE) / s_vol_sqrt_t
        } else {
            0
        };
        
        // Vega = S·φ(d₁)·√T
        let vega_val = (spot * pdf_d1 / SCALE) * sqrt_t / SCALE;
        
        // Theta = -(S·φ(d₁)·σ)/(2√T) - r·K·e^(-rT)·Φ(d₂)
        let term1_num = (spot * pdf_d1 / SCALE) * volatility / SCALE;
        let term1 = if (sqrt_t > 0) {
            (term1_num * SCALE) / (2 * sqrt_t)
        } else {
            0
        };
        let term2 = (rate * strike / SCALE) * discount / SCALE * cdf_d2 / SCALE;
        let theta_val = signed_wad::new(term1 + term2, true);
        
        // Rho = K·T·e^(-rT)·Φ(d₂)
        let rho_val = (strike * time_to_expiry / SCALE) * discount / SCALE * cdf_d2 / SCALE;
        
        (delta, gamma_val, vega_val, theta_val, rho_val)
    }

    /// Compute all Greeks for a put option in one call.
    /// 
    /// # Returns
    /// * `(SignedWad, u256, u256, SignedWad, SignedWad)` - (delta, gamma, vega, theta, rho)
    public fun all_greeks_put(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
    ): (SignedWad, u256, u256, SignedWad, SignedWad) {
        validate_inputs(spot, strike, time_to_expiry, volatility);
        
        // Compute d values once
        let (d1, d2) = d_values::compute_d_values(spot, strike, time_to_expiry, rate, volatility);
        
        // CDF and PDF values
        let cdf_d1 = normal_forward::cdf_standard(&d1);
        let pdf_d1 = normal_forward::pdf_standard(&d1);
        let neg_d2 = signed_wad::negate(&d2);
        let cdf_neg_d2 = normal_forward::cdf_standard(&neg_d2);
        
        // Common terms
        let sqrt_t = transcendental::sqrt_wad(time_to_expiry);
        let r_times_t = (rate * time_to_expiry) / SCALE;
        let neg_r_t = signed_wad::new(r_times_t, true);
        let discount = transcendental::exp_wad(&neg_r_t);
        
        // Delta = Φ(d₁) - 1
        let delta_val = if (cdf_d1 >= SCALE) {
            signed_wad::zero()
        } else {
            signed_wad::new(SCALE - cdf_d1, true)
        };
        
        // Gamma = φ(d₁) / (S·σ·√T) - same as call
        let s_vol_sqrt_t = (spot * volatility / SCALE) * sqrt_t / SCALE;
        let gamma_val = if (s_vol_sqrt_t > 0) {
            (pdf_d1 * SCALE) / s_vol_sqrt_t
        } else {
            0
        };
        
        // Vega = S·φ(d₁)·√T - same as call
        let vega_val = (spot * pdf_d1 / SCALE) * sqrt_t / SCALE;
        
        // Theta = -(S·φ(d₁)·σ)/(2√T) + r·K·e^(-rT)·Φ(-d₂)
        let term1_num = (spot * pdf_d1 / SCALE) * volatility / SCALE;
        let term1 = if (sqrt_t > 0) {
            (term1_num * SCALE) / (2 * sqrt_t)
        } else {
            0
        };
        let term2 = (rate * strike / SCALE) * discount / SCALE * cdf_neg_d2 / SCALE;
        let theta_val = if (term2 >= term1) {
            signed_wad::from_wad(term2 - term1)
        } else {
            signed_wad::new(term1 - term2, true)
        };
        
        // Rho = -K·T·e^(-rT)·Φ(-d₂)
        let rho_mag = (strike * time_to_expiry / SCALE) * discount / SCALE * cdf_neg_d2 / SCALE;
        let rho_val = signed_wad::new(rho_mag, true);
        
        (delta_val, gamma_val, vega_val, theta_val, rho_val)
    }

    // === Internal Functions ===

    fun validate_inputs(spot: u256, strike: u256, time_to_expiry: u256, volatility: u256) {
        assert!(spot > 0, ESpotZero);
        assert!(strike > 0, EStrikeZero);
        assert!(time_to_expiry > 0, ETimeZero);
        assert!(volatility > 0, EVolatilityZero);
    }

    // === Test Functions ===

    #[test]
    fun test_delta_call_atm() {
        // ATM call: delta should be around 0.5
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let delta = delta_call(spot, strike, time, rate, vol);
        
        // ATM delta ≈ 0.5-0.6 (slightly above 0.5 due to drift)
        assert!(delta > SCALE / 2, 0); // > 0.5
        assert!(delta < 7 * SCALE / 10, 1); // < 0.7
    }

    #[test]
    fun test_delta_call_bounds() {
        // Delta should always be in [0, 1]
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let delta = delta_call(spot, strike, time, rate, vol);
        assert!(delta <= SCALE, 0);
    }

    #[test]
    fun test_delta_put_negative() {
        // Put delta should be negative
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let delta = delta_put(spot, strike, time, rate, vol);
        assert!(signed_wad::is_negative(&delta), 0);
    }

    #[test]
    fun test_delta_put_call_relationship() {
        // Put delta = Call delta - 1
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let call_delta = delta_call(spot, strike, time, rate, vol);
        let put_delta = delta_put(spot, strike, time, rate, vol);
        
        // put_delta + 1 should equal call_delta
        let put_delta_abs = signed_wad::abs(&put_delta);
        // call_delta + put_delta_abs should equal SCALE (within tolerance)
        let sum = call_delta + put_delta_abs;
        assert!(sum > SCALE * 99 / 100 && sum < SCALE * 101 / 100, 0);
    }

    #[test]
    fun test_gamma_positive() {
        // Gamma should always be positive
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let gamma_val = gamma(spot, strike, time, rate, vol);
        assert!(gamma_val > 0, 0);
    }

    #[test]
    fun test_gamma_highest_atm() {
        // Gamma should be highest for ATM options
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let gamma_atm = gamma(100 * SCALE, strike, time, rate, vol);
        let gamma_itm = gamma(120 * SCALE, strike, time, rate, vol);
        let gamma_otm = gamma(80 * SCALE, strike, time, rate, vol);

        assert!(gamma_atm > gamma_itm, 0);
        assert!(gamma_atm > gamma_otm, 1);
    }

    #[test]
    fun test_vega_positive() {
        // Vega should always be positive
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let vega_val = vega(spot, strike, time, rate, vol);
        assert!(vega_val > 0, 0);
    }

    #[test]
    fun test_vega_highest_atm() {
        // Vega should be highest for ATM options
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let vega_atm = vega(100 * SCALE, strike, time, rate, vol);
        let vega_itm = vega(120 * SCALE, strike, time, rate, vol);
        let vega_otm = vega(80 * SCALE, strike, time, rate, vol);

        assert!(vega_atm > vega_itm, 0);
        assert!(vega_atm > vega_otm, 1);
    }

    #[test]
    fun test_theta_call_negative() {
        // Theta for ATM call should be negative (time decay)
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let theta = theta_call(spot, strike, time, rate, vol);
        assert!(signed_wad::is_negative(&theta), 0);
    }

    #[test]
    fun test_rho_call_positive() {
        // Rho for call should be positive (higher rates help calls)
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let rho_val = rho_call(spot, strike, time, rate, vol);
        assert!(rho_val > 0, 0);
    }

    #[test]
    fun test_rho_put_negative() {
        // Rho for put should be negative
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let rho_val = rho_put(spot, strike, time, rate, vol);
        assert!(signed_wad::is_negative(&rho_val), 0);
    }

    #[test]
    fun test_all_greeks_call() {
        // Test that all_greeks_call produces same results as individual functions
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let (delta, gamma_val, vega_val, _theta, _rho_val) = 
            all_greeks_call(spot, strike, time, rate, vol);
        
        let delta_individual = delta_call(spot, strike, time, rate, vol);
        let gamma_individual = gamma(spot, strike, time, rate, vol);
        let vega_individual = vega(spot, strike, time, rate, vol);
        
        // Allow small tolerance due to rounding
        let tolerance = SCALE / 1000; // 0.1%
        
        assert!(abs_diff(delta, delta_individual) <= tolerance, 0);
        assert!(abs_diff(gamma_val, gamma_individual) <= tolerance, 1);
        assert!(abs_diff(vega_val, vega_individual) <= tolerance, 2);
    }

    #[test]
    fun test_deep_itm_call_delta() {
        // Deep ITM call: delta should be close to 1
        let spot = 200 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let delta = delta_call(spot, strike, time, rate, vol);
        assert!(delta > SCALE * 95 / 100, 0); // > 0.95
    }

    #[test]
    fun test_deep_otm_call_delta() {
        // Deep OTM call: delta should be close to 0
        let spot = 50 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let delta = delta_call(spot, strike, time, rate, vol);
        assert!(delta < SCALE / 10, 0); // < 0.1
    }

    // Helper function for tests
    #[test_only]
    fun abs_diff(a: u256, b: u256): u256 {
        if (a >= b) { a - b } else { b - a }
    }
}
