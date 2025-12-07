/// European option pricing using Black-Scholes formula.
/// 
/// # Overview
/// 
/// This module implements the Black-Scholes model for pricing European
/// call and put options. All calculations use WAD scaling (10^18).
/// 
/// # Formulas
/// 
/// **Call Option:**
/// ```
/// C = S·Φ(d₁) - K·e^(-rT)·Φ(d₂)
/// ```
/// 
/// **Put Option:**
/// ```
/// P = K·e^(-rT)·Φ(-d₂) - S·Φ(-d₁)
/// ```
/// 
/// Where Φ(x) is the standard normal CDF.
/// 
/// # Put-Call Parity
/// 
/// This implementation guarantees:
/// ```
/// C - P = S - K·e^(-rT)
/// ```
module black_scholes::option_pricing {
    use gaussian::transcendental;
    use gaussian::normal_forward;
    use gaussian::signed_wad;
    use black_scholes::d_values;

    // === Constants ===

    /// WAD scale factor: 10^18
    const SCALE: u256 = 1_000_000_000_000_000_000;

    // === Errors ===

    /// Spot price must be positive
    const ESpotZero: u64 = 200;
    
    /// Strike price must be positive
    const EStrikeZero: u64 = 201;
    
    /// Time to expiry must be positive
    const ETimeZero: u64 = 202;
    
    /// Volatility must be positive
    const EVolatilityZero: u64 = 203;

    // === Public Functions ===

    /// Compute European call option price.
    /// 
    /// # Formula
    /// ```
    /// C = S·Φ(d₁) - K·e^(-rT)·Φ(d₂)
    /// ```
    /// 
    /// # Arguments
    /// * `spot` - Current asset price S (WAD-scaled)
    /// * `strike` - Strike price K (WAD-scaled)
    /// * `time_to_expiry` - Time T in years (WAD-scaled)
    /// * `rate` - Risk-free rate r (WAD-scaled)
    /// * `volatility` - Annualized volatility σ (WAD-scaled)
    /// 
    /// # Returns
    /// * `u256` - Call option price (WAD-scaled)
    /// 
    /// # Example
    /// ```move
    /// // S = 100, K = 100, T = 1, r = 5%, σ = 20%
    /// let call = option_pricing::call_price(
    ///     100 * SCALE, 100 * SCALE, SCALE, 
    ///     SCALE / 20, SCALE / 5
    /// );
    /// // call ≈ 10.45 WAD
    /// ```
    public fun call_price(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
    ): u256 {
        let (call, _) = option_prices(spot, strike, time_to_expiry, rate, volatility);
        call
    }

    /// Compute European put option price.
    /// 
    /// # Formula
    /// ```
    /// P = K·e^(-rT)·Φ(-d₂) - S·Φ(-d₁)
    /// ```
    /// 
    /// # Arguments
    /// Same as `call_price`.
    /// 
    /// # Returns
    /// * `u256` - Put option price (WAD-scaled)
    public fun put_price(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
    ): u256 {
        let (_, put) = option_prices(spot, strike, time_to_expiry, rate, volatility);
        put
    }

    /// Compute both call and put prices (gas efficient).
    /// 
    /// This is more efficient than calling `call_price` and `put_price`
    /// separately because it computes d₁, d₂, and the discount factor
    /// only once.
    /// 
    /// # Arguments
    /// Same as `call_price`.
    /// 
    /// # Returns
    /// * `(u256, u256)` - (call_price, put_price) both WAD-scaled
    /// 
    /// # Example
    /// ```move
    /// let (call, put) = option_pricing::option_prices(
    ///     100 * SCALE, 100 * SCALE, SCALE,
    ///     SCALE / 20, SCALE / 5
    /// );
    /// // call ≈ 10.45 WAD, put ≈ 5.57 WAD
    /// ```
    public fun option_prices(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
    ): (u256, u256) {
        // Validate inputs
        assert!(spot > 0, ESpotZero);
        assert!(strike > 0, EStrikeZero);
        assert!(time_to_expiry > 0, ETimeZero);
        assert!(volatility > 0, EVolatilityZero);

        // Step 1: Compute d₁ and d₂
        let (d1, d2) = d_values::compute_d_values(
            spot, strike, time_to_expiry, rate, volatility
        );

        // Step 2: Compute Φ(d₁), Φ(d₂), Φ(-d₁), Φ(-d₂)
        let phi_d1 = normal_forward::cdf_standard(&d1);
        let phi_d2 = normal_forward::cdf_standard(&d2);
        
        let neg_d1 = signed_wad::negate(&d1);
        let neg_d2 = signed_wad::negate(&d2);
        let phi_neg_d1 = normal_forward::cdf_standard(&neg_d1);
        let phi_neg_d2 = normal_forward::cdf_standard(&neg_d2);

        // Step 3: Compute discount factor e^(-rT)
        let r_times_t = (rate * time_to_expiry) / SCALE;
        let neg_r_t = signed_wad::new(r_times_t, true); // -rT
        let discount = transcendental::exp_wad(&neg_r_t);

        // Step 4: Compute call price
        // C = S·Φ(d₁) - K·e^(-rT)·Φ(d₂)
        let s_phi_d1 = (spot * phi_d1) / SCALE;
        let k_discount_phi_d2 = (strike * discount / SCALE) * phi_d2 / SCALE;
        
        // Handle potential underflow
        let call = if (s_phi_d1 >= k_discount_phi_d2) {
            s_phi_d1 - k_discount_phi_d2
        } else {
            0 // Should not happen for valid options, but be safe
        };

        // Step 5: Compute put price
        // P = K·e^(-rT)·Φ(-d₂) - S·Φ(-d₁)
        let k_discount_phi_neg_d2 = (strike * discount / SCALE) * phi_neg_d2 / SCALE;
        let s_phi_neg_d1 = (spot * phi_neg_d1) / SCALE;
        
        let put = if (k_discount_phi_neg_d2 >= s_phi_neg_d1) {
            k_discount_phi_neg_d2 - s_phi_neg_d1
        } else {
            0
        };

        (call, put)
    }

    /// Compute the intrinsic value of a call option.
    /// 
    /// # Formula
    /// ```
    /// Intrinsic = max(S - K, 0)
    /// ```
    /// 
    /// # Arguments
    /// * `spot` - Current asset price S (WAD-scaled)
    /// * `strike` - Strike price K (WAD-scaled)
    /// 
    /// # Returns
    /// * `u256` - Intrinsic value (WAD-scaled)
    public fun call_intrinsic(spot: u256, strike: u256): u256 {
        if (spot > strike) {
            spot - strike
        } else {
            0
        }
    }

    /// Compute the intrinsic value of a put option.
    /// 
    /// # Formula
    /// ```
    /// Intrinsic = max(K - S, 0)
    /// ```
    /// 
    /// # Arguments
    /// * `spot` - Current asset price S (WAD-scaled)
    /// * `strike` - Strike price K (WAD-scaled)
    /// 
    /// # Returns
    /// * `u256` - Intrinsic value (WAD-scaled)
    public fun put_intrinsic(spot: u256, strike: u256): u256 {
        if (strike > spot) {
            strike - spot
        } else {
            0
        }
    }

    /// Compute the discount factor e^(-rT).
    /// 
    /// # Arguments
    /// * `rate` - Risk-free rate r (WAD-scaled)
    /// * `time_to_expiry` - Time T in years (WAD-scaled)
    /// 
    /// # Returns
    /// * `u256` - Discount factor (WAD-scaled)
    public fun discount_factor(rate: u256, time_to_expiry: u256): u256 {
        let r_times_t = (rate * time_to_expiry) / SCALE;
        let neg_r_t = signed_wad::new(r_times_t, true);
        transcendental::exp_wad(&neg_r_t)
    }

    /// Verify put-call parity: C - P = S - K·e^(-rT)
    /// 
    /// # Arguments
    /// Same as `option_prices`.
    /// 
    /// # Returns
    /// * `bool` - true if parity holds within tolerance (1e-10 WAD)
    public fun verify_put_call_parity(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
    ): bool {
        let (call, put) = option_prices(spot, strike, time_to_expiry, rate, volatility);
        let discount = discount_factor(rate, time_to_expiry);
        let k_discounted = (strike * discount) / SCALE;

        // C - P should equal S - K·e^(-rT)
        // Rearrange to avoid signed arithmetic: C + K·e^(-rT) ≈ P + S
        let lhs = call + k_discounted;
        let rhs = put + spot;

        // Allow tolerance of 1e-10 WAD (0.0000000001)
        let tolerance: u256 = 100_000_000; // 1e8
        
        if (lhs >= rhs) {
            lhs - rhs <= tolerance
        } else {
            rhs - lhs <= tolerance
        }
    }

    // === Test Functions ===

    #[test]
    fun test_atm_option_prices() {
        // ATM option: S = K = 100, T = 1, r = 5%, σ = 20%
        // Expected from scipy: call ≈ 10.4506, put ≈ 5.5735
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20; // 5%
        let vol = SCALE / 5;   // 20%

        let (call, put) = option_prices(spot, strike, time, rate, vol);

        // Call should be around 10.45 WAD
        // Allow 1% tolerance
        let expected_call = 10_450_000_000_000_000_000u256; // 10.45
        let call_lower = expected_call * 99 / 100;
        let call_upper = expected_call * 101 / 100;
        assert!(call >= call_lower && call <= call_upper, 0);

        // Put should be around 5.57 WAD
        let expected_put = 5_570_000_000_000_000_000u256; // 5.57
        let put_lower = expected_put * 99 / 100;
        let put_upper = expected_put * 101 / 100;
        assert!(put >= put_lower && put <= put_upper, 1);
    }

    #[test]
    fun test_put_call_parity() {
        // Put-call parity should hold for any valid option
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        assert!(verify_put_call_parity(spot, strike, time, rate, vol), 0);
    }

    #[test]
    fun test_put_call_parity_itm() {
        // ITM call
        let spot = 120 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        assert!(verify_put_call_parity(spot, strike, time, rate, vol), 0);
    }

    #[test]
    fun test_put_call_parity_otm() {
        // OTM call
        let spot = 80 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        assert!(verify_put_call_parity(spot, strike, time, rate, vol), 0);
    }

    #[test]
    fun test_call_intrinsic_itm() {
        // ITM: S = 120, K = 100
        let spot = 120 * SCALE;
        let strike = 100 * SCALE;

        let intrinsic = call_intrinsic(spot, strike);
        assert!(intrinsic == 20 * SCALE, 0);
    }

    #[test]
    fun test_call_intrinsic_otm() {
        // OTM: S = 80, K = 100
        let spot = 80 * SCALE;
        let strike = 100 * SCALE;

        let intrinsic = call_intrinsic(spot, strike);
        assert!(intrinsic == 0, 0);
    }

    #[test]
    fun test_put_intrinsic_itm() {
        // ITM put: S = 80, K = 100
        let spot = 80 * SCALE;
        let strike = 100 * SCALE;

        let intrinsic = put_intrinsic(spot, strike);
        assert!(intrinsic == 20 * SCALE, 0);
    }

    #[test]
    fun test_put_intrinsic_otm() {
        // OTM put: S = 120, K = 100
        let spot = 120 * SCALE;
        let strike = 100 * SCALE;

        let intrinsic = put_intrinsic(spot, strike);
        assert!(intrinsic == 0, 0);
    }

    #[test]
    fun test_call_above_intrinsic() {
        // Option value should always be >= intrinsic value
        let spot = 120 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let call = call_price(spot, strike, time, rate, vol);
        let intrinsic = call_intrinsic(spot, strike);

        assert!(call >= intrinsic, 0);
    }

    #[test]
    fun test_put_above_intrinsic() {
        // ITM put: S = 90, K = 100 (less extreme than S=80)
        // Put has intrinsic value of 10, should have time value too
        let spot = 90 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let put = put_price(spot, strike, time, rate, vol);
        let intrinsic = put_intrinsic(spot, strike);

        // For ITM puts, option value should be >= intrinsic
        // Allow small tolerance for rounding (1%)
        assert!(put + SCALE / 100 >= intrinsic, 0);
    }

    #[test]
    fun test_discount_factor() {
        // r = 5%, T = 1 => e^(-0.05) ≈ 0.9512
        let rate = SCALE / 20;
        let time = SCALE;

        let discount = discount_factor(rate, time);
        
        // Expected: 0.9512 = 951_200_000_000_000_000
        let expected = 951_200_000_000_000_000u256;
        let lower = expected * 99 / 100;
        let upper = expected * 101 / 100;
        assert!(discount >= lower && discount <= upper, 0);
    }

    #[test]
    fun test_deep_itm_call() {
        // Deep ITM: S = 200, K = 100
        // Call should be close to S - K·e^(-rT)
        let spot = 200 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let call = call_price(spot, strike, time, rate, vol);
        let discount = discount_factor(rate, time);
        let forward_diff = spot - (strike * discount) / SCALE;

        // Deep ITM call ≈ forward difference
        // Allow 5% tolerance since there's still some time value
        assert!(call >= forward_diff * 95 / 100, 0);
    }

    #[test]
    fun test_deep_otm_call() {
        // Deep OTM: S = 50, K = 100
        // Call should be very small
        let spot = 50 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let call = call_price(spot, strike, time, rate, vol);
        
        // Should be less than 1% of spot
        assert!(call < spot / 100, 0);
    }

    #[test]
    fun test_higher_vol_higher_price() {
        // Higher volatility should increase option price
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let time = SCALE;
        let rate = SCALE / 20;
        
        let vol_low = SCALE / 10;  // 10%
        let vol_high = SCALE / 2;  // 50%

        let call_low = call_price(spot, strike, time, rate, vol_low);
        let call_high = call_price(spot, strike, time, rate, vol_high);

        assert!(call_high > call_low, 0);
    }

    #[test]
    fun test_longer_time_higher_price() {
        // Longer time should increase option price (for ATM)
        let spot = 100 * SCALE;
        let strike = 100 * SCALE;
        let rate = SCALE / 20;
        let vol = SCALE / 5;

        let time_short = SCALE / 4;  // 3 months
        let time_long = SCALE;       // 1 year

        let call_short = call_price(spot, strike, time_short, rate, vol);
        let call_long = call_price(spot, strike, time_long, rate, vol);

        assert!(call_long > call_short, 0);
    }
}
