/// On-chain events for Black-Scholes option pricing operations.
/// 
/// # Overview
/// 
/// This module defines events emitted by Black-Scholes pricing and Greeks
/// calculations. Events enable off-chain indexing, analytics, and debugging.
/// 
/// # Events
/// 
/// - `OptionPriceEvent` - Emitted when computing call/put prices
/// - `GreeksEvent` - Emitted when computing option Greeks
/// - `DValuesEvent` - Emitted when computing d₁/d₂ values
/// 
/// # Usage
/// 
/// Events are emitted by entry functions in the `entry` module.
/// Off-chain systems can subscribe to these events for:
/// - Option pricing analytics
/// - Greeks monitoring for risk management
/// - DeFi protocol integration
/// - Transaction debugging
module black_scholes::bs_events {
    use sui::event;

    // === Event Structs ===

    /// Emitted when computing European option prices.
    /// 
    /// Contains all input parameters and the resulting call/put prices.
    /// All values are WAD-scaled (10^18).
    public struct OptionPriceEvent has copy, drop {
        /// Spot price S (WAD-scaled)
        spot: u256,
        /// Strike price K (WAD-scaled)
        strike: u256,
        /// Time to expiry T in years (WAD-scaled)
        time_to_expiry: u256,
        /// Risk-free rate r (WAD-scaled, e.g., 0.05 = 5%)
        rate: u256,
        /// Volatility σ (WAD-scaled, e.g., 0.2 = 20%)
        volatility: u256,
        /// Computed call option price (WAD-scaled)
        call_price: u256,
        /// Computed put option price (WAD-scaled)
        put_price: u256,
        /// Put-call parity holds within tolerance
        parity_verified: bool,
        /// Address that initiated the computation
        caller: address,
    }

    /// Emitted when computing option Greeks.
    /// 
    /// Contains all input parameters and the computed Greeks.
    /// All values are WAD-scaled (10^18).
    public struct GreeksEvent has copy, drop {
        /// Spot price S (WAD-scaled)
        spot: u256,
        /// Strike price K (WAD-scaled)
        strike: u256,
        /// Time to expiry T in years (WAD-scaled)
        time_to_expiry: u256,
        /// Risk-free rate r (WAD-scaled)
        rate: u256,
        /// Volatility σ (WAD-scaled)
        volatility: u256,
        /// Option type: true = call, false = put
        is_call: bool,
        /// Delta: sensitivity to spot price (WAD-scaled, can be negative for puts)
        delta_magnitude: u256,
        delta_negative: bool,
        /// Gamma: rate of change of delta (WAD-scaled, always positive)
        gamma: u256,
        /// Vega: sensitivity to volatility (WAD-scaled, always positive)
        vega: u256,
        /// Theta: time decay (WAD-scaled magnitude, usually negative)
        theta_magnitude: u256,
        theta_negative: bool,
        /// Rho: sensitivity to interest rate (WAD-scaled, can be negative for puts)
        rho_magnitude: u256,
        rho_negative: bool,
        /// Address that initiated the computation
        caller: address,
    }

    /// Emitted when computing d₁ and d₂ values.
    /// 
    /// Useful for debugging and verification of intermediate calculations.
    public struct DValuesEvent has copy, drop {
        /// Spot price S (WAD-scaled)
        spot: u256,
        /// Strike price K (WAD-scaled)
        strike: u256,
        /// Time to expiry T (WAD-scaled)
        time_to_expiry: u256,
        /// Risk-free rate r (WAD-scaled)
        rate: u256,
        /// Volatility σ (WAD-scaled)
        volatility: u256,
        /// d₁ magnitude (WAD-scaled)
        d1_magnitude: u256,
        /// d₁ sign (true = negative)
        d1_negative: bool,
        /// d₂ magnitude (WAD-scaled)
        d2_magnitude: u256,
        /// d₂ sign (true = negative)
        d2_negative: bool,
        /// σ√T term (WAD-scaled)
        vol_sqrt_t: u256,
        /// Address that initiated the computation
        caller: address,
    }

    /// Emitted for intrinsic value calculations.
    public struct IntrinsicValueEvent has copy, drop {
        /// Spot price S (WAD-scaled)
        spot: u256,
        /// Strike price K (WAD-scaled)
        strike: u256,
        /// Call intrinsic: max(S - K, 0)
        call_intrinsic: u256,
        /// Put intrinsic: max(K - S, 0)
        put_intrinsic: u256,
        /// Address that initiated the computation
        caller: address,
    }

    // === Internal Emit Functions ===

    /// Emit an OptionPriceEvent.
    public(package) fun emit_option_price(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
        call_price: u256,
        put_price: u256,
        parity_verified: bool,
        caller: address,
    ) {
        event::emit(OptionPriceEvent {
            spot,
            strike,
            time_to_expiry,
            rate,
            volatility,
            call_price,
            put_price,
            parity_verified,
            caller,
        });
    }

    /// Emit a GreeksEvent.
    public(package) fun emit_greeks(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
        is_call: bool,
        delta_magnitude: u256,
        delta_negative: bool,
        gamma: u256,
        vega: u256,
        theta_magnitude: u256,
        theta_negative: bool,
        rho_magnitude: u256,
        rho_negative: bool,
        caller: address,
    ) {
        event::emit(GreeksEvent {
            spot,
            strike,
            time_to_expiry,
            rate,
            volatility,
            is_call,
            delta_magnitude,
            delta_negative,
            gamma,
            vega,
            theta_magnitude,
            theta_negative,
            rho_magnitude,
            rho_negative,
            caller,
        });
    }

    /// Emit a DValuesEvent.
    public(package) fun emit_d_values(
        spot: u256,
        strike: u256,
        time_to_expiry: u256,
        rate: u256,
        volatility: u256,
        d1_magnitude: u256,
        d1_negative: bool,
        d2_magnitude: u256,
        d2_negative: bool,
        vol_sqrt_t: u256,
        caller: address,
    ) {
        event::emit(DValuesEvent {
            spot,
            strike,
            time_to_expiry,
            rate,
            volatility,
            d1_magnitude,
            d1_negative,
            d2_magnitude,
            d2_negative,
            vol_sqrt_t,
            caller,
        });
    }

    /// Emit an IntrinsicValueEvent.
    public(package) fun emit_intrinsic_value(
        spot: u256,
        strike: u256,
        call_intrinsic: u256,
        put_intrinsic: u256,
        caller: address,
    ) {
        event::emit(IntrinsicValueEvent {
            spot,
            strike,
            call_intrinsic,
            put_intrinsic,
            caller,
        });
    }
}
