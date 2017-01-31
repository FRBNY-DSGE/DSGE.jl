"""
```
compute_forecast(m, system, kal; enforce_zlb = false, shocks = Matrix{S}())

compute_forecast(m, system, z0; enforce_zlb = false, shocks = Matrix{S}())

compute_forecast(system, z0, shocks; enforce_zlb = false)

compute_forecast(T, R, C, Q, Z, D, Z_pseudo, D_pseudo, z0, shocks; enforce_zlb = false)
```

### Inputs

- `m::AbstractModel`: model object. Only needed for the method in which `shocks`
  are not provided.
- `system::System{S}`: state-space system matrices. Alternatively, provide
  transition equation matrices `T`, `R`, `C`; measurement equation matrices `Q`,
  `Z`, `D`; and (possibly empty) pseudo-measurement equation matrices `Z_pseudo`
  and `D_pseudo`.
- `kal::Kalman{S}` or `z0::Vector{S}`: result of running the Kalman filter or
  state vector in the final historical period (aka initial forecast period)

where `S<:AbstractFloat`.

### Keyword Arguments

- `cond_type::Symbol`: one of `:none`, `:semi`, or `:full`, used to determine
  how many periods to forecast ahead. If `cond_type in [:semi, :full]`, the
  forecast horizon is reduced by the number of periods of conditional
  data. Defaults to `:none`.
- `enforce_zlb::Bool`: whether to enforce the zero lower bound. Defaults to
  `false`.
- `shocks::Matrix{S}`: matrix of size `nshocks` x `horizon` of shock innovations
  under which to forecast. If not provided, shocks are drawn according to:

  1. If `forecast_killshocks(m)`, `shocks` is set to a `nshocks` x `horizon`
     matrix of zeros
  2. Otherwise, if `forecast_tdist_shocks(m)`, draw `horizons` many shocks from a
     `Distributions.TDist(forecast_tdist_df_val(m))`
  3. Otherwise, draw `horizons` many shocks from a
     `DegenerateMvNormal(zeros(nshocks), sqrt(system[:QQ]))`

### Outputs

- `states::Matrix{S}`: matrix of size `nstates` x `horizon` of forecasted states
- `obs::Matrix{S}`: matrix of size `nobs` x `horizon` of forecasted observables
- `pseudo::Matrix{S}`: matrix of size `npseudo` x `horizon` of forecasted
  pseudo-observables. If `!forecast_pseudoobservables(m)` or the provided
  `Z_pseudo` and `D_pseudo` matrices are empty, then `pseudo` will be empty.
- `shocks::Matrix{S}`: matrix of size `nshocks` x `horizon` of shock innovations
"""
function compute_forecast{S<:AbstractFloat}(m::AbstractModel, system::System{S},
    kal::Kalman{S}; cond_type::Symbol = :none, enforce_zlb::Bool = false,
    shocks::Matrix{S} = Matrix{S}())

    draw_z0(kal::Kalman) = rand(DegenerateMvNormal(kal[:zend], kal[:Pend]))
    z0 = if forecast_draw_z0(m)
        draw_z0(kal)
    else
        kal[:zend]
    end

    compute_forecast(m, system, z0; cond_type = cond_type, enforce_zlb = enforce_zlb,
                     shocks = shocks)
end

function compute_forecast{S<:AbstractFloat}(m::AbstractModel, system::System{S},
    z0::Vector{S}; cond_type::Symbol = :none, enforce_zlb::Bool = false,
    shocks::Matrix{S} = Matrix{S}())

    # Numbers of things
    nshocks = n_shocks_exogenous(m)
    horizon = forecast_horizons(m; cond_type = cond_type)

    # Populate shocks matrix
    if isempty(shocks)
        if forecast_kill_shocks(m)
            shocks = zeros(S, nshocks, horizon)
        else
            μ = zeros(S, nshocks)
            σ = sqrt(system[:QQ])
            dist = if forecast_tdist_shocks(m)
                # Use t-distributed shocks
                ν = forecast_tdist_df_val(m)
                DegenerateDiagMvTDist(μ, σ, ν)
            else
                # Use normally distributed shocks
                DegenerateMvNormal(μ, σ)
            end

            shocks = rand(dist, horizon)

            # Forecast without anticipated shocks
            if n_anticipated_shocks(m) > 0
                ind_ant1 = m.exogenous_shocks[:rm_shl1]
                ind_antn = m.exogenous_shocks[symbol("rm_shl$(n_anticipated_shocks(m))")]
                ant_shock_inds = ind_ant1:ind_antn
                shocks[ant_shock_inds, :] = 0
            end
        end
    end

    # Get variables necessary to enforce the zero lower bound in the forecast
    ind_r = m.observables[:obs_nominalrate]
    ind_r_sh = m.exogenous_shocks[:rm_sh]
    zlb_value = forecast_zlb_value(m)

    compute_forecast(system, z0, shocks; enforce_zlb = enforce_zlb,
        ind_r = ind_r, ind_r_sh = ind_r_sh, zlb_value = zlb_value)
end


function compute_forecast{S<:AbstractFloat}(system::System{S}, z0::Vector{S},
    shocks::Matrix{S}; enforce_zlb::Bool = false, ind_r::Int = -1,
    ind_r_sh::Int = -1, zlb_value::S = 0.13/4)

    # Unpack system
    T, R, C = system[:TTT], system[:RRR], system[:CCC]
    Q, Z, D = system[:QQ], system[:ZZ], system[:DD]

    Z_pseudo, D_pseudo = if !isnull(system.pseudo_measurement)
        system[:ZZ_pseudo], system[:DD_pseudo]
    else
        Matrix{S}(), Vector{S}()
    end

    compute_forecast(T, R, C, Q, Z, D, Z_pseudo, D_pseudo, z0, shocks;
        enforce_zlb = enforce_zlb, ind_r = ind_r, ind_r_sh = ind_r_sh,
        zlb_value = zlb_value)
end

function compute_forecast{S<:AbstractFloat}(T::Matrix{S}, R::Matrix{S},
    C::Vector{S}, Q::Matrix{S}, Z::Matrix{S}, D::Vector{S}, Z_pseudo::Matrix{S},
    D_pseudo::Vector{S}, z0::Vector{S}, shocks::Matrix{S}; enforce_zlb::Bool = false,
    ind_r::Int = -1, ind_r_sh::Int = -1, zlb_value::S = 0.13/4)

    # Setup
    nshocks = size(R, 2)
    nstates = size(T, 2)
    nobs    = size(Z, 1)
    npseudo = size(Z_pseudo, 1)
    horizon = size(shocks, 2)

    # Define our iteration function
    function iterate(z_t1, ϵ_t)
        z_t = C + T*z_t1 + R*ϵ_t

        # Change monetary policy shock to account for 0.13 interest rate bound
        if enforce_zlb
            interest_rate_forecast = getindex(D + Z*z_t, ind_r)
            if interest_rate_forecast < zlb_value
                # Solve for interest rate shock causing interest rate forecast to be exactly ZLB
                ϵ_t[ind_r_sh] = 0.
                z_t = C + T*z_t1 + R*ϵ_t
                ϵ_t[ind_r_sh] = getindex((zlb_value - D[ind_r] - Z[ind_r, :]*z_t) / (Z[ind_r, :]*R[:, ind_r_sh]), 1)

                # Forecast again with new shocks
                z_t = C + T*z_t1 + R*ϵ_t

                # Confirm procedure worked
                interest_rate_forecast = getindex(D + Z*z_t, ind_r)
                @assert interest_rate_forecast >= zlb_value - 0.01
            end
        end
        return z_t, ϵ_t
    end

    # Iterate state space forward
    states = zeros(S, nstates, horizon)
    states[:, 1], shocks[:, 1] = iterate(z0, shocks[:, 1])
    for t in 2:horizon
        states[:, t], shocks[:, t] = iterate(states[:, t-1], shocks[:, t])
    end

    # Apply measurement and pseudo-measurement equations
    obs    = D .+ Z*states
    pseudo = if !isempty(Z_pseudo) && !isempty(D_pseudo)
        D_pseudo .+ Z_pseudo * states
    else
        Matrix{S}()
    end

    # Return forecasts
    return states, obs, pseudo, shocks
end
