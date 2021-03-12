# Some utility functions to keep code clean in statespace_functions.jl

function compute_gensys_gensys2_regimes(m::AbstractDSGEModel)

    regime_switching = haskey(get_settings(m), :regime_switching) && get_setting(m, :regime_switching)
    n_regimes        = (regime_switching && haskey(get_settings(m), :n_regimes)) ? get_setting(m, :n_regimes) : 1

    if haskey(get_settings(m), :gensys2) && get_setting(m, :gensys2)
        if haskey(get_settings(m), :regime_eqcond_info)
            sorted_eqcond = sort!(collect(get_setting(m, :regime_eqcond_info)), by = x -> x[1]) # x -> x[1] sorts by keys
            first_gensys2_ind = findfirst(x -> x[2].alternative_policy.key in # finds the first value (EqcondEntry) whose key
                                          get_setting(m, :temporary_altpolicy_names), sorted_eqcond) # is a temporary altpolicy
            first_gensys2_regime = !isnothing(first_gensys2_ind) ? sorted_eqcond[first_gensys2_ind][1] : nothing
        elseif haskey(get_settings(m), :gensys_regimes) && haskey(get_settings(m), :gensys2_regimes)
            return get_setting(m, :gensys_regimes), get_setting(m, :gensys2_regimes)
        else
            first_gensys2_regime = nothing
        end

        if isnothing(first_gensys2_regime)
            # throw(GensysError("No equilibrium conditions in any regime are being temporarily replaced, " *
            #                   "but the setting :gensys2 is true."))
            gensys2_regimes = Vector{UnitRange{Int}}(undef, 0)
            gensys_regimes  = UnitRange{Int}[1:n_regimes]
        else
            last_gensys2_regime = haskey(get_settings(m), :temporary_altpolicy_length) ?
                min(first_gensys2_regime + get_setting(m, :temporary_altpolicy_length), n_regimes) :
                n_regimes # NOTE removed a +1 here--if tests start failing, check here first

            gensys_regimes = UnitRange{Int}[1:(first_gensys2_regime - 1)]
            if last_gensys2_regime != n_regimes
                append!(gensys_regimes, [(last_gensys2_regime + 1):n_regimes])
            end
            gensys2_regimes = [first_gensys2_regime-1:last_gensys2_regime]
        end
    else
        gensys2_regimes = Vector{UnitRange{Int}}(undef, 0)
        gensys_regimes  = UnitRange{Int}[1:n_regimes]
    end

    return gensys_regimes, gensys2_regimes
end

"""
```
compute_system_function(system::System{S}) where S<:AbstractFloat
```

### Inputs

- `system::System`: The output of compute_system(m), i.e. the matrix outputs from solving a given model, m.

### Outputs

- `Φ::Function`: transition equation
- `Ψ::Function`: measurement equation
- `F_ϵ::Distributions.MvNormal`: shock distribution
- `F_u::Distributions.MvNormal`: measurement error distribution
"""
function compute_system_function(system::System{S}) where S<:AbstractFloat
    # Unpack system
    TTT    = system[:TTT]
    RRR    = system[:RRR]
    CCC    = system[:CCC]
    QQ     = system[:QQ]
    ZZ     = system[:ZZ]
    DD     = system[:DD]
    EE     = system[:EE]

    # Define transition and measurement functions
    @inline Φ(s_t1::Vector{S}, ϵ_t::Vector{S}) = TTT*s_t1 + RRR*ϵ_t + CCC
    @inline Ψ(s_t::Vector{S}) = ZZ*s_t + DD

    # Define shock and measurement error distributions
    nshocks = size(QQ, 1)
    nobs    = size(EE, 1)
    F_ϵ = Distributions.MvNormal(zeros(nshocks), QQ)
    F_u = Distributions.MvNormal(zeros(nobs),    EE)

    return Φ, Ψ, F_ϵ, F_u
end

function zero_system_constants(system::System{S}) where S<:AbstractFloat
    system = copy(system)

    system.transition.CCC = zeros(size(system[:CCC]))
    system.measurement.DD = zeros(size(system[:DD]))
    system.pseudo_measurement.DD_pseudo = zeros(size(system[:DD_pseudo]))

    return system
end

function zero_system_constants(system::RegimeSwitchingSystem{S}) where S<:AbstractFloat
    system = copy(system)

    for i in 1:n_regimes(system)
        system.transitions[i].CCC = zeros(size(system[i, :CCC]))
        system.measurements[i].DD = zeros(size(system[i, :DD]))
        system.pseudo_measurements[i].DD_pseudo = zeros(size(system[i, :DD_pseudo]))
    end

    return system
end

mutable struct ForwardExpectationsMemo{T}
    time_varying_memo::Dict{Int64, T} # for products of time-varying TTT matrices
    permanent_memo::Dict{Int64, T} # for products of non time-varying TTT matrices
end

"""
```
ForwardExpectationsMemo(TTTs::Vector{<: AbstractMatrix{S}},
                        first_tv_period::Int64, last_tv_period::Int64,
                        first_perm_period::Int64, min_perm_power::Int64 = 0,
                        max_perm_power::Int64 = 0) where {S <: Real}
```
computes the memo dictionaries of the necessary products/powers of TTTs
for computing forward expectations of states and sums of states.

### Inputs
- `TTTs`: complete sequence of time-varying transition matrices from regime 1 to the final regime
- `first_tv_period`: the first period with time-varying matrices
- `last_tv_period`: the last period in which the transition matrix is believed to have changed relative to the previous period.
- `first_perm_period`: the first period in which the transition matrix in `TTTs` is permanently imposed. This period
    should be from the perspective of an omniscient econometrician rather than the agents' perspective.
- `min_perm_power`: minimum power of the permanent matrix to be calculated, i.e. we calculate at least `TTTs[first_perm_period] ^ min_perm_power`
- `max_perm_power`: maximum power of the permanent matrix to be calculated, i.e. we calculate at most `TTTs[first_perm_period] ^ max_perm_power`

### Notes
- To clarify what `last_tv_period` should be, note that
if every matrix in `TTTs` is time-varying, AND agents believe every matrix is time-varying, then `last_tv_period = length(TTTs)`.
However, if agents are myopic and believe that the last time-varying matrix is the second to last one, then
`last_tv_period = length(TTTs) - 1`. This flexibility allows the user to specify different degrees of awareness about `TTTs`.

- Note that when `last_tv_period == first_perm_period`, we do not calculate

```
time_varying_memo[last_tv_period] = TTTs[first_perm_period] * TTTs[first_perm_period - 1] * ... TTTs[t + 1]
```

Instead, we calculate

```
time_varying_memo[last_tv_period] = TTTs[first_perm_period - 1] * ... TTTs[t + 1]
```

The reason is that it is more efficient to include that first `TTTs[first_perm_period]` in the powers of `TTT[first_perm_period]`
computed for `permanent_memo`. So if `length(TTTs) == 7`, to get the correct `k`-periods ahead forward
expectations from some period `t < 7`, you need to run

```
memo = ForwardExpectationsMemo(TTTs, t, 7, 7, 1, t + k + 1 - 7)
# or equivalently . . .
# memo = ForwardExpectationsMemo(TTTs, t, length(TTTs), length(TTTs), 1, t + k + 1 - length(TTTs))
```
"""
function ForwardExpectationsMemo(TTTs::Vector{<: AbstractMatrix{S}},
                                 first_tv_period::Int64, last_tv_period::Int64,
                                 first_perm_period::Int64, min_perm_power::Int64 = 0,
                                 max_perm_power::Int64 = 0) where {S <: Real}

    @assert first_tv_period <= last_tv_period
    @assert min_perm_power <= max_perm_power
    @assert last_tv_period <= first_perm_period

    tv_memo   = Dict{Int64, eltype(TTTs)}()
    perm_memo = Dict{Int64, eltype(TTTs)}()

    # When the last time-varying period equals the first period in which the
    # transition matrix is permanently imposed, we calculate the products
    # of the time-varying matrices
    if last_tv_period == first_perm_period
        last_tv_period -= 1
    end
    tv_memo[last_tv_period] = TTTs[last_tv_period]
    for k in (last_tv_period - 1):-1:(first_tv_period + 1)
        tv_memo[k] = tv_memo[k + 1] * TTTs[k]
    end

    if !(min_perm_power == max_perm_power == 0)
        TTT_perm                  = TTTs[first_perm_period]
        perm_memo[min_perm_power] = TTT_perm ^ min_perm_power
        for k in (min_perm_power + 1):max_perm_power
            perm_memo[k] = perm_memo[k - 1] * TTT_perm
        end
    end

    return ForwardExpectationsMemo(tv_memo, perm_memo)
end


mutable struct ForwardExpectedSumMemo{T}
    time_varying_memo::Dict{Int64, Dict{Int64, T}} # for products of time-varying TTT matrices
    permanent_memo::Dict{Int64, Dict{Int64, T}} # for products of non time-varying TTT matrices
end

function ForwardExpectedSum(TTTs::Vector{<: AbstractMatrix{S}},
                            first_tv_period::Int64, min_last_tv_period::Int64, max_last_tv_period::Int64,
                            first_perm_period::Int64, min_perm_power::Int64 = 0,
                            max_perm_power::Int64 = 0) where {S <: Real}

    @assert first_tv_period <= last_tv_period
    @assert min_perm_power <= max_perm_power
    @assert last_tv_period <= first_perm_period

    tv_memo   = Dict{Int64, eltype(TTTs)}()
    perm_memo = Dict{Int64, eltype(TTTs)}()

    # When the last time-varying period equals the first period in which the
    # transition matrix is permanently imposed, we calculate the products
    # of the time-varying matrices
    if last_tv_period == first_perm_period
        last_tv_period -= 1
    end
    tv_memo[last_tv_period] = TTTs[last_tv_period]
    for k in (last_tv_period - 1):-1:(first_tv_period + 1)
        tv_memo[k] = tv_memo[k + 1] * TTTs[k]
    end

    if !(min_perm_power == max_perm_power == 0)
        TTT_perm                  = TTTs[first_perm_period]
        perm_memo[min_perm_power] = TTT_perm ^ min_perm_power
        for k in (min_perm_power + 1):max_perm_power
            perm_memo[k] = perm_memo[k - 1] * TTT_perm
        end
    end

    return ForwardExpectationsMemo(tv_memo, perm_memo)
end


"""
```
k_periods_ahead_expectations(TTT, CCC, TTTs, CCCs, t, k; permanent_t = length(TTTs),
                             integ_series = false, memo = nothing)
```

calculates the matrices associated with the expected state `k` periods ahead from `t`.
This function should NOT be used with linear state space system matrices with any unit roots.

The `TTT` and `CCC` inputs are the transition matrix and constant vector associated with
the current period `t`, while the `TTTs` and `CCCs` are vectors containing the time-varying
transition matrices and constant vectors, such that `TTTs[t]` retrieves the time-varying
transition matrix associated with period `t` and `TTTs[t + k]` retrieves the time-varying
transition matrix associated with period `t + k`. The optional argument `permanent_t`
indicates the period for which the matrices/vectors are no longer time-varying, i.e.
if `t >= permanent_t`, then `TTTs[permanent_t]` is the transition matrix.

The formula implemented by this function is
```
𝔼ₜ[sₜ₊ₖ] = (∏ⱼ=₁ᵏ Tₜ₊ⱼ) sₜ + (∑ₘ₌₁ᵏ⁻¹ (∏ⱼ₌ₘ₊₁ᵏ Tₜ₊ⱼ) Cₜ₊ₘ) + Cₜ₊ₖ.
```
Additional simplifications are made if it is known that `t + k > permanent_t`
since this implies some matrices are the same. This recognition reduces
unnecessary computations.

### Keyword Arguments
- `integ_series::Bool`: set to true if there are some transition matries in `TTT`
    that result in integrated series, in which case we cannot speed up
    computations by using left-divides.

- `memo::Union{ForwardExpectationsMemo, Nothing}`: pass a properly formed
    `ForwardExpectationsMemo` to avoid calculating unnecessary products and powers of
    the matrices in `TTTs`. Typically, the memo you want to compute is

```
# min_t is minimum t you will use, maximum_t is maximum t you will use, and
# max_k is the maximum window for forward expectations.
memo = ForwardExpectationsMemo(TTTs, min_t, length(TTTs), length(TTTs), min_t + max_k - length(TTTs),
                               max_t + max_k + 1 - length(TTTs))
```
"""
function k_periods_ahead_expectations(TTT::AbstractMatrix, CCC::AbstractVector,
                                      TTTs::Vector{<: AbstractMatrix}, CCCs::Vector{<: AbstractVector},
                                      t::Int, k::Int, permanent_t::Int = length(TTTs);
                                      integ_series::Bool = false,
                                      memo::Union{ForwardExpectationsMemo, Nothing} = nothing)

    if isempty(TTTs) || isempty(CCCs)
        Tᵏ = TTT^k
        if all(CCC .≈ 0.)
            return Tᵏ, CCC
        else
            if integ_series # are there integrated series? If so, we have to accumulate instead of doing a left-division
                T_memo = Dict{Int, typeof(TTT)}()
                T_memo[1] = TTT
                for i in 2:k
                    T_memo[i] = T_memo[i - 1] * TTT
                end
                Tᵏsum = sum([T_memo[j] for j in 1:k])
            else
                Tᵏsum = (I - TTT) \ (I - Tᵏ)
            end
            return Tᵏ, Tᵏsum * CCC
        end
    else
        if t + k <= permanent_t
            # Cannot save computation speed by not calculating further times b/c always time-varying
            t_plus_k = t + k
            if isnothing(memo)
                T_memo = Dict{Int, eltype(TTTs)}()
                T_memo[k] = TTTs[t_plus_k]
                for i in (k-1):-1:1
                    T_memo[i] = T_memo[i + 1] * TTTs[t + i]
                end
            end

            C_accum = deepcopy(CCCs[t_plus_k])
            for i in 1:(k - 1)
                if isnothing(memo)
                    C_accum .+= T_memo[i + 1] * CCCs[t + i]
                elseif t_plus_k == permanent_t
                    # This is an edge case b/c the memo computation will not add
                    # the permanent matrix to the product (and instead stores it
                    # in the memo.permanent_memo dict)
                    if i == k - 1
                        C_accum .+= TTTs[permanent_t] * CCCs[t + i]
                    else
                        C_accum .+= TTTs[permanent_t] * memo.time_varying_memo[t + (i + 1)] * CCCs[t + i]
                    end
                else
                    C_accum .+= memo.time_varying_memo[t + (i + 1)] * CCCs[t + i]
                end
            end

            if isnothing(memo)
                return T_memo[1], C_accum
            elseif t_plus_k == permanent_t
                return TTTs[permanent_t] * memo.time_varying_memo[t + 1], C_accum
            else
                return memo.time_varying_memo[t + 1], C_accum
            end
        elseif t == permanent_t
            # None of the matrices are time-varying anymore
            if all(CCCs[permanent_t] .≈ 0.)
                # TODO: we may want to add a way to figure out if the power requested here
                # could be pre-computed with the memo.
                return TTTs[permanent_t]^k, CCCs[permanent_t]
            else
                if isnothing(memo)
                    Tᵏₜ₊₁ = TTTs[permanent_t]^k
                else
                    Tᵏₜ₊₁ = memo.permanent_memo[k]
                end
                if integ_series
                    T_memo = Dict{Int, eltype(TTTs)}()
                    T_memo[1] = TTTs[permanent_t]
                    for i in 2:k
                        T_memo[i] = T_memo[i - 1] * TTTs[permanent_t]
                    end
                    Tᵏsum = sum([T_memo[j] for j in 1:k])
                else
                    Tᵏsum = (I - TTTs[permanent_t]) \ (I - Tᵏₜ₊₁)
                end
                return Tᵏₜ₊₁, Tᵏsum * CCCs[permanent_t]
            end
        else
            # Computation time can be saved by realizing some matrices are not time-varying
            h = (permanent_t - 1) - t # last time of time-variation is permanent_t - 1
            if isnothing(memo)
                Tᵏ⁻ʰₜ₊ₕ₊₁ = (k == h) ? Diagonal(ones(length(CCCs[permanent_t]))) : TTTs[permanent_t]^(k - h) # why using Diagonal instead of I?
            else
                # Tᵏ⁻ʰₜ₊ₕ₊₁ = (k == h) ? I : memo.permanent_memo[k - h] # why using Diagonal instead of I?
                Tᵏ⁻ʰₜ₊ₕ₊₁ = (k == h) ? Diagonal(ones(length(CCCs[permanent_t]))) : memo.permanent_memo[k - h] # why using Diagonal instead of I?
            end

            if isnothing(memo)
                T_memo = Dict{Int, eltype(TTTs)}()
            end

            if h > 0
                if isnothing(memo)
                    T_memo[h] = TTTs[t + h] # maps i to ∏ⱼ₌ᵢʰ Tₜ₊ⱼ, so T_memo[h] = Tₜ₊ₕ, T_memo[h-1] = Tₜ₊ₕ * Tₜ₊ₕ₋₁, ...
                    for i in (h-1):-1:1
                        T_memo[i] = T_memo[i + 1] * TTTs[t + i]
                    end
                    T_accum = Tᵏ⁻ʰₜ₊ₕ₊₁ * T_memo[1]
                else
                    T_accum = Tᵏ⁻ʰₜ₊ₕ₊₁ * memo.time_varying_memo[t + 1]
                end
            else
                T_accum = Tᵏ⁻ʰₜ₊ₕ₊₁
            end

            if h == 0 # Nothing to accumulate from the past
                C_accum = zeros(size(T_accum, 1))
            else
                C_accum = deepcopy(CCCs[t + h])
                for i in 1:(h - 1)
                    if isnothing(memo)
                        C_accum .+= T_memo[i + 1] * CCCs[t + i]
                    else
                        C_accum .+= memo.time_varying_memo[t + i + 1] * CCCs[t + i]
                    end
                end
                C_accum .= Tᵏ⁻ʰₜ₊ₕ₊₁ * C_accum
            end

            if all(CCCs[permanent_t] .≈ 0.)
                return T_accum, C_accum
            else
                if integ_series
                    T_memo = Dict{Int, eltype(TTTs)}()
                    T_memo[1] = TTTs[permanent_t]
                    for i in 2:(k - h - 1)
                        T_memo[i] = T_memo[i - 1] * TTTs[permanent_t]
                    end
                    Tᵏ⁻ᵐsum = sum([T_memo[j] for j in 1:(k - h - 1)]) + I
                else
                    Tᵏ⁻ᵐsum = (I - TTTs[permanent_t]) \ (I - Tᵏ⁻ʰₜ₊ₕ₊₁)
                end
                C_accum .+= Tᵏ⁻ᵐsum * CCCs[permanent_t]
                return T_accum, C_accum
            end
        end
    end
end

"""
```
k_periods_ahead_expected_sums(TTT, CCC, TTTs, CCCs, t, k; permanent_t = length(TTTs),
                              integ_series = false, memo = nothing)
```

calculates the matrices associated with the sum of the expected states between periods
`t + 1` and  `t + k`. This function should NOT be used with
linear state space system matrices with any unit roots.

The `TTT` and `CCC` inputs are the transition matrix and constant vector associated with
the current period `t`, while the `TTTs` and `CCCs` are vectors containing the time-varying
transition matrices and constant vectors, such that `TTTs[t]` retrieves the time-varying
transition matrix associated with period `t` and `TTTs[t + k]` retrieves the time-varying
transition matrix associated with period `t + k`. The optional argument `permanent_t`
indicates the period for which the matrices/vectors are no longer time-varying, i.e.
if `t >= permanent_t`, then `TTTs[permanent_t]` is the transition matrix.

The formula implemented by this function is
```
∑ⱼ₌₁ᵏ 𝔼ₜ[sₜ₊ⱼ] = ∑ⱼ₌₁ᵏ(∏ⱼ=₁ᵏ Tₜ₊ⱼ) sₜ + ∑ᵣ₌₁ᵏ⁻¹(I + ∑ⱼ₌ᵣ₊₁ᵏ (∏ₘ₌ᵣ₊₁ʲ Tₜ₊ₘ))Cₜ₊ᵣ + Cₜ₊ₖ.
```
Additional simplifications are made if it is known that `t + k > permanent_t`
since this implies some matrices are the same. This recognition reduces
unnecessary computations.

### Keyword Arguments
- `integ_series::Bool`: set to true if there are some transition matries in `TTT`
    that result in integrated series, in which case we cannot speed up
    computations by using left-divides.

- `memo::Union{ForwardExpectedSumMemo, Nothing}`: pass a properly formed
    `ForwardExpectationsMemo` to avoid calculating unnecessary products and powers of
    the matrices in `TTTs`. Typically, the memo you want to compute is

```
# min_t is minimum t you will use, maximum_t is maximum t you will use, and
# max_k is the maximum window for forward expectations.
memo = ForwardExpectedSumMemo(TTTs, min_t, length(TTTs), length(TTTs), min_t + max_k - length(TTTs),
                               max_t + max_k + 1 - length(TTTs))
```
"""
function k_periods_ahead_expected_sums(TTT::AbstractMatrix, CCC::AbstractVector,
                                       TTTs::Vector{<: AbstractMatrix}, CCCs::Vector{<: AbstractVector},
                                       t::Int, k::Int, permanent_t::Int = length(TTTs);
                                       integ_series::Bool = false,
                                       memo::Union{ForwardExpectedSumMemo, Nothing} = nothing)

    if integ_series # Do this by directly summing the k-periods ahead expectations. Not fully efficient but also not the typical case
        T_accum = Vector{eltype(TTTs)}(undef, k)
        C_accum = Vector{eltype(CCCs)}(undef, k)
        for i in 1:k
            T_accum[i], C_accum[i] = k_periods_ahead_expectations(TTT, CCC, TTTs, CCCs, t, i,
                                                                  integ_series = integ_series, memo = memo)
        end
        return sum(T_accum), sum(C_accum)
    end

    if isempty(TTTs) || isempty(CCCs)
        Tᵏsum = (I - TTT) \ (TTT - TTT^(k + 1))
        if all(CCC .≈ 0.)
            return Tᵏsum, CCC
        else
            TTTʲmemo = Vector{typeof(TTT)}(undef, k)
            TTTʲmemo[1] = TTT
            for q in 2:k
                TTTʲmemo[q] = TTTʲmemo[q - 1] * TTT
            end
            return Tᵏsum, (I + sum([(I - TTT) \ (I - TTTʲmemo[k - q + 1]) for q in 1:(k - 1)])) * CCC
        end
    else
        if t + k <= permanent_t
            # Cannot save computation speed by not calculating further times b/c always time-varying
            total_Tsum = zeros(eltype(TTTs[t]), size(TTTs[t]))
            total_Csum = zeros(eltype(CCCs[t]), size(CCCs[t]))

            # T_accum and C_accum store expected values in period i for some i
            T_accum = TTTs[t+k]
            C_accum = CCCs[t+k]
            total_Tsum .+= T_accum
            total_Csum .+= C_accum

            # Time-varying matrices so need to calculate in loop
            for i in (t+2):(t+k)
                T_accum = TTTs[i] * T_accum
                total_Tsum .+= T_accum

                C_accum = TTTs[i] * C_accum + CCCs[i]
                total_Csum .+= C_accum
            end
#=            for i in 1:k
                tmp1, tmp2 = k_periods_ahead_expectations(TTT, CCC, TTTs, CCCs, t, i, permanent_t,
                                                          integ_series = integ_series, memo = memo)
                total_Tsum .+= tmp1
                total_Csum .+= tmp2
            end=#

            return total_Tsum, total_Csum

        elseif t == permanent_t
            # None of the matrices are time-varying anymore
            Tᵏsum = (I - TTTs[permanent_t]) \ (TTTs[permanent_t] - TTTs[permanent_t]^(k + 1))
            if all(CCC .≈ 0.)
                return Tᵏsum, CCCs[permanent_t]
            else
                TTTʲmemo = Vector{eltype(TTTs)}(undef, k)
                TTTʲmemo[1] = TTTs[permanent_t]
                for q in 2:k
                    TTTʲmemo[q] = TTTʲmemo[q - 1] * TTTs[permanent_t]
                end
                return Tᵏsum, (I + sum([(I - TTTs[permanent_t]) \ (I - TTTʲmemo[k - q + 1]) for q in 1:(k - 1)])) * CCCs[permanent_t]
            end
        else
            total_Tsum = zeros(eltype(TTTs[t]), size(TTTs[t]))
            total_Csum = zeros(eltype(CCCs[t]), size(CCCs[t]))

            # Do we need to adjust for C after permanent_t?
            do_C = !(all(CCCs[permanent_t] .≈ 0.0))

            # T_accum and C_accum store the expectation of T and C at time i
            ## So, T_accum at time i will store ∏_{j=t+1}^i T_j
            ## C_memo at time i stores ∑_{q=t+1}^i C_q + ∑_{j=t+1}^q T_j C_q
            # total_sum will stores sum over the expectations

            T_accum = TTTs[t+1]
            C_accum = CCCs[t+1]
            total_Tsum .+= T_accum
            total_Csum .+= C_accum

            # Time-varying matrices so need to calculate in loop
            for i in (t+2):permanent_t
                T_accum = TTTs[i] * T_accum
                total_Tsum .+= T_accum

                C_accum = TTTs[i] * C_accum + CCCs[i]
                total_Csum .+= C_accum
            end

            # Post permanent_t, fixed matrices, so can compute here
            T_sums = (I - TTTs[permanent_t]) \ (TTTs[permanent_t] - TTTs[permanent_t]^(t+k-permanent_t+1))
            total_Tsum .+= T_sums * T_accum

            # If C unchanging after permanent_t, can compute with sum expression
            ## else do for loop
            if !do_C
                C_tmp = T_sums * C_accum
                total_Csum .+= C_tmp
            else
                for i in (permanent_t+1):(t+k)
                    C_accum = TTTs[permanent_t] * C_accum + CCCs[permanent_t]
                    total_Csum .+= C_accum
                end
#=            # Computation time can be saved by realizing some matrices are not time-varying
            T_accum = Vector{eltype(TTTs)}(undef, k)
            C_accum = Vector{eltype(CCCs)}(undef, k)
            for i in 1:k
                T_accum[i], C_accum[i] = k_periods_ahead_expectations(TTT, CCC, TTTs, CCCs, t, i,
                                                                      integ_series = integ_series, memo = memo)=#
            end

            return total_Tsum, total_Csum

            # This code seems to work (match directly adding each k_period_ahead_expectations)
            # but sometimes causes mild numerical differences,
            # so we calculate the expected sum by adding each k_period_ahead_expectations.
            #=            h = (permanent_t - 1) - t # last time of time-variation is permanent_t - 1

            Tₜ₊ₕ₊₁_memo = Vector{eltype(TTTs)}(undef, k - h)
            Tₜ₊ₕ₊₁_memo[1] = TTTs[permanent_t] # maps j to Tₜ₊ₕ₊₁ʲ⁻ʰ for j in (h + 1):k or Tₜ₊ₕ₊₁ʲ for j in 1:(k-h)
            for j in 2:(k - h)
            Tₜ₊ₕ₊₁_memo[j] = TTTs[permanent_t] * Tₜ₊ₕ₊₁_memo[j - 1]
            end

            Tₜ₊ₕ₊₁ʲsum = (I - TTTs[permanent_t]) \ (TTTs[permanent_t] - TTTs[permanent_t] ^ (k - h + 1)) # ∑ⱼ₌₁ᵏ⁻ʰ Tₜ₊ₕ₊₁ʲ

            TC_memo = Vector{eltype(TTTs)}(undef, h) # matrices to be used for both accumulated TTT and CCC
            TC_memo[h] = TTTs[t + h] # maps i to ∏ⱼ₌ᵢʰ Tₜ₊ⱼ, so T_memo[h] = Tₜ₊ₕ, T_memo[h-1] = Tₜ₊ₕ * Tₜ₊ₕ₋₁, ...
            for i in (h-1):-1:1
            TC_memo[i] = TC_memo[i + 1] * TTTs[t + i]
            end

            T1_term = Vector{eltype(TTTs)}(undef, h) # first part of the accumulated TTT matrix
            T1_term[1] = TTTs[t + 1] # maps i to ∏ⱼ₌₁ⁱ Tₜ₊ⱼ, so T_memo[h] = Tₜ₊ₕ * ⋯  * Tₜ₊₁
            for i in 2:(h - 1)
            T1_term[i] = TTTs[t + i] * T1_term[i - 1]
            end
            T1_term[h] = TC_memo[1] # This one was calculated already

            # second part of the accumulated TTT matrix
            # maps j to Tₜ₊ₕ₊₁ʲ ∏ᵣ₌₁ʰ Tₜ₊ᵣ, so T_memo[h] = Tₜ₊ₕ * ⋯  * Tₜ₊₁
            T2_term = [Tₜ₊ₕ₊₁_memo[j] * T1_term[h] for j in 1:(k - h)]

            T_accum = sum(T1_term) + sum(T2_term) # calculated accumulated matrix

            # Calculate final portions of accumulated CCC vector first
            I₊Tʲ = (I + Tₜ₊ₕ₊₁ʲsum)
            C_accum   = I₊Tʲ * CCCs[t + h]
            if any(.!(CCCs[permanent_t] .≈ 0.))
            C_accum .+= (I + (k - h - 1) .* I₊Tʲ) * CCCs[permanent_t]
            end
            if h > 1
            C1_term = Vector{eltype(TTTs)}(undef, h - 1)
            C2_term = Vector{eltype(TTTs)}(undef, h - 1)

            for q in 1:(h - 1)
            C1_term[q] = sum([prod([TTTs[t + m] for m in (q + 1):j]) for j in (q + 1):h])
            C2_term[q] = sum([Tₜ₊ₕ₊₁_memo[j] * TC_memo[q + 1] for j in 1:(k - h)])
            end
            C_accum .+= sum([(I + C1_term[q] + C2_term[q]) * CCCs[t + q] for q in 1:(h - 1)])
            end=#
        end
    end
end
