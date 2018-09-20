"""
`Transition{T<:AbstractFloat}`

The transition equation of the state-space model takes the form

    `s_t = TTT*s_{t-1} + RRR*ϵ_t + CCC`

The `Transition` type stores the coefficient `Matrix{T}`s (`TTT`, `RRR`) and constant `Vector{T} CCC`.
"""
mutable struct Transition{T<:AbstractFloat}
    TTT::Matrix{T}
    RRR::Matrix{T}
    CCC::Vector{T}
end
function Transition(TTT::Matrix{T}, RRR::Matrix{T}) where T <: AbstractFloat
    CCC = zeros(eltype(TTT), size(TTT, 1))
    Transition{T}(TTT, RRR, CCC)
end
function Transition(TTT::Matrix{T}, RRR::Matrix{T}, CCC::Matrix{T}) where T <: AbstractFloat
    Transition{T}(TTT, RRR, collect(CCC))
end
function Base.getindex(eq::Transition, d::Symbol)
    if d in (:TTT, :RRR, :CCC)
        return getfield(eq, d)
    else
        throw(KeyError(d))
    end
end

"""
`Measurement{T<:AbstractFloat}`

The measurement equation of the state-space model takes the form

   `y_t = ZZ*s_t + DD + u_t`

where the error `u_t` is the measurement error, which is uncorrelated with the
shocks in the transition equation `ϵ_t`.

### Fields

If `Ns` is the number of states `s_t`, `Ny` is the number of
observables `y_t`, and `Ne` is the number of shocks `ϵ_t`:

- `ZZ`: the `Ny` x `Ns`  measurement matrix
- `DD`: the `Ny` x 1 constant vector
- `QQ`: the `Ne` x `Ne` covariance matrix for the shocks `ϵ_t`
- `EE`: the `Ny` x `Ny` covariance matrix for the measurement error `η_t`
"""
mutable struct Measurement{T<:AbstractFloat}
    ZZ::Matrix{T}
    DD::Vector{T}
    QQ::Matrix{T}
    EE::Matrix{T}
end

function Base.getindex(M::Measurement, d::Symbol)
    if d in (:ZZ, :DD, :QQ, :EE)
        return getfield(M, d)
    else
        throw(KeyError(d))
    end
end

function measurement(m::AbstractModel, trans::Transition; shocks::Bool=true)
    TTT = trans[:TTT]
    RRR = trans[:RRR]
    CCC = trans[:CCC]
    measurement(m, TTT, RRR, CCC; shocks=shocks)
end

"""
```
PseudoMeasurement{T<:AbstractFloat}
```

The pseudo-measurement equation of the state-space model takes the form

   `x_t = ZZ_pseudo*s_t + DD_pseudo`

### Fields

Let `Ns` be the number of states `s_t` and `Nx` be the number of
pseudo-observables `x_t`:

- `ZZ_pseudo`: the `Nx` x `Ns` pseudo-measurement matrix
- `DD_pseudo`: the `Nx` x 1 constant vector
"""
mutable struct PseudoMeasurement{T<:AbstractFloat}
    ZZ_pseudo::Matrix{T}
    DD_pseudo::Vector{T}
end

function Base.getindex(M::PseudoMeasurement, d::Symbol)
    if d in (:ZZ_pseudo, :DD_pseudo)
        return getfield(M, d)
    else
        throw(KeyError(d))
    end
end

"""
`System{T<:AbstractFloat}`

A type containing the transition and measurement equations for a
state-space model. The matrices may be directly indexed: `sys[:TTT]`
returns `sys.transition.TTT`, etc.
"""
mutable struct System{T<:AbstractFloat}
    transition::Transition{T}
    measurement::Measurement{T}
    pseudo_measurement::PseudoMeasurement{T}
end

function System(transition::Transition{T}, measurement::Measurement{T}) where T <: AbstractFloat
    # Initialize empty pseudo-measurement equation
    _n_states = size(transition.TTT, 1)
    _n_pseudo = 0
    ZZ_pseudo = zeros(_n_pseudo, _n_states)
    DD_pseudo = zeros(_n_pseudo)
    pseudo_measurement = PseudoMeasurement(ZZ_pseudo, DD_pseudo)

    return System(transition, measurement, pseudo_measurement)
end

function Base.getindex(system::System, d::Symbol)
    if d in (:transition, :measurement, :pseudo_measurement)
        return getfield(system, d)
    elseif d in fieldnames(typeof(system.transition))
        return getfield(system.transition, d)
    elseif d in fieldnames(typeof(system.measurement))
        return getfield(system.measurement, d)
    elseif d in fieldnames(typeof(system.pseudo_measurement))
        return getfield(system.pseudo_measurement, d)
    else
        throw(KeyError(d))
    end
end

function Base.copy(system::System)
    trans = Transition(system[:TTT], system[:RRR], system[:CCC])
    meas  = Measurement(system[:ZZ], system[:DD], system[:QQ], system[:EE])
    pseudo_meas = PseudoMeasurement(system[:ZZ_pseudo], system[:DD_pseudo])
    return System(trans, meas, pseudo_meas)
end

"""
```
compute_system(m; apply_altpolicy = false)
```

Given the current model parameters, compute the state-space system
corresponding to model `m`. Returns a `System` object.
"""
function compute_system(m::AbstractModel{T}; apply_altpolicy = false) where T <: AbstractFloat
    # Solve model
    TTT, RRR, CCC = solve(m; apply_altpolicy = apply_altpolicy)
    transition_equation = Transition(TTT, RRR, CCC)

    # Solve measurement equation
    measurement_equation = measurement(m, TTT, RRR, CCC)

    type_tuple = (typeof(m), Matrix{T}, Matrix{T}, Vector{T})
    if hasmethod(pseudo_measurement, type_tuple)
        # Solve pseudo-measurement equation
        pseudo_measurement_equation = pseudo_measurement(m, TTT, RRR, CCC)
        return System(transition_equation, measurement_equation, pseudo_measurement_equation)
    else
        return System(transition_equation, measurement_equation)
    end
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
function compute_system_function(system::System{S}) where S <: AbstractFloat
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

function zero_system_constants(system::System{S}) where S <: AbstractFloat
    system = copy(system)

    system.transition.CCC = zeros(size(system[:CCC]))
    system.measurement.DD = zeros(size(system[:DD]))
    system.pseudo_measurement.DD_pseudo = zeros(size(system[:DD_pseudo]))

    return system
end
