"""
`Transition{T<:Real}`

The transition equation of the state-space model takes the form

    `s_t = TTT*s_{t-1} + RRR*ϵ_t + CCC`

The `Transition` type stores the coefficient `Matrix{T}`s (`TTT`, `RRR`) and constant `Vector{T} CCC`.
"""
mutable struct Transition{T<:Real}
    TTT::Matrix{T}
    RRR::Matrix{T}
    CCC::Vector{T}
end
function Transition(TTT::Matrix{T}, RRR::Matrix{T}) where T<:Real
    CCC = zeros(eltype(TTT), size(TTT, 1))
    Transition{T}(TTT, RRR, CCC)
end
function Transition(TTT::Matrix{T}, RRR::Matrix{T}, CCC::Matrix{T}) where T<:Real
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
`Measurement{T<:Real}`

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
mutable struct Measurement{T<:Real}
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

function measurement(m::AbstractDSGEModel, trans::Transition; shocks::Bool=true)
    TTT = trans[:TTT]
    RRR = trans[:RRR]
    CCC = trans[:CCC]
    measurement(m, TTT, RRR, CCC; shocks=shocks)
end

"""
```
PseudoMeasurement{T<:Real}
```

The pseudo-measurement equation of the state-space model takes the form

   `x_t = ZZ_pseudo*s_t + DD_pseudo`

### Fields

Let `Ns` be the number of states `s_t` and `Nx` be the number of
pseudo-observables `x_t`:

- `ZZ_pseudo`: the `Nx` x `Ns` pseudo-measurement matrix
- `DD_pseudo`: the `Nx` x 1 constant vector
"""
mutable struct PseudoMeasurement{T<:Real}
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

abstract type AbstractSystem{T} end

"""
`System{T <: Real}`

A mutable struct containing the transition and measurement equations for a
state-space model. The matrices may be directly indexed: `sys[:TTT]`
returns `sys.transition.TTT`, `sys[:ZZ]` returns `sys.measurement.ZZ`, etc.
"""
mutable struct System{T <: Real} <: AbstractSystem{T}
    transition::Transition{T}
    measurement::Measurement{T}
    pseudo_measurement::PseudoMeasurement{T}
end

function System(transition::Transition{T}, measurement::Measurement{T}) where T<:Real
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
`RegimeSwitchingSystem{T <: Real}`

A mutable struct containing the transition and measurement equations for a
state-space model with regime-switching.
The matrices may be directly indexed: `sys[1, :TTT]`
returns `sys.regime[1].transition.TTT`, etc.
"""
mutable struct RegimeSwitchingSystem{T <: Real} <: AbstractSystem{T}
    transitions::Vector{Transition{T}}
    measurements::Vector{Measurement{T}}
    pseudo_measurements::Vector{PseudoMeasurement{T}}

    RegimeSwitchingSystem{T}(transitions, measurements, pseudo_measurements) where {T <: Real} = length(transitions) == length(measurements) ?
        new(transitions, measurements, pseudo_measurements) :
        error("The number of Transition (n = $(length(transitions)))" *
              " and Measurement (n = $(length(measurements))) objects must match.")
end

RegimeSwitchingSystem(transitions::Vector{Transition{T}}, measurements::Vector{Measurement{T}},
                      pseudo_measurements::Vector{PseudoMeasurement{T}}) where {T <: Real} =
                          RegimeSwitchingSystem{T}(transitions, measurements, pseudo_measurements)

function RegimeSwitchingSystem(transitions::Vector{Transition{T}},
                               measurements::Vector{Measurement{T}}) where {T<:Real}

    # Initialize empty pseudo-measurement equation
    _n_states = size(transitions[1].TTT, 1)
    _n_pseudo = 0
    ZZ_pseudo = zeros(_n_pseudo, _n_states)
    DD_pseudo = zeros(_n_pseudo)
    pseudo_measurement = PseudoMeasurement(ZZ_pseudo, DD_pseudo)
    pseudo_measurements = [pseudo_measurement for i in 1:length(transitions)]

    return RegimeSwitchingSystem(transitions, measurements, pseudo_measurements)
end

function RegimeSwitchingSystem(systems::Vector{System{T}}) where {T <: Real}
    n_regimes           = length(systems)
    transitions         = Vector{Transition{T}}(undef, n_regimes)
    measurements        = Vector{Measurement{T}}(undef, n_regimes)
    pseudo_measurements = Vector{PseudoMeasurement{T}}(undef, n_regimes)
    for i in 1:n_regimes
        transitions[i]         = systems[i].transition
        measurements[i]        = systems[i].measurement
        pseudo_measurements[i] = systems[i].pseudo_measurement
    end
    return RegimeSwitchingSystem(transitions, measurements, pseudo_measurements)
end

# Retrieve System correspond to a regime
function System(system::RegimeSwitchingSystem, regime::Int)
    return System(system.transitions[regime],
                  system.measurements[regime],
                  system.pseudo_measurements[regime])
end

# Get fields of the RegimeSwitchingSystem
function Base.getindex(system::RegimeSwitchingSystem{T},
                       d::Symbol) where {T<:Real}
    if d in (:transitions, :measurements, :pseudo_measurements)
        return getfield(system, d)
    elseif d == :regimes
        return 1:length(system.transitions)
    else
        throw(KeyError(d))
    end
end

# Get a specific regime
function Base.getindex(system::RegimeSwitchingSystem{T},
                       d::Int) where {T<:Real}
    if d < 1 || d > length(system.transitions)
        throw(BoundsError(system.transitions,d))
    else
        return System(system, d)
    end
end

# Get specific matrix or type from specific regime
function Base.getindex(system::RegimeSwitchingSystem{T},
                       regime::Int, d::Symbol) where {T <: Real}
    if d == :transition
        system.transitions[regime]
    elseif d == :measurement
        system.measurements[regime]
    elseif d == :pseudo_measurement
        system.pseudo_measurements[regime]
    elseif d in (:TTT, :RRR, :CCC)
        system.transitions[regime][d]
    elseif d in (:ZZ, :DD, :QQ, :EE)
        system.measurements[regime][d]
    elseif d in (:ZZ_pseudo, :DD_pseudo)
        system.pseudo_measurements[regime][d]
    else
        throw(KeyError(d))
    end
end

function n_regimes(system::RegimeSwitchingSystem{T}) where {T<:Real}
    return length(system.transitions)
end

function Base.copy(system::RegimeSwitchingSystem{T}) where {T<:Real}
    n_reg = n_regimes(system)
    transitions         = Vector{Transition{T}}(undef, n_reg)
    measurements        = Vector{Measurement{T}}(undef, n_reg)
    pseudo_measurements = Vector{PseudoMeasurement{T}}(undef, n_reg)

    for i in 1:length(system[:transitions])
        trans = Transition(system[i,:TTT], system[i,:RRR], system[i,:CCC])
        meas  = Measurement(system[i,:ZZ], system[i,:DD],
                            system[i,:QQ], system[i,:EE])
        pseudo_meas = PseudoMeasurement(system[i,:ZZ_pseudo], system[i,:DD_pseudo])

        transitions[i]         = trans
        measurements[i]        = meas
        pseudo_measurements[i] = pseudo_meas
    end
    return RegimeSwitchingSystem(transitions, measurements, pseudo_measurements)
end

function Base.deepcopy(system::RegimeSwitchingSystem{T}) where {T<:Real}
    n_reg = n_regimes(system)
    transitions         = Vector{Transition{T}}(undef, n_reg)
    measurements        = Vector{Measurement{T}}(undef, n_reg)
    pseudo_measurements = Vector{PseudoMeasurement{T}}(undef, n_reg)

    for i in 1:length(system[:transitions])
        trans = Transition(deepcopy(system[i,:TTT]), deepcopy(system[i,:RRR]), deepcopy(system[i,:CCC]))
        meas  = Measurement(deepcopy(system[i,:ZZ]), deepcopy(system[i,:DD]),
                            deepcopy(system[i,:QQ]), deepcopy(system[i,:EE]))
        pseudo_meas = PseudoMeasurement(deepcopy(system[i,:ZZ_pseudo]), deepcopy(system[i,:DD_pseudo]))

        transitions[i]         = trans
        measurements[i]        = meas
        pseudo_measurements[i] = pseudo_meas
    end
    return RegimeSwitchingSystem(transitions, measurements, pseudo_measurements)
end

"""
`TimeVaryingInformationSetSystem`

A mutable struct containing the transition and measurement equations for a
state-space model with regime-switching and time-varying information sets.
"""
mutable struct TimeVaryingInformationSetSystem{T <: Real} <: AbstractSystem{T}
    transitions::Vector{Vector{Transition{T}}}
    measurements::Vector{Measurement{T}}
    pseudo_measurements::Vector{PseudoMeasurement{T}}
    information_set::Vector{UnitRange{Int}}
    select::Vector{Int}
end

function TimeVaryingInformationSetSystem(transitions::Vector{Vector{Transition{T}}},
                                         measurements::Vector{Measurement{T}},
                                         information_set::Vector{UnitRange{Int}},
                                         select::Vector{Int}) where {T <: Real}

    # Initialize empty pseudo-measurement equation
    _n_states = size(transitions[select[1]][1].TTT, 1)
    _n_pseudo = 0
    ZZ_pseudo = zeros(_n_pseudo, _n_states)
    DD_pseudo = zeros(_n_pseudo)
    pseudo_measurement = PseudoMeasurement(ZZ_pseudo, DD_pseudo)
    pseudo_measurements = [pseudo_measurement for i in 1:length(measurements)]

    return TimeVaryingInformationSetSystem(transitions, measurements, pseudo_measurements, information_set, select)
end

function TimeVaryingInformationSetSystem(transitions::Vector{Vector{Transition{T}}},
                                         measurements::Vector{Measurement{T}},
                                         pseudo_measurements::Vector{PseudoMeasurement{T}},
                                         information_set::Vector{UnitRange{Int}},
                                         select::Vector{Int}) where {T <: Real}

    # Initialize empty pseudo-measurement equation
    return TimeVaryingInformationSetSystem{T}(transitions, measurements, pseudo_measurements, information_set, select)
end


function System(system::TimeVaryingInformationSetSystem, select::Int, regime::Int)
    return System(system.transitions[select][regime],
                  system.measurements[regime],
                  system.pseudo_measurements[regime])
end

function Base.getindex(system::TimeVaryingInformationSetSystem, d::Symbol)
    if d in (:transitions, :measurements, :pseudo_measurements, :information_set, :select)
        return getfield(system, d)
    elseif d == :regimes
        return 1:length(system.measurements)
    else
        throw(KeyError(d))
    end
end

function Base.getindex(system::TimeVaryingInformationSetSystem, regime::Int,
                       d::Symbol)
    if d == :measurement
        system.measurements[regime]
    elseif d == :pseudo_measurement
        system.pseudo_measurements[regime]
    elseif d in (:ZZ, :DD, :QQ, :EE)
        system.measurements[regime][d]
    elseif d in (:ZZ_pseudo, :DD_pseudo)
        system.pseudo_measurements[regime][d]
    else
        throw(KeyError(d))
    end
end

# Get a specific regime
function Base.getindex(system::TimeVaryingInformationSetSystem{T},
                       select::Int, d::Int) where {T <: Real}
    if d < 1 || d > length(system.measurements)
        throw(BoundsError(system.measurements, d))
    else
        return System(system, select, d)
    end
end

# Get specific matrix or type from specific regime
function Base.getindex(system::TimeVaryingInformationSetSystem{T},
                       select::Int, regime::Int, d::Symbol) where {T <: Real}
    if d == :transition
        system.transitions[select][regime]
    elseif d in (:TTT, :RRR, :CCC)
        system.transitions[select][regime][d]
    else
        throw(KeyError(d))
    end
end

function n_regimes(system::TimeVaryingInformationSetSystem)
    return length(system.measurements)
end
