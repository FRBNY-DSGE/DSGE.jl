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
function Transition(TTT::Matrix{T}, RRR::Matrix{T}) where T<:AbstractFloat
    CCC = zeros(eltype(TTT), size(TTT, 1))
    Transition{T}(TTT, RRR, CCC)
end
function Transition(TTT::Matrix{T}, RRR::Matrix{T}, CCC::Matrix{T}) where T<:AbstractFloat
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

function measurement(m::AbstractDSGEModel, trans::Transition; shocks::Bool=true)
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

A mutable struct containing the transition and measurement equations for a
state-space model. The matrices may be directly indexed: `sys[:TTT]`
returns `sys.transition.TTT`, etc.
"""
mutable struct System{T<:AbstractFloat}
    transition::Transition{T}
    measurement::Measurement{T}
    pseudo_measurement::PseudoMeasurement{T}
end

function System(transition::Transition{T}, measurement::Measurement{T}) where T<:AbstractFloat
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
`RegimeSwitchingSystem{N<:Int,T<:AbstractFloat}`

A mutable struct containing the transition and measurement equations for a
state-space model with regime-switching.
The matrices may be directly indexed: `sys[1][:TTT]`
returns `sys.regime[1].transition.TTT`, etc.
"""
mutable struct RegimeSwitchingSystem{T<:AbstractFloat}
    transitions::Vector{Transition{T}}
    measurements::Vector{Measurement{T}}
    pseudo_measurements::Vector{PseudoMeasurement{T}}
end

function RegimeSwitchingSystem(transitions::Vector{Transition{T}},
                               measurements::Vector{Measurement{T}}) where {T<:AbstractFloat}
    if length(transitions) != length(measurements)
        error("The number of Transition (n = $(length(transitions)))" *
              " and Measurement (n = $(length(measurements))) objects must match.")
    end

    # Initialize empty pseudo-measurement equation
    _n_states = size(transitions[1].TTT, 1)
    _n_pseudo = 0
    ZZ_pseudo = zeros(_n_pseudo, _n_states)
    DD_pseudo = zeros(_n_pseudo)
    pseudo_measurement = PseudoMeasurement(ZZ_pseudo, DD_pseudo)
    pseudo_measurements = [pseudo_measurement for i in 1:length(transitions)]

    return RegimeSwitchingSystem(transitions, measurements, pseudo_measurements)
end

function RegimeSwitchingSystem(transitions::Vector{Transition{T}},
                               measurement::Measurement{T},
                               pseudo_measurements::Vector{PseudoMeasurement{T}}) where {T<:AbstractFloat}

   #= if length(transitions) != length(measurements)
        error("The number of Transition (n = $(length(transitions)))" *
              " and Measurement (n = $(length(measurements))) objects must match.")
    end=#
    return RegimeSwitchingSystem(transitions, [measurement; measurement], pseudo_measurements)
end


# Retrieve System correspond to a regime
function System(system::RegimeSwitchingSystem, regime::Int)
    return System(system.transitions[regime],
                  system.measurements[regime],
                  system.pseudo_measurements[regime])
end

# Get dictionary of different regime matrices
function Base.getindex(system::RegimeSwitchingSystem{T},
                       d::Symbol) where {T<:AbstractFloat}
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
                       d::Int) where {T<:AbstractFloat}
    if d < 1 || d > length(system.transitions)
        throw(BoundsError(system.transitions,d))
    else
        return System(system, d)
    end
end

# Get specific matrix or type from specific regime
function Base.getindex(system::RegimeSwitchingSystem{T},
                       d::Tuple{N,Symbol}) where {N<:Int,T<:AbstractFloat}
    if d[2] == :transition
        system.transitions[d[1]]
    elseif d[2] == :measurement
        system.measurements[d[1]]
    elseif d[2] == :pseudo_measurement
        system.pseudo_measurements[d[1]]
    elseif d in (:TTT, :RRR, :CCC)
        system.transitions[d[1]][d[2]]
    elseif d in (:ZZ, :DD, :QQ, :EE)
        system.measurements[d[1]][d[2]]
    elseif d in (:ZZ_pseudo, :DD_pseudo)
        system.pseudo_measurements[d[1]][d[2]]
    else
        throw(KeyError(d))
    end
end

function n_regimes(system::RegimeSwitchingSystem{T}) where {T<:AbstractFloat}
    return length(system.transitions)
end

function Base.copy(system::RegimeSwitchingSystem{T}) where {T<:AbstractFloat}
    transitions         = OrderedDict{N,T}()
    measurements        = OrderedDict{N,T}()
    pseudo_measurements = OrderedDict{N,T}()

    for i in 1:length(system[:transitions])
        trans = Transition(system[(i,:TTT)], system[(i,:RRR)], system[(i,:CCC)])
        meas  = Measurement(system[(i,:ZZ)], system[(i,:DD)],
                            system[(i,:QQ)], system[(i,:EE)])
        pseudo_meas = PseudoMeasurement(system[(i,:ZZ_pseudo)], system[(i,:DD_pseudo)])

        transitions[i]         = trans
        measurements[i]        = meas
        pseudo_measurements[i] = pseudo_meas
    end
    return RegimeSwitchingSystem(transitions, measurements, pseudo_measurements)
end

"""
```
compute_system(m; apply_altpolicy = false)
```

Given the current model parameters, compute the state-space system
corresponding to model `m`. Returns a `System` object.
"""
function compute_system(m::AbstractDSGEModel{T}; apply_altpolicy::Bool = false,
                        regime_switching::Bool = false, num_regimes::Int = 2,
                        verbose::Symbol = :high) where T<:AbstractFloat

    solution_method = get_setting(m, :solution_method)

    # Solve model
    if regime_switching
        if solution_method == :gensys
            TTTs = Vector{Matrix{T}}(undef,num_regimes)
            RRRs = Vector{Matrix{T}}(undef,num_regimes)
            CCCs = Vector{Vector{T}}(undef,num_regimes)
            transition_equations = Vector{Transition{T}}(undef,num_regimes)

            for i = 1:num_regimes
                TTTs[i], RRRs[i], CCCs[i] = solve(m; apply_altpolicy = apply_altpolicy,
                                                  regime_switching = true, regime = i, verbose = verbose)
                transition_equations[i] = Transition(TTTs[i], RRRs[i], CCCs[i])
            end
            # Measurement eqn doesn't depend on the time-varying parameters (so just pick 1st)
            measurement_equation = measurement(m, TTTs[1], RRRs[1], CCCs[1])

            pseudo_measurement_equations = pseudo_measurement(m, TTTs, RRRs, CCCs)
            return RegimeSwitchingSystem(transition_equations,
                                         measurement_equation,
                                         pseudo_measurement_equations)


            # Infer which measurement and pseudo-measurement equations to use
        #=    type_tuple = (typeof(m), Vector{Matrix{T}}, Vector{Matrix{T}}, Vector{Vector{T}})
            if hasmethod(measurement, type_tuple)
                measurement_equations = measurement(m, TTTs, RRRs, CCCs)
            else
                measurement_equation = measurement(m, TTTs[1], RRRs[1], CCCs[1])
            end

            if hasmethod(pseudo_measurement, type_tuple)
                pseudo_measurement_equations = pseudo_measurement(m, TTTs, RRRs, CCCs)
                return RegimeSwitchingSystem(transition_equations,
                                             measurement_equations,
                                             pseudo_measurement_equations)
            else
                return RegimeSwitchingSystem(transition_equations, measurement_equations)
            end=#
        else
            throw("solution_method $solution_method has not been implemented.")
        end
    else
        if solution_method == :gensys
            TTT, RRR, CCC = solve(m; apply_altpolicy = apply_altpolicy, verbose = verbose)
            transition_equation = Transition(TTT, RRR, CCC)

            # Solve measurement equation
            measurement_equation = measurement(m, TTT, RRR, CCC)
        elseif solution_method == :klein
            # Unpacking the method from solve to hang on to TTT_jump
            if m.spec == "het_dsge"
                TTT_jump, TTT_state, eu = klein(m)
            else
                TTT_jump, TTT_state, eu = klein(m)
            end
            if eu==-1
                throw(KleinError())
            end

            TTT, RRR = klein_transition_matrices(m, TTT_state, TTT_jump)
            CCC = zeros(n_model_states(m))

            if m.spec == "real_bond_mkup"
                GDPeqn = construct_GDPeqn(m, TTT_jump)
                TTT, RRR, CCC = augment_states(m, TTT, TTT_jump, RRR, CCC, GDPeqn)
                # Measurement (needs the additional TTT_jump argument)
                measurement_equation = measurement(m, TTT, TTT_jump, RRR, CCC, GDPeqn)
            elseif m.spec == "het_dsge"
                TTT, RRR, CCC = augment_states(m, TTT, TTT_jump, RRR, CCC)
                measurement_equation = measurement(m, TTT, RRR, CCC)
            else
                TTT, RRR, CCC        = augment_states(m, TTT, RRR, CCC)
                measurement_equation = measurement(m, TTT, RRR, CCC)
            end

            transition_equation = Transition(TTT, RRR, CCC)
        else
            throw("solution_method $solution_method has not been implemented.")
        end

        type_tuple = (typeof(m), Matrix{T}, Matrix{T}, Vector{T})
        if hasmethod(pseudo_measurement, type_tuple)
            # Solve pseudo-measurement equation
            pseudo_measurement_equation = pseudo_measurement(m, TTT, RRR, CCC)
            return System(transition_equation, measurement_equation, pseudo_measurement_equation)
        else
            return System(transition_equation, measurement_equation)
        end
    end
end

"""
```
compute_system(m::PoolModel{T})
```

Given the current model parameters, compute the state-space system
corresponding to the PoolModel model `m`.

Outputs

```
Φ: state transition function
Ψ: likelihood function, given weights on underlying models (the states) and predictive densities
F_ϵ: structural shock distribution
F_u: likelihood function measurement error distribution
F_λ: initial distribution of λ for state transition function
"""
function compute_system(m::PoolModel{T};
                        verbose::Symbol = :high) where T<:AbstractFloat
    Φ, F_ϵ, F_λ = transition(m)
    Ψ, F_u = measurement(m)
    return Φ, Ψ, F_ϵ, F_u, F_λ
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
