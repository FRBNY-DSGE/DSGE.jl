"""
`Transition{T<:AbstractFloat}`

The transition equation of the state-space model takes the form

   `s_{t} = TTT*s_{t-1} + RRR*ϵ_t + CCC`

The `Transition` type stores the coefficient `Matrix{T}`s (`TTT`, `RRR`) and constant `Vector{T} CCC`.
"""
type Transition{T<:AbstractFloat}
    TTT::Matrix{T}
    RRR::Matrix{T}
    CCC::Vector{T}
end
function Transition{T<:AbstractFloat}(TTT::Matrix{T}, RRR::Matrix{T})
    CCC = zeros(eltype(TTT), size(TTT, 1))
    Transition{T}(TTT, RRR, CCC)
end
function Transition{T<:AbstractFloat}(TTT::Matrix{T}, RRR::Matrix{T}, CCC::Matrix{T})
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

where the error `η_t` is the measurement error, which is uncorrelated with the
shocks in the transition equation `ϵ_t`.

### Fields

If `Nz` is the number of states `s_t`, `Ny` is the number of
observables `y_t`, and `Ne` is the number of shocks `ϵ_t`:

- `ZZ`: the `Ny` x `Nz`  measurement matrix
- `DD`: the `Ny` x 1 constant vector
- `QQ`: the `Ne` x `Ne` covariance matrix for the shocks `ϵ_t`
- `EE`: the `Ny` x `Ny` covariance matrix for the measurement error `η_t`
"""
type Measurement{T<:AbstractFloat}
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
`System{T<:AbstractFloat}`

A type containing the transition and measurement equations for a
state-space model. The matrices may be directly indexed: `sys[:TTT]`
returns `sys.transition.TTT`, etc.
"""
type System{T<:AbstractFloat}
    transition::Transition{T}
    measurement::Measurement{T}
    pseudo_measurement::Nullable{PseudoObservableMapping{T}}
end

function System{T<:AbstractFloat}(transition::Transition{T}, measurement::Measurement{T})
    return System(transition, measurement, Nullable{PseudoObservableMapping{T}}())
end

function Base.getindex(system::System, d::Symbol)
    if d in (:transition, :measurement, :pseudo_measurement)
        return getfield(system, d)
    elseif d in fieldnames(system.transition)
        return getfield(system.transition, d)
    elseif d in fieldnames(system.measurement)
        return getfield(system.measurement, d)
    elseif !isnull(system.pseudo_measurement) && d in [:ZZ_pseudo, :DD_pseudo]
        return getfield(get(system.pseudo_measurement), d)
    elseif isnull(system.pseudo_measurement) && d in [:ZZ_pseudo, :DD_pseudo]
        throw(PseudoMeasurementUndefError())
    else
        throw(KeyError(d))
    end
end

function Base.copy(system::System)
    trans = Transition(system[:TTT], system[:RRR], system[:CCC])
    meas  = Measurement(system[:ZZ], system[:DD], system[:QQ], system[:EE])
    pseudo_measurement = if isnull(system.pseudo_measurement)
        Nullable{PseudoObservableMapping}()
    else
        pseudo_inds = get(system.pseudo_measurement).inds
        pseudo_measurement = Nullable(PseudoObservableMapping(pseudo_inds,
                                 system[:ZZ_pseudo], system[:DD_pseudo]))
    end
    return System(trans, meas, pseudo_measurement)
end

type PseudoMeasurementUndefError <: Exception
end

Base.showerror(io::IO, e::PseudoMeasurementUndefError) = print(io, "A pseudo-measurement
    equation was not defined for the model from which this state-space system was computed.")