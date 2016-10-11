"""
`Transition{T<:AbstractFloat}`

The transition equation of the state-space model takes the form

   `s_{t} = TTT*s_{t-1} + RRR*ε_{t} + CCC`
   
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
function Transition{T<:AbstractFloat}(TTT::Matrix{T}, RRR::Matrix{T}, CCC::Vector{T})
    Transition{T}(TTT, RRR, CCC)
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

   `y_t = ZZ*s_t +  DD + u_t`

where the error `u_t = η_t + MM*ϵ_t` is the measurement error (`η_t`)
plus the error induced by the shocks in the transition equation (`MM*ϵ_t`).

### Fields

If `Nz` is the number of states `s_t`, `Ny` is the number of
observables `y_t`, and `Ne` is the number of shocks `ϵ_t`:

- `ZZ`: the `Ny` x `Nz`  measurement matrix
- `DD`: the `Ny` x 1 constant vector
- `QQ`: the `Ne` x `Ne` covariance matrix for the shocks ϵ_t
- `EE`: the `Ny` x `Ny` variance of measurement error `η_t`
- `MM`: an `Ny` x `Ne` matrix mapping shocks ϵ_t to error in measurement equation
- `VVall`: an `Nz+Ny` x `Nz+Ny` matrix for a time-invariant variance
  matrix for the error in the transition equation and the error in the
  measurement equation (`[RRR*ϵ_t, u_t]'`). Thus:
  - `VVall[1:Nz,1:Nz]` is the covariance of RRR*ϵ_t (`RRR * QQ * RRR'`)
  - `VVall[Nz+1:Nz+Ny, Nz+1:Nz+Ny] = EE + MM*QQ*MM'`
"""
type Measurement{T<:AbstractFloat}
    ZZ::Matrix{T}
    DD::Vector{T}
    QQ::Matrix{T}
    EE::Matrix{T}
    MM::Matrix{T}
    VVall::Matrix{T}
end
function Base.getindex(M::Measurement, d::Symbol)
    if d in (:ZZ, :DD, :QQ, :EE, :MM, :VVall)
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
end
function Base.getindex(system::System, d::Symbol)
    if d in (:transition, :measurement)
        return getfield(system, d)
    elseif d in fieldnames(system.transition)
        return getfield(system.transition, d)
    elseif d in fieldnames(system.measurement)
        return getfield(system.measurement, d)
    else
        throw(KeyError(d))
    end
end
