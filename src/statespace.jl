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
