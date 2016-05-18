type Transition{T<:AbstractFloat}
    TTT::Matrix{T}
    RRR::Matrix{T}
    CCC::Matrix{T}
end
function Transition{T<:AbstractFloat}(TTT::Matrix{T}, RRR::Matrix{T})
    CCC = zeros(eltype(TTT), size(TTT, 1), 1)
    Transition{T}(TTT, RRR, CCC)
end
function Transition{T<:AbstractFloat}(TTT::Matrix{T}, RRR::Matrix{T}, CCC::Matrix{T})
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
    DD::Matrix{T}
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
function Base.getindex(sys::System, d::Symbol)
    if d in (:transition, :measurement)
        return getfield(sys, d)
    elseif d in fieldnames(sys.transition)
        return getfield(sys.transition, d)
    elseif d in fieldnames(sys.measurement)
        return getfield(sys.measurement, d)
    else
        throw(KeyError(d))
    end
end
