type Transition{T<:AbstractFloat}
    TTT::Matrix{T}
    RRR::Matrix{T}
    CCC::Matrix{T}
end
function Base.getindex(eq::Transition, d::Symbol)
    if d in (:TTT, :RRR, :CCC)
        return getfield(eq, d)
    else
        throw(KeyError(d))
    end
end

type System{T<:AbstractFloat}
    transition::Transition{T}
    measurement::Measurement{T}
end
function Base.getindex(sys::System, d::Symbol)
    if d in (:transition, :measurement)
        return getfield(sys, d)
    else
        throw(KeyError(d))
    end
end
