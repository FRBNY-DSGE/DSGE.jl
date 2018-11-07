"""
```
transform_transition_matrices(m::AbstractModel, T::Matrix{T},
                              R::Matrix{T}, C::Matrix{T}) where {T<:AbstractFloat}
```

### Arguments

-`m`: the m object
-`T`, `R`, and `C`: matrices of the state transition equation

### Return values

- `TTT`, `RRR`, and `CCC`: transform transition equation matrices to accomodate states simulated
    between observations

### Description
"""
function transform_transition_matrices(m::AbstractModel, TT::Matrix{T},
                                       R::Matrix{T}, C::Vector{T};
                                       track_lag::Bool = true) where {T<:AbstractFloat}

    freq = get_setting(m, :state_simulation_freq)
    if track_lag
        TTT = zeros(size(TT, 1) * (freq+1), size(TT, 2) * (freq+1)) # pre-allocate
        RRR = zeros(size(R, 1) * (freq+1), size(R, 2) * freq) # pre-allocate
    else
        TTT = zeros(size(TT, 1) * freq, size(TT, 2) * freq) # pre-allocate
        RRR = zeros(size(R, 1) * freq, size(R, 2) * freq) # pre-allocate
    end

    # Construct powers of TT iteratively
    TT_powers = Dict{Int64, Matrix{T}}()
    for i = 1:freq
        if i == 1
            TT_powers[i] = TT
        else
            TT_powers[i] = TT_powers[i - 1] * TT
        end
    end

    # Insert directly into TTT
    if track_lag
        for i = 1:freq+1
            if i == 1
                TTT[1:size(TT,1), 1 + size(TT,2) * freq:end] = eye(TT) # track last lag
            else
                TTT[1 + size(TT,1)*(i-1):size(TT,1)*i, 1 + size(TT,2)*freq:end] = TT_powers[i-1]
            end
        end
    else
        for i = 1:freq
            TTT[1 + size(TT,1)*(i-1):size(TT,1)*i, 1 + size(TT,2)*(freq-1):end] = TT_powers[i]
        end
    end

    # Create RRR
    TR_powers = Dict{Int64, Matrix{T}}()
    TR_powers[0] = R

    # Construct T^i * R iteratively
    for i = 1:(freq - 1)
        TR_powers[i] = TT * TR_powers[i - 1]
    end

    # Fill in RRR and CCC
    if track_lag
        RRR[1+size(R,1):size(R, 1)*2, 1:size(R,2)] = TR_powers[0]
        for i in 3:freq+1
            for j = (i-2):-1:0
                RRR[1+size(R,1)*(i-1):size(R,1)*i, 1+size(R,2)*j:size(R,2)*(j+1)] = TR_powers[i-j-2]
            end
        end
        CCC = vcat(zeros(C), repeat(C, freq))
    else
        RRR[1:size(R, 1), 1:size(R,2)] = TR_powers[0]
        for i in 2:freq
            for j = (i-1):-1:0
                # inserts values from right nonzero blocks to left
                RRR[1+size(R,1)*(i-1):size(R,1)*i, 1+size(R,2)*j:size(R,2)*(j+1)] = TR_powers[i-j-1]
            end
        end
        CCC = repeat(C, freq) #vec(repeat(C, 1, freq))
    end

    return TTT, RRR, CCC
end
