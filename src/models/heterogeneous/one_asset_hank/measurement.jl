"""
```
measurement(m::KrusellSmith{T}, TTT::Matrix{T},
            RRR::Matrix{T}, CCC::Vector{T}) where {T<:AbstractFloat}
```

Assign measurement equation

```
y_t = ZZ*s_t + DD + u_t
```

where

```
Var(ϵ_t) = QQ
Var(u_t) = EE
Cov(ϵ_t, u_t) = 0
```
"""
function measurement(m::OneAssetHANK{T},
                     TTT::Matrix{T},
                     RRR::Matrix{T},
                     CCC::Vector{T},
                     inverse_basis::Matrix{T}) where {T<:AbstractFloat}
    endo = m.endogenous_states

    track_lag = get_setting(m, :track_lag)
    freq = get_setting(m, :state_simulation_freq)

    QQ = eye(freq) * 1/freq # freq = no. of partitions of quarter
    EE = zeros(1,1)
    DD = [0.]
    if track_lag
        ZZ = zeros(1, n_states(m))
        ZZ[endo[:output][1]] = 1
        augmented_inverse_basis = zeros(n_states(m), (freq+1)*size(TTT, 1))
        augmented_inverse_basis[:,end - size(inverse_basis,2)+1:end] = inverse_basis
        ZZ = ZZ*augmented_inverse_basis
    else
        ZZ = zeros(1, n_states(m))
        ZZ[endo[:output][1]] = 1
        augmented_inverse_basis = zeros(n_states(m), freq*size(TTT, 1))
        augmented_inverse_basis[:,end - size(inverse_basis,2)+1:end] = inverse_basis
        ZZ = ZZ*augmented_inverse_basis
    end
    return Measurement(ZZ, DD, QQ, EE)

end

function measurement(m::OneAssetHANK{T},
                     TTT::Matrix{T},
                     RRR::Matrix{T},
                     CCC::Vector{T},
                     inverse_basis::SparseMatrixCSC{T,Int64}) where {T<:AbstractFloat}
    endo = m.endogenous_states

    track_lag = get_setting(m, :track_lag)
    freq = get_setting(m, :state_simulation_freq)

    QQ = eye(freq) * 1/freq # freq = no. of partitions of quarter
    EE = zeros(1,1)
    DD = [0.]
    if track_lag
        ZZ = zeros(1, n_states(m))
        ZZ[endo[:output][1]] = 1
        augmented_inverse_basis = zeros(n_states(m), (freq+1)*size(TTT, 1))
        augmented_inverse_basis[:, end - size(inverse_basis,2) + 1 : end] = inverse_basis
        ZZ = ZZ * augmented_inverse_basis
    else
        ZZ = zeros(1, n_states(m))
        ZZ[endo[:output][1]] = 1
        augmented_inverse_basis = zeros(n_states(m), freq*size(TTT, 1))
        augmented_inverse_basis[:,end - size(inverse_basis,2)+1:end] = inverse_basis
        ZZ = ZZ*augmented_inverse_basis
    end
    return Measurement(ZZ, DD, QQ, EE)

end

function measurement(m::OneAssetHANK{T},
                     TT::Matrix{T},
                     RR::Matrix{T},
                     CC::Vector{T},
                     inverse_basis::Matrix{T},
                     dt::T,
                     track_lag::Bool = true;
                     abs_err::Float64 = 1e-8,
                     max_iter::Int64 = size(TT, 1)*2) where {T<:AbstractFloat}

    endo = m.endogenous_states
    freq = get_setting(m, :state_simulation_freq)
    QQ   = eye(freq, freq) * dt # shocks are Brownian motion -> variance is dt
    EE   = zeros(1,1) * dt
    DD   = [0.]

    if track_lag
          ZZ = zeros(1, bas_sz2 * (freq + 1))

        # Numerically approximate integral of expm(T*dt) * s_0 using the formula
        # dt * (I + (TT * dt)/2! + (TT * dt)^2/3! + ...)
        TT_prev_power = TT * dt^2 # save previous power so can recursively multiply
        polynomial = eye(TT) * dt + TT_prev_power/2 # initialize polynomial
        for i = 2:max_iter+1
            TT_prev_power = TT_prev_power * TT * dt # compute next power

            # Check if error is small enough
            if maximum(abs.(TT_prev_power./factorial(i+1))) < abs_err
                break
            end

            # Otherwise add onto polynomial
            polynomial += TT_prev_power/factorial(i+1)
        end

        # Change basis from reduced basis to full basis on left side
        integ_output = inverse_basis * polynomial

        # Grab row corresponding to output, fill in ZZ
        integ_output = integ_output[endo[:output], :]
        for i = 1:freq
            ZZ[1, 1 + size(inverse_basis, 2) * (i-1) : size(inverse_basis, 2) * i] = integ_output
        end

        # Change of basis into reduced basis
        # augment_inverse_basis = zeros(size(inverse_basis,1), size(inverse_basis,2) * (freq+1))
        # for i = 1:(freq+1)
        #     augment_inverse_basis[:, 1+size(inverse_basis,2)*(i-1):size(inverse_basis,2)*i] = inverse_basis
        # end
    else
        ZZ = zeros(1, size(inverse_basis,2) * freq)

        # Numerically approximate integral of expm(T*(-dt)) * s_t using the formula
        # -dt * (I + (TT * -dt)/2! + (TT * -dt)^2/3! + ...)
        # Note the negative dt is to conduct a backward motion
        TT_prev_power = TT * dt ^ 2 # save previous power so can recursively multiply
        polynomial    = eye(TT) * (-dt) + TT_prev_power / 2

        for i = 2:max_iter+1
            TT_prev_power = TT_prev_power * TT * (-dt) # compute next power

            # Check if error is small enough
            if maximum(abs.(TT_prev_power ./ factorial(i + 1))) < abs_err
                break
            end

            # Otherwise add onto polynomial
            polynomial += TT_prev_power / factorial(i + 1)
        end

        # Change basis from reduced basis to full basis on left side
        integ_output = inverse_basis * polynomial

        # Grab row corresponding to output, fill in ZZ
        integ_output = integ_output[endo[:output], :]
        for i = 1:freq
            ZZ[1, 1+size(inverse_basis,1)*(i-1):size(inverse_basis,1)*i] = integ_output
        end

        # Change of basis into reduced basis
        # augment_inverse_basis = zeros(size(inverse_basis,1), size(inverse_basis,2) * freq)
        # for i = 1:freq
        #     augment_inverse_basis[:, 1+size(inverse_basis,2)*(i-1):size(inverse_basis,2)*i] = inverse_basis
        # end
    end
    return Measurement(ZZ, DD, QQ, EE)
end
