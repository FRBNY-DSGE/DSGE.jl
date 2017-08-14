#This function should be identical to kalman_filter here: https://github.com/FRBNY-DSGE/StateSpaceRoutines.jl/blob/master/src/kalman_filter.jl
#However, this one saves a vector of likelihoods for each time step

function kalman_filter_likvec{S<:AbstractFloat}(regime_indices::Vector{Range{Int64}},
    data::Matrix{S}, TTTs::Vector{Matrix{S}}, RRRs::Vector{Matrix{S}}, CCCs::Vector{Vector{S}},
    QQs::Vector{Matrix{S}}, ZZs::Vector{Matrix{S}}, DDs::Vector{Vector{S}}, EEs::Vector{Matrix{S}},log_lik_vec::Vector{Float64},
    z0::Vector{S} = Vector{S}(), P0::Matrix{S} = Matrix{S}();
    allout::Bool = true, n_presample_periods::Int = 0)

    # Dimensions
    T  = size(data,    2) # number of periods of data
    Nz = size(TTTs[1], 1) # number of states
    Ny = size(ZZs[1],  1) # number of observables

    # Initialize outputs
    z = z0
    P = P0
    log_likelihood = zero(S)

    if allout
        pred  = zeros(S, Nz, T)
        vpred = zeros(S, Nz, Nz, T)
        filt  = zeros(S, Nz, T)
        vfilt = zeros(S, Nz, Nz, T)
        yprederror    = zeros(S, Ny, T)
        ystdprederror = zeros(S, Ny, T)
        rmse  = zeros(S, 1, Ny)
        rmsd  = zeros(S, 1, Ny)
    end
 # Iterate through regimes
    for i = 1:length(regime_indices)
        ts = regime_indices[i]
        regime_data = data[:, ts]
        T0 = i == 1 ? n_presample_periods : 0

        if allout
      
            L, z, P, pred[:,ts], vpred[:,:,ts], filt[:,ts], vfilt[:,:,ts],
            yprederror[:,ts], ystdprederror[:,ts], _, _, z0_, P0_        =
                kalman_filter_likvec(regime_data, TTTs[i], RRRs[i], CCCs[i], QQs[i], ZZs[i], DDs[i], EEs[i],log_lik_vec, z, P;
                              allout = true, n_presample_periods = T0)

            # If `n_presample_periods > 0`, then `z0_` and `P0_` are returned as
            # the filtered values at the end of the presample/beginning of the
            # main sample (i.e. not the same the `z0` and `P0` passed into this
            # method, which are from the beginning of the presample). If we are
            # in the first regime, we want to reassign `z0` and `P0`
            # accordingly.
            if i == 1
                z0, P0 = z0_, P0_
            end
        else
              
            L, z, P = kalman_filter_likvec(regime_data, TTTs[i], RRRs[i], CCCs[i], QQs[i], ZZs[i], DDs[i], EEs[i], log_lik_vec, z, P,
                                    allout = false, n_presample_periods = T0)
        end
        log_likelihood += L
      end

    if allout
        rmse = sqrt(mean((yprederror.^2)', 1))
        rmsd = sqrt(mean((ystdprederror.^2)', 1))

        return log_likelihood, z, P, pred, vpred, filt, vfilt, yprederror, ystdprederror, rmse, rmsd, z0, P0, log_lik_vec
    else
        return log_likelihood, z, P
    end
end
function kalman_filter_likvec{S<:AbstractFloat}(data::Matrix{S},
    TTT::Matrix{S}, RRR::Matrix{S}, CCC::Vector{S},
    QQ::Matrix{S}, ZZ::Matrix{S}, DD::Vector{S}, EE::Matrix{S}, log_lik_vec::Vector{Float64},
    z0::Vector{S} = Vector{S}(), P0::Matrix{S} = Matrix{S}();
    allout::Bool = true, n_presample_periods::Int = 0)
  
    # Dimensions
    T  = size(data, 2) # number of periods of data
    Nz = size(TTT,  1) # number of states
    Ne = size(RRR,  2) # number of shocks
    Ny = size(ZZ,   1) # number of observables

    # Populate initial conditions if they are empty
    if isempty(z0) || isempty(P0)
        e, _ = eig(TTT)
        if all(abs(e) .< 1.)
            z0 = (UniformScaling(1) - TTT)\CCC
            P0 = solve_discrete_lyapunov(TTT, RRR*QQ*RRR')
        else
            z0 = CCC
            P0 = 1e6 * eye(Nz)
        end
    end

    z = z0
    P = P0

    # Initialize outputs
    log_likelihood = zero(S)
   
    if allout
        pred          = zeros(S, Nz, T)
        vpred         = zeros(S, Nz, Nz, T)
        filt          = zeros(S, Nz, T)
        vfilt         = zeros(S, Nz, Nz, T)
        yprederror    = NaN*zeros(S, Ny, T)
        ystdprederror = NaN*zeros(S, Ny, T)
    end

    V = RRR*QQ*RRR' # V = Var(z_t) = Var(Rϵ_t)
for t = 1:T
        # Index out rows of the measurement equation for which we have
        # nonmissing data in period t
        nonmissing = !isnan(data[:, t])
        y_t  = data[nonmissing, t]
        ZZ_t = ZZ[nonmissing, :]
        DD_t = DD[nonmissing]
        EE_t = EE[nonmissing, nonmissing]
        Ny_t = length(y_t)

        ## Forecast
        z = TTT*z + CCC                 # z_{t|t-1} = TTT*z_{t-1|t-1} + CCC
        P = TTT*P*TTT' + RRR*QQ*RRR'    # P_{t|t-1} = Var s_{t|t-1} = TTT*P_{t-1|t-1}*TTT' + RRR*QQ*RRR'
        V = ZZ_t*P*ZZ_t' + EE_t         # V_{t|t-1} = Var y_{t|t-1} = ZZ*P_{t|t-1}*ZZ' + EE
        V = (V+V')/2

        dy = y_t - ZZ_t*z - DD_t        # dy  = y_t - y_{t|t-1} = prediction error
        ddy = V\dy                      # ddy = (1/V_{t|t-1})dy = weighted prediction error

        if allout
            pred[:, t]                   = z
            vpred[:, :, t]               = P
            yprederror[nonmissing, t]    = dy
            ystdprederror[nonmissing, t] = dy ./ sqrt(diag(V))
        end

        ## Compute marginal log-likelihood, log P(y_t|y_1,...y_{t-1},θ)
        ## log P(y_1,...,y_T|θ) ∝ log P(y_1|θ) + log P(y_2|y_1,θ) + ... + P(y_T|y_1,...,y_{T-1},θ)
        if t > n_presample_periods
            log_likelihood += -log(det(V))/2 - first(dy'*ddy/2) - Ny_t*log(2*pi)/2
            push!(log_lik_vec, -log(det(V))/2 - first(dy'*ddy/2) - Ny_t*log(2*pi)/2) 
        end

        ## Update
        z = z + P'*ZZ_t'*ddy            # z_{t|t} = z_{t|t-1} + P_{t|t-1}'*ZZ'*(1/V_{t|t-1})dy
        P = P - P'*ZZ_t'/V*ZZ_t*P       # P_{t|t} = P_{t|t-1} - P_{t|t-1}'*ZZ'*(1/V_{t|t-1})*ZZ*P_{t|t-1}

        if allout
            filt[:, t]     = z
            vfilt[:, :, t] = P
        end

    end # of loop through periods
if allout && n_presample_periods > 0
        mainsample_periods = n_presample_periods+1:T

        # If we choose to discard presample periods, then we reassign `z0`
        # and `P0` to be their values at the end of the presample/beginning
        # of the main sample
        z0 = filt[:,     n_presample_periods]
        P0 = vfilt[:, :, n_presample_periods]

        pred          = pred[:,     mainsample_periods]
        vpred         = vpred[:, :, mainsample_periods]
        filt          = filt[:,     mainsample_periods]
        vfilt         = vfilt[:, :, mainsample_periods]
        yprederror    = yprederror[:,  mainsample_periods]
        ystdprederror = ystdprederror[:, mainsample_periods]
    end

    if allout
        rmse = sqrt(mean((yprederror.^2)', 1))
        rmsd = sqrt(mean((ystdprederror.^2)', 1))

        return log_likelihood, z, P, pred, vpred, filt, vfilt, yprederror, ystdprederror, rmse, rmsd, z0, P0, log_lik_vec
    else
        return log_likelihood, z, P
    end
end