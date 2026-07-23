

function compute_system!(m::kCSC;
                              tvis::Bool = false,
                              set_regime_eqcond_info!::Union{Function, Nothing} = nothing,
                              verbose::Symbol = :high)
    if !isnothing(set_regime_eqcond_info!)
        set_regime_eqcond_info!(m)
    end
    # t0 = time()

    #update_ss_v5
    DSGE.update_ss_csc!(m)
    # t1 = time()
    # println("updatess: $(t1 - t0) seconds")
    
    #check certain values >0

    #state_reduc_tv_copula 
    DSGE.state_reduc_csc(m)
    # t2 = time()
    # println("statereduc: $(t2 - t1) seconds")


    #maybe square some values?

    #compute jacob TESTING ONLY, read csv jacob instead of actually compute jacob for testing
    F21_ad, F22_ad, F23_ad, F24_ad, F41_ad, F42_ad, F43_ad, F44_ad  = DSGE.jacobian!(m)
    # _jcsv = joinpath(@__DIR__, "tests", "csv")
    # F21_ad = Matrix(CSV.read(joinpath(_jcsv, "F21_ad.csv"), DataFrame))
    # F22_ad = Matrix(CSV.read(joinpath(_jcsv, "F22_ad.csv"), DataFrame))
    # F23_ad = Matrix(CSV.read(joinpath(_jcsv, "F23_ad.csv"), DataFrame))
    # F24_ad = Matrix(CSV.read(joinpath(_jcsv, "F24_ad.csv"), DataFrame))
    # F41_ad = Matrix(CSV.read(joinpath(_jcsv, "F41_ad.csv"), DataFrame))
    # F42_ad = Matrix(CSV.read(joinpath(_jcsv, "F42_ad.csv"), DataFrame))
    # F43_ad = Matrix(CSV.read(joinpath(_jcsv, "F43_ad.csv"), DataFrame))
    # F44_ad = Matrix(CSV.read(joinpath(_jcsv, "F44_ad.csv"), DataFrame))


    # t3 = time()
    # println("jacob: $(t3 - t2) seconds")
    # SGU solver
    hx, gx, F1_aux, F2_aux, F3_aux, F4_aux, indicator_1 = DSGE.SGU_solver(m, F21_ad, F22_ad, F23_ad, F24_ad, F41_ad, F42_ad, F43_ad, F44_ad)

    # t4 = time()
    # println("sgu: $(t4 - t3) seconds")

    if indicator_1 == 0
        return nothing
    end

    df = DSGE.load_data_bbq(m)

    # t5 = time()
    # println("load data: $(t5 - t4) seconds")

    H_aux, P_ref, Q_ref, HQ, SIGMA_full, SIGMA = DSGE.pq(m, hx, gx, F1_aux, F2_aux, F3_aux, F4_aux)

    # t6 = time()
    # println("pq1: $(t6 - t5) seconds")

    regime_inds, y, Ts, Rs, Cs, Qs, Zs, Ds, Es = DSGE.prepare_PQ_regimes(m, df, H_aux, P_ref, Q_ref, HQ, SIGMA_full, SIGMA)

    # t7 = time()
    # println("pq1: $(t7 - t6) seconds")

    transition = Transition(Ts[1], Rs[1], Cs[1])
    measurement = Measurement(Zs[1], Ds[1], Qs[1], Es[1])

    return System(transition, measurement)
end

# function compute_system(m::kCSC;
#                         tvis::Bool = false,
#                         set_regime_eqcond_info!::Union{Function, Nothing} = nothing,
#                         verbose::Symbol = :high)
#     return compute_system_nozlb!(m; tvis = tvis,
#                                 set_regime_eqcond_info! = set_regime_eqcond_info!,
#                                 verbose = verbose)
# end
