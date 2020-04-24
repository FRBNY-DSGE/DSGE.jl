using DSGE, ModelConstructors, Dates, Test

@testset "Calculating indicies and matrices for pre- and post-ZLB regimes" begin
    m  = Model1002("ss10")
    m <= Setting(:regime_switching, false)
    m <= Setting(:date_presample_start, Date(1959, 9, 30))
    T = subtract_quarters(Date(2015, 9, 30), Date(1959, 9, 30)) + 1
    izlb = subtract_quarters(date_zlb_start(m), Date(1959, 9, 30)) + 1
    data = zeros(n_observables(m), T) # pseudo-data matrix, just needed for dimensions
    m <= Setting(:n_mon_anticipated_shocks, 6)
    system = compute_system(m)
    @test_throws ErrorException DSGE.zlb_regime_indices(m, data, Date(1900, 1, 1))
    @test_throws ErrorException DSGE.zlb_regime_matrices(m, system, Date(1900, 1, 1))
    @test DSGE.zlb_regime_indices(m, data) == [1:(izlb - 1), izlb:T]
    out = DSGE.zlb_regime_matrices(m, system)
    for (i, x) in enumerate([:TTT, :RRR, :CCC, :QQ, :ZZ, :DD, :EE])
        @test length(out[i]) == 2
        if x == :QQ
            @test @test_matrix_approx_eq out[i][2] system[x]
            QQtmp = copy(system[x])
            inds = DSGE.inds_shocks_no_ant(m)
            QQtmp[inds[end] + 1:end, inds[end] + 1:end] .= 0.
            @test @test_matrix_approx_eq out[i][1] QQtmp
        else
            @test @test_matrix_approx_eq out[i][1] system[x]
            @test @test_matrix_approx_eq out[i][2] system[x]
        end
    end
    m <= Setting(:date_zlb_start, Date(1959, 3, 31))
    @test DSGE.zlb_regime_indices(m, data) == [1:T]
    out = DSGE.zlb_regime_matrices(m, system)
    for (i, x) in enumerate([:TTT, :RRR, :CCC, :QQ, :ZZ, :DD, :EE])
        @test length(out[i]) == 1
        @test @test_matrix_approx_eq out[i][1] system[x]
    end
    @test DSGE.zlb_regime_indices(m, Matrix{Float64}(undef, 0, 0)) == [1:0]
    m <= Setting(:n_mon_anticipated_shocks, 0)
    @test DSGE.zlb_regime_indices(m, data) == [1:T]
    for (i, x) in enumerate([:TTT, :RRR, :CCC, :QQ, :ZZ, :DD, :EE])
        @test length(out[i]) == 1
        @test @test_matrix_approx_eq out[i][1] system[x]
    end
end

@testset "Calculate the regime indices and matrices with regime switching (incl. pre- and post-ZLB regimes)" begin
    m  = AnSchorfheide(testing = true)
    m <= Setting(:regime_switching, true)
    m <= Setting(:n_regimes, 3)
    df = load_data(m)
    m <= Setting(:date_presample_start, df[1, :date])
    T = length(df[!, :date])
    izlb = findfirst(df[!, :date] .== date_zlb_start(m))

    regime_dates_dicts = [Dict{Int, Date}(1 => date_presample_start(m),
                                          2 => DSGE.quartertodate("2010-Q1"),
                                          3 => DSGE.quartertodate("2012-Q4")),
                          Dict{Int, Date}(1 => date_presample_start(m),
                                          2 => DSGE.quartertodate("1980-Q2"),
                                          3 => DSGE.quartertodate("2012-Q4")),
                          Dict{Int, Date}(1 => date_presample_start(m),
                                          2 => DSGE.quartertodate("1980-Q2"),
                                          3 => DSGE.quartertodate("2003-Q4")),
                          Dict{Int, Date}(1 => date_presample_start(m),
                                          2 => DSGE.quartertodate("2008-Q4"),
                                          3 => DSGE.quartertodate("2012-Q4")),
                          Dict{Int, Date}(1 => date_presample_start(m),
                                          2 => DSGE.quartertodate("2000-Q2"),
                                          3 => DSGE.quartertodate("2008-Q4")),
                          Dict{Int, Date}(1 => date_zlb_start(m),
                                          2 => DSGE.quartertodate("2010-Q4"),
                                          3 => DSGE.quartertodate("2013-Q1")),
                          Dict{Int, Date}(1 => DSGE.next_quarter(date_zlb_start(m)),
                                          2 => DSGE.quartertodate("2010-Q4"),
                                          3 => DSGE.quartertodate("2013-Q1"))]
    regime_dates_ans = []
    for i in 1:length(regime_dates_dicts)
        dic = regime_dates_dicts[i]
        ind1 = findfirst(df[!, :date] .== dic[1])
        ind2 = findfirst(df[!, :date] .== dic[2])
        ind3 = findfirst(df[!, :date] .== dic[3])

        push!(regime_dates_ans, [ind1:(ind2 - 1), ind2:(ind3 - 1), ind3:T])
    end
    zlb_plus_regime_dates_ans = []
    push!(zlb_plus_regime_dates_ans, ([1:(izlb - 1), izlb:regime_dates_ans[1][1][end], regime_dates_ans[1][2], regime_dates_ans[1][3]], 2, true))
    push!(zlb_plus_regime_dates_ans, ([regime_dates_ans[2][1], regime_dates_ans[2][2][1]:(izlb - 1),
                                       izlb:regime_dates_ans[2][2][end], regime_dates_ans[2][3]], 3, true))
    push!(zlb_plus_regime_dates_ans, ([regime_dates_ans[3][1], regime_dates_ans[3][2], regime_dates_ans[3][3][1]:(izlb - 1),
                                       izlb:T], 4, true))
    push!(zlb_plus_regime_dates_ans, ([regime_dates_ans[4][1], izlb:regime_dates_ans[4][2][end], regime_dates_ans[4][3]], 2, false))
    push!(zlb_plus_regime_dates_ans, ([regime_dates_ans[5][1], regime_dates_ans[5][2], izlb:T], 3, false))
    push!(zlb_plus_regime_dates_ans, ([1:(length(regime_dates_ans[6][1])),
                                       (length(regime_dates_ans[6][1]) + 1):(length(regime_dates_ans[6][2]) + length(regime_dates_ans[6][1])),
                                       (length(regime_dates_ans[6][2]) +
                                        length(regime_dates_ans[6][1]) + 1):regime_dates_ans[6][end][end]], 1, false))
    push!(zlb_plus_regime_dates_ans, ([1:(length(regime_dates_ans[7][1])),
                                       (length(regime_dates_ans[7][1]) + 1):(length(regime_dates_ans[7][2]) + length(regime_dates_ans[7][1])),
                                       (length(regime_dates_ans[7][2]) +
                                        length(regime_dates_ans[7][1]) + 1):regime_dates_ans[7][end][end]], 1, false))

    out_regime_dates = []
    out_zlb_plus_regime_dates = []
    data = df_to_matrix(m, df)
    for i in 1:length(regime_dates_dicts)
        m <= Setting(:regime_dates, regime_dates_dicts[i])
        push!(out_regime_dates, DSGE.regime_indices(m))
        if i == 7
            m <= Setting(:date_zlb_start, DSGE.next_quarter(date_zlb_start(m)))
            m <= Setting(:date_presample_start, date_zlb_start(m))
            push!(out_zlb_plus_regime_dates, DSGE.zlb_plus_regime_indices(m, data))
            m <= Setting(:date_zlb_start, DSGE.prev_quarter(date_zlb_start(m)))
            m <= Setting(:date_presample_start, df[1, :date])
        elseif i == 6
            m <= Setting(:date_presample_start, date_zlb_start(m))
            push!(out_zlb_plus_regime_dates, DSGE.zlb_plus_regime_indices(m, data))
            m <= Setting(:date_presample_start, df[1, :date])
        else
            push!(out_zlb_plus_regime_dates, DSGE.zlb_plus_regime_indices(m, data))
        end
        @test_throws ErrorException DSGE.zlb_plus_regime_indices(m, data, Date(1900, 1, 1))
    end

    for (a, b) in zip(regime_dates_ans, out_regime_dates)
        @test a == b
    end
    for (a, b) in zip(zlb_plus_regime_dates_ans, out_zlb_plus_regime_dates)
        @test a == b
    end

    m = Model1002("ss10")
    m <= Setting(:regime_switching, true)
    m <= Setting(:n_regimes, 3)
    m <= Setting(:date_presample_start, df[1, :date]) # this df is from AnSchorfheide; we assume we use a data set of the same number of periods
    T = length(df[!, :date])
    izlb = findfirst(df[!, :date] .== date_zlb_start(m))

    # Prepare the expected answers
    zlb_plus_regime_mats_ans = []
    push!(zlb_plus_regime_mats_ans, [true, false, false, false]) # Booleans for whether we have a pre-ZLB matrix or not
    push!(zlb_plus_regime_mats_ans, [true, true, false, false])
    push!(zlb_plus_regime_mats_ans, [true, true, true, false])
    push!(zlb_plus_regime_mats_ans, [true, false, false])
    push!(zlb_plus_regime_mats_ans, [true, true, false])
    push!(zlb_plus_regime_mats_ans, [false, false, false])
    push!(zlb_plus_regime_mats_ans, [false, false, false])


    system = compute_system(m)
    # system = RegimeSwitchingSystem([system, system])
    prezlb = deepcopy(system[1])
    inds = DSGE.inds_shocks_no_ant(m)
    prezlb[:QQ][inds[end] + 1:end, inds[end] + 1:end] .= 0.
    postzlb = deepcopy(system[2])

    out_zlb_plus_regime_mats = []
    for i = 1:length(regime_dates_dicts)
        if i == 7
            m <= Setting(:date_zlb_start, DSGE.next_quarter(date_zlb_start(m)))
            m <= Setting(:date_presample_start, date_zlb_start(m))
            push!(out_zlb_plus_regime_mats, DSGE.zlb_plus_regime_matrices(m, system, length(out_zlb_plus_regime_dates[i][1]),
                                                                          ind_zlb_start = out_zlb_plus_regime_dates[i][2],
                                                                          splice_zlb_regime = out_zlb_plus_regime_dates[i][3]))
            m <= Setting(:date_zlb_start, DSGE.prev_quarter(date_zlb_start(m)))
            m <= Setting(:date_presample_start, df[1, :date])
        elseif i == 6
            m <= Setting(:date_presample_start, date_zlb_start(m))
            push!(out_zlb_plus_regime_mats, DSGE.zlb_plus_regime_matrices(m, system, length(out_zlb_plus_regime_dates[i][1]),
                                                                          ind_zlb_start = out_zlb_plus_regime_dates[i][2],
                                                                          splice_zlb_regime = out_zlb_plus_regime_dates[i][3]))
            m <= Setting(:date_presample_start, df[1, :date])
        else
            push!(out_zlb_plus_regime_mats, DSGE.zlb_plus_regime_matrices(m, system, length(out_zlb_plus_regime_dates[i][1]),
                                                                          ind_zlb_start = out_zlb_plus_regime_dates[i][2],
                                                                          splice_zlb_regime = out_zlb_plus_regime_dates[i][3]))
        end

        @test_throws ErrorException DSGE.zlb_plus_regime_matrices(m, system, length(out_zlb_plus_regime_dates[i][1]),
                                                                  Date(1900, 1, 1);
                                                                  ind_zlb_start = out_zlb_plus_regime_dates[i][2],
                                                                  splice_zlb_regime = out_zlb_plus_regime_dates[i][3])
        for j in 1:length(out_zlb_plus_regime_mats[i])
            @test length(out_zlb_plus_regime_mats[i][j]) == length(out_zlb_plus_regime_dates[i][1])
        end
        for (j, is_prezlb) in enumerate(zlb_plus_regime_mats_ans[i])
            if is_prezlb
                for (k, x) in enumerate([:TTT, :RRR, :CCC, :QQ, :ZZ, :DD, :EE])
                    @test @test_matrix_approx_eq out_zlb_plus_regime_mats[i][k][j] prezlb[x]
                end
            else
                for (k, x) in enumerate([:TTT, :RRR, :CCC, :QQ, :ZZ, :DD, :EE])
                    @test @test_matrix_approx_eq out_zlb_plus_regime_mats[i][k][j] postzlb[x]
                end
            end
        end
    end
    # Now consider a case where the regimes have different values
    system[2, :DD][1] = 1e3
    system[3, :DD][2] = 1e3
    prezlb = deepcopy(system[1])
    inds = DSGE.inds_shocks_no_ant(m)
    prezlb[:QQ][inds[end] + 1:end, inds[end] + 1:end] .= 0.
    postzlb1 = deepcopy(system[1])
    postzlb2 = deepcopy(system[2])
    postzlb3 = deepcopy(system[3])
    out = DSGE.zlb_plus_regime_matrices(m, system, length(out_zlb_plus_regime_dates[1][1]),
                                        ind_zlb_start = out_zlb_plus_regime_dates[1][2],
                                        splice_zlb_regime = out_zlb_plus_regime_dates[1][3])
    for (j, x) in enumerate([:TTT, :RRR, :CCC, :QQ, :ZZ, :DD, :EE])
        @test @test_matrix_approx_eq out[j][1] prezlb[x]
        @test @test_matrix_approx_eq out[j][2] postzlb1[x]
        @test @test_matrix_approx_eq out[j][3] postzlb2[x]
        @test @test_matrix_approx_eq out[j][4] postzlb3[x]
    end
end
