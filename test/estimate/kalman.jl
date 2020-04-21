# Test regime indices and matrix computations for integration with kalman_likelihood
using DSGE, ModelConstructors, Dates, Test

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
push!(zlb_plus_regime_dates_ans, ([1:(izlb - 1), izlb:regime_dates_ans[1][1][end], regime_dates_ans[1][2], regime_dates_ans[1][3]], false))
push!(zlb_plus_regime_dates_ans, ([regime_dates_ans[2][1], regime_dates_ans[2][2][1]:(izlb - 1),
                                izlb:regime_dates_ans[2][2][end], regime_dates_ans[2][3]], false))
push!(zlb_plus_regime_dates_ans, ([regime_dates_ans[3][1], regime_dates_ans[3][2], regime_dates_ans[2][3][1]:(izlb - 1),
                                izlb:T], false))
push!(zlb_plus_regime_dates_ans, ([regime_dates_ans[4][1], izlb:regime_dates_ans[4][2][end], regime_dates_ans[4][3]], true))
push!(zlb_plus_regime_dates_ans, ([regime_dates_ans[5][1], regime_dates_ans[5][2], izlb:T], true))
push!(zlb_plus_regime_dates_ans, (regime_dates_ans[6], true))
push!(zlb_plus_regime_dates_ans, ([1:(length(regime_dates_ans[7][1])),
                                   (length(regime_dates_ans[7][1]) + 1):(length(regime_dates_ans[7][2]) + length(regime_dates_ans[7][1])),
                                   (length(regime_dates_ans[7][2]) + length(regime_dates_ans[7][1]) + 1):regime_dates_ans[7][end][end]], true))

out_regime_dates = []
out_zlb_plus_regime_dates = []
data = df_to_matrix(m, df)
for i in 1:length(regime_dates_dicts)
    m <= Setting(:regime_dates, regime_dates_dicts[i])
    push!(out_regime_dates, DSGE.regime_indices(m))
    if i != 7
        push!(out_zlb_plus_regime_dates, DSGE.zlb_plus_regime_indices(m, data))
    else
        push!(out_zlb_plus_regime_dates, DSGE.zlb_plus_regime_indices(m, data, regime_dates_dicts[7][1]))
    end
end

@testset "Verify regime indices" begin
    for (a, b) in zip(regime_dates_ans, out_regime_dates)
        @test a == b
    end
    for (a, b) in zip(zlb_plus_regime_dates_ans, out_zlb_plus_regime_dates)
        @test a[1] == b[1]
        @test a[2] == b[2]
    end
end

# Test on ZLB regime matrices using
