using DSGE, Plots

# This script loads in matrices from the DSGE in Del Negro, Schorfheide, Smets, and Wouters (2007)
# and computes the DSGE impulse responses to exogenous structural shocks and the
# corresponding DSGE-VECM impulse responses using the rotation matrix from the DSGE
# to identify the mapping from the innovations in the VECM to the structural shocks.
# The script creates a dictionary called `plots_dict`, whose keys are the names of
# different shocks and values are the plots of the impulse responses of different
# observables to the shock.
@time begin
    fp = dirname(@__FILE__)

    matdata     = load(joinpath(fp, "../../test/reference/dsgevecm_lambda_irfs.jld2"))
    horizon     = 20
    nshocks     = size(matdata["RRR"], 2)
    obs_dict    = Dict{Symbol, Int}(:output_growth => 1, :consumption_growth => 2,
                                 :investment_growth => 3, :hours => 4,
                                 :wage_growth => 5, :inflation => 6, :nominalrate => 7)
    # The shocks are technology, labor preference, investment price, preference,
    # government spending, price markup, and monetary policy
    shocks_dict = Dict{Symbol, Int}(:z => 1, :φ => 2, :μ => 3, :b => 4, :g => 5,
                                    :λ_f => 6, :r => 7)

    # DSGEVECM IRFs in deviations
    ŷ_dsgevecm = DSGE.impulse_responses(matdata["TTT"], matdata["RRR"], matdata["ZZ"],
                                        vec(matdata["DD"]), matdata["MM"],
                                        matdata["QQ"], Int(matdata["k"]), Int(matdata["nvar"]),
                                        Int(matdata["coint"]), matdata["cct_sim"], matdata["sig_sim"],
                                        matdata["cointvec"], horizon, deviations = true)

    # DSGE IRFs in deviations
    dsge_states = zeros(size(matdata["TTT"], 1), horizon, nshocks)
    for i in 1:nshocks
        shock = zeros(nshocks)
        shock[i] = -sqrt(matdata["QQ"][i, i])
        dsge_states[:, 1, i] = matdata["RRR"] * shock
        for t in 2:horizon
            dsge_states[:, t, i] = matdata["TTT"] * dsge_states[:, t - 1, i]
        end
    end
    ŷ_dsge = matdata["ZZ"] * dsge_states[:, :, i]

    # Plot!
    plots_dict = Dict{Symbol, Dict{Symbol, Plots.Plot}}()
    for (shock, j) in  shocks_dict
        for (obs, i) in obs_dict
            plots_dict[shock][obs] = plot(1:horizon, ŷ_dsge[i, :, j], label = "DSGE",
                                           color = :black, linestyle = :dash, left_margin = 20px)
            plot!(1:horizon, ŷ_dsgevecm[i, :, j], label = "DSGE-VECM", color = :black,
                  left_margin = 20px)
        end
    end
end
