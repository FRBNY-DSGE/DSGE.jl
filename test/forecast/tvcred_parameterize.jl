function tvcred_parameterize!(m::AbstractDSGEModel)
    # Parameterize
    for i in 1:5
        adj = (i == 2 || i == 3) ? 0.25 : 1.
        ModelConstructors.set_regime_val!(m[:σ_g], i, adj * m[:σ_g].value)
        ModelConstructors.set_regime_val!(m[:σ_b], i, adj * m[:σ_b].value)
        ModelConstructors.set_regime_val!(m[:σ_μ], i, adj * m[:σ_μ].value)
        ModelConstructors.set_regime_val!(m[:σ_ztil], i, adj * m[:σ_ztil].value)
        ModelConstructors.set_regime_val!(m[:σ_λ_f], i, adj * m[:σ_λ_f].value)
        ModelConstructors.set_regime_val!(m[:σ_λ_w], i, adj * m[:σ_λ_w].value)
        ModelConstructors.set_regime_val!(m[:σ_r_m], i, adj * m[:σ_r_m].value)
        ModelConstructors.set_regime_val!(m[:σ_σ_ω], i, adj * m[:σ_σ_ω].value)
        ModelConstructors.set_regime_val!(m[:σ_μ_e], i, adj * m[:σ_μ_e].value)
        ModelConstructors.set_regime_val!(m[:σ_γ], i, adj * m[:σ_γ].value)
        ModelConstructors.set_regime_val!(m[:σ_π_star], i, adj * m[:σ_π_star].value)
        ModelConstructors.set_regime_val!(m[:σ_lr], i, adj * m[:σ_lr].value)
        ModelConstructors.set_regime_val!(m[:σ_z_p], i, adj * m[:σ_z_p].value)
        ModelConstructors.set_regime_val!(m[:σ_tfp], i, adj * m[:σ_tfp].value)
        err_adj = (i in 2:4) ? 10. : adj
        ModelConstructors.set_regime_val!(m[:σ_gdpdef], i, err_adj * m[:σ_gdpdef].value)
        ModelConstructors.set_regime_val!(m[:σ_corepce], i, err_adj * m[:σ_corepce].value)
        ModelConstructors.set_regime_val!(m[:σ_gdp], i, adj * m[:σ_gdp].value)
        ModelConstructors.set_regime_val!(m[:σ_gdi], i, adj * m[:σ_gdi].value)

        for j = 1:DSGE.n_mon_anticipated_shocks(m)
            ModelConstructors.set_regime_val!(m[Symbol("σ_r_m$(j)")], i, adj * m[Symbol("σ_r_m$(j)")])
        end
    end

    ModelConstructors.set_regime_val!(m[:σ_z_p], 2, 0.)
    ModelConstructors.set_regime_val!(m[:σ_z_p], 3, 0.)
    ModelConstructors.set_regime_val!(m[:σ_r_m], 2, 10. * m[:σ_r_m].value)
    ModelConstructors.set_regime_val!(m[:σ_r_m], 3, 10. * m[:σ_r_m].value)
    for i = 2:3
        for j = 1:DSGE.n_mon_anticipated_shocks(m)
            ModelConstructors.set_regime_val!(m[Symbol("σ_r_m$(j)")], i, m[Symbol("σ_r_m$(j)")].value)
        end
    end

    # Set non-standard shocks
    ModelConstructors.set_regime_val!(m[:σ_φ], 1, 0.)
    ModelConstructors.set_regime_val!(m[:σ_φ], 2, 400.)
    ModelConstructors.set_regime_val!(m[:σ_φ], 3, 400.)
    ModelConstructors.set_regime_val!(m[:σ_φ], 4, 400.)
    ModelConstructors.set_regime_val!(m[:σ_φ], 5, 0.)
    ModelConstructors.set_regime_val!(m[:σ_φ_prop], 1, 0.)
    ModelConstructors.set_regime_val!(m[:σ_φ_prop], 2, 2.)
    ModelConstructors.set_regime_val!(m[:σ_φ_prop], 3, 0.)
    ModelConstructors.set_regime_val!(m[:σ_φ_prop], 4, 0.)
    ModelConstructors.set_regime_val!(m[:σ_φ_prop], 5, 0.)
    ModelConstructors.set_regime_val!(m[:σ_biidc], 1, 0.)
    ModelConstructors.set_regime_val!(m[:σ_biidc], 2, 4.)
    ModelConstructors.set_regime_val!(m[:σ_biidc], 3, 4.)
    ModelConstructors.set_regime_val!(m[:σ_biidc], 4, 4.)
    ModelConstructors.set_regime_val!(m[:σ_biidc], 5, 0.)
    ModelConstructors.set_regime_val!(m[:σ_biidc_prop], 1, 0.)
    ModelConstructors.set_regime_val!(m[:σ_biidc_prop], 2, 2.)
    ModelConstructors.set_regime_val!(m[:σ_biidc_prop], 3, 0.)
    ModelConstructors.set_regime_val!(m[:σ_biidc_prop], 4, 0.)
    ModelConstructors.set_regime_val!(m[:σ_biidc_prop], 5, 0.)
    ModelConstructors.set_regime_val!(m[:σ_ygap], 1, 0.)
    ModelConstructors.set_regime_val!(m[:σ_ygap], 2, 0.)
    ModelConstructors.set_regime_val!(m[:σ_ygap], 3, 20.)
    ModelConstructors.set_regime_val!(m[:σ_ygap], 4, 0.)
    ModelConstructors.set_regime_val!(m[:σ_ygap], 5, 0.)
    ModelConstructors.set_regime_val!(m[:σ_pgap], 1, 0.)
    ModelConstructors.set_regime_val!(m[:σ_pgap], 2, 0.)
    ModelConstructors.set_regime_val!(m[:σ_pgap], 3, 20.)
    ModelConstructors.set_regime_val!(m[:σ_pgap], 4, 0.)
    ModelConstructors.set_regime_val!(m[:σ_pgap], 5, 0.)
    ModelConstructors.set_regime_val!(m[:σ_ziid], 1, 0.)
    ModelConstructors.set_regime_val!(m[:σ_ziid], 2, 5.)
    ModelConstructors.set_regime_val!(m[:σ_ziid], 3, 5.)
    ModelConstructors.set_regime_val!(m[:σ_ziid], 4, 5.)
    ModelConstructors.set_regime_val!(m[:σ_ziid], 5, 0.)
    ModelConstructors.set_regime_val!(m[:σ_ziid_prop], 1, 0.)
    ModelConstructors.set_regime_val!(m[:σ_ziid_prop], 2, 2.)
    ModelConstructors.set_regime_val!(m[:σ_ziid_prop], 3, 0.)
    ModelConstructors.set_regime_val!(m[:σ_ziid_prop], 4, 0.)
    ModelConstructors.set_regime_val!(m[:σ_ziid_prop], 5, 0.)
    m[:σ_biidc1].prior = m[:σ_biidc].prior
    m[:σ_biidc1].valuebounds = m[:σ_biidc].valuebounds
    ModelConstructors.set_regime_val!(m[:σ_biidc1], 1, 0.)
    ModelConstructors.set_regime_val!(m[:σ_biidc1], 2, 0.)
    ModelConstructors.set_regime_val!(m[:σ_biidc1], 3, 0.)
    ModelConstructors.set_regime_val!(m[:σ_biidc1], 4, ModelConstructors.regime_val(m[:σ_biidc], 2))
    ModelConstructors.set_regime_val!(m[:σ_biidc1], 5, 0.)
end

function set_regime_vals_fnct(m, n)
    if n > 6
        for p in m.parameters
            regime_switching = !isempty(p.regimes)
            if regime_switching
                if haskey(p.regimes, :value)
                    for i in 6:n
                        ModelConstructors.set_regime_val!(p, i, ModelConstructors.regime_val(p, 1))
                    end
                end
            end
        end
    end
end
