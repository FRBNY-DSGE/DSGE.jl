function jacobian(m::HetDSGESimpleTaylor)

    # Load in endogenous state and eq cond indices
    endo = m.endogenous_states #=augment_model_states(m.endogenous_states_unnormalized,
                                n_model_states_unnormalized(m)) =#
    eq   = m.equilibrium_conditions

    # Load in parameters, steady-state parameters, and grids
    r::Float64     = m[:r].value
    α::Float64     = m[:α].value
    H::Float64     = m[:H].value
    δ::Float64     = m[:δ].value
    μ_sp::Float64  = m[:μ_sp].value
    ρ_sp::Float64  = m[:ρ_sp].value
    σ_sp::Float64  = m[:σ_sp].value
    γ::Float64     = m[:γ].value
    g::Float64     = m[:g].value
    η::Float64     = m[:η].value
    ρB::Float64    = m[:ρB].value
    ρG::Float64    = m[:ρG].value
    ρZ::Float64    = m[:ρZ].value
    ρμ::Float64    = m[:ρμ].value
    ρlamw::Float64 = m[:ρlamw].value
    ρlamf::Float64 = m[:ρlamf].value
    ρmon::Float64  = m[:ρmon].value
    spp::Float64   = m[:spp].value
    lamw::Float64  = m[:lamw].value
    ϕh::Float64    = m[:ϕh].value
    Φw::Float64    = m[:Φw].value
    lamf::Float64  = m[:lamf].value
    Φp::Float64    = m[:Φp].value
    ρR::Float64    = m[:ρR].value
    ψπ::Float64    = m[:ψπ].value
    ψy::Float64    = m[:ψy].value

    R = 1 + r

    ell::Vector{Float64}  = m[:lstar].value
    c::Vector{Float64}    = m[:cstar].value
    μ::Vector{Float64}    = m[:μstar].value
    β::Float64            = m[:βstar].value

    T::Float64     = m[:Tstar].value
    ω::Float64     = m[:ωstar].value
    xstar::Float64 = m[:xstar].value
    ystar::Float64 = m[:ystar].value
    Rk::Float64    = m[:Rkstar].value
    kstar::Float64 = m[:kstar].value


    xgrid::Vector{Float64} = m.grids[:xgrid].points
    xwts::Vector{Float64}  = m.grids[:xgrid].weights
    sgrid::Vector{Float64} = m.grids[:sgrid].points
    swts::Vector{Float64}  = m.grids[:sgrid].weights
    fgrid::Matrix{Float64} = m.grids[:fgrid]

    xswts = kron(swts,xwts)

    zlo::Float64 = get_setting(m, :zlo)
    zhi::Float64 = get_setting(m, :zhi)

    nx::Int = get_setting(m, :nx)
    ns::Int = get_setting(m, :ns)

    xlo = get_setting(m, :xlo)
#####
    mindens = 1e-8
    rescale_xwts = true
    truncate = true
    if truncate == true
        oldnx = nx
        nx = maximum(find(μ[1:nx]+μ[nx+1:2*nx].>mindens)) # used to be 1e-8
        μ = μ[[1:nx;oldnx+1:oldnx+nx]]
        ell = ell[[1:nx;oldnx+1:oldnx+nx]]
        # experiment with this:
        if rescale_xwts == true
            xgrid   = xgrid[1:nx] #Evenly spaced grid
            xhi = xgrid[nx]
            xscale = xhi-xlo
            xwts     = (xscale/nx)*ones(nx)          #quadrature weights
            xswts = kron(swts,xwts)
        end
    end

    m <= Setting(:nx, nx)
    m <= Setting(:xhi, xhi)
    m <= Setting(:xscale, xscale)
    m.grids[:xgrid] =  Grid(uniform_quadrature(xscale), xlo, xhi, nx, scale = xscale)
    xswts = kron(swts,xwts)

    nxns = nx*ns
    nvars = 4*nxns+35
    nscalars = 35 # number of eqs which output scalars
    nyscalars = 13 # number of scalar jumps
    nxscalars = nscalars - nyscalars # number of scalar states
    m <= Setting(:nxns, nxns)
    m <= Setting(:nvars, nvars)
    m <= Setting(:nscalars, nscalars)
    m <= Setting(:nyscalars, nyscalars)
    m <= Setting(:nxscalars, nxscalars)

    endo[:l′_t1]    = 1:nxns            # lagged ell function
    endo[:μ′_t1]    = nxns+1:2*nxns     # lagged wealth distribution
    #endogenous scalar-valued states
    endo[:k′_t]  = 2*nxns+1:2*nxns+1             # capital –dont get confused with steadystate object K
    endo[:R′_t1] = 2*nxns+2:2*nxns+2              # lagged real interest rate
    endo[:i′_t1] = 2*nxns+3:2*nxns+3             # lagged nominal interest rate
    endo[:π′_t1] = 2*nxns+4:2*nxns+4             # lagged inflation
    endo[:π′_t2] = 2*nxns+5:2*nxns+5             # double-lagged inflation
    endo[:π′_t3] = 2*nxns+6:2*nxns+6             # triple-lagged inflation
    endo[:y′_t1] = 2*nxns+7:2*nxns+7             # lagged gdp
    endo[:y′_t2] = 2*nxns+8:2*nxns+8             # double-lagged gdp
    endo[:y′_t3] = 2*nxns+9:2*nxns+9             # triple-lagged gdp
    endo[:y′_t4] = 2*nxns+10:2*nxns+10           # quad-lagged gdp
    endo[:z′_t1] = 2*nxns+11:2*nxns+11           # lag tfp growth
    endo[:z′_t2] = 2*nxns+12:2*nxns+12           # double-lag tfp growth
    endo[:z′_t3] = 2*nxns+13:2*nxns+13           # triple-lag tfp growth
    endo[:w′_t1] = 2*nxns+14:2*nxns+14           # lag real wages
    endo[:I′_t1] = 2*nxns+15:2*nxns+15           # lag investment–don't get this confused with i, the nominal interest rate
    # exogenous scalar-valued states:
    endo[:BP]    = 2*nxns+16:2*nxns+16        # discount factor shock
    endo[:GP]    = 2*nxns+17:2*nxns+17        # govt spending
    endo[:ZP]    = 2*nxns+18:2*nxns+18        # tfp growth
    endo[:MUP]   = 2*nxns+19:2*nxns+19        # investment shock
    endo[:LAMWP] = 2*nxns+20:2*nxns+20        # wage markup
    endo[:LAMFP] = 2*nxns+21:2*nxns+21        # price markup
    endo[:MONP]  = 2*nxns+22:2*nxns+22        # monetary policy shock

 # function-valued jumps
    endo[:l′_t]  = 2*nxns+23:3*nxns+22 # ell function
    endo[:μ′_t]  = 3*nxns+23:4*nxns+22 # wealth distribution
    #scalar-valued jumps
    endo[:R′_t]  = 4*nxns+23:4*nxns+23        # real interest rate
    endo[:i′_t]  = 4*nxns+24:4*nxns+24        # nominal interest rate
    endo[:t′_t]  = 4*nxns+25:4*nxns+25        # transfers + dividends
    endo[:w′_t]  = 4*nxns+26:4*nxns+26        # real wage
    endo[:L′_t]  = 4*nxns+27:4*nxns+27        # hours worked
    endo[:π′_t]  = 4*nxns+28:4*nxns+28        # inflation
    endo[:wageinflation′_t] = 4*nxns+29:4*nxns+29        # nominal wage inflation
    endo[:mu′_t] =  4*nxns+30:4*nxns+30        # average marginal utility
    endo[:y′_t]  = 4*nxns+31:4*nxns+31       # gdp
    endo[:I′_t]  = 4*nxns+32:4*nxns+32        # investment
    endo[:mc′_t] = 4*nxns+33:4*nxns+33        # marginal cost - this is ζ in HetDSGEₖd.pdf
    endo[:Q′_t]  = 4*nxns+34:4*nxns+34        # Tobin's qfunction
    endo[:capreturn′_t] = 4*nxns+35:4*nxns+35        # return on capital

endo[:l_t1] = endo[:l′_t1] .+ nvars #LELL
endo[:μ_t1] = endo[:μ′_t1] .+ nvars #LM
endo[:k_t] = endo[:k′_t] .+ nvars #KK
endo[:R_t1] = endo[:R′_t1] .+ nvars #LRR
endo[:i_t1] = endo[:i′_t1] .+ nvars #LII
endo[:π_t1] = endo[:π′_t1]  .+ nvars #LPI
endo[:π_t2] = endo[:π′_t2]  .+ nvars #L2PI
endo[:π_t3] = endo[:π′_t3]  .+ nvars #L3PI
endo[:y_t1] = endo[:y′_t1] .+ nvars #LY
endo[:y_t2] = endo[:y′_t2] .+ nvars #L2Y
endo[:y_t3] = endo[:y′_t3] .+ nvars #l3Y
endo[:y_t4] = endo[:y′_t4] .+ nvars #L4Y
endo[:z_t1] = endo[:z′_t1] .+ nvars #LZ
endo[:z_t2] = endo[:z′_t2] .+ nvars #L2Z
endo[:z_t3] = endo[:z′_t3] .+ nvars #L3Z
endo[:w_t1] = endo[:w′_t1] .+ nvars #LW
endo[:I_t1] = endo[:I′_t1] .+ nvars #LX
endo[:B] = endo[:BP] .+ nvars #B
endo[:G] = endo[:GP] .+ nvars #G
endo[:Z] = endo[:ZP] .+ nvars #Z
endo[:MU] = endo[:MUP] .+ nvars #MU
endo[:LAMW] = endo[:LAMWP] .+ nvars #LAMW
endo[:LAMF] = endo[:LAMFP] .+ nvars #LAMF
endo[:MON] = endo[:MONP] .+ nvars #MON
endo[:ELL] = endo[:l′_t] .+ nvars #ELL
endo[:μ_t] = endo[:μ′_t] .+ nvars #M
endo[:R_t] = endo[:R′_t] .+ nvars #RR
endo[:i_t] = endo[:i′_t] .+ nvars #II
endo[:t_t] = endo[:t′_t]  .+ nvars #TT
endo[:w_t] = endo[:w′_t] .+ nvars #W
endo[:L_t] = endo[:L′_t] .+ nvars #HH
endo[:π_t] = endo[:π′_t] .+ nvars #PI
endo[:wageinflation_t] = endo[:wageinflation′_t] .+ nvars #PIW
endo[:mu_t] = endo[:mu′_t] .+ nvars #LAM
endo[:y_t] = endo[:y′_t] .+ nvars #Y
endo[:I_t] = endo[:I′_t]  .+ nvars #X
endo[:mc_t] = endo[:mc′_t]  .+ nvars #MC
endo[:Q_t] = endo[:Q′_t]  .+ nvars #Q
endo[:capreturn_t] = endo[:capreturn′_t] .+ nvars #RK

eqconds = m.equilibrium_conditions
    # function blocks which output a function
    eqconds[:eq_euler]              = 1:nxns
    eqconds[:eq_kolmogorov_fwd]     = nx*ns+1:2*nx*ns
    eqconds[:eq_lag_ell]            = 2*nxns+1:3*nxns
    eqconds[:eq_lag_wealth]         = 3*nxns+1:4*nxns
    # function blocks which map functions to scalars
    eqconds[:eq_market_clearing]    = 4*nxns+1:4*nxns+1
    eqconds[:eq_lambda]             = 4*nxns+2:4*nxns+2
    #scalar blocks involving endogenous variables
    eqconds[:eq_transfers]            = 4*nxns+3:4*nxns+3 # transfers
    eqconds[:eq_investment]          = 4*nxns+4:4*nxns+4 # investment
    eqconds[:eq_tobin_q]              = 4*nxns+5:4*nxns+5 # tobin's q
    eqconds[:eq_capital_accumulation] = 4*nxns+6:4*nxns+6 # capital accumulation
    eqconds[:eq_wage_phillips] = 4*nxns+7:4*nxns+7 # wage phillips curve
    eqconds[:eq_price_phillips] = 4*nxns+8:4*nxns+8 # price phillips curve
    eqconds[:eq_marginal_cost]  = 4*nxns+9:4*nxns+9 # marginal cost
    eqconds[:eq_gdp]  = 4*nxns+10:4*nxns+10 # gdp
    eqconds[:eq_optimal_kl] = 4*nxns+11:4*nxns+11 # optimal K/L ratio
    eqconds[:eq_taylor] = 4*nxns+12:4*nxns+12 # taylor rule
    eqconds[:eq_fisher] = 4*nxns+13:4*nxns+13 # fisher eqn
    eqconds[:eq_nominal_wage_inflation] = 4*nxns+14:4*nxns+14 # nominal wage inflation
    # lagged variables
    eqconds[:LR] = 4*nxns+15:4*nxns+15 # LR
    eqconds[:LI] = 4*nxns+16:4*nxns+16 # LI
    eqconds[:LPI] = 4*nxns+17:4*nxns+17 # LPI
    eqconds[:L2PI] = 4*nxns+18:4*nxns+18 # L2PI
    eqconds[:L3PI] = 4*nxns+19:4*nxns+19 # L3PI
    eqconds[:LY] = 4*nxns+20:4*nxns+20 # LY
    eqconds[:L2Y] = 4*nxns+21:4*nxns+21 # L2Y
    eqconds[:L3Y] = 4*nxns+22:4*nxns+22 # L3Y
    eqconds[:L4Y] = 4*nxns+23:4*nxns+23 # L4Y
    eqconds[:LZ]  = 4*nxns+24:4*nxns+24 # LZ
    eqconds[:L2Z] = 4*nxns+25:4*nxns+25 # L2Z
    eqconds[:L3Z] = 4*nxns+26:4*nxns+26 # L3Z
    eqconds[:LW]  = 4*nxns+27:4*nxns+27 # LW
    eqconds[:LX] = 4*nxns+28:4*nxns+28 # LX
    # shocks
    eqconds[:F33] = 4*nxns+29:4*nxns+29 # discount factor B
    eqconds[:F34] = 4*nxns+30:4*nxns+30 # govt spending G
    eqconds[:F35] = 4*nxns+31:4*nxns+31 # tfp growth Z
    eqconds[:F36] = 4*nxns+32:4*nxns+32 # investment MU
    eqconds[:F37] = 4*nxns+33:4*nxns+33 # wage mkup LAMW
    eqconds[:F38] = 4*nxns+34:4*nxns+34 # price mkup LAMF
    eqconds[:F39] = 4*nxns+35:4*nxns+35 # monetary policy MON

####


    zgrid = range(zlo,stop = zhi,length = nx)
    zwts = (zhi-zlo)/nx
    sumz=0.
    for i=1:nx
        sumz += mollifier_hetdsge(zgrid[i],zhi,zlo)*zwts
    end

    qp(z) = dmollifier_hetdsge(z, zhi, zlo)
    qfunction(x) = mollifier_hetdsge(x, zhi, zlo) #/sumz

    unc = 1 ./ ell .<= repeat(xgrid,ns) .+ η

    # Euler Equation
    dF1_dELL, dF1_dRZ, dF1_dELLP, dF1_dWHP, dF1_dTTP, ee =
        euler_equation_hetdsge(nx, ns, qp, qfunction, xgrid, sgrid, fgrid, unc, xswts,
                               R, γ, β, η, ell, T, ω, H)

    # KF Equation
    dF2_dWH, dF2_dRZ, dF2_dTT,dF2_dELL, bigΨ, dF2_dM =
        kolmogorov_fwd_hetdsge(nx, ns, qfunction, qp, xgrid, sgrid, fgrid, unc, xswts,
                               R, γ, ell, μ, η, T, ω, H, ee)

    # Market clearing, lambda function
    c = min.(1 ./ ell,repeat(xgrid,ns).+η)
    lam = (xswts.*μ)'*(1 ./ c) # average marginal utility which the union uses to set wages
    ϕ = lam*ω/(H^ϕh) # now that we know lam in steady state, choose disutility to target hours H


    nvars = get_setting(m, :nvars)
    # Make the Jacobian
    JJ = zeros(nvars, 2*nvars)

    # Euler equation
    JJ[eq[:eq_euler],endo[:l′_t]] = dF1_dELLP
    JJ[eq[:eq_euler],endo[:ZP]]   = -dF1_dRZ
    JJ[eq[:eq_euler],endo[:w′_t]]   = dF1_dWHP
    JJ[eq[:eq_euler],endo[:L′_t]]  = dF1_dWHP
    JJ[eq[:eq_euler],endo[:t′_t]]  = dF1_dTTP
    JJ[eq[:eq_euler],endo[:B]]    = ell
    JJ[eq[:eq_euler],endo[:ELL]]  = dF1_dELL
    JJ[eq[:eq_euler],endo[:R_t]]   = dF1_dRZ

    # KF eqn
    JJ[eq[:eq_kolmogorov_fwd],endo[:μ_t1]]   = dF2_dM
    JJ[eq[:eq_kolmogorov_fwd],endo[:R_t1]]  = dF2_dRZ
    JJ[eq[:eq_kolmogorov_fwd],endo[:l_t1]] = dF2_dELL
    JJ[eq[:eq_kolmogorov_fwd],endo[:Z]]    = -dF2_dRZ
    JJ[eq[:eq_kolmogorov_fwd],endo[:μ_t]]    = Matrix{Float64}(I, nxns, nxns)
    JJ[eq[:eq_kolmogorov_fwd],endo[:w_t]]    = dF2_dWH
    JJ[eq[:eq_kolmogorov_fwd],endo[:L_t]]   = dF2_dWH
    JJ[eq[:eq_kolmogorov_fwd],endo[:t_t]]   = dF2_dTT

    # update lagged ELL
    JJ[eq[:eq_lag_ell],endo[:l′_t1]] = Matrix{Float64}(I, nxns, nxns)
    JJ[eq[:eq_lag_ell],endo[:ELL]]   = -Matrix{Float64}(I, nxns, nxns)

    # update lagged M
    JJ[eq[:eq_lag_wealth],endo[:μ′_t1]] = Matrix{Float64}(I, nxns, nxns)
    JJ[eq[:eq_lag_wealth],endo[:μ_t]]   = -Matrix{Float64}(I, nxns, nxns)

    # mkt clearing
    JJ[first(eq[:eq_market_clearing]),first(endo[:y_t])]   = ystar/g
    JJ[first(eq[:eq_market_clearing]),first(endo[:G])]   = -ystar/g
    JJ[first(eq[:eq_market_clearing]),first(endo[:I_t])]   = -xstar
    JJ[eq[:eq_market_clearing],endo[:ELL]] = (μ.*unc.*xswts.*c)'
    JJ[eq[:eq_market_clearing],endo[:μ_t]]   = -(xswts.*c)' # note, now we linearize

    # lambda = average marginal utility
    JJ[first(eq[:eq_lambda]),first(endo[:mu_t])] = lam
    JJ[first(eq[:eq_lambda]),endo[:μ_t]]   = -(xswts./c)' # note, now we linearize
    JJ[first(eq[:eq_lambda]),endo[:ELL]] = -(xswts.*unc.*μ./c)'

    # transfer
    JJ[first(eq[:eq_transfers]),first(endo[:t_t])]  = T
    JJ[first(eq[:eq_transfers]),first(endo[:capreturn_t])]  = -Rk*kstar
    JJ[first(eq[:eq_transfers]),first(endo[:k_t])]  = -Rk*kstar
    JJ[first(eq[:eq_transfers]),first(endo[:Z])]  = Rk*kstar
    JJ[first(eq[:eq_transfers]),first(endo[:I_t])]   = xstar
    JJ[first(eq[:eq_transfers]),first(endo[:mc_t])]  = ystar
    JJ[first(eq[:eq_transfers]),first(endo[:y_t])]   = (1-1/g)
    JJ[first(eq[:eq_transfers]),first(endo[:G])]   = (ystar/g)

    # investment
    JJ[first(eq[:eq_investment]),first(endo[:Q_t])]  = 1.
    JJ[first(eq[:eq_investment]),first(endo[:MU])] = 1.
    JJ[first(eq[:eq_investment]),first(endo[:I′_t])] = spp*(exp(3*γ))/R
    JJ[first(eq[:eq_investment]),first(endo[:ZP])] = spp*(exp(3*γ))/R
    JJ[first(eq[:eq_investment]),first(endo[:I_t])]  = -spp*(exp(3*γ))/R - spp*exp(2*γ)
    JJ[first(eq[:eq_investment]),first(endo[:I_t1])] = spp*exp(2*γ)
    JJ[first(eq[:eq_investment]),first(endo[:Z])]  = -spp*exp(2*γ)

    # tobin's q
    JJ[first(eq[:eq_tobin_q]),first(endo[:R_t])]  = R
    JJ[first(eq[:eq_tobin_q]),first(endo[:Q_t])]   = R
    JJ[first(eq[:eq_tobin_q]),first(endo[:capreturn′_t])] = -Rk
    JJ[first(eq[:eq_tobin_q]),first(endo[:Q′_t])]  = -(1-δ)

    # capital accumulation
    JJ[first(eq[:eq_capital_accumulation]),first(endo[:k′_t])] = 1.
    JJ[first(eq[:eq_capital_accumulation]),first(endo[:k_t])]  = -(1-δ)
    JJ[first(eq[:eq_capital_accumulation]),first(endo[:Z])]   = (1-δ)
    JJ[first(eq[:eq_capital_accumulation]),first(endo[:MU])]  = -xstar/kstar
    JJ[first(eq[:eq_capital_accumulation]),first(endo[:I_t])]   = -xstar/kstar

    # wage phillips curve
    JJ[first(eq[:eq_wage_phillips]),first(endo[:wageinflation_t])]  = -1.
    JJ[first(eq[:eq_wage_phillips]),first(endo[:LAMW])] = (ϕ*H^ϕh)/Φw
    JJ[first(eq[:eq_wage_phillips]),first(endo[:L_t])]   = (ϕ*(H^ϕh)*(1+lamw)/lamw*Φw)*ϕh
    JJ[first(eq[:eq_wage_phillips]),first(endo[:mu_t])]  = -(ϕ*(H^ϕh)*(1+lamw)/lamw*Φw)
    JJ[first(eq[:eq_wage_phillips]),first(endo[:w_t])]    = -(ϕ*(H^ϕh)*(1+lamw)/lamw*Φw)
    JJ[first(eq[:eq_wage_phillips]),first(endo[:wageinflation′_t])]  = β

    # price phillips curve
    JJ[first(eq[:eq_price_phillips]),first(endo[:π_t])]   = -1.
    JJ[first(eq[:eq_price_phillips]),first(endo[:mc_t])]  = (1+lamf)/(lamf*Φp)
    JJ[first(eq[:eq_price_phillips]),first(endo[:LAMF])] = 1/Φp
    JJ[first(eq[:eq_price_phillips]),first(endo[:π′_t])]  = 1/R

    # marginal cost
    JJ[first(eq[:eq_marginal_cost]),first(endo[:mc_t])] = 1.
    JJ[first(eq[:eq_marginal_cost]),first(endo[:w_t])]  = -(1-α)
    JJ[first(eq[:eq_marginal_cost]),first(endo[:capreturn_t])] = -α

    # gdp
    JJ[first(eq[:eq_gdp]),first(endo[:y_t])]  = 1.
    JJ[first(eq[:eq_gdp]),first(endo[:Z])]  = α
    JJ[first(eq[:eq_gdp]),first(endo[:k_t])] = -α
    JJ[first(eq[:eq_gdp]),first(endo[:L_t])] = -(1-α)

    # optimal k/l ratio
    JJ[first(eq[:eq_optimal_kl]),first(endo[:capreturn_t])] = 1.
    JJ[first(eq[:eq_optimal_kl]),first(endo[:w_t])]  = -1.
    JJ[first(eq[:eq_optimal_kl]),first(endo[:L_t])] = -1.
    JJ[first(eq[:eq_optimal_kl]),first(endo[:k_t])] = 1.
    JJ[first(eq[:eq_optimal_kl]),first(endo[:Z])]  = -1.

    # taylor rule
    JJ[first(eq[:eq_taylor]),first(endo[:i_t])]   = -1.
    JJ[first(eq[:eq_taylor]),first(endo[:π_t])]   = ψπ
    JJ[first(eq[:eq_taylor]),first(endo[:MONP])] = 1.

    # fisher eqn
    JJ[first(eq[:eq_fisher]),first(endo[:R_t])]  = 1.
    JJ[first(eq[:eq_fisher]),first(endo[:π′_t])] = 1.
    JJ[first(eq[:eq_fisher]),first(endo[:i_t])]  = -1.

    # wage inflation
    JJ[first(eq[:eq_nominal_wage_inflation]),first(endo[:wageinflation_t])] = 1.
    JJ[first(eq[:eq_nominal_wage_inflation]),first(endo[:π_t])]  = -1.
    JJ[first(eq[:eq_nominal_wage_inflation]),first(endo[:Z])]   = -1.
    JJ[first(eq[:eq_nominal_wage_inflation]),first(endo[:w_t])]   = -1.
    JJ[first(eq[:eq_nominal_wage_inflation]),first(endo[:w_t1])]  = 1.

    # update lagged variables
    JJ[first(eq[:LR]),first(endo[:R′_t1])] = 1.
    JJ[first(eq[:LR]),first(endo[:R_t])]   = -1.

    JJ[first(eq[:LI]),first(endo[:i′_t1])] = 1.
    JJ[first(eq[:LI]),first(endo[:i_t])]   = -1.

    JJ[first(eq[:LPI]),first(endo[:π′_t1])] = 1.
    JJ[first(eq[:LPI]),first(endo[:π_t])]   = -1.

    JJ[first(eq[:L2PI]),first(endo[:π′_t2])] = 1.
    JJ[first(eq[:L2PI]),first(endo[:π_t1])]   = -1.

    JJ[first(eq[:L3PI]),first(endo[:π′_t3])] = 1.
    JJ[first(eq[:L3PI]),first(endo[:π_t2])]  = -1.

    JJ[first(eq[:LY]),first(endo[:y′_t1])] = 1.
    JJ[first(eq[:LY]),first(endo[:y_t])]   = -1.

    JJ[first(eq[:L2Y]),first(endo[:y′_t2])] = 1.
    JJ[first(eq[:L2Y]),first(endo[:y_t1])]   = -1.

    JJ[first(eq[:L3Y]),first(endo[:y′_t3])] = 1.
    JJ[first(eq[:L3Y]),first(endo[:y_t2])]  = -1.

    JJ[first(eq[:L4Y]),first(endo[:y′_t4])] = 1.
    JJ[first(eq[:L4Y]),first(endo[:y_t3])]  = -1.

    JJ[first(eq[:LZ]),first(endo[:z′_t1])] = 1.
    JJ[first(eq[:LZ]),first(endo[:Z])]   = -1.

    JJ[first(eq[:L2Z]),first(endo[:z′_t2])] = 1.
    JJ[first(eq[:L2Z]),first(endo[:z_t1])]   = -1.

    JJ[first(eq[:L3Z]),first(endo[:z′_t3])] = 1.
    JJ[first(eq[:L3Z]),first(endo[:z_t2])]  = -1.

    JJ[first(eq[:LW]),first(endo[:w′_t1])] = 1.
    JJ[first(eq[:LW]),first(endo[:w_t])]   = -1.

    JJ[first(eq[:LX]),first(endo[:I′_t1])] = 1.
    JJ[first(eq[:LX]),first(endo[:I_t])]   = -1.

    # discount factor shock
    JJ[first(eq[:F33]),first(endo[:BP])] = 1.
    JJ[first(eq[:F33]),first(endo[:B])]  = -ρB

    # g/y shock
    JJ[first(eq[:F34]),first(endo[:GP])] = 1.
    JJ[first(eq[:F34]),first(endo[:G])]  = -ρG

    # tfp growth shock
    JJ[first(eq[:F35]),first(endo[:ZP])] = 1.
    JJ[first(eq[:F35]),first(endo[:Z])]  = -ρZ

    # investment shock
    JJ[first(eq[:F36]),first(endo[:MUP])] = 1.
    JJ[first(eq[:F36]),first(endo[:MU])]  = -ρμ

    # wage mkup shock
    JJ[first(eq[:F37]),first(endo[:LAMWP])] = 1.
    JJ[first(eq[:F37]),first(endo[:LAMW])]  = -ρlamw

    # price mkup shock
    JJ[first(eq[:F38]),first(endo[:LAMFP])] = 1.
    JJ[first(eq[:F38]),first(endo[:LAMF])]  = -ρlamf

    # monetary policy shock
    JJ[first(eq[:F39]),first(endo[:MONP])] = 1.
    JJ[first(eq[:F39]),first(endo[:MON])]  = -ρmon

    if !m.testing && get_setting(m, :normalize_distr_variables)
        JJ = normalize(m, JJ)
    end

    return JJ
end

function euler_equation_hetdsge_simple_taylor(nx::Int, ns::Int,
                                qp::Function, qfunction::Function,
                                xgrid::Vector{Float64}, sgrid::Vector{Float64},
                                fgrid::Matrix{Float64},
                                unc::BitArray,
                                xswts::Vector{Float64},
                                R::Float64, γ::Float64, β::Float64,
                                η::Float64, ell::Vector{Float64}, T::Float64,
                                ω::Float64, H::Float64)
    nxns = nx*ns
    ee  = zeros(nxns,nxns) # ee[i,j] takes you from i to j
    ξ   = zeros(nxns,nxns)
    Ξ   = zeros(nxns,nxns)
    dF1_dELL = zeros(nxns,nxns)
    dF1_dRZ  = zeros(nxns)
    dF1_dELLP = zeros(nxns,nxns)
    dF1_dWHP = zeros(nxns)
    dF1_dTTP = zeros(nxns)
    ellRHS = zeros(nxns)
    for iss=1:ns
        for ia=1:nx
            i  = nx*(iss-1)+ia
            sumELL = 0.
            sumRZ  = 0.
            sumWH  = 0.
            sumTT  = 0.
            sum_ellRHS = 0.
            for isp=1:ns
                sum_ellRHS = 0.
                for iap=1:nx
                    ip = nx*(isp-1)+iap
                    ee[i,ip] = (xgrid[iap] - R*(exp(-γ))*max(xgrid[ia]-1/ell[i], -η) - T)/(ω*H*sgrid[isp])
                    ξ[i,ip] = ((β*R*xswts[i]*exp(-γ))/(ω*H*sgrid[isp])^2)*max(ell[ip],1/(xgrid[iap]+η))*qp(ee[i,ip])*fgrid[iss,isp]
                    Ξ[i,ip] = ((β*R*xswts[i]*exp(-γ))/(ω*H*sgrid[isp]))*max(ell[ip],1/(xgrid[iap]+η))*qfunction(ee[i,ip])*fgrid[iss,isp]
                    sumELL += ξ[i,ip]*R*(exp(-γ))*unc[i]/ell[i]
                    sumRZ  += ξ[i,ip]*R*(exp(-γ))*max(xgrid[ia] - 1/ell[i],-η)
                    dF1_dELLP[i,ip] = Ξ[i,ip]*unc[ip]
                    sumWH  += Ξ[i,ip] + ξ[i,ip]*ee[i,ip]/(ω*H*sgrid[isp])
                    sumTT  += ξ[i,ip]*T
                    sum_ellRHS += Ξ[i,ip]
                end
            end
            ellRHS[i] = sum_ellRHS
            dF1_dELL[i,i] = -ell[i] - sumELL
            dF1_dRZ[i] = ellRHS[i] - sumRZ
            dF1_dWHP[i] = -sumWH
            dF1_dTTP[i] = -sumTT
        end
    end
    return dF1_dELL, dF1_dRZ, dF1_dELLP, dF1_dWHP, dF1_dTTP, ee
end

function kolmogorov_fwd_hetdsge_simple_taylor(nx::Int, ns::Int,
                                qfunction::Function, qp::Function,
                                xgrid::Vector{Float64}, sgrid::Vector{Float64},
                                fgrid::Matrix{Float64}, unc::BitArray,
                                xswts::Vector{Float64},
                                R::Float64, γ::Float64,
                                ell::Vector{Float64}, μ::Vector{Float64},
                                η::Float64, T::Float64, ω::Float64, H::Float64, ee::Matrix{Float64})
    nxns = nx*ns
    bigΨ    = zeros(nxns,nxns) # this is also dF2_dM
    smallψ = zeros(nxns,nxns)
    dF2_dRZ = zeros(nxns)
    dF2_dELL = zeros(nxns,nxns)
    dF2_dM = zeros(nxns,nxns)
    dF2_dWH = zeros(nxns)
    dF2_dTT = zeros(nxns)

    for isp=1:ns
        for iap=1:nx
            ip  = nx*(isp-1)+iap
            sumWH = 0.
            sumRZ = 0.
            sumTT = 0.
            for iss=1:ns
                for ia=1:nx
                    i = nx*(iss-1)+ia
                    bigΨ[ip,i] = xswts[i]*μ[i]*qfunction(ee[i,ip])*fgrid[iss,isp]/(ω*H*sgrid[isp])
                    smallψ[ip,i] = xswts[i]*μ[i]*qp(ee[i,ip])*fgrid[iss,isp]/((ω*H*sgrid[isp])^2)
                    sumWH += bigΨ[ip,i] + smallψ[ip,i]*ee[i,ip]/(ω*H*sgrid[isp])
                    sumRZ += smallψ[ip,i]*(R*exp(-γ))*max(xgrid[ia] - 1/ell[i],-η)
                    dF2_dELL[ip,i] = smallψ[ip,i]*(R*exp(-γ))*(unc[i]/ell[i])
                    sumTT += smallψ[ip,i]*T
                    dF2_dM[ip,i] = xswts[i]*qfunction(ee[i,ip])*fgrid[iss,isp]/(ω*H*sgrid[isp]) # note, now we linearize
                end
            end
            dF2_dWH[ip] = -sumWH
            dF2_dRZ[ip] = -sumRZ
            dF2_dTT[ip] = -sumTT
        end
    end
    return dF2_dWH, dF2_dRZ, dF2_dTT, dF2_dELL, bigΨ, dF2_dM
end

function normalize(m::HetDSGESimpleTaylor, JJ::Matrix{Float64})

    Qx, _, Qleft, Qright = compose_normalization_matrices(m)

    m <= Setting(:n_predetermined_variables, size(Qx, 1))

	Jac1 = Qleft*sparse(JJ)*Qright

    return Jac1
end

function compose_normalization_matrices(m::HetDSGESimpleTaylor)
    nx::Int = get_setting(m, :nx)
    ns::Int = get_setting(m, :ns)
    nscalars::Int = get_setting(m, :nscalars)
    nyscalars::Int = get_setting(m, :nyscalars)
    nxscalars::Int = get_setting(m, :nxscalars)

    # Create PPP matrix
    P1 = kron(Matrix{Float64}(I, ns,ns),ones(nx,1))
    Ptemp = Matrix{Float64}(I, nx, nx)
    Ptemp = Ptemp[:, 2:end]
    P2 = kron(Matrix{Float64}(I, ns, ns), Ptemp)
    P  = hcat(P1, P2)

    Q,R = qr(P)
    Q = Array(Q)
    S         = Q[:, ns+1:end]'

    nxns = nx*ns

    Qleft     = sparse(cat(Matrix{Float64}(I, nxns, nxns),S,
                    Matrix{Float64}(I, nxns, nxns),S,
                    Matrix{Float64}(I, nscalars, nscalars), dims = [1 2]))
    Qx        = sparse(cat(Matrix{Float64}(I, nxns, nxns),S,
                    Matrix{Float64}(I, nxscalars, nxscalars), dims = [1 2]))
    Qy        = sparse(cat(Matrix{Float64}(I, nxns, nxns),S,
                    Matrix{Float64}(I, nyscalars, nyscalars), dims = [1 2]))
    Qright    = sparse(cat(Qx',Qy',Qx',Qy', dims = [1,2]))

    return Qx, Qy, Qleft, Qright
end
