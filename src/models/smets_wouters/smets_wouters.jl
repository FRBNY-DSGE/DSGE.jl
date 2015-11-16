"""
TODO: Decide whether this is the right place for this documentation...

SmetsWouters{T} <: AbstractDSGEModel{T}

The SmetsWouters type defines the structure of the FRBNY DSGE
model. We can then concisely pass around a Model object to the
remaining steps of the model (solve, estimate, and forecast).

### Fields

#### Parameters and Steady-States
* `parameters::Vector{AbstractParameter}`: Vector of all time-invariant model parameters.

* `steady_state::Vector`: Model steady-state values, computed as a
function of elements of `parameters`.

* `keys::Dict{Symbol,Int}`: Maps human-readable names for all model
parameters and steady-states to their indices in `parameters` and
`steady_state`.

#### Inputs to Measurement and Equilibrium Condition Equations

The following fields are dictionaries that map human-readible names to
row and column indices in the matrix representations of of the
measurement equation and equilibrium conditions.

* `endogenous_states::Dict{Symbol,Int}`: Maps each state to a column
in the measurement and equilibrium condition matrices.

* `exogenous_shocks::Dict{Symbol,Int}`: Maps each shock to a column in
the measurement and equilibrium condition matrices.

* `expected_shocks::Dict{Symbol,Int}`: Maps each expected shock to a
column in the measurement and equilibrium condition matrices.

* `equilibrium_conditions::Dict{Symbol,Int}`: Maps each equlibrium
condition to a row in the model's equilibrium condition matrices.

* `endogenous_states_postgensys::Dict{Symbol,Int}`: Maps lagged states
to their columns in the measurement and equilibrium condition
equations. These are added after Gensys solves the model.

* `observables::Dict{Symbol,Int}`: Maps each observable to a row in
the model's measurement equation matrices.

#### Model Specifications and Settings

* `spec::AbstractString`: The model specification identifier, "m990",
cached here for filepath computation.

* `subspec::AbstractString`: The model subspecification number,
indicating that some parameters from the original model spec ("ss0")
are initialized differently. Cached here for filepath computation. 



* `settings::Dict{Symbol,Setting}`: Settings/flags that affect
computation without changing the economic or mathematical setup of
the model.

* `test_settings::Dict{Symbol,Setting}`: Settings/flags for testing mode

#### Other Fields

* `rng::MersenneTwister`: Random number generator. By default, it is
seeded to ensure replicability in algorithms that involve randomness
(such as Metropolis-Hastings).

* `testing::Bool`: Indicates whether the model is in testing mode. If
`true`, settings from `m.test_settings` are used in place of those in
`m.settings`.

* `_filestrings::SortedDict{Symbol,AbstractString,ForwardOrdering}`:
An alphabetized list of setting identifier strings. These are
concatenated and appended to the filenames of all output files to
avoid overwriting the output of previous estimations/forecasts that
differ only in their settings, but not in their underlying
mathematical structure. 
"""
type SmetsWouters{T} <: AbstractDSGEModel{T}
    parameters::ParameterVector{T}                  # vector of all time-invariant model parameters
    steady_state::ParameterVector{T}                # model steady-state values
    keys::Dict{Symbol,Int}                          # human-readable names for all the model
                                                    # parameters and steady-num_states

    endogenous_states::Dict{Symbol,Int}             # these fields used to create matrices in the
    exogenous_shocks::Dict{Symbol,Int}              # measurement and equilibrium condition equations.
    expected_shocks::Dict{Symbol,Int}               #
    equilibrium_conditions::Dict{Symbol,Int}        #
    endogenous_states_postgensys::Dict{Symbol,Int}  #
    observables::Dict{Symbol,Int}                   #

    spec::ASCIIString                               # Model specification number 
    subspec::ASCIIString                            # Model subspecification 
    settings::Dict{Symbol,Setting}                  # Settings/flags for computation
    test_settings::Dict{Symbol,Setting}             # Settings/flags for testing mode
    rng::MersenneTwister                            # Random number generator
    testing::Bool                                   # Whether we are in testing mode or not
    _filestrings::SortedDict{Symbol,AbstractString, ForwardOrdering} # The strings we will print to a filename
end

description(m::SmetsWouters) = "Smets-Wouters Model"

#=
doc"""
Inputs: `m:: SmetsWouters`

Description:
Initializes indices for all of `m`'s states, shocks, and equilibrium conditions.
"""
=#
function initialize_model_indices!(m::SmetsWouters)
    # Endogenous states
    endogenous_states = [[
        :y_t, :c_t, :i_t, :qk_t, :k_t, :kbar_t, :u_t, :rk_t, :mc_t,
        :π_t, :μ_ω_t, :w_t, :L_t, :R_t, :g_t, :b_t, :μ_t, :z_t,
        :λ_f_t, :λ_f_t1, :λ_w_t, :λ_w_t1, :rm_t, :Ec_t, :Eqk_t, :Ei_t,
        :Eπ_t, :EL_t, :Erk_t, :Ew_t, :y_f_t, :c_f_t,
        :i_f_t, :qk_f_t, :k_f_t, :kbar_f_t, :u_f_t, :rk_f_t, :w_f_t,
        :L_f_t, :r_f_t, :Ec_f_t, :Eqk_f_t, :Ei_f_t, :EL_f_t, :Erk_f_t, :ztil_t];
        [symbol("rm_tl$i") for i = 1:num_anticipated_shocks(m)]]

    # Exogenous shocks
    exogenous_shocks = [[
        :g_sh, :b_sh, :μ_sh, :z_sh, :λ_f_sh, :λ_w_sh, :rm_sh];
        [symbol("rm_shl$i") for i = 1:num_anticipated_shocks(m)]]

    # Expectations shocks
    expected_shocks = [
        :Ec_sh, :Eqk_sh, :Ei_sh, :Eπ_sh, :EL_sh, :Erk_sh, :Ew_sh, :Ec_f_sh,
        :Eqk_f_sh, :Ei_f_sh, :EL_f_sh, :Erk_f_sh]

    # Equilibrium conditions
    equilibrium_conditions = [[
        :euler, :inv, :capval, :output, :caputl, :capsrv, :capev,
        :mkupp, :phlps, :caprnt, :msub, :wage, :mp, :res, :eq_g, :eq_b, :eq_μ, :eq_z,
        :eq_λ_f, :eq_λ_w, :eq_rm, :eq_λ_f1, :eq_λ_w1, :eq_Ec,
        :eq_Eqk, :eq_Ei, :eq_Eπ, :eq_EL, :eq_Erk, :eq_Ew, :euler_f, :inv_f,
        :capval_f, :output_f, :caputl_f, :capsrv_f, :capev_f, :mkupp_f, :caprnt_f, :msub_f,
        :res_f, :eq_Ec_f, :eq_Eqk_f, :eq_Ei_f, :eq_EL_f, :eq_Erk_f, :eq_ztil];
        [symbol("eq_rml$i") for i=1:num_anticipated_shocks(m)]]

    # Additional states added after solving model
    # Lagged states and observables measurement error
    endogenous_states_postgensys = [
        :y_t1, :c_t1, :i_t1, :w_t1, :π_t1, :L_t1, :Et_π_t]

    # Measurement equation observables
    observables = [[
        :g_y,         # quarterly output growth
        :g_hours,     # aggregate hours growth
        :g_w,         # real wage growth
        :π_gdpdef,    # inflation (GDP deflator)
        :R_n,         # nominal interest rate
        :g_c,         # consumption growth
        :g_i];         # investment growth
        [symbol("R_n$i") for i=1:num_anticipated_shocks(m)]] # compounded nominal rates

    for (i,k) in enumerate(endogenous_states);            m.endogenous_states[k]            = i end
    for (i,k) in enumerate(exogenous_shocks);             m.exogenous_shocks[k]             = i end
    for (i,k) in enumerate(expected_shocks);              m.expected_shocks[k]              = i end
    for (i,k) in enumerate(equilibrium_conditions);       m.equilibrium_conditions[k]       = i end
    for (i,k) in enumerate(endogenous_states);            m.endogenous_states[k]            = i end
    for (i,k) in enumerate(endogenous_states_postgensys); m.endogenous_states_postgensys[k] = i+length(endogenous_states) end
    for (i,k) in enumerate(observables);                  m.observables[k]                  = i end
end


function SmetsWouters(subspec::AbstractString="ss0")

    # Model-specific specifications
    spec               = split(basename(@__FILE__),'.')[1]   
    subspec            = subspec
    settings           = Dict{Symbol,Setting}()
    test_settings      = Dict{Symbol,Setting}()
    rng                = MersenneTwister()        # Random Number Generator
    testing            = false                       
    _filestrings       = SortedDict{Symbol,AbstractString, ForwardOrdering}()
    
    # initialize empty model
    m = SmetsWouters{Float64}(
            # model parameters and steady state values
            @compat(Vector{AbstractParameter{Float64}}()), @compat(Vector{Float64}()), Dict{Symbol,Int}(),

            # model indices
            Dict{Symbol,Int}(), Dict{Symbol,Int}(), Dict{Symbol,Int}(), Dict{Symbol,Int}(), Dict{Symbol,Int}(), Dict{Symbol,Int}(),
                          
            spec,
            subspec,
            settings,
            test_settings,
            rng,
            testing,
            _filestrings)

    # Set settings
    settings_smets_wouters(m)
    default_test_settings(m)
    
    # Initialize parameters
    m <= parameter(:α,      0.24, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     Normal(0.30, 0.05),         fixed=false,
                   description="α: Capital elasticity in the intermediate goods sector's Cobb-Douglas production function.",
                   texLabel="\\alpha")

    m <= parameter(:ζ_p,   0.7813, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     BetaAlt(0.5, 0.1),          fixed=false,
                   description="ζ_p: The Calvo parameter. In every period, (1-ζ_p) of the intermediate goods producers optimize prices. ζ_p of them adjust prices according to steady-state inflation, π_star.",
                   texLabel="\\zeta_p")

    m <= parameter(:ι_p,   0.3291, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     BetaAlt(0.5, 0.15),         fixed=false,
                   description="ι_p: The persistence of last period's inflation in the equation that describes the intertemporal change in prices for intermediate goods producers who cannot adjust prices. The change in prices is a geometric average of steady-state inflation (π_star, with weight (1-ι_p)) and last period's inflation (π_{t-1})).",
                   texLabel="\\iota_p")

    m <= parameter(:δ,      0.025,  fixed=true,
                   description="δ: The capital depreciation rate.", texLabel="\\delta" )     

    ## TODO
    m <= parameter(:Upsilon,  1.000,  (0., 10.),     (1e-5, 0.),      DSGE.Exponential(),    GammaAlt(1., 0.5),          fixed=true,
                   description="Υ: This is the something something.",
                   texLabel="\\mathcal{\\Upsilon}") 

    ## TODO - ask Marc and Marco
    m <= parameter(:Φ,   1.4672, (1., 10.),     (1.00, 10.00),   DSGE.Exponential(),    Normal(1.25, 0.12),         fixed=false,
                   description="Φ: This is the something something.",
                   texLabel="\\Phi")

    m <= parameter(:S′′,       6.3325, (-15., 15.),   (-15., 15.),     Untransformed(),  Normal(4., 1.5),            fixed=false,
                   description="S'': The second derivative of households' cost of adjusting investment.", texLabel="S\\prime\\prime")

    m <= parameter(:h,        0.7205, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     BetaAlt(0.7, 0.1),          fixed=false,
                   description="h: Consumption habit persistence.", texLabel="h")

    ## TODO - ask Marc and Marco
    m <= parameter(:ppsi,     0.2648, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     BetaAlt(0.5, 0.15),         fixed=false,
                   description="ppsi: This is the something something.", texLabel="ppsi")

    ## TODO - ask Marc and Marco
    m <= parameter(:ν_l,     2.8401, (1e-5, 10.),   (1e-5, 10.),     DSGE.Exponential(),    Normal(2, 0.75),            fixed=false,
                   description="ν_l: The coefficient of relative risk aversion on the labor term of households' utility function.", texLabel="\nu_l")
    
    m <= parameter(:ζ_w,   0.7937, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     BetaAlt(0.5, 0.1),          fixed=false,
                   description="ζ_w: (1-ζ_w) is the probability with which households can freely choose wages in each period. With probability ζ_w, wages increase at a geometrically weighted average of the steady state rate of wage increases and last period's productivity times last period's inflation.",
                   texLabel="\\zeta_w")

    #TODO: Something to do with intertemporal changes in wages? – Check in tex file
    m <= parameter(:ι_w,   0.4425, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     BetaAlt(0.5, 0.15),         fixed=false,
                   description="ι_w: This is the something something.",
                   texLabel="\\iota_w")

    m <= parameter(:λ_w,      1.5000,                                                                               fixed=true,
                   description="λ_w: The wage markup, which affects the elasticity of substitution between differentiated labor services.",
                   texLabel="\\lambda_w")     

    m <= parameter(:β, 0.7420, (1e-5, 10.),   (1e-5, 10.),     DSGE.Exponential(),    GammaAlt(0.25, 0.1),        fixed=false,  scaling = x -> 1/(1 + x/100),
                   description="β: Discount rate.",       
                   texLabel="\\beta ")

    m <= parameter(:ψ1,  1.7985, (1e-5, 10.),   (1e-5, 10.00),   DSGE.Exponential(),    Normal(1.5, 0.25),          fixed=false,
                   description="ψ₁: Weight on inflation gap in monetary policy rule.",
                   texLabel="\\psi_1")

    m <= parameter(:ψ2,  0.0893, (-0.5, 0.5),   (-0.5, 0.5),     Untransformed(),  Normal(0.12, 0.05),         fixed=false,
                   description="ψ₂: Weight on output gap in monetary policy rule.",
                   texLabel="\\psi_2")

    m <= parameter(:ψ3, 0.2239, (-0.5, 0.5),   (-0.5, 0.5),     Untransformed(),  Normal(0.12, 0.05),         fixed=false,
                   description="ψ₃: Weight on rate of change of output gap in the monetary policy rule.",
                   texLabel="\\psi_3")

    m <= parameter(:π_star,   0.7000, (1e-5, 10.),   (1e-5, 10.),     DSGE.Exponential(),    GammaAlt(0.62, 0.1),        fixed=true,  scaling = x -> 1 + x/100,
                   description="π_star: The steady-state rate of inflation.",  
                   texLabel="\\pi_*")   

    m <= parameter(:σ_c, 1.2312, (1e-5, 10.),   (1e-5, 10.),     DSGE.Exponential(),    Normal(1.5, 0.37),          fixed=false,
                   description="σ_c: This is the something something.",
                   texLabel="\\sigma_{c}")
    
    m <= parameter(:ρ,      .8258, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     BetaAlt(0.75, 0.10),        fixed=false,
                   description="ρ: The degree of inertia in the monetary policy rule.",
                   texLabel="\\rho")

    m <= parameter(:ϵ_p,     10.000,                                                                               fixed=true,
                   description="ϵ_p: This is the something something.",
                   texLabel="\\varepsilon_{p}")     

    m <= parameter(:ϵ_w,     10.000,                                                                               fixed=true,
                   description="ϵ_w: This is the something something.",
                   texLabel="\\varepsilon_{w}")     
    

    # exogenous processes - level
    m <= parameter(:γ,      0.3982, (-5.0, 5.0),     (-5., 5.),     Untransformed(), Normal(0.4, 0.1),            fixed=false, scaling = x -> x/100, 
                   description="γ: The log of the steady-state growth rate of technology.", 
                   texLabel="\\gamma")

    m <= parameter(:Lmean,  875., (-1000., 1000.), (-1e3, 1e3),   Untransformed(), Normal(-45, 5),   fixed=false,
                   description="Lmean: This is the something something.",
                   texLabel="Lmean")

    m <= parameter(:g_star,    0.1800,                                                                               fixed=true,
                   description="g_star: This is the something something.",
                   texLabel="g_*")  

    # exogenous processes - autocorrelation 
    m <= parameter(:ρ_g,      0.9930, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(),    BetaAlt(0.5, 0.2),           fixed=false,
                   description="ρ_g: AR(1) coefficient in the government spending process.",
                   texLabel="\\rho_g")

    m <= parameter(:ρ_b,      0.2703, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(),    BetaAlt(0.5, 0.2),           fixed=false,
                   description="ρ_b: AR(1) coefficient in the intertemporal preference shifter process.",
                   texLabel="\\rho_b")

    m <= parameter(:ρ_μ,     0.5724, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(),    BetaAlt(0.5, 0.2),           fixed=false,
                   description="ρ_μ: AR(1) coefficient in capital adjustment cost process.",
                   texLabel="\\rho_{\\mu}")

    m <= parameter(:ρ_z,      0.9676, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(),    BetaAlt(0.5, 0.2),           fixed=false,
                   description="ρ_z: AR(1) coefficient in the technology process.",
                   texLabel="\\rho_z")

    m <= parameter(:ρ_λ_f,    0.8692, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(),    BetaAlt(0.5, 0.2),           fixed=false,
                   description="ρ_λ_f: AR(1) coefficient in the price mark-up shock process.",
                   texLabel="\\rho_{\\lambda_f}") 

    m <= parameter(:ρ_λ_w,    0.9546, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(),    BetaAlt(0.5, 0.2),           fixed=false,
                   description="ρ_λ_w: AR(1) coefficient in the wage mark-up shock process.", # CHECK THIS
                   texLabel="\\rho_{\\lambda_f}")

    #monetary policy shock - see eqcond
    m <= parameter(:ρ_rm,     0.3000, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(),    BetaAlt(0.5, 0.2),           fixed=false,
                   description="ρ_rm: AR(1) coefficient in the monetary policy shock process.", # CHECK THIS
                   texLabel="\\rho_{rm}")


    # exogenous processes - standard deviation
    m <= parameter(:σ_g,      0.6090, (1e-8, 5.),      (1e-8, 5.),    DSGE.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
                   description="σ_g: The standard deviation of the government spending process.",
                   texLabel="\\sigma_{g}")
    
    m <= parameter(:σ_b,      0.1818, (1e-8, 5.),      (1e-8, 5.),    DSGE.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
                   description="σ_b: This is the something something.",
                   texLabel="\\sigma_{b}")

    m <= parameter(:σ_μ,     0.4601, (1e-8, 5.),      (1e-8, 5.),    DSGE.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
                   description="σ_μ: The standard deviation of the exogenous marginal efficiency of investment shock process.",
                   texLabel="\\sigma_{\\mu}")

    ## TODO
    m <= parameter(:σ_z,      0.4618, (1e-8, 5.),      (1e-8, 5.),    DSGE.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
                   description="σ_z: This is the something something.",
                   texLabel="\\sigma_{z}")

    m <= parameter(:σ_λ_f,    0.1455, (1e-8, 5.),      (1e-8, 5.),    DSGE.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
                   description="σ_λ_f: The mean of the process that generates the price elasticity of the composite good.  Specifically, the elasticity is (1+λ_{f,t})/(λ_{f_t}).",
                   texLabel="\\sigma_{\\lambda_f}")

    ## TODO
    m <= parameter(:σ_λ_w,    0.2089, (1e-8, 5.),      (1e-8, 5.),    DSGE.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
                   description="σ_λ_w: This is the something something.",
                   texLabel="\\sigma_{\\lambda_w}")

    ## TODO
    m <= parameter(:σ_rm,     0.2397, (1e-8, 5.),      (1e-8, 5.),    DSGE.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
                   description="σ_rm: This is the something something.",
                   texLabel="\\sigma_{rm}")
    
    
    ## TODO - look at etas in eqcond and figure out what they are, then confirm with Marc and Marco
    m <= parameter(:η_gz,       0.0500, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(),
                   BetaAlt(0.50, 0.20), fixed=false,
                   description="η_gz: Has to do with correlation of g and z shocks.",
                   texLabel="\\eta_{gz}")        

    m <= parameter(:η_λ_f,      0.7652, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(),
                   BetaAlt(0.50, 0.20),         fixed=false,
                   description="η_λ_f: This is the something something.",
                   texLabel="\\eta_{\\lambda_f}")

    m <= parameter(:η_λ_w,      0.8936, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(),
                   BetaAlt(0.50, 0.20),         fixed=false,
                   description="η_λ_w: AR(2) coefficient on wage markup shock process.",
                   texLabel="\\eta_{\\lambda_w}")
    

    # steady states
    m <= SteadyStateParameter(:zstar,  NaN, description="steady-state growth rate of productivity", texLabel="\\z_*")
    m <= SteadyStateParameter(:rstar,   NaN, description="steady-state something something", texLabel="\\r_*")
    m <= SteadyStateParameter(:Rstarn,  NaN, description="steady-state something something", texLabel="\\R_*_n")
    m <= SteadyStateParameter(:rkstar,  NaN, description="steady-state something something", texLabel="\\BLAH")
    m <= SteadyStateParameter(:wstar,   NaN, description="steady-state something something", texLabel="\\w_*")
    m <= SteadyStateParameter(:Lstar,   NaN, description="steady-state something something", texLabel="\\L_*") 
    m <= SteadyStateParameter(:kstar,   NaN, description="Effective capital that households rent to firms in the steady state.", texLabel="\\k_*")
    m <= SteadyStateParameter(:kbarstar, NaN, description="Total capital owned by households in the steady state.", texLabel="\\bar{k}_*")
    m <= SteadyStateParameter(:istar,  NaN, description="Detrended steady-state investment", texLabel="\\i_*")
    m <= SteadyStateParameter(:ystar,  NaN, description="steady-state something something", texLabel="\\y_*")
    m <= SteadyStateParameter(:cstar,  NaN, description="steady-state something something", texLabel="\\c_*")
    m <= SteadyStateParameter(:wl_c,   NaN, description="steady-state something something", texLabel="\\wl_c")

    initialize_model_indices!(m)
    initialize_subspec(m)
    steadystate!(m)
    return m
end


#=
doc"""
Inputs: `m::SmetsWouters`

Description: (Re)calculates the model's steady-state values. `steadystate!(m)` must be called whenever the parameters of `m` are updated. 
"""
=#
# (Re)calculates steady-state values
function steadystate!(m::SmetsWouters)
    m[:zstar]    = log(1+m[:γ]) + m[:α]/(1-m[:α])*log(m[:Upsilon])
    m[:rstar]    = exp(m[:σ_c]*m[:zstar]) / m[:β]
    m[:Rstarn]   = 100*(m[:rstar]*m[:π_star] - 1)
    m[:rkstar]   = m[:rstar]*m[:Upsilon] - (1-m[:δ])
    m[:wstar]    = (m[:α]^m[:α] * (1-m[:α])^(1-m[:α]) * m[:rkstar]^(-m[:α]) / m[:Φ])^(1/(1-m[:α]))
    m[:Lstar]    = 1.
    m[:kstar]    = (m[:α]/(1-m[:α])) * m[:wstar] * m[:Lstar] / m[:rkstar]
    m[:kbarstar] = m[:kstar] * (1+m[:γ]) * m[:Upsilon]^(1 / (1-m[:α]))
    m[:istar]    = m[:kbarstar] * (1-((1-m[:δ])/((1+m[:γ]) * m[:Upsilon]^(1/(1-m[:α])))))
    m[:ystar]    = (m[:kstar]^m[:α]) * (m[:Lstar]^(1-m[:α])) / m[:Φ]
    m[:cstar]    = (1-m[:g_star])*m[:ystar] - m[:istar]
    m[:wl_c]     = (m[:wstar]*m[:Lstar])/(m[:cstar]*m[:λ_w])

    return m
end


function settings_smets_wouters(m::SmetsWouters)

    default_settings(m)
    
    # Anticipated shocks
    m <= Setting(:num_anticipated_shocks,         0, "Number of anticipated policy shocks")
    m <= Setting(:num_anticipated_shocks_padding, 20, "Padding for anticipated policy shocks")
    m <= Setting(:num_anticipated_lags,  26, "Number of periods back to incorporate zero bound expectations")

    # Estimation
    m <= Setting(:reoptimize, true, true, "reop", "whether to re-find mode")
    m <= Setting(:recalculate_hessian, true, true, "ch", "whether to calculate the hessian")

    # Data vintage
    m <= Setting(:data_vintage, "150827", true, "vint", "Date of data")

    m.settings
end

