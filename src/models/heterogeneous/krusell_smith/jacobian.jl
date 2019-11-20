function jacobian(m::KrusellSmith)

    # Load in endogenous state and eq cond indices
    endo = augment_model_states(m.endogenous_states_unnormalized,
                                n_model_states_unnormalized(m))
    eq   = m.equilibrium_conditions

    # Load in grid settings
    nw::Int = get_setting(m, :nw)
    ns::Int = get_setting(m, :ns)

    zhi::Float64 = get_setting(m, :zhi)
    zlo::Float64 = get_setting(m, :zlo)
    smoother::Float64 = get_setting(m, :smoother)

    wgrid = get_grid(m, :wgrid)
    sgrid = get_grid(m, :sgrid)

    # Make the Jacobian
    JJ = zeros(4*nw+2,8*nw+4)

    # Create auxiliary variables
    g::Array{Float64,1} = trunc_lognpdf(sgrid, m[:μ_s].value, m[:σ_s].value)                   # density g (of skill, s) evaluated on sgrid

    m[:Wstar] = MPL(1.0, m[:Lstar].value, m[:Kstar].value, m[:α].value)            # Steady state wages
    m[:Rstar] = MPK(1.0, m[:Lstar].value, m[:Kstar].value, m[:α].value, m[:δ].value)     # Steady state net return on capital

    # KF equation
    KFmollificand::Array{Float64, 2} = (repeat(wgrid.points, 1, nw*ns) - m[:Rstar]*repeat(wgrid.points'-m[:cstar].value', nw, ns))/m[:Wstar] - kron(sgrid.points',ones(nw, nw))
    μRHS::Array{Float64,1} = mollifier.(KFmollificand, zhi, zlo, smoother)*kron(sgrid.weights.*(g/m[:Wstar]), wgrid.weights.*m[:μstar].value)

    # These matrices (from the Euler equation) correspond to the matrices in the PDF documentation on Dropbox
    # (note that Γ is Λ removing the Σ_i^nx, and the ω^x_i terms from the summation.
    ξ = Array{Float64, 2}(undef, nw, nw) #zeros(nw, nw)
    Γ = Array{Float64, 2}(undef, nw, nw) #zeros(nw, nw)

    # These matrices (from the KF equation) corresp. to the lyx document on dropbox, script A and ζ
    AA = Array{Float64,2}(undef, nw, nw)
    ζ = Array{Float64, 2}(undef, nw, nw)

    # i_w is the ith entry in the nw discretization of x
    # i_wp is i_w' the grid of potential future x's
    for i_w in 1:nw
        for i_wp in 1:nw
            # The following hold the sums over i=1̣̣̣̣, ..., ns
            sum_ns_ξ::Float64 = 0.0
            sum_ns_Γ::Float64 = 0.0
            sum_ns_AA::Float64 = 0.0
            sum_ns_ζ::Float64 = 0.0
            # The first part of ξ and Γ which only relies on i_w' so we calculate outside of 1:ns loop
            front = ((m[:cstar].value[i_wp].^(-m[:γ]))/m[:Wstar])
            first_part_mollificand = (wgrid.points[i_wp] - m[:Rstar] * (wgrid.points[i_w] -
                                                                        m[:cstar].value[i_w]))/m[:Wstar]
            for iss in 1:ns
                mollificand::Float64 = first_part_mollificand - sgrid.points[iss]
                # dm and mm are the mollifiers which correspond to q(mollificand) in the paper
                dm = dmollifier(mollificand, zhi, zlo, smoother)
                mm = mollifier(mollificand, zhi, zlo, smoother)
                skill_distrib_weights = g[iss] * sgrid.weights[iss]

                sum_ns_ξ += front * dm * skill_distrib_weights
                sum_ns_Γ += front * mm * skill_distrib_weights
                sum_ns_AA+= (1.0/m[:Wstar]) * mm * skill_distrib_weights
                sum_ns_ζ += (1.0/m[:Wstar]) * dm * skill_distrib_weights
            end
            ξ[i_w,i_wp] = sum_ns_ξ
            Γ[i_w,i_wp] = sum_ns_Γ
            AA[i_w,i_wp] = sum_ns_AA
            ζ[i_w,i_wp] = sum_ns_ζ
        end
    end

    dRdK::Float64 = -(1-m[:α])*((m[:Rstar]+m[:δ]-1)/m[:Kstar])
    dWdK::Float64 = (m[:α]*m[:Wstar]/m[:Kstar])
    dRdZ::Float64 = m[:Rstar]+m[:δ]-1.0
    dWdZ::Float64 = m[:Wstar].value
    lRHS = m[:β]*m[:Rstar]*Γ*wgrid.weights

    # Fill in Jacobian

    # Euler equation (EE)
    term1_EE::Array{Float64, 1} = m[:β]*Γ*wgrid.weights - m[:β]*(m[:Rstar]/m[:Wstar]) * ((wgrid.points-m[:cstar].value) .* (ξ*wgrid.weights))
    term2_EE::Array{Float64, 1} = -(1.0/m[:Wstar]) * (lRHS+((m[:β]*m[:Rstar])/m[:Wstar]) * ξ.*(repeat(wgrid.points', nw, 1) - m[:Rstar]*repeat(wgrid.points-m[:cstar].value,1,nw))*wgrid.weights)

    JJ[eq[:eq_euler], endo[:K′_t]]  = term1_EE*dRdK + term2_EE*dWdK
    JJ[eq[:eq_euler], endo[:z′_t]] = term1_EE*dRdZ + term2_EE*dWdZ

    JJ[eq[:eq_euler], endo[:l′_t]] = m[:β]*m[:Rstar]*Γ*Matrix(Diagonal(((wgrid.weights.*(m[:lstar].value.^(-(1.0+m[:γ])/m[:γ]))
                                                                         .*(m[:lstar].value.^(-1.0/m[:γ]) .<= wgrid.points))./m[:cstar].value)))

    JJ[eq[:eq_euler], endo[:l_t]]  = -(eye(nw) + Matrix(Diagonal((m[:β]/m[:γ])*(m[:Rstar]*m[:Rstar]/m[:Wstar])*(ξ*wgrid.weights)
                                                                 .*(m[:lstar].value.^(-(1.0+m[:γ])/m[:γ])).*(m[:lstar].value.^(-1.0/m[:γ]) .<= wgrid.points))))

    # KF Equation
    JJ[eq[:eq_kolmogorov_fwd], endo[:μ_t1]]   = AA'*Matrix(Diagonal(wgrid.weights))

    JJ[eq[:eq_kolmogorov_fwd], endo[:l_t1]] = (-(m[:Rstar]/m[:Wstar])*(1.0/m[:γ])*ζ') * Matrix(Diagonal( m[:μstar].value.*wgrid.weights.*(m[:lstar].value.^(-(1.0+m[:γ])/m[:γ])) .*(m[:lstar].value.^(-1.0/m[:γ]) .<= wgrid.points)))

    term1_KF::Array{Float64,1} = -(1.0/m[:Wstar])*(ζ'.*repeat((wgrid.points-m[:cstar].value)',nw,1))*(m[:μstar].value.*wgrid.weights)
    term2_KF::Array{Float64,1} = -(μRHS/m[:Wstar] + (1.0/(m[:Wstar]*m[:Wstar]))*(ζ'.*(repeat(wgrid.points,1,nw) - m[:Rstar]*repeat( (wgrid.points-m[:cstar].value)',nw,1)))*(m[:μstar].value.*wgrid.weights))

    JJ[eq[:eq_kolmogorov_fwd], endo[:K_t]]  = term1_KF*dRdK + term2_KF*dWdK
    JJ[eq[:eq_kolmogorov_fwd], endo[:z_t]]   = term1_KF*dRdZ + term2_KF*dWdZ

    JJ[eq[:eq_kolmogorov_fwd], endo[:μ_t]]    = -eye(nw)

    # DEFN of LM(t+1) = M(t)
    JJ[eq[:eq_μ], endo[:μ′_t1]]  = eye(nw)

    JJ[eq[:eq_μ], endo[:μ_t]]    = -eye(nw)

    # DEFN of LELL(t+1) = ELL(t)
    JJ[eq[:eq_l], endo[:l′_t1]]= eye(nw)

    JJ[eq[:eq_l], endo[:l_t]]  = -eye(nw)

    # LOM K
    JJ[eq[:eq_K_law_of_motion], endo[:l_t]]  = -(1.0/m[:γ])*m[:μstar].value'*(Matrix(Diagonal(wgrid.weights.*(m[:lstar].value.^(-(1.0+m[:γ])/m[:γ]))
                                                                                              .*(m[:lstar].value.^(-1.0/m[:γ]) .<= wgrid.points))))

    JJ[eq[:eq_K_law_of_motion], endo[:μ_t]]  .= -(wgrid.points-m[:cstar].value)'*Matrix(Diagonal(wgrid.weights))

    JJ[eq[:eq_K_law_of_motion], endo[:K′_t]]  .= 1.0

    # TFP
    JJ[eq[:eq_TFP], endo[:z′_t]]   .= 1.0

    JJ[eq[:eq_TFP], endo[:z_t]]    .= -m[:ρ_z]

    if !m.testing && get_setting(m, :normalize_distr_variables)
        JJ = normalize(m, JJ)
    end

    return JJ
end

function normalize(m::KrusellSmith, JJ::Matrix{Float64})
    nw = get_setting(m, :nw)

    P = eye(nw)
    P[:,1] = ones(nw)

    Q = Matrix{Float64}(undef, nw, nw)
    Q::Matrix{Float64}, _ = qr(P)

    S = Matrix{Float64}(undef, nw, nw-1)
    S::Matrix{Float64} = Q[:, 2:end]'

    Qf = zeros(4*nw, 4*nw+2)
    Qf[1:nw, 1:nw] = eye(nw)
    Qf[nw+1:2*nw-1, nw+1:2*nw] = S
    Qf[2*nw:3*nw-2, 2*nw+1:3*nw] = S
    Qf[3*nw-1:4*nw, 3*nw+1:4*nw+2] = eye(nw+2)

    # Qx and Qy are for normalizing any variables that represent distributions
    # the S component is for normalizing a distribution, the other identity portions
    # are for non-distributional variables
    # The distribution in the states, x, is the lagged μ
    # The distribution in the jumps, y, is the current μ
    Qx = cat(S, eye(nw), [1], [1], dims = [1 2])   # 161 x 162
    Qy = cat(S, eye(nw), dims = [1 2])             # 159 x 160

    m <= Setting(:n_predetermined_variables, size(Qx, 1))

    # whenever outputs and/or inputs are densities,
    # pre/postmultiply the jacobians from step 2 by an appropriate
    # matrix so things integrate to 1

    Qright = cat(Qx',Qy',Qx',Qy', dims = [1 2])
    Jac1 = Qf*JJ*Qright

    return Jac1
end
