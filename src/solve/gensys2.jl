function gensys_cplus(m::AbstractDSGEModel, Γ0::Matrix{Float64}, Γ1::Matrix{Float64}, C::Vector{Float64},
                      Ψ::Matrix{Float64}, Π::Matrix{Float64}, TTT::Matrix{Float64}, RRR::Matrix{Float64},
                      CCC::Vector{Float64})
    T_switch = get_setting(m, :n_rule_periods) + 1 #get_setting(m, :n_regimes) - 1
    exp_eq_ind = sum(Π, dims = 2)
    Γ0_til = zeros(size(Γ0))
    Γ1_til = zeros(size(Γ1))
    Γ2_til = zeros(size(Γ0))
    C_til = C
    Ψ_til = Ψ

    for row in 1:length(exp_eq_ind)
        if exp_eq_ind[row] == 0
            # Not expectational equation
            Γ1_til[row, :] = Γ1[row, :]
            Γ0_til[row, :] = Γ0[row, :]
            Γ2_til[row, :] .= 0.
        else
            # expectational equation
            Γ0_til[row, findfirst(Γ1[row, :] .> 0)] = -1
            Γ2_til[row, findfirst(Γ0[row, :] .> 0)] = 1
        end
    end

    Tcal = Vector{Matrix{Float64}}(undef, T_switch)
    Rcal = Vector{Matrix{Float64}}(undef, T_switch)
    Ccal = Vector{Vector{Float64}}(undef, T_switch)

    Tcal[end] = TTT
    Rcal[end] = RRR
    Ccal[end] = CCC

    for t = 1:(T_switch-1)
        Tcal[end-t] = (Γ2_til*TTT + Γ0_til)\Γ1_til
        Rcal[end-t] = (Γ2_til*TTT + Γ0_til)\Ψ_til
        Ccal[end-t] = (Γ2_til*TTT + Γ0_til)\(C_til - Γ2_til*CCC)

        TTT = Tcal[end-t]
        CCC = Ccal[end-t]
        RRR = Rcal[end-t]
    end
    return Tcal, Rcal, Ccal
end


function gensys_cplus(m::AbstractDSGEModel, Γ0s::Vector{Matrix{Float64}}, Γ1s::Vector{Matrix{Float64}},
                      Cs::Vector{Vector{Float64}},
                      Ψs::Vector{Matrix{Float64}}, Πs::Vector{Matrix{Float64}},
                      TTT::Matrix{Float64}, RRR::Matrix{Float64},
                      CCC::Vector{Float64})
    T_switch = get_setting(m, :n_rule_periods) + 1
    Γ0_tils = Vector{Matrix{Float64}}(undef, length(Γ0s))
    Γ1_tils = Vector{Matrix{Float64}}(undef, length(Γ0s))
    Γ2_tils = Vector{Matrix{Float64}}(undef, length(Γ0s))
    C_tils = Vector{Vector{Float64}}(undef, length(Γ0s))
    Ψ_tils = Vector{Matrix{Float64}}(undef, length(Γ0s))
    for i in 1:length(Γ0s)
        exp_eq_ind = sum(Πs[i], dims = 2)
        Γ0_tils[i] = zeros(size(Γ0s[i]))
        Γ1_tils[i] = zeros(size(Γ1s[i]))
        Γ2_tils[i] = zeros(size(Γ0s[i]))
        C_tils[i] = Cs[i]
        Ψ_tils[i] = Ψs[i]

        for row in 1:length(exp_eq_ind)
            if exp_eq_ind[row] == 0
                # Not expectational equation
                Γ1_tils[i][row, :] = Γ1s[i][row, :]
                Γ0_tils[i][row, :] = Γ0s[i][row, :]
                Γ2_tils[i][row, :] .= 0.
            else
                # expectational equation
                Γ0_tils[i][row, findfirst(Γ1s[i][row, :] .> 0)] = -1
                Γ2_tils[i][row, findfirst(Γ0s[i][row, :] .> 0)] = 1
            end
        end
    end

    Tcal = Vector{Matrix{Float64}}(undef, T_switch)
    Rcal = Vector{Matrix{Float64}}(undef, T_switch)
    Ccal = Vector{Vector{Float64}}(undef, T_switch)

    Tcal[end] = TTT
    Rcal[end] = RRR
    Ccal[end] = CCC

    for t = 1:(T_switch-1)
        Tcal[end-t] = (Γ2_tils[end-t]*TTT + Γ0_tils[end-t])\Γ1_tils[end-t]
        Rcal[end-t] = (Γ2_tils[end-t]*TTT + Γ0_tils[end-t])\Ψ_tils[end-t]
        Ccal[end-t] = (Γ2_tils[end-t]*TTT + Γ0_tils[end-t])\(C_tils[end-t] - Γ2_tils[end-t]*CCC)

        TTT = Tcal[end-t]
        CCC = Ccal[end-t]
        RRR = Rcal[end-t]
    end
    return Tcal, Rcal, Ccal
end



#=

function gensys_cplus_old(m::AbstractDSGEModel, Γ0::Matrix{Float64}, Γ1::Matrix{Float64}, C::Vector{Float64},
                      Ψ::Matrix{Float64}, TTT::Matrix{Float64}, RRR::Matrix{Float64},
                      CCC::Vector{Float64})
    T_switch = 2 #get_setting(m, :T_switch)
    nstates = n_states(m)
    nant = n_mon_anticipated_shocks(m)
    n_exo = length(m.exogenous_shocks)
    n_shocks = n_exo - nant
    ind_states_exp = findall(x -> occursin("E", string(x)), collect(keys(m.endogenous_states)))
    ind_eq_exp = findall(x -> occursin("E", string(x)), collect(keys(m.equilibrium_conditions)))

    non_exp_exp_keys = map(x-> Symbol(string(x)[2:end]), collect(keys(m.endogenous_states))[ind_states_exp])
    non_exp_exp_inds = findall(x->x in non_exp_exp_keys, collect(keys(m.endogenous_states)))

    n_expected = length(ind_states_exp)

    num_z = nstates - n_expected #n_end+n_exo
    # Initializing the Γ2 matrix
    Γ2_tilde = zeros(nstates, nstates)

    # taking the Expectations
    Γ0_tilde = copy(Γ0)
    Γ0_tilde[ind_eq_exp,:] = Γ0_tilde[ind_eq_exp,:]*0.
    Γ0_tilde[:,ind_states_exp] = Γ0_tilde[:,ind_states_exp]*0.

    # moving the expectation from Γ0 to Γ2
  #  for i_ind = 1:n_expected
    Γ2_tilde[:, non_exp_exp_inds] = Γ0[:, ind_eq_exp]

    # identity map of expectational states
    Γ2_tilde[ind_eq_exp, non_exp_exp_inds] .= 1
    Γ0_tilde[ind_eq_exp, ind_states_exp] .= -1
  #  end

    #taking expectations
    #Γ2_tilde(num_z+1:end-1,:) = []
    #Γ2_tilde(:,num_z+1:end-1) = []

    #taking the Expectations
    Γ1_tilde = Γ1
    Γ1_tilde[ind_eq_exp,:] .= 0 #[]
    Γ1_tilde[:,ind_states_exp] .= 0 #[]
    # Taking expectation
    Ψ_tilde = Ψ
    Ψ_tilde[ind_eq_exp, :] =. 0 #[]
    Ψ_tilde[:,n_shocks+1:end] =. 0 #[] #taking out anticipated policy shocks
    # taking expectations out of C
    C_tilde = C #[C(1:num_z)C(end)]


    ## This solves for the TTT and RRR Matrix using the methods in "Solving
    ## Linear Rational Expectations Methods with Anticipated Policy Changes
    # taking out expectation equations
    T11 = TTT #[Not(ind_eq_exp), Not(ind_states_exp)]
    R11 = RRR #[Not(ind_eq_exp), 1:n_shocks]
    C11 = CCC #[Not(ind_eq_exp)]
#=    T11 = TTT
    T11[ind_eq_exp,:] .= NaN
    T11[:,ind_state_exp] .= NaN

    R11 = RRR
    R11[ind_eq,:] .= NaN
    R11[:,nshocks+1:end] .= NaN  #taking out anticipated policy shocks

    C11 = CCC
    C11[ind_eq_exp] .= NaN =#



    #initialize Tcal and Rcal matrices
    Tcal = zeros(size(T11, 1), size(T11, 2), T_switch+1) #zeros(num_z+1, num_z+1, T_switch+1)
    Rcal = zeros(size(R11, 1), size(R11, 2), T_switch+1) #num_z+1, nshocks, T_switch+1)
    Ccal = zeros(length(C11), T_switch+1) #num_z+1, T_switch+1)

    Tcal[:,:,end] = T11
    Rcal[:,:,end] = R11
    Ccal[:,end] = C11

    for t = 1:T_switch
        Tcal[:,:,end-t] = (Γ2_tilde*T11 + Γ0_tilde)\Γ1_tilde
        Rcal[:,:,end-t] = (Γ2_tilde*T11 + Γ0_tilde)\Ψ_tilde
        Ccal[:,end-t] = (Γ2_tilde*T11 + Γ0_tilde)\(C_tilde - Γ2_tilde*C11)


        # Tcal(:,:,end-t) = pinv(Γ2_tilde*T11 + Γ0_tilde)*Γ1_tilde
        # Rcal(:,:,end-t) = pinv(Γ2_tilde*T11 + Γ0_tilde)*PSI_tilde
        # Ccal(:,end-t) = pinv(Γ2_tilde*T11 + Γ0_tilde)*(C_tilde - Γ2_tilde*C11)

        T11 = Tcal[:,:,end-t]
    end
        C11 = Ccal[:,end-t]
        R11= Rcal[:,:,end-t]
    end
end
=#
