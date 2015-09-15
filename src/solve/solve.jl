using Debug
# using MATLAB
include("../../test/util.jl")


# Outputs TTT, RRR, CCC - matrices of the state transition equation:
#   S_t = TTT*S_{t-1} + RRR*ε_t + CCC
@debug function solve(model::AbstractDSGEModel; verbose::Bool=false)

    # Get equilibrium condition matrices
    Γ0, Γ1, C, Ψ, Π  = eqcond(model) #(, augmented=false)

    ## These big, commented-out "if verbose" blocks  are all code that checks the output of eqcond
    ## against the identical point in the Matlab code. It will not work as written because
    ## username is undefined, and unnecessary binary files are not included in the commit.
    ## The code will be deleted in a commit immediately following this one; just want to keep it
    ## here in the comments in case we want to refer back to it at some point.
    
    ## path =
    ## "/home/$username/.julia/v0.3/DSGE/test/estimate/metropolis_hastings/m990-no_reoptimize_no_recalc_hessian/"

    ## solvepath = "/home/$username/.julia/v0.3/DSGE/test/estimate/solve/"
    
    ## if verbose
    
    ##     println("\ntesting model parameters (from inside valid0 loop) against matlab")
    ##     mf = MatFile("$solvepath/para4.mat")
    ##     para_mat = get_variable(mf, "para") jl_para = similar(para_mat)
    ##     for (i,p) in enumerate(model.parameters) if i <=
    ##         length(para_mat) jl_para[i] = p.value end end
    ##     test_matrix_eq(para_mat, jl_para, ε=1e-12, noisy=true)
    ##     close(mf)
        
    ##     println("\ntesting steady state parameters (from inside valid0 loop) against matlab")
    ##     mf = MatFile("$solvepath/steadystate.mat")
    ##     ss_para_mat = get_variable(mf, "steadyStateParams")
    ##     test_matrix_eq(ss_para_mat, float64(model.steady_state), ε=1e-12, noisy=true)
    ##     close(mf)

    ##     ss_diffs = ss_para_mat - float64(model.steady_state)
    ##     ss_names = [:zstar, :rstar, :Rstarn, :rkstar, :wstar, :Lstar, :kstar, :kbarstar, :istar, :ystar,
    ##         :cstar, :wl_c, :nstar, :vstar, :zeta_spsigw, :zeta_spmue, :zeta_nRk, :zeta_nR, :zeta_nqk,
    ##         :zeta_nn, :zeta_nmue, :zeta_nsigw]

    ##     for i in 1:length(ss_diffs)
    ##         @printf "steadystate = %s, diff = %f, mat = %f, jl = %f\n" ss_names[i] ss_diffs[i] ss_para_mat[i] model.steady_state[i]
    ##     end
                    
    ##     println("\ntesting eqcond_mat (from inside valid0 loop) against results of eqcond")
    ##     mf = MatFile("$path/eqcond_mat.mat")

    ##     G0 = get_variable(mf, "G0")
    ##     println(typeof(Γ0),typeof(G0))
        
    ##     println("testing G0")    
    ##     test_matrix_eq(G0, Γ0; ε=1e-12, noisy=true)
    ##     #abs_diff_entries, opp_sign_entries = find_matrix_diffs(G0, Γ0; ε=1e-12)
    ##     get_diff_symbols(model, G0, Γ0, ε=1e-12)
        
    ##     println("testing G1")
    ##     G1 = get_variable(mf, "G1")
    ##     test_matrix_eq(G1, Γ1; ε=1e-12, noisy=true)

    ##     println("testing C")
    ##     C_mat = get_variable(mf, "C")
    ##     test_matrix_eq(C_mat, C, ε=1e-12, noisy=true)

    ##     println("testing PSI")
    ##     PSI = get_variable(mf, "PSI")
    ##     test_matrix_eq(PSI, Ψ, ε=1e-12, noisy=true)

    ##     println("testing PIE")
    ##     PIE = get_variable(mf, "PIE")
    ##     test_matrix_eq(PIE, Π, ε=1e-12, noisy=true)

    ##     close(mf)
        
    ##     @bp

    ## end # of verbose 
    
    # Solve model
    TTT_gensys, CCC_gensys, RRR_gensys = gensys(Γ0, Γ1, C, Ψ, Π, 1+1e-6)
    TTT_gensys = real(TTT_gensys)
    RRR_gensys = real(RRR_gensys)
    CCC_gensys = reshape(CCC_gensys, length(CCC_gensys), 1)

    ## if verbose
    ##     mf = MatFile("$path/TTT_gensys.mat")
    ##     TTT_gensys_mat = get_variable(mf, "TTT")
    ##     RRR_gensys_mat = get_variable(mf, "RRR")
    ##     close(mf)

    ##     println("Testing TTT and RRR when they're returned from gensys")
    ##     println("TTT_gensys:")
    ##     test_matrix_eq(TTT_gensys_mat, TTT_gensys; ε=1e-12, noisy=true)
    ##     println("RRR_gensys:")
    ##     test_matrix_eq(RRR_gensys_mat, RRR_gensys; ε=1e-12, noisy=true)

    ##     @bp
    ## end
    
    # Augment states
    TTT, RRR, CCC = augment_states(model, TTT_gensys, RRR_gensys, CCC_gensys)

    ## if verbose
    ##     mf = MatFile("$path/TTT_solve.mat")
    ##     TTT_mat = get_variable(mf, "TTT")
    ##     RRR_mat = get_variable(mf, "RRR")
    ##     close(mf)

    ##     println("Testing the TTT matrix returned from solve against that returned from dsgesolv before valid0 loop:")
    ##     test_matrix_eq(TTT_mat, TTT; ε=1e-12, noisy=true);
        

    ##     println("Testing the RRR matrix returned from solve against that returned from dsgesolv before valid0 loop:")
    ##     test_matrix_eq(RRR_mat, RRR; ε=1e-12, noisy=true);


    ##     mf = MatFile("$path/TTT_solve_valid0.mat")
    ##     TTT_mat_valid0 = get_variable(mf, "TTT_old")
    ##     RRR_mat_valid0 = get_variable(mf, "RRR_old")
    ##     close(mf)

    ##     println("Testing the TTT matrix returned from solve against that returned from dsgesolv in the valid0 loop:")
    ##     test_matrix_eq(TTT_mat_valid0, TTT; ε=1e-12, noisy=true);

    ##     println("Testing the RRR matrix returned from solve against that returned from dsgesolv in the valid0 loop:")
    ##     test_matrix_eq(RRR_mat_valid0, RRR; ε=1e-12, noisy=true);
        
    ##     @bp
    ## end # of verbose
    
    return TTT, RRR, CCC
end

# Some of our observables are growth rates, which is calculated as a
# linear combination of a present and lagged state. To capture the lagged state,
# we assign to it an index. In addition, we also need to expand the
# matrices of the state transition equation to accommodate the extra state.
# In dsgesolv.m, AddTTT is appended to TTT to capture the lagged state
# value in the current state vector, while AddRRR and AddCCC augment
# RRR and CCC with the appropriate number of zeros.

# These additional states are added after the model is solved to reduce the load on gensys
function augment_states{T<:FloatingPoint}(m::AbstractDSGEModel, TTT::Matrix{T}, RRR::Matrix{T}, CCC::Matrix{T})
    endo = m.endogenous_states
    endo_addl = m.endogenous_states_postgensys
    exo = m.exogenous_shocks

    n_endo = num_states(m)
    n_exo = num_shocks_exogenous(m)
    @assert (n_endo, n_endo) == size(TTT)
    @assert (n_endo, n_exo) == size(RRR)
    @assert (n_endo, 1) == size(CCC)

    # Initialize augmented matrices
    numAdd = 12
    TTT_aug = zeros(n_endo + numAdd, n_endo + numAdd)
    TTT_aug[1:n_endo, 1:n_endo] = TTT
    RRR_aug = [RRR; zeros(numAdd, n_exo)]
    CCC_aug = [CCC; zeros(numAdd, 1)]



    ### TTT modifications

    # Track Lags
    TTT_aug[endo_addl[:y_t1], endo[:y_t]] = 1.0
    TTT_aug[endo_addl[:c_t1], endo[:c_t]] = 1.0
    TTT_aug[endo_addl[:i_t1], endo[:i_t]] = 1.0
    TTT_aug[endo_addl[:w_t1], endo[:w_t]] = 1.0
    TTT_aug[endo_addl[:pi_t1], endo[:pi_t]] = 1.0
    TTT_aug[endo_addl[:L_t1], endo[:L_t]]  = 1.0
    TTT_aug[endo_addl[:u_t1], endo[:u_t]] = 1.0

    # Expected inflation
    TTT_aug[endo_addl[:Et_pi_t], 1:n_endo] = (TTT^2)[endo[:pi_t], :]

    # The 8th column of AddTTT corresponds to "v_lr" which is set equal to
    # e_lr –measurements errors for the two real wage observables built in
    # as exogenous structural shocks.
    TTT_aug[endo_addl[:lr_t], endo_addl[:lr_t]] = m[:ρ_lr].scaledvalue
    TTT_aug[endo_addl[:tfp_t], endo_addl[:tfp_t]] = m[:ρ_tfp].scaledvalue
    TTT_aug[endo_addl[:e_gdpdef], endo_addl[:e_gdpdef]] = m[:ρ_gdpdef].scaledvalue
    TTT_aug[endo_addl[:e_pce], endo_addl[:e_pce]] = m[:ρ_pce].scaledvalue



    ### RRR modfications

    # Expected inflation
    RRR_aug[endo_addl[:Et_pi_t], :] = (TTT*RRR)[endo[:pi_t], :]

    # Measurement Error on long rate
    RRR_aug[endo_addl[:lr_t], exo[:lr_sh]] = 1.0

    # Measurement Error on TFP
    RRR_aug[endo_addl[:tfp_t], exo[:tfp_sh]] = 1.0

    # Measurement Error on GDP Deflator
    RRR_aug[endo_addl[:e_gdpdef], exo[:gdpdef_sh]] = 1.0

    # Measurement Error on Core PCE
    RRR_aug[endo_addl[:e_pce], exo[:pce_sh]] = 1.0



    ### CCC Modifications

    # Expected inflation
    CCC_aug[endo_addl[:Et_pi_t], :] = (CCC + TTT*CCC)[endo[:pi_t], :]



    return TTT_aug, RRR_aug, CCC_aug
end
