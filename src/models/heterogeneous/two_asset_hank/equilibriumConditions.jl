using SparseArrays
function equilibriumConditions(vars)

    #----------------------------------------------------------------
    # Housekeeping
    #----------------------------------------------------------------
    # Unpack variables
    V						= reshape(vars[1:n_v] .+ vars_SS[1:n_v],I,J,N)  # value function
    g						= vars[n_v+1:n_v+n_g] .+ vars_SS[n_v+1:n_v+n_g] # distribution
    dab_aux					= reshape(dab,I*J*N,1)
    dab_g_aux               = reshape(dab_g,I_g*J_g*N,1)
    g_end					= (1 - sum(g .* dab_g_aux[1:end-1])) / dab_g[I_g,J_g,N]
    gg						= [g;g_end]
    K						= vars[n_v+n_g+1] + vars_SS[n_v+n_g+1]     # aggregate capital
    r_b                     = vars[n_v+n_g+2] + vars_SS[n_v+n_g+2]

    if aggregate_variables == 1
        aggY                    = vars[n_v+n_g+3] + vars_SS[n_v+n_g+3] # aggregate output
        aggC                    = vars[n_v+n_g+4] + vars_SS[n_v+n_g+4] # aggregate consumption
    elseif distributional_variables == 1
        C_Var                   = vars[n_v+n_g+3] + vars_SS[n_v+n_g+3] # consumption inequality
	    earn_Var                = vars[n_v+n_g+4] + vars_SS[n_v+n_g+4] # earnings inequality
    elseif distributional_variables_1 == 1
        C_WHTM                  = vars[n_v+n_g+3] + vars_SS[n_v+n_g+3] # consumption of wealthy hand-to-mouth
	    C_PHTM                  = vars[n_v+n_g+4] + vars_SS[n_v+n_g+4] # consumption of poor hand-to-mouth
        #C_NHTM                  = vars[n_v+n_g+3] + vars_SS[n_v+n_g+3]
    end

    aggZ      				= vars[n_v+n_g+n_p+1] + vars_SS[n_v+n_g+n_p+1]      # aggregate Z

    V_Dot					= vars[nVars+1:nVars+n_v]
    g_Dot					= vars[nVars+n_v+1:nVars+n_v+n_g]
    aggZ_Dot     			= vars[nVars+n_v+n_g+n_p+1]

    VEErrors				= vars[2*nVars+1:2*nVars+n_v]
    aggZ_Shock   			= vars[2*nVars+nEErrors+1]

    # Prices
    if permanent == 0
        w                   = (1 - aalpha) * (K ^ aalpha) * exp(aggZ) ^ (1-aalpha) * n_SS ^ (-aalpha)
        r_a                 = aalpha * (K ^ (aalpha - 1)) * (exp(aggZ) * n_SS) ^ (1 - aalpha) - ddelta
    elseif permanent == 1
        w                   = (1 - aalpha) * (K ^ aalpha) * (n_SS) ^ (-aalpha)
        r_a                 = aalpha * (K ^ (aalpha - 1)) * (n_SS) ^ (1 - aalpha) - ddelta
    end

    # Auxiliary variables
    r_b_borr				= r_b .+ borrwedge_SS

    #----------------------------------------------------------------
    # HJB Equation
    #----------------------------------------------------------------

    # Set of grids (calls external script)
    include("Subroutines/set_grids.jl")

    # Other necessary objects
    r_a_grid	 = repeat([r_a],I,J,N)
    l_grid		 = permutedims(repeat(ones(N,1),1,I,J), [2 3 1])
    y_shock      = y .* exp.(kappa * aggZ * (y .- y_mean) ./ std(y))
    y_shock_mean = dot(y_shock, y_dist) #y_shock * y_dist
    y_shock      = y_shock ./ y_shock_mean .* y_mean
    y_grid       = permutedims(repeat(y_shock',1,I,J), [2 3 1])

    ###
    # Compute derivatives w.r.t. illiquid assets a
    ###

    # Preparations
    VaF 			= 0*V
    VaB 			= 0*V
    Vamin 			= 0

    # Forward difference
    VaF[:,1:J-1,:] 	= (V[:,2:J,:]-V[:,1:J-1,:])./daf_grid[:,1:J-1,:]
    VaF[:,1:J-1,:] 	= max.(VaF[:,1:J-1,:],Vamin)

    # Backward difference
    VaB[:,2:J,:] 	= (V[:,2:J,:]-V[:,1:J-1,:])./dab_grid[:,2:J,:]
    VaB[:,2:J,:] 	= max.(VaB[:,2:J,:],Vamin)

    ###
    # Compute derivatives w.r.t. liquid assets b
    ###

    # Preparations (necessary to ensure that everything is a dual number, required by derivative software)
    VbF 			= 0*V
    VbB 			= 0*V
    Vbmin 			= 1e-8

    # Forward difference
    VbF[1:I-1,:,:]	= (V[2:I,:,:]-V[1:I-1,:,:])./dbf_grid[1:I-1,:,:]
    VbF[1:I-1,:,:] 	= max.(VbF[1:I-1,:,:],Vbmin)

    # Backward difference
    VbB[2:I,:,:] 	= (V[2:I,:,:]-V[1:I-1,:,:])./dbb_grid[2:I,:,:]
    VbB[2:I,:,:] 	= max.(VbB[2:I,:,:],Vbmin)

    ###
    # Consumption decision
    ###

    # Preparations
    cF 				= 0 * V
    sF 				= 0 * V
    HcF 			= 0 * V

    cB 				= 0 * V
    sB 				= 0 * V
    HcB 			= 0 * V

    c0 				= 0 * V
    s0 				= 0 * V
    Hc0 			= 0 * V

    # Decisions conditional on a particular direction of derivative
    cF[1:I-1,:,:] 	= VbF[1:I-1,:,:].^(-1/ggamma)
    cF[I,:,:] 		= zeros(1,J,N)
    if permanent == 1
        sF[1:I-1,:,:] 	= ((1-xxi)-tau_I) * w * l_grid[1:I-1,:,:] .* y_grid[1:I-1,:,:] .+ b_grid[1:I-1,:,:] .* (r_b_grid[1:I-1,:,:] .+ ddeath*pam - aggZ) .+ trans_grid[1:I-1,:,:] - cF[1:I-1,:,:]
    elseif permanent == 0
        sF[1:I-1,:,:] 	= ((1-xxi)-tau_I) * w * l_grid[1:I-1,:,:] .* y_grid[1:I-1,:,:] .+ b_grid[1:I-1,:,:] .* (r_b_grid[1:I-1,:,:] .+ ddeath*pam) .+ trans_grid[1:I-1,:,:] - cF[1:I-1,:,:]
    end
    sF[I,:,:] 		= zeros(1,J,N)
    HcF[1:I-1,:,:] 	= u_fn(cF[1:I-1,:,:], ggamma) .+ VbF[1:I-1,:,:] .* sF[1:I-1,:,:]
    HcF[I,:,:] 		= -1e12*ones(1,J,N)
    validF 			= (sF .> 0)

    cB[2:I,:,:] 	= VbB[2:I,:,:].^(-1/ggamma)
    if permanent == 1
        cB[1,:,:] 		= ((1-xxi)-tau_I) * w * l_grid[1,:,:] .* y_grid[1,:,:] .+ b_grid[1,:,:] .* (r_b_grid[1,:,:] .+ ddeath*pam - aggZ) .+ trans_grid[1,:,:]
        sB[2:I,:,:] 	= ((1-xxi)-tau_I) * w * l_grid[2:I,:,:] .* y_grid[2:I,:,:] .+ b_grid[2:I,:,:] .* (r_b_grid[2:I,:,:] .+ ddeath*pam - aggZ) .+ trans_grid[2:I,:,:] - cB[2:I,:,:]
    elseif permanent == 0
        cB[1,:,:] 	= ((1-xxi)-tau_I) * w * l_grid[1,:,:] .* y_grid[1,:,:] .+
            b_grid[1,:,:] .* (r_b_grid[1,:,:] .+ ddeath*pam) .+ trans_grid[1,:,:]
        sB[2:I,:,:] = ((1-xxi)-tau_I) * w * l_grid[2:I,:,:] .*
            y_grid[2:I,:,:] .+ b_grid[2:I,:,:] .* (r_b_grid[2:I,:,:] .+ ddeath*pam) .+
            trans_grid[2:I,:,:] - cB[2:I,:,:]
    end
    sB[1,:,:] 		= zeros(1,J,N)
    HcB[:,:,:] 		= u_fn(cB,ggamma) .+ VbB .* sB
    validB 			= (sB .< 0)

if permanent == 1
    c0[:,:,:] 		= ((1-xxi)-tau_I) * w * l_grid[:,:,:] .* y_grid[:,:,:] .+ b_grid[:,:,:] .* (r_b_grid[:,:,:] .+ ddeath*pam - aggZ) .+ trans_grid[:,:,:]
elseif permanent == 0
    c0[:,:,:] 		= ((1-xxi)-tau_I) * w * l_grid[:,:,:] .* y_grid[:,:,:] .+ b_grid[:,:,:] .* (r_b_grid[:,:,:] .+ ddeath*pam) .+ trans_grid[:,:,:]
end
s0[:,:,:] 		= zeros(I,J,N)
Hc0[:,:,:] 		= u_fn(c0,ggamma)

# Which direction to use
IcF				= IcF_SS
IcB				= IcB_SS
Ic0				= Ic0_SS

# Optimal consumption and liquid savings
c 				= IcF .* cF .+ IcB .* cB .+ Ic0 .* c0
s 				= IcF .* sF .+ IcB .* sB .+ Ic0 .* s0
u 				= u_fn(c,ggamma)

###
# Deposit decision
###

# Preparations
dFB 			= 0*V
HdFB 			= 0*V
dBF 			= 0*V
HdBF 			= 0*V
dBB 			= 0*V
HdBB 			= 0*V

# Decisions conditional on a particular direction of derivative
dFB[2:I,1:J-1,:]	= opt_deposits(VaF[2:I,1:J-1,:], VbB[2:I,1:J-1,:], a_grid[2:I,1:J-1,:], chi0, chi1, chi2, a_lb)
dFB[:,J,:]			= zeros(I,1,N)
dFB[1,1:J-1,:] 		= zeros(1,J-1,N)
HdFB[2:I,1:J-1,:] 	= VaF[2:I,1:J-1,:] .* dFB[2:I,1:J-1,:] - VbB[2:I,1:J-1,:] .* (dFB[2:I,1:J-1,:] .+ adj_cost_fn(dFB[2:I,1:J-1,:],a_grid[2:I,1:J-1,:], chi0, chi1, chi2, a_lb))
HdFB[:,J,:] 		= -1.0e12 * ones(I,1,N)
HdFB[1,1:J-1,:] 	= -1.0e12 * ones(1,J-1,N)
validFB 			= (dFB .> 0) .* (HdFB .> 0)

dBF[1:I-1,2:J,:] 	= opt_deposits(VaB[1:I-1,2:J,:],VbF[1:I-1,2:J,:],a_grid[1:I-1,2:J,:], chi0, chi1, chi2, a_lb)
dBF[:,1,:] 			= zeros(I,1,N)
dBF[I,2:J,:]		= zeros(1,J-1,N)
HdBF[1:I-1,2:J,:] 	= VaB[1:I-1,2:J,:] .* dBF[1:I-1,2:J,:] - VbF[1:I-1,2:J,:] .* (dBF[1:I-1,2:J,:] .+ adj_cost_fn(dBF[1:I-1,2:J,:],a_grid[1:I-1,2:J,:], chi0, chi1, chi2, a_lb))
HdBF[:,1,:] 		= -1.0e12 * ones(I,1,N)
HdBF[I,2:J,:] 		= -1.0e12 * ones(1,J-1,N)
validBF 			= (dBF .<= -adj_cost_fn(dBF,a_grid, chi0, chi1, chi2, a_lb)) .* (HdBF .> 0)

VbB[1,2:J,:] 		= u_fn(cB[1,2:J,:],ggamma)
dBB[:,2:J,:]		= opt_deposits(VaB[:,2:J,:],VbB[:,2:J,:],a_grid[:,2:J,:], chi0, chi1, chi2, a_lb)
dBB[:,1,:] 			= zeros(I,1,N)
HdBB[:,2:J,:] 		= VaB[:,2:J,:] .* dBB[:,2:J,:] - VbB[:,2:J,:] .* (dBB[:,2:J,:] .+ adj_cost_fn(dBB[:,2:J,:],a_grid[:,2:J,:], chi0, chi1, chi2, a_lb))
HdBB[:,1,:] 		= -1.0e12 * ones(I,1,N)
validBB 			= (dBB .> -adj_cost_fn(dBB,a_grid, chi0, chi1, chi2, a_lb)) .*
    (dBB .<= 0) .* (HdBB .> 0)

# Which direction to use
IcFB				= IcFB_SS
IcBF				= IcBF_SS
IcBB				= IcBB_SS
Ic00				= Ic00_SS

# Optimal deposit decision
d = IcFB .* dFB .+ IcBF .* dBF .+ IcBB .* dBB .+ Ic00 .* zeros(I,J,N)

## Interpolate
d_g = reshape(interp_decision*d[:],I_g,J_g,N)
s_g = reshape(interp_decision*s[:],I_g,J_g,N)
c_g = reshape(interp_decision*c[:],I_g,J_g,N)

###
# Construct transition matrix
###

tmp = reshape(gg,I_g,J_g,N)
audriftB = 0*tmp
audriftF = 0*tmp
budriftB = 0*tmp
budriftF = 0*tmp

chi = 0*V
yy = 0*V
zeta = 0*V

X = 0*V
Y = 0*V
Z = 0*V

chiu = 0*tmp
yyu = 0*tmp
zetau = 0*tmp

Xu = 0*tmp
Yu = 0*tmp
Zu = 0*tmp

# Compute drifts for HJB
if permanent == 1
    adriftB = min.(d,0) .+ min.(a_grid .* (r_a_grid .+ ddeath*pam - aggZ) .+
                                xxi * w * l_grid .* y_grid,0)
    adriftF = max.(d,0) .+ max.(a_grid .* (r_a_grid .+ ddeath*pam - aggZ) .+
                                xxi * w * l_grid .* y_grid,0)
elseif permanent == 0
    adriftB = min.(d,0) .+ min.(a_grid .* (r_a_grid .+ ddeath*pam) .+
                                xxi * w * l_grid .* y_grid,0)
    adriftF = max.(d,0) .+ max.(a_grid .* (r_a_grid .+ ddeath*pam) .+
                                xxi * w * l_grid .* y_grid,0)
end

bdriftB = min.(-d - adj_cost_fn(d,a_grid, chi0, chi1, chi2, a_lb),0) .+ min.(s,0)
bdriftF = max.(-d - adj_cost_fn(d,a_grid, chi0, chi1, chi2, a_lb),0) .+ max.(s,0)

# Compute drifts for KFE
if permanent == 0
    audriftB[1:I_g-1,:,:] = min.(d_g[1:I_g-1,:,:] .+ a_g_grid[1:I_g-1,:,:] .*
                                 (r_a_g_grid[1:I_g-1,:,:] .+ ddeath*pam) .+ xxi *
                                 w * l_g_grid[1:I_g-1,:,:] .* y_g_grid[1:I_g-1,:,:],0)
    audriftB[I_g,:,:] = min.(d_g[I_g,:,:] .+ a_g_grid[I_g,:,:] .*
                             (r_a_g_grid[I_g,:,:] .+ ddeath*pam) .+ xxi * w *
                             l_g_grid[I_g,:,:] .* y_g_grid[I_g,:,:],0)
    audriftF[1:I_g-1,:,:] = max.(d_g[1:I_g-1,:,:] .+ a_g_grid[1:I_g-1,:,:] .*
                                 (r_a_g_grid[1:I_g-1,:,:] .+ ddeath*pam) .+
                                 xxi * w * l_g_grid[1:I_g-1,:,:] .*
                                 y_g_grid[1:I_g-1,:,:],0)
    audriftF[I_g,:,:] = max.(d_g[I_g,:,:] .+ a_g_grid[I_g,:,:] .*
                             (r_a_g_grid[I_g,:,:] .+ ddeath*pam) .+
                             xxi * w * l_g_grid[I_g,:,:] .* y_g_grid[I_g,:,:],0)

    budriftB[1:I_g-1,:,:] = min.(s_g[1:I_g-1,:,:] - d_g[1:I_g-1,:,:] -
                                 adj_cost_fn(d_g[1:I_g-1,:,:], a_g_grid[1:I_g-1,:,:],
                                             chi0, chi1, chi2, a_lb),0)
    budriftB[I_g,:,:] = min.(s_g[I_g,:,:] - d_g[I_g,:,:] -
                             adj_cost_fn(d_g[I_g,:,:],a_g_grid[I_g,:,:],
                                         chi0, chi1, chi2, a_lb),0)
    budriftF[1:I_g-1,:,:] = max.(s_g[1:I_g-1,:,:] - d_g[1:I_g-1,:,:] -
                                 adj_cost_fn(d_g[1:I_g-1,:,:], a_g_grid[1:I_g-1,:,:],
                                             chi0, chi1, chi2, a_lb),0)
    budriftF[I_g,:,:] = max.(s_g[I_g,:,:] - d_g[I_g,:,:] -
                             adj_cost_fn(d_g[I_g,:,:],a_g_grid[I_g,:,:],
                                         chi0, chi1, chi2, a_lb),0)
elseif permanent == 1
    audriftB[1:I_g-1,:,:] = min.(d_g[1:I_g-1,:,:] .+ a_g_grid[1:I_g-1,:,:] .*
                                 (r_a_g_grid[1:I_g-1,:,:] .+ ddeath*pam) .+ xxi *
                                 w * l_g_grid[1:I_g-1,:,:] .* y_g_grid[1:I_g-1,:,:]
                                 - aggZ * a_g_grid[1:I_g-1,:,:], 0)
    audriftB[I_g,:,:] = min.(d_g[I_g,:,:] .+ a_g_grid[I_g,:,:] .* (r_a_g_grid[I_g,:,:] .+ ddeath*pam) .+ xxi * w * l_g_grid[I_g,:,:] .* y_g_grid[I_g,:,:] - aggZ * a_g_grid[I_g,:,:],0)
    audriftF[1:I_g-1,:,:] = max.(d_g[1:I_g-1,:,:] .+ a_g_grid[1:I_g-1,:,:] .* (r_a_g_grid[1:I_g-1,:,:] .+ ddeath*pam) .+ xxi * w * l_g_grid[1:I_g-1,:,:] .* y_g_grid[1:I_g-1,:,:] - aggZ * a_g_grid[1:I_g-1,:,:],0)
    audriftF[I_g,:,:] = max.(d_g[I_g,:,:] .+ a_g_grid[I_g,:,:] .* (r_a_g_grid[I_g,:,:] .+ ddeath*pam) .+ xxi * w * l_g_grid[I_g,:,:] .* y_g_grid[I_g,:,:] - aggZ * a_g_grid[I_g,:,:],0)

    budriftB[1:I_g-1,:,:] = min.(s_g[1:I_g-1,:,:] - d_g[1:I_g-1,:,:] - adj_cost_fn(d_g[1:I_g-1,:,:],a_g_grid[1:I_g-1,:,:], chi0, chi1, chi2, a_lb) - aggZ * b_g_grid[1:I_g-1,:,:],0)
    budriftB[I_g,:,:] = min.(s_g[I_g,:,:] - d_g[I_g,:,:] - adj_cost_fn(d_g[I_g,:,:],a_g_grid[I_g,:,:], chi0, chi1, chi2, a_lb) - aggZ * b_g_grid[I_g,:,:],0)
    budriftF[1:I_g-1,:,:] = max.(s_g[1:I_g-1,:,:] - d_g[1:I_g-1,:,:] - adj_cost_fn(d_g[1:I_g-1,:,:],a_g_grid[1:I_g-1,:,:], chi0, chi1, chi2, a_lb) - aggZ * b_g_grid[1:I_g-1,:,:],0)
    budriftF[I_g,:,:] = max.(s_g[I_g,:,:] - d_g[I_g,:,:] -
                             adj_cost_fn(d_g[I_g,:,:],a_g_grid[I_g,:,:],
                                         chi0, chi1, chi2, a_lb) - aggZ *
                             b_g_grid[I_g,:,:],0)
end

# Derive transition matrices
include("Subroutines/transition_a_deriva.jl")
include("Subroutines/transition_b_deriva.jl")
include("Subroutines/transition_c_deriva.jl")

# full transition matrix
A = aa + bb + cc

#----------------------------------------------------------------
# KFE
#----------------------------------------------------------------

# KFE
AT = aau + bbu + ccu
AT = AT'
gg_tilde            = dab_g_tilde_mat * gg
gIntermediate       = AT * gg_tilde
dab_g_tilde_mat_inv = spdiagm(0 => vec(repeat(1.0 ./ dab_g_tilde,N,1)))
gIntermediate       = dab_g_tilde_mat_inv * gIntermediate

dab_g_small = reshape(dab_g[:,:,1],I_g*J_g,1)

loc = findall(b .== 0)
dab_g_small       = dab_g_small ./ dab_g_small[loc]*ddeath
dab_g_small[loc] .= 0.0
death_process         = -ddeath * my_speye(I_g*J_g)
death_process[loc,:]  = vec(dab_g_small)
death_process         = kron(my_speye(N), death_process)

gIntermediate = gIntermediate + death_process * gg

# find death-corrected savings

a_g_grid_aux = reshape(a_g_grid,I_g*J_g*N,1)
a_save = sum(gIntermediate .* dab_g_aux .* a_g_grid_aux)

b_g_grid_aux = reshape(b_g_grid,I_g*J_g*N,1)
b_save = sum(gIntermediate .* dab_g_aux .* b_g_grid_aux)

# consumption for low types and high types

c_low     = c[:,:,1]
dab_low   = dab[:,:,1]
c_high    = c[:,:,2]
dab_high  = dab[:,:,2]

#----------------------------------------------------------------
# Compute equilibrium conditions
#----------------------------------------------------------------

# HJB equation
if permanent == 0
    hjbResidual = reshape(u,I*J*N,1) + A * reshape(V,I*J*N,1) + V_Dot + VEErrors - (rrho + ddeath) * reshape(V,I*J*N,1)
elseif permanent == 1
    hjbResidual = reshape(u,I*J*N,1) + A * reshape(V,I*J*N,1) + V_Dot + VEErrors - (rrho + ddeath - (1 - ggamma) * aggZ) * reshape(V,I*J*N,1)
end

# KFE
gResidual 			= g_Dot - gIntermediate[1:n_g,1]

# Aggregate conditions
if permanent == 0

	if K_liquid == 1
        K_out           = sum((a_g_grid[:] .+ b_g_grid[:]) .* gg .* dab_g[:])
    else
        K_out           = sum(a_g_grid[:] .* gg .* dab_g[:])
    end
    if r_b_fix == 1
        r_b_out         = r_b_SS
    elseif r_b_phi == 1
        r_b_out			= sum(b_g_grid[:] .* gg .* dab_g[:])
    elseif B_fix == 1
        r_b_out         = 0
	elseif K_liquid == 1
		r_b_out			= r_a_out - illiquid_wedge
    end

    if aggregate_variables == 1
        aggY_out			= (K ^ aalpha) * (exp(aggZ) * n_SS) ^ (1 - aalpha)
        aggC_out			= sum(c_g[:] .* gg .* dab_g[:])
    elseif distributional_variables == 1
        C_Var_out           = sum(log(c_g[:]).^2 .* gg .* dab_g[:]) - sum(log(c_g[:]) .* gg .* dab_g[:])^2
        earn                = log((1-tau_I) * w * l_g_grid .* y_g_grid + b_g_grid .* (r_b_g_grid + ddeath*pam) + trans_grid + a_g_grid .* (r_a_g_grid + ddeath*pam))
        earn_Var_out        = sum(earn[:].^2 .* gg .* dab_g[:]) - sum(earn[:] .* gg .* dab_g[:])^2
    elseif distributional_variables_1 == 1
        WHTM_indicator      = zeros(I_g,J_g,N)
        WHTM_indicator[b_g_0pos:b_g_0pos+1,a_g_0pos+2:end,:] = 1
        WHTM_out            = sum(WHTM_indicator[:] .* gg .* dab_g[:])
        #C_WHTM_out          = sum(WHTM_indicator[:] .* c_g[:] .* gg .* dab_g[:])/WHTM_out
        C_WHTM_out          = sum(WHTM_indicator[:] .* c_g[:] .* gg .* dab_g[:])

        PHTM_indicator      = zeros(I_g,J_g,N)
        PHTM_indicator[b_g_0pos:b_g_0pos+1,a_g_0pos:a_g_0pos+2:end,:] = 1
        PHTM_out            = sum(PHTM_indicator[:] .* gg .* dab_g[:])
        #C_PHTM_out          = sum(PHTM_indicator[:] .* c_g[:] .* gg .* dab_g[:])/PHTM_out
        C_PHTM_out          = sum(PHTM_indicator[:] .* c_g[:] .* gg .* dab_g[:])

        NHTM_indicator      = zeros(I_g,J_g,N)
        NHTM_indicator[b_g_0pos+2:end,a_g_0pos+2:end,:] = 1
        #NHTM_indicator      = zeros(I_g,J_g,N) NHTM_indicator[b_g_0pos+2:end,:,:] = 1
        NHTM_indicator      = zeros(I_g,J_g,N)
        NHTM_indicator[b_g_0pos+3:end,2:end,:] = 1
        #NHTM_indicator(3:b_g_0pos-1,:,:) = 1
        NHTM_out            = sum(NHTM_indicator[:] .* gg .* dab_g[:])
        C_NHTM_out          = sum(NHTM_indicator[:] .* c_g[:] .* gg .* dab_g[:])/NHTM_out
    end

elseif permanent == 1

	if K_liquid == 1
	    K_out               = sum((a_g_grid[:] + b_g_grid[:]) .* gg .* dab_g[:])
	elseif K_liquid == 0
		K_out               = sum(a_g_grid[:] .* gg .* dab_g[:])
	end
    if r_b_fix == 1
        r_b_out         = r_b_SS
    elseif r_b_phi == 1
        r_b_out			= sum(b_g_grid[:] .* gg .* dab_g[:])
    elseif B_fix == 1
        r_b_out         = 0
    elseif K_liquid == 1
        r_b_out         = r_a_out - illiquid_wedge
    end

    if aggregate_variables == 1
        aggY_out			= (K ^ aalpha) * (n_SS ^ (1 - aalpha))
        aggC_out			= sum(c_g[:] .* gg .* dab_g[:])
    elseif distributional_variables == 1
        C_Var_out           = sum(log(c_g[:]).^2 .* gg .* dab_g[:]) - sum(log(c_g[:]) .* gg .* dab_g[:])^2
        earn                = log((1-tau_I) * w * l_g_grid .* y_g_grid + b_g_grid .* (r_b_g_grid + ddeath*pam) + trans_grid + a_g_grid .* (r_a_g_grid + ddeath*pam))
        earn_Var_out        = sum(earn[:].^2 .* gg .* dab_g[:]) - sum(earn[:] .* gg .* dab_g[:])^2
    elseif distributional_variables_1 == 1
        WHTM_indicator      = zeros(I_g,J_g,N)
        WHTM_indicator[b_g_0pos:b_g_0pos+1,a_g_0pos+2:end,:] = 1
        WHTM_out            = sum(WHTM_indicator[:] .* gg .* dab_g[:])
        #C_WHTM_out          = sum(WHTM_indicator[:] .* c_g[:] .* gg .* dab_g[:])/WHTM_out
        C_WHTM_out          = sum(WHTM_indicator[:] .* c_g[:] .* gg .* dab_g[:])

        PHTM_indicator      = zeros(I_g,J_g,N)
        PHTM_indicator[b_g_0pos:b_g_0pos+1,a_g_0pos:a_g_0pos+2:end,:] = 1
        PHTM_out            = sum(PHTM_indicator[:] .* gg .* dab_g[:])
        #C_PHTM_out          = sum(PHTM_indicator(:] .* c_g(:] .* gg .* dab_g(:])/PHTM_out
        C_PHTM_out          = sum(PHTM_indicator[:] .* c_g[:] .* gg .* dab_g[:])

        NHTM_indicator      = zeros(I_g,J_g,N)
        NHTM_indicator[b_g_0pos+2:end,:,:] = 1
        #NHTM_indicator(3:b_g_0pos-1,:,:] = 1
        NHTM_out            = sum(NHTM_indicator[:] .* gg .* dab_g[:])
        C_NHTM_out          = sum(NHTM_indicator[:] .* c_g[:] .* gg .* dab_g[:])/NHTM_out
    end
end

K_Residual			= K_out - K
if r_b_fix == 1
    r_b_Residual    = r_b_out - r_b
elseif r_b_phi == 1
    r_b_Residual    = r_b_out - B_SS * exp(1/pphi * (r_b - r_b_SS))
elseif B_fix == 1
    r_b_Residual    = r_b_out - b_save
elseif K_liquid == 1
	r_b_Residual	= r_b_out - r_b
end

if aggregate_variables == 1
    Y_Residual			= aggY_out - aggY
    C_Residual			= aggC_out - aggC
elseif distributional_variables == 1
    C_Var_Residual      = C_Var_out - C_Var
    earn_Var_Residual   = earn_Var_out - earn_Var
elseif distributional_variables_1 == 1
    C_WHTM_Residual     = C_WHTM_out - C_WHTM
	C_PHTM_Residual     = C_PHTM_out - C_PHTM
    #C_NHTM_Residual     = C_NHTM_out - C_NHTM
end

# Law of motion for aggregate tfp shock
aggZ_Residual 	= aggZ_Dot - (-nnu_aggZ * aggZ + ssigma_aggZ * aggZ_Shock)

    #  Package results
    if aggregate_variables == 1
        vResidual = [hjbResidual; gResidual; K_Residual; r_b_Residual; Y_Residual; C_Residual; aggZ_Residual]
    elseif distributional_variables == 1
        vResidual = [hjbResidual; gResidual; K_Residual; r_b_Residual; C_Var_Residual; earn_Var_Residual; aggZ_Residual]
    elseif distributional_variables_1 == 1
        vResidual = [hjbResidual; gResidual; K_Residual; r_b_Residual; C_WHTM_Residual; C_PHTM_Residual; aggZ_Residual]
    else
        vResidual = [hjbResidual; gResidual; K_Residual; r_b_Residual; aggZ_Residual]
    end
    return vResidual
end
