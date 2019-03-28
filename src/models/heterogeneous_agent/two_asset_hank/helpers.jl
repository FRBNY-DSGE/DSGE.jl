"""
THE GUTS OF EQCOND.

"""
#=
@inline function eqcond_helper(V::Array{T,3}, I::S, J::S, I_g::S, J_g::S,
                               N::S, chi0::T, chi1::T, chi2::T, a_lb::T,
                               ggamma::T, permanent::Bool,
                               interp_decision::SparseMatrixCSC{T, S},
                               ddeath::T, pam::T, aggZ::T, xxi::T,
                               tau_I::T, w::T, trans::T,
                               y_grid::Array{Complex{T},3}, b_grid::Array{T,3},
                               r_b_grid::Array{T,3}, alb_grid::Array{T,3},
                               daf_grid::Array{T,3}, dab_grid::Array{T,3},
                               dbf_grid::Array{T,3}, dbb_grid::Array{T,3},
                               IcF::BitArray{3}, IcB::BitArray{3}, Ic0::Array{S,3},
                               IcFB::BitArray{3}, IcBF::BitArray{3}, IcBB::BitArray{3},
                               Ic00::BitArray{3}) where {S<:Int64, T<:Float64}
=#
@inline function eqcond_helper(V, I_g, J_g, chi0, chi1, chi2,
                               a_lb, ggamma, permanent, interp_decision,
                               ddeath, pam, aggZ, xxi, tau_I, w, trans,
                               r_b_vec, alb_vec,
                               daf_grid, dab_grid, dbf_grid, dbb_grid,
                               IcF, IcB, Ic0, IcFB, IcBF, IcBB, Ic00,
                               y, b)
    #----------------------------------------------------------------
    # HJB Equation
    #----------------------------------------------------------------
    I, J, N = size(V)

    # ---- Liquid ----
    VaF   = similar(V)
    VaB   = similar(V)
    Vamin = 0.0
    Va_dif = V[:,2:J,:] - V[:,1:J-1,:]

    # Forward difference
    VaF[:,1:J-1,:] = max.((Va_dif) ./ daf_grid[:,1:J-1,:], Vamin)
    VaF[:,J,:]    .= 0.0

    # Backward difference
    VaB[:,2:J,:] = max.((Va_dif) ./ dab_grid[:,2:J,:], Vamin)
    VaB[:,1,:]  .= 0.0

    # ---- Illiquid ----
    VbF   = similar(V)
    VbB   = similar(V)
    Vbmin = 1e-8
    Vb_dif = V[2:I,:,:] - V[1:I-1,:,:]

    # Forward difference
    VbF[1:I-1,:,:] = max.((Vb_dif) ./ dbf_grid[1:I-1,:,:], Vbmin)
    VbF[I,:,:]    .= 0.0

    # Backward difference
    VbB[2:I,:,:] = max.((Vb_dif) ./ dbb_grid[2:I,:,:], Vbmin)
    VbB[1,:,:]  .= 0.0

    #----------------------------------------------------------------
    # Consumption decision
    #----------------------------------------------------------------
    perm_const = (permanent == 1) ? ddeath * pam - aggZ : ddeath * pam
    c0_c = ((1-xxi) - tau_I) * w

    # Optimal consumption and liquid savings
    c, s = similar(V), similar(V)
    for i=1:I, j=1:J, n=1:N
        c0 = c0_c * real(y[1,n]) + b[i] * (r_b_vec[i] + perm_const) + trans

        # Decisions conditional on a particular direction of derivative
        cF = (i==I) ? 0.0 : VbF[i,j,n] ^ (-1 / ggamma)
        cB = (i==1) ? c0  : VbB[i,j,n] ^ (-1 / ggamma)

        sF = (i==I) ? 0.0 : c0 - cF
        sB = (i==1) ? 0.0 : c0 - cB

        c[i,j,n] = IcF[i,j,n] * cF + IcB[i,j,n] * cB + Ic0[i,j,n] * c0
        s[i,j,n] = IcF[i,j,n] * sF + IcB[i,j,n] * sB
    end
    u = u_fn.(c, ggamma)

    #----------------------------------------------------------------
    # Deposit decision
    #----------------------------------------------------------------
    d = similar(V)
    for i=1:I, j=1:J, n=1:N
        dFB = (i == 1 || j == J) ? 0.0 :
            opt_deposits(VaF[i,j,n], VbB[i,j,n], alb_vec[j], chi0, chi1, chi2)

        dBF = (i == I || j == 1) ? 0.0 :
            opt_deposits(VaB[i,j,n], VbF[i,j,n], alb_vec[j], chi0, chi1, chi2)

        dBB = (j == 1) ? 0.0 :
            opt_deposits(VaB[i,j,n], VbB[i,j,n], alb_vec[j], chi0, chi1, chi2)

        d[i,j,n] = IcFB[i,j,n] * dFB + IcBF[i,j,n] * dBF + dBB * IcBB[i,j,n]
    end

    # Interpolate
    d_g = reshape(interp_decision * vec(d), I_g, J_g, N)
    s_g = reshape(interp_decision * vec(s), I_g, J_g, N)
    c_g = reshape(interp_decision * vec(c), I_g, J_g, N)

    return c, s, u, d, d_g, s_g, c_g
end

"""
```
@inline function create_y_grid(y_size::Int64, ygrid_new::Int64)
```
Function initializes income grid.
"""
@inline function create_y_grid(y_size::Int64, ygrid_new::Int64, dataroot::String)

    y      = Array{Float64}(undef, 0)
    y_dist = Array{Float64}(undef, 0)
    λ      = Array{Float64}(undef, 0)

    if y_size == 2
        y  = [0.1, 1]

        # Transition probabilities
        λ1 = 0.5 # expected duration of unemployment is 2 quarters
        λ2 = (λ1 / (y[2] * .93 - y[1])) * (y[2] - y[2] * .93) # unemployment rate 7
        λ  = [-λ1 λ1; λ2 -λ2]

    elseif y_size == 30
        # Load grid and transition probabilities from external calibration
        y      = load(dataroot * "income_grid_30.jld2", "y")
        y_dist = load(dataroot * "income_grid_30.jld2", "y_dist")
        λ      = load(dataroot * "income_transition_30.jld2", "lambda")

        # Create grid of income shocks
        y = exp.(y)

    elseif y_size == 33
        # Load grid and transition probabilities from external calibration
        y      = load(dataroot * "income_grid.jld2", "y")
        y_dist = load(dataroot * "income_grid.jld2", "y_dist")
        λ      = load(dataroot * "income_transition.jld2", "lambda")

        # Create grid of income shocks
        y = exp.(y)
    end

    # Compute stationary income distribution
    y_dist = stat_dist(λ')   # stationary distribution
    y_mean = y * y_dist

    if     ygrid_new == 0
        y = y ./ y_mean
    elseif ygrid_new == 1
        y = 0.1907858 .* y ./ y_mean
    end

    y_mean = dot(y, y_dist)
    N = size(y, 2)
    return λ, y, y_mean, y_dist, N
end

"""
```
@inline function create_a_grid(agrid_new::Int64, J::Int64, J_g::Int64, amin::Int64, amax::Int64)
```
Function initializes illiquid asset grids.
"""
@inline function create_a_grid(agrid_new::Int64, J::Int64, J_g::Int64, amin::Int64, amax::Int64)
    a   = Array{Float64}(undef, 0)
    a_g = Array{Float64}(undef, 0)

    if agrid_new == 0

        a 		    = range(0, stop = 1, length = J) # set the grid
        coeff_power = 0.9
        power       = 8
        a 			= amax*((1 - coeff_power) * a + coeff_power * (a.^power))

        a_g         = range(0, stop = 1, length = J_g)
        amax_g      = 150
        coeff_power = 0.9
        power       = 8
        a_g         = amax_g*((1-coeff_power) * a_g + coeff_power * (a_g.^power))

    elseif agrid_new == 1

        agridparam = 0.15
        a_raw = range(0, stop = 1, length = J)
        a_raw = a_raw.^(1/agridparam)

        a = amin .+ (amax - amin) * a_raw

        for i = 1:9
            a[i] = (i-1) * a[10]/(10 - 1)
        end

        ## Define new grid points for g here. For now, it is created the same.
        agridparam = 0.15

        a_raw = range(0, stop = 1, length = J_g)
        a_raw = a_raw.^(1/agridparam)
        a_g = amin .+ (amax - amin) * a_raw

        for i = 1:9
            a_g[i] = (i-1) * a_g[10]/(10 - 1)
        end

    elseif agrid_new == 2

        agridparam = 0.15

        a_raw = range(0, stop = 1, length = J)
        a_raw = a_raw.^(1/agridparam)
        a = amin .+ (amax - amin) * a_raw

        for i = 1:9
            a[i] = (i-1) * a[10]/(10 - 1)
        end

        ## Define new grid points for g here. For now, it is created the same.
        agridparam = 0.21

        a_raw = range(0, stop = 1, length = J_g)
        a_raw = a_raw.^(1/agridparam)
        a_g = amin .+ (amax - amin) * a_raw

        for i = 1:9
            a_g[i] = (i-1) * a_g[10]/(10. - 1.)
        end

    end

    a_g_0pos = findall(!iszero, a_g .== 0)

    return a, a_g, a_g_0pos
end

"""
```
@inline function create_b_grid(bgrid_new::Int64, I::Int64)
```
Function to initialize liquid asset grid.
"""
@inline function create_b_grid(bgrid_new::Int64, I::Int64, I_g::Int64)
    I = 50 # TEMP
    I_neg = 10
    I_pos = 40

    bmin = Array{Float64}(undef, 0)
    bmax = Array{Float64}(undef, 0)
    b    = Array{Float64}(undef, 0)
    b_g  = Array{Float64}(undef, 0)

    if bgrid_new == 0

        # Lower bound
        b_    = -2		# lower bound on liquid assets

        # Preparations
        I_neg = Int(I / 10) # number of grid points on the negative part of the real line
        bmin  = b_ 		# lower bound
        bmax  = 80

        # Grid for positive part
        b_pos 			= range(0, stop = 1, length = Int(I-I_neg))
        coeff_power 	= 0.9
        power = 8
        b_pos 			= bmax*((1 - coeff_power) * b_pos + coeff_power * (b_pos.^power))

        # Grid for negative part
        b_neg 			= range(0, stop = 1, length = Int(I_neg+1))
        coeff_power 	= 0
        power           = 1
        b_neg 			= -bmin*((1 - coeff_power) * b_neg + coeff_power * (b_neg.^power)) .+ bmin
        b_neg           = b_neg[1:end-1]

        # Put together grid on liquid assets
        b = [b_neg;b_pos]

        ## Define new grid points for g here. For now, it is created the same.
        I_g_neg = 5
        bmax_g  = 70

        # Grid for positive part
        b_pos 			= range(0, stop = 1, length = Int(I_g-I_g_neg))
        coeff_power 	= 0.9
        power = 8
        b_pos 			= bmax_g*((1 - coeff_power) * b_pos + coeff_power * (b_pos.^power))

        # Grid for negative part
        b_neg 			= range(0, stop = 1, length = Int(I_g_neg+1))
        coeff_power 	= 0
        power = 1
        b_neg = -bmin*((1 - coeff_power) * b_neg + coeff_power * (b_neg.^power)) .+ bmin
        b_neg = b_neg[1:end-1]

        # Put together grid on liquid assets
        b_g = [b_neg;b_pos]

    elseif bgrid_new == 1

        # set parameters
        I_neg = I - I_pos
        bmin  = -1
        b0    = 0
        bmax  = 40
        bgridparam = 0.35
        bgridparam_neg = 0.4
        b_ = bmin

        # positive part
        bp_raw = range(0, stop = 1, length = I_pos)
        bp_raw = bp_raw .^ (1 / bgridparam)
        bp     = b0 .+ (bmax - b0) * bp_raw

        # negative part
        bn_raw = range(0, stop = 1, length = Int64(I_neg/2+1))
        bn_raw = bn_raw.^(1/bgridparam_neg)
        bn = bmin .+ (bmin/2 - bmin) * bn_raw
        bn = [bn; Vector{Float64}(undef, Int64(I_neg/2-1))]
        for i = Int64(I_neg/2+2) : I_neg
            bn[i] = b0 - (bn[Int64(I_neg)+2-i] - bn[1])
        end

        # Put everything together
        b = [bn; bp]

        ## Define new grid points for g here. For now, it is created the same.
        # set parameters
        I_neg = I_g - I_pos
        bmin  = -1
        b0    = 0
        bmax  = 40
        bgridparam = 0.35
        bgridparam_neg = 0.4
        b_ = bmin

        # positive part
        bp_raw = range(0, stop = 1, length = I_pos)
        bp_raw = bp_raw.^(1/bgridparam)
        bp = b0 .+ (bmax - b0) * bp_raw

        # negative part
        bn_raw = range(0, stop = 1, length = Int64(I_neg/2+1))
        bn_raw = bn_raw.^(1/bgridparam_neg)
        bn = bmin .+ (bmin/2 - bmin) * bn_raw
        bn = [bn;Vector{Float64}(undef,Int64(I_neg/2-1))]
        for i = Int64(I_neg/2+2):I_neg
            bn[i] = b0 - (bn[Int64(I_neg+2-i)] - bn[1])
        end

        # Put everything together
        b_g = [bn; bp]

    elseif bgrid_new == 2

        # set parameters
        bmin = -1
        b0   = 0
        bmax = 40
        bgridparam = 0.25
        bgridparam_neg = 0.4
        b_ = bmin

        # positive part
        bp_raw = range(0, stop = 1, length = I_pos)
        bp_raw = bp_raw.^(1/bgridparam)
        bp = b0 .+ (bmax - b0) * bp_raw

        # negative part
        bn_raw = range(0, stop = 1, length = Int(I_neg/2+1))
        bn_raw = bn_raw.^(1/bgridparam_neg)
        bn = bmin .+ (bmin/2 - bmin) * bn_raw
        bn = [bn;Vector{Float64}(undef,Int(I_neg/2-1))]

        for i = I_neg/2+2:I_neg
            bn[i] = b0 - (bn[I_neg+2-i] - bn[1])
        end

        # put everything together
        b = [bn;bp]

        ## Define new grid points for g here. For now, it is created the same.
        # set parameters
        I_pos = 40
        I_neg = I_g - I_pos
        bmin = -1
        b0 = 0
        bmax = 40
        bgridparam = 0.6
        bgridparam_neg = 0.4
        b_ = bmin

        # positive part
        bp_raw = range(0, stop = 1, length = I_pos)
        bp_raw = bp_raw.^(1 / bgridparam)
        bp = b0 .+ (bmax - b0) * bp_raw

        # negative part
        bn_raw = range(0, stop = 1, length = Int(I_neg/2+1))
        bn_raw = bn_raw.^(1/bgridparam_neg)
        bn = bmin .+ (bmin/2 - bmin) * bn_raw
        bn = [bn;Vector{Float64}(undef,Int(I_neg/2-1))]
        for i = Int(I_neg/2+2):I_neg
            bn[i] = b0 - (bn[Int(I_neg+2-i)] - bn[1])
        end

        # put everything together
        b_g = [bn;bp]
    end
    b_g_0pos = findall(b_g .== 0)

    return I_neg, bmin, bmax, b, b_g, b_g_0pos
end

@inline function set_grids(a, b, a_g, b_g, y, I, J, I_g, J_g, N, w, r_a, r_b, r_b_borr, trans)
    b_grid     = permutedims(repeat(b , 1, J, N), [1 2 3])
    a_grid     = permutedims(repeat(a , 1, I, N), [2 1 3])
    y_grid     = permutedims(repeat(y', 1, I, J), [2 3 1])
    r_b_grid   = r_b .* (b_grid .>= 0) + r_b_borr .* (b_grid .< 0)
    trans_grid = trans * ones(I,J,N)

    dbf_grid            = Array{Float64}(undef, I,J,N)
    dbb_grid            = Array{Float64}(undef, I,J,N)
    dbf_grid[1:I-1,:,:] = b_grid[2:I,:,:] - b_grid[1:I-1,:,:]
    dbf_grid[I,:,:]     = dbf_grid[I-1,:,:]
    dbb_grid[2:I,:,:]   = b_grid[2:I,:,:] - b_grid[1:I-1,:,:]
    dbb_grid[1,:,:]     = dbb_grid[2,:,:]

    daf_grid            = Array{Float64}(undef, I, J, N)
    dab_grid            = Array{Float64}(undef, I, J, N)
    daf_grid[:,1:J-1,:] = a_grid[:,2:J,:] - a_grid[:,1:J-1,:]
    daf_grid[:,J,:]     = daf_grid[:,J-1,:]
    dab_grid[:,2:J,:]   = a_grid[:,2:J,:] - a_grid[:,1:J-1,:]
    dab_grid[:,1,:]     = dab_grid[:,2,:]

    db_tilde      = 0.5*(dbb_grid[:,1,1] + dbf_grid[:,1,1])
    db_tilde[1]   = 0.5*dbf_grid[1,1,1]
    db_tilde[end] = 0.5*dbb_grid[end,1,1]
    da_tilde      = 0.5*(dab_grid[1,:,1] + daf_grid[1,:,1])
    da_tilde[1]   = 0.5 * daf_grid[1,1,1]
    da_tilde[end] = 0.5*dab_grid[1,end,1]

    dab_tilde      = kron(da_tilde, db_tilde)
    dab_tilde_grid = reshape(repeat(dab_tilde, N, 1), I, J, N)
    dab_tilde_mat  = spdiagm(0 => vec(repeat(dab_tilde, N, 1)))

    r_a_grid = repeat([r_a], I, J, N)
    w_grid	 = repeat([w], I, J, N)
    l_grid   = permutedims(repeat(ones(N,1),1,I,J),[2 3 1])

    b_g_grid     = permutedims(repeat(b_g, 1, J_g, N),  [1 2 3])
    a_g_grid     = permutedims(repeat(a_g, 1, I_g, N),  [2 1 3])
    y_g_grid     = permutedims(repeat(y', 1, I_g, J_g), [2 3 1])
    r_b_g_grid   = r_b .* (b_g_grid .>= 0) + r_b_borr .* (b_g_grid .< 0)
    trans_g_grid = trans * ones(I_g,J_g,N)

    dbf_g_grid = Array{Float64}(undef, I_g, J_g, N)
    dbf_g_grid[1:I_g-1,:,:] = b_g_grid[2:I_g,:,:] - b_g_grid[1:I_g-1,:,:]
    dbf_g_grid[I_g,:,:] = dbf_g_grid[I_g-1,:,:]
    dbb_g_grid = Array{Float64}(undef, I_g,J_g,N)
    dbb_g_grid[2:I_g,:,:] = b_g_grid[2:I_g,:,:] - b_g_grid[1:I_g-1,:,:]
    dbb_g_grid[1,:,:] = dbb_g_grid[2,:,:]

    daf_g_grid = Array{Float64}(undef, I_g, J_g, N)
    daf_g_grid[:,1:J_g-1,:] = a_g_grid[:,2:J_g,:] - a_g_grid[:,1:J_g-1,:]
    daf_g_grid[:,J_g,:] = daf_g_grid[:,J_g-1,:]
    dab_g_grid = Array{Float64}(undef, I_g, J_g, N)
    dab_g_grid[:,2:J_g,:] = a_g_grid[:,2:J_g,:] - a_g_grid[:,1:J_g-1,:]
    dab_g_grid[:,1,:] = dab_g_grid[:,2,:]

    db_g_tilde       = 0.5*(dbb_g_grid[:,1,1] + dbf_g_grid[:,1,1])
    db_g_tilde[1]    = 0.5*dbf_g_grid[1,1,1]
    db_g_tilde[end]  = 0.5*dbb_g_grid[end,1,1]
    da_g_tilde       = 0.5*(dab_g_grid[1,:,1] + daf_g_grid[1,:,1])
    da_g_tilde[1]    = 0.5*daf_g_grid[1,1,1]
    da_g_tilde[end]  = 0.5*dab_g_grid[1,end,1]
    dab_g_tilde      = kron(da_g_tilde, db_g_tilde)
    dab_g_tilde_grid = reshape(repeat(dab_g_tilde,N,1),I_g,J_g,N)
    dab_g_tilde_mat  = spdiagm(0 => vec(repeat(dab_g_tilde,N,1)))

    r_a_g_grid	     = repeat([r_a], I_g, J_g, N)
    w_g_grid	     = repeat([w], I_g, J_g,N)
    l_g_grid		 = permutedims(repeat(ones(N,1),1,I_g,J_g),[2 3 1])

    return a_grid, a_g_grid, b_grid, b_g_grid, y_grid, y_g_grid, r_a_grid, r_b_grid, r_a_g_grid, r_b_g_grid, daf_grid, daf_g_grid, dab_grid, dab_g_grid, dab_tilde_grid, dab_g_tilde_grid, dab_g_tilde_mat, dab_g_tilde, dbf_grid, dbf_g_grid, dbb_grid, dbb_g_grid, trans_grid, trans_g_grid, l_grid, l_g_grid, w_grid
end

@inline function transition(I_g, J_g, N, I, J, ddeath, pam, xxi, w, chi0, chi1, chi2, a_lb,
                            l_grid, l_g_grid, y_grid, y_g_grid, d, dab_grid, daf_grid,
                            dab_g_grid, daf_g_grid, dbb_grid, dbf_grid, dbb_g_grid, dbf_g_grid,
                            d_g, a_grid, a_g_grid, s, s_g, r_a_grid, r_a_g_grid)

    audriftB = Array{Float64,3}(undef, I_g,J_g,N)
    audriftF = Array{Float64,3}(undef, I_g,J_g,N)
    budriftB = Array{Float64,3}(undef, I_g,J_g,N)
    budriftF = Array{Float64,3}(undef, I_g,J_g,N)

    chi  = Array{Float64,3}(undef, I,J,N)
    yy   = Array{Float64,3}(undef, I,J,N)
    zeta = Array{Float64,3}(undef, I,J,N)

    X = Array{Float64,3}(undef, I,J,N)
    Y = Array{Float64,3}(undef, I,J,N)
    Z = Array{Float64,3}(undef, I,J,N)

    chiu  = Array{Float64,3}(undef, I_g,J_g,N)
    yyu   = Array{Float64,3}(undef, I_g,J_g,N)
    zetau = Array{Float64,3}(undef, I_g,J_g,N)

    Xu = Array{Float64,3}(undef, I_g,J_g,N)
    Yu = Array{Float64,3}(undef, I_g,J_g,N)
    Zu = Array{Float64,3}(undef, I_g,J_g,N)

    # compute all drifts
    adriftB = min.(d,0) .+ min.(a_grid .* (r_a_grid .+ ddeath*pam) .+ xxi * w * l_grid .* y_grid, 0)
    adriftF = max.(d,0) .+ max.(a_grid .* (r_a_grid .+ ddeath*pam) .+ xxi * w * l_grid .* y_grid, 0)

    bdriftB = min.(-d .- adj_cost_fn(d, a_grid, chi0, chi1, chi2, a_lb),0) .+ min.(s,0)
    bdriftF = max.(-d .- adj_cost_fn(d, a_grid, chi0, chi1, chi2, a_lb),0) .+ max.(s,0)

    audriftB[1:I_g-1,:,:] = min.(d_g[1:I_g-1,:,:] .+ a_g_grid[1:I_g-1,:,:] .*
                                (r_a_g_grid[1:I_g-1,:,:] .+ ddeath*pam) .+ xxi * w *
                                l_g_grid[1:I_g-1,:,:] .* y_g_grid[1:I_g-1,:,:], 0)
    audriftB[I_g,:,:]     = min.(d_g[I_g,:,:] .+ a_g_grid[I_g,:,:] .*
                                 ((r_a_g_grid[I_g,:,:] .+ ddeath*pam) .+ xxi * w *
                                  l_g_grid[I_g,:,:] .* y_g_grid[I_g,:,:]), 0)

    audriftF[1:I_g-1,:,:] = max.(d_g[1:I_g-1,:,:] .+ a_g_grid[1:I_g-1,:,:] .*
                                (r_a_g_grid[1:I_g-1,:,:] .+ ddeath*pam) + xxi * w *
                                l_g_grid[1:I_g-1,:,:] .* y_g_grid[1:I_g-1,:,:], 0)
    audriftF[I_g,:,:]     = max.(d_g[I_g,:,:] .+ a_g_grid[I_g,:,:] .*
                                (r_a_g_grid[I_g,:,:] .+ ddeath*pam) .+ xxi * w *
                                l_g_grid[I_g,:,:] .* y_g_grid[I_g,:,:], 0)

    budriftB[1:I_g-1,:,:] = min.(s_g[1:I_g-1,:,:] - d_g[1:I_g-1,:,:] -
                                adj_cost_fn(d_g[1:I_g-1,:,:], a_g_grid[1:I_g-1,:,:],
                                            chi0, chi1, chi2, a_lb), 0)
    budriftB[I_g,:,:]     = min.(s_g[I_g,:,:] - d_g[I_g,:,:] -
                                adj_cost_fn(d_g[I_g,:,:], a_g_grid[I_g,:,:],
                                            chi0, chi1, chi2, a_lb), 0)

    budriftF[1:I_g-1,:,:] = max.(s_g[1:I_g-1,:,:] - d_g[1:I_g-1,:,:] -
                                adj_cost_fn(d_g[1:I_g-1,:,:],a_g_grid[1:I_g-1,:,:],
                                            chi0, chi1, chi2, a_lb), 0)
    budriftF[I_g,:,:]     = max.(s_g[I_g,:,:] - d_g[I_g,:,:] -
                                adj_cost_fn(d_g[I_g,:,:], a_g_grid[I_g,:,:],
                                            chi0, chi1, chi2, a_lb),0)

    # Transition a
    chi[:,2:J,:] = -adriftB[:,2:J,:] ./ dab_grid[:,2:J,:]
    chi[:,1,:] = zeros(I,1,N)

    yy[:,2:J-1,:] = adriftB[:,2:J-1,:] ./ dab_grid[:,2:J-1,:] -
                    adriftF[:,2:J-1,:] ./ daf_grid[:,2:J-1,:]
    yy[:,1,:] = -adriftF[:,1,:] ./ daf_grid[:,1,:]
    yy[:,J,:] =  adriftB[:,J,:] ./ dab_grid[:,J,:]

    zeta[:,1:J-1,:] = adriftF[:,1:J-1,:] ./ daf_grid[:,1:J-1,:]
    zeta[:,J,:]     = zeros(I,1,N)

    centdiag = reshape(yy,I*J,N)
    lowdiag  = reshape(chi,I*J,N)
    lowdiag  = circshift(lowdiag,-I)
    updiag   = reshape(zeta,I*J,N)
    updiag   = circshift(updiag,I)

    centdiag = reshape(centdiag,I*J*N,1)
    updiag   = reshape(updiag,I*J*N,1)
    lowdiag  = reshape(lowdiag,I*J*N,1)

    aa = spdiagm(-I => vec(lowdiag)[1:end-I], 0 => vec(centdiag), I => vec(updiag)[I+1:end])

    chiu[:,2:J_g,:] = -audriftB[:,2:J_g,:] ./ dab_g_grid[:,2:J_g,:]
    chiu[:,1,:]     = zeros(I_g,1,N)

    yyu[:,2:J_g-1,:] = audriftB[:,2:J_g-1,:] ./ dab_g_grid[:,2:J_g-1,:] -
        audriftF[:,2:J_g-1,:] ./ daf_g_grid[:,2:J_g-1,:]
    yyu[:,1,:]   = -audriftF[:,1,:]   ./ daf_g_grid[:,1,:]
    yyu[:,J_g,:] =  audriftB[:,J_g,:] ./ dab_g_grid[:,J_g,:]

    zetau[:,1:J_g-1,:] = audriftF[:,1:J_g-1,:] ./ daf_g_grid[:,1:J_g-1,:]
    zetau[:,J_g,:]     = zeros(I_g,1,N)

    centdiagu = reshape(yyu, I_g * J_g, N)
    lowdiagu  = reshape(chiu,I_g * J_g, N)
    lowdiagu  = circshift(lowdiagu, -I_g)
    updiagu   = reshape(zetau, I_g * J_g, N)
    updiagu   = circshift(updiagu, I_g)

    centdiagu = reshape(centdiagu,I_g * J_g * N, 1)
    updiagu   = reshape(updiagu,  I_g * J_g * N, 1)
    lowdiagu  = reshape(lowdiagu, I_g * J_g * N, 1)

    aau = spdiagm(0 => vec(centdiagu), -I_g => vec(lowdiagu)[1:end-I_g],
                     I_g => vec(updiagu)[I_g+1:end])

    X[2:I,:,:]   = -bdriftB[2:I,:,:] ./ dbb_grid[2:I,:,:]
    X[1,:,:]     = zeros(1,J,N)
    Y[2:I-1,:,:] = bdriftB[2:I-1,:,:] ./ dbb_grid[2:I-1,:,:] -
                   bdriftF[2:I-1,:,:] ./ dbf_grid[2:I-1,:,:]
    Y[1,:,:]     = -bdriftF[1,:,:] ./ dbf_grid[1,:,:]
    Y[I,:,:]     = bdriftB[I,:,:] ./ dbb_grid[I,:,:]
    Z[1:I-1,:,:] = bdriftF[1:I-1,:,:] ./ dbf_grid[1:I-1,:,:]
    Z[I,:,:]     = zeros(1,J,N)

    centdiag = reshape(Y,I*J,N)
    lowdiag  = reshape(X,I*J,N)
    lowdiag  = circshift(lowdiag,-1)
    updiag   = reshape(Z,I*J,N)
    updiag   = circshift(updiag,1)

    centdiag = reshape(centdiag,I*J*N,1)
    updiag   = reshape(updiag,I*J*N,1)
    lowdiag  = reshape(lowdiag,I*J*N,1)

    bb = spdiagm(0 => vec(centdiag), 1 => vec(updiag)[2:end], -1 => vec(lowdiag)[1:end-1])

    Xu[2:I_g,:,:]   = -budriftB[2:I_g,:,:] ./ dbb_g_grid[2:I_g,:,:]
    Xu[1,:,:]       = zeros(1,J_g,N)
    Yu[2:I_g-1,:,:] = budriftB[2:I_g-1,:,:] ./ dbb_g_grid[2:I_g-1,:,:] -
        budriftF[2:I_g-1,:,:] ./ dbf_g_grid[2:I_g-1,:,:]
    Yu[1,:,:]       = -budriftF[1,:,:] ./ dbf_g_grid[1,:,:]
    Yu[I_g,:,:]     = budriftB[I_g,:,:] ./ dbb_g_grid[I_g,:,:]
    Zu[1:I_g-1,:,:] = budriftF[1:I_g-1,:,:] ./ dbf_g_grid[1:I_g-1,:,:]
    Zu[I_g,:,:]     = zeros(1,J_g,N)

    centdiagu = reshape(Yu,I_g*J_g,N)
    lowdiagu = reshape(Xu,I_g*J_g,N)
    lowdiagu = circshift(lowdiagu,-1)
    updiagu = reshape(Zu,I_g*J_g,N)
    updiagu = circshift(updiagu,1)

    centdiagu = reshape(centdiagu,I_g*J_g*N,1)
    updiagu   = reshape(updiagu,I_g*J_g*N,1)
    lowdiagu  = reshape(lowdiagu,I_g*J_g*N,1)

    bbu = spdiagm(0 => vec(centdiagu), 1 => vec(updiagu)[2:end], -1 => vec(lowdiagu)[1:end-1])

    return aa, bb, aau, bbu
end

@inline function transition_deriva(permanent, ddeath, pam, xxi, w, chi0, chi1,
                                   chi2, a_lb, dab_grid, daf_grid, dab_g_grid, daf_g_grid,
                                   dbb_grid, dbf_grid, dbb_g_grid, dbf_g_grid,
                                   d, d_g, s, s_g, r_a, aggZ,
                                   a, a_g, b, b_g, y)
    I, J, N  = size(d)
    I_g, J_g = size(d_g)

    # Compute drifts for HJB
    perm_const = (permanent == 1) ? ddeath * pam - aggZ : ddeath * pam

    chi, zeta = similar(d), similar(d)
    X, Z      = similar(d), similar(d)
    for i=1:I, j=1:J, n=1:N
        α = a[j] * (r_a + perm_const) + (xxi * w * real(y[1,n]))

        chi[i,j,n]  = (j==1) ? 0.0 : -(min(norm(d[i,j,n]), 0) + min(α, 0)) / dab_grid[i,j,n]
        zeta[i,j,n] = (j==J) ? 0.0 :  (max(norm(d[i,j,n]), 0) + max(α, 0)) / daf_grid[i,j,n]

        X[i,j,n] = (i==1) ? 0.0 : -(min(-d[i,j,n] -
                                  adj_cost_fn(d[i,j,n], a[j], chi0, chi1, chi2, a_lb), 0) +
                              min(s[i,j,n], 0)) / dbb_grid[i,j,n]
        Z[i,j,n] = (i==I) ? 0.0 : (max(-d[i,j,n] -
                                 adj_cost_fn(d[i,j,n], a[j], chi0, chi1, chi2, a_lb), 0) +
                             max(s[i,j,n], 0)) / dbf_grid[i,j,n]
    end
    yy = -(chi .+ zeta)

    chi  = reshape(chi, I*J, N)
    chi  = circshift(chi, -I)
    zeta = reshape(zeta, I*J, N)
    zeta = circshift(zeta, I)

    aa = spdiagm(0 => vec(yy), I => vec(zeta)[I+1:end], -I => vec(chi)[1:end-I])

    Y = -(X .+ Z)

    X  = reshape(X, I*J, N)
    X  = circshift(X, -1)
    Z  = reshape(Z, I*J, N)
    Z  = circshift(Z, 1)

    bb = spdiagm(0 => vec(Y), 1 => vec(Z)[2:end], -1 => vec(X)[1:end-1])

    chiu, zetau = similar(d_g), similar(d_g)
    Xu, Zu      = similar(d_g), similar(d_g)
    for i=1:I_g, j=1:J_g, n=1:N

        # Compute drifts for KFE -- is it possible that below (ddeathpam) is incorrect?
        audrift = d_g[i,j,n] + a_g[j] * (r_a + ddeath*pam) + xxi * w * y[1,n] -
            (permanent == 1 ? aggZ * a_g[j] : 0.0)

        budrift = s_g[i,j,n] - d_g[i,j,n] -
            adj_cost_fn(d_g[i,j,n], a_g[j], chi0, chi1, chi2, a_lb) -
            (permanent == 1 ? aggZ * b_g[i] : 0.0)

        chiu[i,j,n]  = (j==1)   ? 0.0 : -min(audrift, 0) / dab_g_grid[i,j,n]
        zetau[i,j,n] = (j==J_g) ? 0.0 :  max(audrift, 0) / daf_g_grid[i,j,n]

        Xu[i,j,n] = (i==1)   ? 0.0 : -min(budrift, 0) / dbb_g_grid[i,j,n]
        Zu[i,j,n] = (i==I_g) ? 0.0 :  max(budrift, 0) / dbf_g_grid[i,j,n]
    end

    yyu   = -(chiu .+ zetau)
    chiu  = reshape(chiu,I_g*J_g,N)
    chiu  = circshift(chiu,-I_g)
    zetau = reshape(zetau,I_g*J_g,N)
    zetau = circshift(zetau,I_g)

    aau = spdiagm(0 => vec(yyu), I_g => vec(zetau)[I_g+1:end], -I_g => vec(chiu)[1:end-I_g])

    Yu = -(Xu .+ Zu)
    Xu = reshape(Xu, I_g*J_g, N)
    Xu = circshift(Xu, -1)
    Zu = reshape(Zu, I_g*J_g, N)
    Zu = circshift(Zu, 1)

    bbu = spdiagm(0 => vec(Yu), 1 => vec(Zu)[2:end], -1 => vec(Xu)[1:end-1])

    return aa, bb, aau, bbu
end

# This file contains additional helper functions for computing the steady state.

# Discretize a state space dimension of I gridpoints
# Where agridparam is the "bending" coefficient of the grid
# i.e. grid_param = 1 implies a uniform grid
# And grid_min/grid_max are the start and endpoints of the grid
@inline function construct_asset_grid(I::Int64, grid_param::Int64, grid_min::Float64, grid_max::Float64)
    a  = collect(range(0, stop=1, length=I))
    a  = a .^ (1 / grid_param)
    a  = grid_min .+ (grid_max - grid_min) * a
    return a
end

# Provide the income and asset grids
@inline function initialize_diff_grids(a::Vector{Float64}, I::Int64, J::Int64)
    daf    = similar(a) # forward difference for a
    dab    = similar(a) # backward difference for a
    adelta = similar(a) # size of differences

    for i=1:I
        # Create a grid of lengths of overlapping intervals in a dimension.
        # The purpose is generally to compute Riemann integrals
        # by taking midpoint Riemann sums and dividing by two to adjust for
        # or average out the added Lebesgue measure given by using overlapping intervals.
        daf[i] = (i==I) ? a[end] - a[end-1] : a[i+1] - a[i]
        dab[i] = (i==1) ? a[2] - a[1]       : a[i] - a[i-1]

        # Ex:    adelta[2] =    (a_3 - a_1) / 2; if we integrate f, then
        # f[2] * adelta[2] = f_2(a_3 - a_1) / 2
        if i==1
            adelta[1]   = 0.5 * daf[1]
        elseif i==I
            adelta[end] = 0.5 * daf[end-1]
        else
            adelta[i]   = 0.5 * (daf[i-1] + daf[i])
        end
    end
    azdelta = repeat(adelta, J)

    return daf, dab, azdelta
end

# P is the Markov transition matrix
# n_income_states is the number of income states (the number of discrete states in the distribution)
# iter_num is the number of iterations willing to be accepted for convergence
@inline function compute_stationary_income_distribution(P::Matrix{Float64}, n_income_states::Int64;
                                                iter_num::Int64 = 50)
    Pt = P'
    g_z = fill(1/n_income_states, n_income_states)
    for n = 1:iter_num
        g_z_new = (speye(n_income_states) - Pt * 1000)\g_z
        diff    = maximum(abs.(g_z_new - g_z))
        if diff < 1e-5
            break
        end
        g_z = g_z_new
    end
    return g_z
end

# initial_ygrid is the unscaled income values for each income state
# income_distr is the stationary income distribution
# meanlabeff is the mean labor efficiency value that scales z such that
# The expected value of income with respect to income_distr is meanlabeff
# n_gridpoints should be the same as the number of grid points in the asset process
@inline function construct_labor_income_grid(initial_ygrid::Vector{Float64},
                                             income_distr::Vector{Float64},
                                             meanlabeff::Float64, n_gridpoints::Int64)
    z       = exp.(initial_ygrid)
    z_bar   = dot(z, income_distr)
    z  	    = (meanlabeff/z_bar) .* z
    zz      = ones(n_gridpoints, 1) * z'
    return zz
end

@inline function construct_household_problem_functions(V::Matrix{S}, w::T,
                                         coefrra::R, frisch::R, labtax::R,
                                         labdisutil::R) where {R<:AbstractFloat,T<:Real,S<:Number}

    @inline function util(c::U, h::U) where {U<:Number}
        f(x::U) = coefrra == 1.0 ? log(x) : x^(1-coefrra) / (1-coefrra)
        return f(c) - labdisutil * (h ^ (1 + 1/frisch)/(1 + 1/frisch))
    end

    @inline income(h::U, z::Float64, profshare::V, lumptransfer::V, r::T,
                   a::Float64) where {T<:Number,U<:Number,V<:Number} = h * z * w * (1 - labtax) + lumptransfer + profshare + r * a

    @inline labor(z::U, val::V) where {U<:Number,V<:Number} = (z * w * (1 - labtax) * val / labdisutil) ^ frisch

    return util, income, labor
end

@inline function construct_initial_diff_matrices(V::Matrix{T},
                                                 Vaf::Matrix{T}, Vab::Matrix{T},
                                                 income::Function, labor::Function,
                                                 h::Matrix{U}, h0::Matrix{U},
                                                 zz::Matrix{S}, profshare::Matrix{T},
                                                 lumptransfer::T,
                                                 amax::S, amin::S, coefrra::S, r::S,
                                                 daf::Vector{S}, dab::Vector{S},
                                                 maxhours::S) where {S<:Number,T<:Number,U<:Number}
    I,J = size(V)
    cf  = similar(V)
    hf  = similar(V)
    cb  = similar(V)
    hb  = similar(V)

    for j=1:J, i=1:I
        if i==I
            Vaf[end, j] = income(h[end, j], zz[end, j], profshare[end, j], lumptransfer, r, amax) ^ (-coefrra)
            Vab[i, j]   = (V[i, j] - V[i-1, j]) / dab[i]
        elseif i==1
            Vaf[i, j]   = (V[i+1, j] - V[i, j]) / daf[i]
            Vab[1, j]   = income(h0[1, j], zz[1, j], profshare[1, j], lumptransfer, r, amin) ^ (-coefrra)
        else
            Vaf[i, j]   = (V[i+1, j] - V[i, j]) / daf[i]
            Vab[i, j]   = (V[i, j] - V[i-1, j]) / dab[i]
        end
        cf[i,j] = Vaf[i,j] ^ (-1 / coefrra)
        cb[i,j] = Vab[i,j] ^ (-1 / coefrra)

        hf[i,j] = min(norm(labor(zz[i,j], Vaf[i,j])), maxhours)
        hb[i,j] = min(norm(labor(zz[i,j], Vab[i,j])), maxhours)
    end

    return Vaf, Vab, cf, hf, cb, hb
end

# For initialization
@inline function calculate_ss_equil_vars(zz::Matrix{Float64}, m_ss::Float64, meanlabeff::Float64,
                                 lumptransferpc::Float64, govbondtarget::Float64)

    N_ss         = complex(1/3) # steady state hours: so that quarterly GDP = 1 in s.s
    Y_ss         = complex(1.)
    B_ss         = govbondtarget * Y_ss
    profit_ss    = complex((1 - m_ss) * Y_ss)
    profshare    = zz / meanlabeff * profit_ss
    lumptransfer = complex(lumptransferpc * Y_ss)

    return N_ss, Y_ss, B_ss, profit_ss, profshare, lumptransfer
end

@inline function calculate_ss_equil_vars(zz::Matrix{Float64},
                                 h::Matrix{ComplexF64},
                                 g::Vector{ComplexF64}, azdelta::Vector{Float64},
                                 aa::Matrix{Float64}, m_ss::Float64,
                                 meanlabeff::Float64, lumptransferpc::Float64,
                                 govbondtarget::Float64)

    # equilibrium objects
    Y_ss = N_ss  = sum(vec(zz) .* vec(h) .* g .* azdelta)
    Y_ss         = N_ss
    B_ss         = sum(g .* vec(aa) .* azdelta)
    profit_ss    = (1 - m_ss) * Y_ss
    profshare    = zz / meanlabeff * profit_ss
    lumptransfer = lumptransferpc * Y_ss
    bond_err     = B_ss / Y_ss - govbondtarget

    return N_ss, Y_ss, B_ss, profit_ss, profshare, lumptransfer, bond_err
end

@inline function hours_iteration(income::Function, labor::Function,
                                 zz::Matrix{Float64},
                                 profshare::Matrix{T},
                                 lumptransfer::T,
                                 aa::Matrix{Float64},
                                 coefrra::Float64, r::S,
                                 cf::Matrix{T}, hf::Matrix{T},
                                 cb::Matrix{T}, hb::Matrix{T},
                                 c0::Matrix{T}, h0::Matrix{U},
                                 maxhours::Float64,
                                 niter_hours::Int64) where {S<:Number,T<:Number,U<:Number}
    I, J = size(zz)
    for ih = 1:niter_hours
        for j=1:J, i=1:I
            if i==I
                cf[end, j] = income(hf[end, j], zz[end, j], profshare[end, j], lumptransfer, r, aa[end, j])
                hf[end, j] = labor(zz[end, j], cf[end, j] ^ (-coefrra))
                hf[end, j] = min(norm(hf[end, j]), maxhours)
            elseif i==1
                cb[1, j] = income(hb[1, j], zz[1, j], profshare[1, j], lumptransfer, r, aa[1, j])
                hb[1, j] = labor(zz[1, j], cb[1, j] ^ (-coefrra))
                hb[1, j] = min(norm(hb[1, j]), maxhours)
            end
            c0[i, j] = income(h0[i, j], zz[i, j], profshare[i, j], lumptransfer, r, aa[i, j])
            h0[i, j] = labor(zz[i, j], c0[i, j]^(-coefrra))
            h0[i, j] = min(norm(h0[i, j]), maxhours)
        end
    end
    return cf, hf, cb, hb, c0, h0
end

# For compute_steady_state; choose upwinding direction
@inline function upwind(ρ::Float64, V::Matrix{T}, args...;
                        Δ_HJB::Float64 = 1e6) where {T<:Number}
    A, u, h, c, s = upwind(args...)
    I, J = size(u)
    B    = (1 / Δ_HJB + ρ) * speye(T, I*J) - A
    b    = reshape(u, I*J) + reshape(V, I*J) / Δ_HJB
    V    = reshape(B \ b, I, J)
    return V, A, u, h, c, s
end

# For equilibrium_conditions; choose upwinding direction.
@inline function upwind(util::Function,
                        A_switch::SparseMatrixCSC{S, Int64},
                        cf::Matrix{T}, cb::Matrix{T},
                        c0::Matrix{T}, hf::Matrix{T},
                        hb::Matrix{T}, h0::Matrix{T},
                        sf::Matrix{T}, sb::Matrix{T},
                        Vaf::Matrix{T}, Vab::Matrix{T},
                        daf::Vector{Float64},
                        dab::Vector{Float64}) where {T<:Number, S<:Number}
    I,J = size(sb)
    h   = similar(sb)
    c   = similar(sb)
    s   = similar(sb)
    u   = similar(sb)
    X   = similar(sb)
    Z   = similar(sb)
    Y   = similar(sb)

    for i in eachindex(sb)

        Vf = (cf[i] > 0) * (util(cf[i], hf[i]) + sf[i] * Vaf[i]) + (cf[i] <= 0) * (-1e12)
        Vb = (cb[i] > 0) * (util(cb[i], hb[i]) + sb[i] * Vab[i]) + (cb[i] <= 0) * (-1e12)
        V0 = (c0[i] > 0) * util(c0[i], h0[i]) + (c0[i] <= 0) * (-1e12)

        Iunique = (sb[i] < 0) * (1 - (sf[i] > 0)) + (1 - (sb[i] < 0)) * (sf[i] > 0)
        Iboth = (sb[i] < 0) * (sf[i] > 0)

        Ib = Iunique * (sb[i] < 0) * (Vb > V0) + Iboth * (Vb == max(max(Vb, Vf), V0))
        If = Iunique * (sf[i] > 0) * (Vf > V0) + Iboth * (Vf == max(max(Vb, Vf), V0))
        I0 = 1 - Ib - If

        h[i] = hf[i] * If + hb[i] * Ib + h0[i] * I0
        c[i] = cf[i] * If + cb[i] * Ib + c0[i] * I0
        s[i] = sf[i] * If + sb[i] * Ib
        u[i] = util(c[i], h[i])

        # Construct A matrix
        X[i] = -Ib * sb[i] / dab[((i-1) % I) + 1]
        Z[i] =  If * sf[i] / daf[((i-1) % I) + 1]
        Y[i] = -Z[i] - X[i]
    end

    # R: Pretty sure the indexing for this ought be the same as in the KrusellSmith model
    X[1,:] .= T == ComplexF64 ? complex(0.) : 0.
    Z[I,:] .= T == ComplexF64 ? complex(0.) : 0.

    A = spdiagm(-1 => reshape(X,I*J)[2:I*J], 0 => reshape(Y,I*J), 1 => reshape(Z,I*J)[1:I*J-1]) + A_switch

    return A, u, h, c, s
end


# Using the market clearing condition on bonds to determine whether or not
# an equilibrium has been reached
@inline function check_bond_market_clearing(bond_err::ComplexF64, crit_S::Float64,
                                            r::Float64, r_min::Float64, r_max::Float64,
                                            r_ρ::Float64, ρ_min::Float64, ρ_max::Float64,
                                            iter_r::Bool, iter_ρ::Bool)
    clearing_condition = false
    # Using the market clearing condition on bonds to determine whether or not
    # an equilibrium has been reached
    if abs(bond_err) > crit_S
        if bond_err > 0
            if iter_r
                r_max  = r
                r      = 0.5 * (r + r_min)
            elseif iter_ρ
                ρ_min = r_ρ
                r_ρ   = 0.5 * (r_ρ + ρ_max)
            end
        else
            if iter_r
                r_min  = r
                r      = 0.5 * (r + r_max)
            elseif iter_ρ
                ρ_max = r_ρ
                r_ρ   = 0.5 * (r_ρ + ρ_min)
            end
        end
    else
        clearing_condition = true
    end
    return r, r_min, r_max, r_ρ, ρ_min, ρ_max, clearing_condition
end
