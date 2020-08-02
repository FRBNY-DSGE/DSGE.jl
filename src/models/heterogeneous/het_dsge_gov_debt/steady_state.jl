function steadystate!(m::HetDSGEGovDebt;
                      βlo::S = 0.5*exp(m[:γ].scaledvalue)/(1 + m[:r].scaledvalue),
                      βhi::S = exp(m[:γ].scaledvalue)/(1 + m[:r].scaledvalue),
                      excess::S = 5000.,
                      tol::S = 1e-4,
                      maxit::Int64 = 20,
                      βband::S = 1e-2) where {S<:AbstractFloat}
    # If we have already solved for βstar (i.e. it's not NaN) and we only want to
    # estimate the non steady state parameters, there's no need to recompute
    # zlo/zhi, etc.
    if !isnan(m[:βstar].value) && get_setting(m, :estimate_only_non_steady_state_parameters)
        return
    else
        reset_grids!(m)
        nx = get_setting(m, :nx)
        ns = get_setting(m, :ns)
        ni = get_setting(m, :ni)
        nz = get_setting(m, :nz)

        # Determines whether random or not
        us, zs = if get_setting(m, :fix_random_matrices)
            get_setting(m, :us), get_setting(m, :zs)
        else
            generate_us_and_zs(ni, nz)
        end

        # This is a test setting; TODO: Remove because now we can just fix the randomness
        if get_setting(m, :steady_state_only)
            find_steadystate!(m; βlo = βlo, βhi = βhi, excess = excess, tol = tol, maxit = maxit,
                              βband = βband)
            return
        else
            # Do you want to calibrate for matching the income moments?
            if get_setting(m, :calibrate_income_targets)

                target = get_setting(m, :calibration_targets)
                lower  = get_setting(m, :calibration_targets_lb)
                upper  = get_setting(m, :calibration_targets_ub)
                m[:sH_over_sL], m[:zlo], _, _ = best_fit(m[:pLH].value, m[:pHL].value,
                                                         target, lower, upper, us, zs)
            else
                m[:varlinc], m[:vardlinc] = skill_moments(m[:sH_over_sL].value, m[:zlo].value,
                                                          m[:pLH].value, m[:pHL].value,
                                                          us, zs, ni)
            end
            # Whether one calibrates or not, zhi needs to be updated here
            m[:zhi] = 2.0 - m[:zlo].value

            # Construct Markov transition matrix for skill
            f, sgrid, swts, sscale = persistent_skill_process(m[:sH_over_sL].value, m[:pLH].value,
                                                              m[:pHL].value, get_setting(m, :ns))
            m.grids[:sgrid] = Grid(sgrid, swts, sscale)
            m.grids[:fgrid] = f

            xgrid, xwts, xlo, xhi, xscale = cash_grid(sgrid, m[:ωstar].value, m[:H].value,
                                                      m[:r].scaledvalue, m[:η].value,
                                                      m[:γ].scaledvalue,
                                                      m[:Tstar].value, m[:zlo].value, nx)

            m.grids[:xgrid] = Grid(uniform_quadrature(xscale), xlo, xhi, nx, scale = xscale)
            m <= Setting(:xlo, xlo)
            m <= Setting(:xhi, xhi)
            m <= Setting(:xscale, xscale)
            xswts = kron(swts, xwts)
            m.grids[:weights_total] = xswts

            # Once have updated grids, can call steady state and compute other two moments
            find_steadystate!(m; βlo = βlo, βhi = βhi,
                              excess = excess, tol = tol, maxit = maxit,
                              βband = βband)
            m[:mpc] = ave_mpc(m[:μstar].value,   m[:cstar].value, xgrid, xswts, nx, ns)
            m[:pc0] = frac_zero(m[:μstar].value, m[:cstar].value, xgrid, xswts, ns)
        end
    end
end

function find_steadystate!(m::HetDSGEGovDebt;
                           βlo::S = 0.5*exp(m[:γ].scaledvalue)/(1 + m[:r].scaledvalue),
                           βhi::S = exp(m[:γ].scaledvalue)/(1 + m[:r].scaledvalue),
                           excess::S = 5000.,
                           tol::S = 1e-4,
                           maxit::Int64 = 20,
                           βband::S = 1e-2) where {S<:AbstractFloat}
    # Load settings
    nx = get_setting(m, :nx)
    ns = get_setting(m, :ns)

    xlo = get_setting(m, :xlo)
    xhi = get_setting(m, :xhi)

    zlo = m[:zlo].value
    zhi = m[:zhi].value

    # Load parameters
    R  = 1 + m[:r].scaledvalue
    H  = m[:H].value
    η  = m[:η].value
    γ  = m[:γ].scaledvalue
    ω  = m[:ωstar].value
    T  = m[:Tstar].value
    bg = m[:bg].value

    # Load grids
    xgrid = m.grids[:xgrid].points
    sgrid = m.grids[:sgrid].points
    xswts = m.grids[:weights_total]
    f     = m.grids[:fgrid]

    βlo_temp = βlo
    βhi_temp = βhi

    reject  = false
    counter = 1

    n  = ns*nx
    c  = zeros(n)
    bp = zeros(n)
    KF = zeros(n,n)
    μ  = zeros(n)

    # Initial guess
    β   = 1.0
    Win = Vector{Float64}(undef, n)
    Win_guess = ones(n)

    if get_setting(m, :use_last_βstar) && !isnan(m[:βstar].value)

        # This short-circuits computation of the policy function
        βlo = βhi = m[:βstar].value
    elseif !isnan(m[:βstar].value)

        # If one has computed β* before, we first bisect into a neighborhood around it
        βlo_temp = m[:βstar].value - βband
        βhi_temp = m[:βstar].value + βband

        c, bp, Win, KF, reject = policy_hetdsgegovdebt(nx, ns, βlo_temp, R, ω, H, η, T, γ,
                                                       zhi, zlo, xgrid, sgrid, xswts, Win_guess,
                                                       f, damp = get_setting(m, :policy_damp), maxit = get_setting(m, :policy_maxit))
        excess_lo, μ = compute_excess(xswts, KF, bp, bg)

        if excess_lo < 0 && abs(excess_lo) > tol
            βlo = βlo_temp
            c, bp, Win, KF, reject = policy_hetdsgegovdebt(nx, ns, βhi_temp, R, ω, H, η, T, γ, zhi,
                                                           zlo, xgrid, sgrid, xswts, Win_guess, f,
                                                           damp = get_setting(m, :policy_damp), maxit = get_setting(m, :policy_maxit))
            excess_hi, μ = compute_excess(xswts, KF, bp, bg)

            if excess_hi > 0
                βhi = βhi_temp
            end
        elseif excess_lo >= 0
            βhi = βlo_temp
        end
    end

    while abs(excess) > tol && counter < maxit # clearing markets
        β = (βlo + βhi) / 2.0
        c, bp, Win, KF, reject = policy_hetdsgegovdebt(nx, ns, β, R, ω, H, η, T, γ, zhi, zlo,
                                                       xgrid, sgrid, xswts, Win_guess, f,
                                                       damp = get_setting(m, :policy_damp), maxit = get_setting(m, :policy_maxit))
        excess, μ = compute_excess(xswts, KF, bp, bg)
        # bisection
        if excess > 0
            βhi = β
        elseif excess < 0
            βlo = β
        end
        counter += 1
        if get_setting(m, :use_last_βstar) && !isnan(m[:βstar].value)
            break
        end
    end

  #=  @show β
    This all was for testing when we were trying to compare the old steadystate to the new steadystate
    GR.inline("png")
    p = Plots.plot(xgrid, μ[1:300], label = "lowskill")
    Plots.plot!(p, xgrid, μ[301:600], label = "highskill")
    Plots.savefig(p, "agrid_vs_D_together_old.png")
    p = Plots.plot(xgrid, μ[1:300], label = "lowskill")
    Plots.savefig(p, "agrid_vs_D_lowskill_old.png")
    p = Plots.plot(xgrid, μ[301:600], label = "highskill")
    Plots.savefig(p, "agrid_vs_D_highskill_old.png")

    Plots.plot(μ)
    Plots.savefig("mu_old.png")

    p = Plots.plot(xgrid, c[1:300], label = "lowskill")
    Plots.plot!(p, xgrid, c[301:600], label = "highskill")
    Plots.savefig(p, "agrid_vs_cpol_old.png")

    p = Plots.plot(xgrid, xgrid - c[1:300], label = "lowskill")
    Plots.plot!(p, xgrid, xgrid - c[301:600], label = "highskill")
    Plots.savefig(p, "agrid_minus_cpol_old.png")

    p = Plots.plot()
    for mult = 1:6
        for i = 1:mult*100
            Plots.plot!(p, KF[:, i])
        end
        Plots.savefig(p,"KF_old_mult=$mult.png")
    end
    for i = 1:600
        p = Plots.plot(KF[:, i])
        Plots.savefig(p, "KF_dir/KF_old_$i.png")
    end=#


    # If policy function does not converge, we signal to likelihood that should reject
    m <= Setting(:auto_reject, reject)
    m[:lstar]  = Win
    m[:cstar]  = c
    m[:μstar]  = μ
    m[:βstar]  = β
    m[:β_save] = β
    nothing
end

# These are for computing moments to compare against targets (the next few functions)
function ave_mpc(m::AbstractArray, c::AbstractArray, agrid::AbstractArray,
                 aswts::AbstractArray, na::Int, ns::Int)
	mpc = 0.
	for iss=1:ns
		i = na*(iss-1)+1
		mpc += aswts[i]*m[i]*(c[i+1] - c[i])/(agrid[2] - agrid[1])
		for ia=2:na
			i = na*(iss-1)+ia
			mpc += aswts[i]*m[i]*(c[i] - c[i-1])/(agrid[ia] - agrid[ia-1])
		end
	end
	return mpc
end

function frac_zero(m::AbstractArray, c::AbstractArray, agrid::AbstractArray,
                   aswts::AbstractArray, ns::Int)
	return sum(aswts .* m .* (c .== repeat(agrid, ns)))
end

function ssample(us::Matrix{S}, P::Matrix{S}, πss::AbstractArray,
                 ni::Int) where {S<:AbstractFloat}
	shist = ones(Int,ni,8)
	for i=1:ni
		if us[i,1] > πss[1]
			shist[i,1] = 2
		end
		for t=2:8
			if us[i,t] > P[shist[i,t-1],1]
				shist[i,t] = 2
			end
		end
	end
	return shist
end

function ln_annual_inc(zhist::Matrix{S}, us::Matrix{S}, zlo::S, P::Matrix{S},
                       πss::AbstractArray, sgrid::AbstractArray,
                       ni::Int) where {S<:AbstractFloat}
  	s_inds = ssample(us,P,πss,ni)
	linc1 = zeros(ni)
	linc2 = zeros(ni)
	for i=1:ni
		inc1 = 0.
		for t=1:4
			zshock = 1. + (1. - zlo)*(zhist[i,t]-1.)
			sshock = (sgrid[1] + (sgrid[2] - sgrid[1])*(s_inds[i,t] - 1.))
			inc1 += zshock*sshock
		end
		inc2 = 0.
		for t=5:8
			zshock = 1. + (1. - zlo)*(zhist[i,t]-1.)
			sshock = (sgrid[1] + (sgrid[2] - sgrid[1])*(s_inds[i,t] - 1.))
			inc2 += zshock*sshock
		end
		linc1[i] += log(inc1)
     	linc2[i] += log(inc2)
	end
	return linc1, linc2
end

function skill_moments(sH_over_sL::Real, zlo::Real, pLH::S, pHL::S, us::Matrix{S},
                       zs::Matrix{S}, ni::Int = 10000) where {S<:AbstractFloat}
	πL    = pHL / (pLH + pHL)
	πss   = [πL; 1.0-πL]
	P     = [[1.0-pLH pLH]; [pHL 1.0-pHL]]
	slo   = 1.0 / (πL + (1-πL) * sH_over_sL)
	shi   = sH_over_sL * slo
	sgrid = [slo; shi]
	linc1, linc2 = ln_annual_inc(zs, us, zlo, P, πss, sgrid, ni)
	return var(linc1), var(linc2 - linc1)
end

"""
```
function loss(x::Vector{S}, target::Vector{S}) where {S<:AbstractFloat}
```

Computes for given vector and target vector, the sum of absolute loss.
"""
loss(x::Vector{S}, target::Vector{S}) where {S<:AbstractFloat} = sum(abs.(x-target))

"""
```
function best_fit(pLH::S, pHL::S, target::Vector{S}, lower::Vector{S}, upper::Vector{S},
                  us::Matrix{S}, zs::Matrix{S}, max_iter::Int = 20,
                  initial_guess::Vector{S} = [6.3, 0.03]) where {S<:AbstractFloat}
```

Uses Nelder-Mead (gradient descent algorithm) to optimize for income moments.
"""
function best_fit(pLH::S, pHL::S, target::Vector{S}, lower::Vector{S}, upper::Vector{S},
                  us::Matrix{S}, zs::Matrix{S}, max_iter::Int = 20,
                  initial_guess::Vector{S} = [6.3, 0.03]) where {S<:AbstractFloat}

    skill_moments_f(x) = loss(collect(skill_moments(x[1], x[2], pLH, pHL, us, zs)), target)

    res = optimize(skill_moments_f, lower, upper, initial_guess, Fminbox(NelderMead()),
                   Optim.Options(f_calls_limit = max_iter))

    sH_over_sL_argmin, zlo_argmin = Optim.minimizer(res)
    min_varlinc, min_vardlinc = skill_moments(sH_over_sL_argmin, zlo_argmin, pLH,
                                              pHL, us, zs)
    return sH_over_sL_argmin, zlo_argmin, min_varlinc, min_vardlinc
end

function zsample(uz::Matrix{S}, zgrid::AbstractArray, zcdf::AbstractArray,
                 ni::Int, nz::Int) where {S<:AbstractFloat}
    zave = 0.5*zgrid[1:nz-1]+0.5*zgrid[2:nz]
	zs = zeros(ni,8)
	for i=1:ni
		for t=1:8
			for iz=1:nz-1
				if zcdf[iz] < uz[i,t] <= zcdf[iz+1]
					zs[i,t] = zave[iz]
				end
			end
		end
	end
	return zs
end

# This function computes excess savings. The KF that pops out of steadystate is nx*ns by nx*ns but we want to reduce that to a vector of length nx*ns which is the stationary distribution (this is what the eigenstuff is doing). Stationary because no matter where you are now, you want to have the same probability of going to all the different gridpoints tomorrow.
@inline function compute_excess(xswts::Vector{S}, KF::Matrix{S}, bp::Vector{S},
                                bg::S; print_warning::Bool = false,
                                tol::S = 2e-1) where {S<:Float64}
    LPMKF = xswts[1] * KF

    # Find eigenvalue closest to 1
    (D,V) = (eigen(LPMKF)...,)
    max_D = argmax(abs.(D))
    D = D[max_D]

    if abs(D - 1) > tol && print_warning
        @warn "Your eigenvalue is too far from 1, something is wrong."
    end

    # Pick eigenvector associated w/ largest eigenvalue and moving it back to values
    μ = real(V[:,max_D])
    # @show μ
    μ = μ ./ dot(xswts, μ) # Scale of eigenvectors not determinate: rescale to integrate to 1
    excess = dot(xswts, (μ .* bp)) - bg # Compute excess supply of savings, which is a fn of w
    return excess, μ
end

# policy function (consumption and ell = 1/consumption are the "policies" being decied on here)
function policy_hetdsgegovdebt(nx::Int, ns::Int, β::S, R::S, ω::S, H::S, η::S,
                               T::S, γ::S, zhi::S, zlo::S, xgrid::Vector{S},
                               sgrid::Vector{S}, xswts::Vector{S}, Win::Vector{S},
                               f::Matrix{S}, dist::S = 1., tol::S = 1e-4;
                               maxit::Int64 = 500, damp::S = 0.5) where {S<:AbstractFloat}
    n    = nx*ns
    c    = zeros(n)                  # consumption
    bp   = Vector{Float64}(undef, n) # savings
    Wout = Vector{Float64}(undef, length(Win))
    counter = 1
    reject = false
    qfunction(x::Float64) = mollifier_hetdsgegovdebt(x, zhi, zlo)

    while dist>tol && counter<maxit
        # compute c(w) given guess for Win = β*R*E[u'(c_{t+1})]
        for iss in 1:ns
            for ia in 1:nx
                c[nx*(iss-1)+ia] = min(1/Win[nx*(iss-1)+ia],xgrid[ia]+η)
            end
        end

        bp = R*(exp(-γ))*(repeat(xgrid, ns) - c)  # compute bp(w) given guess for Win
        # Wout and Win are the ell function (in the paper!). The point is we're trying to get the policy function (ell) to converge (which is the business with Win and Wout because you need to compare the difference between Win and Wout to see if it converged)
        Wout = parameterized_expectations_hetdsgegovdebt(nx, ns, β, R, ω, H, T, γ,
                                                         qfunction, xgrid,
                                                         sgrid, xswts, c, bp, f)
        dist = maximum(abs.(Wout - Win))
        Win  = damp*Wout + (1.0-damp) * Win
        counter += 1
    end
    if counter == maxit
        @warn "Euler iteration did not converge"
        reject = true
        return c, bp, Wout, zeros(n,n), reject
    end
    tr = kolmogorov_fwd_hetdsgegovdebt(nx, ns, ω, H, T, R, γ, qfunction, xgrid, sgrid, bp, f)
    return c, bp, Wout, tr, reject
end

# This function is the formula for ell_t in the middle top of page 7 of the paper
@inline function parameterized_expectations_hetdsgegovdebt(nx::Int, ns::Int, β::S,
                                                           R::S, ω::S, H::S, T::S, γ::S,
                                                           qfunc::Function, xgrid::Vector{S},
                                                           sgrid::Vector{S}, xswts::Vector{S},
                                                           c::Vector{S}, bp::Vector{S},
                                                           f::Matrix{S}) where {S<:AbstractFloat}
    l_out = zeros(nx*ns)
    for iss=1:ns     # loop off s (only 2 values)
        for ia=1:nx # loop over a
            sumn = 0.0
            for isp=1:ns # loop over s' (only 2 values)
                for iap=1:nx # loop over a'
                    # Name in code                                                 : notation in paper
                    #xswts                                                         : weights for the integral approximation
                    #c[nx*(isp-1)+iap]                                             : c_{t+1}(a', s')
                    # qfunc                                                        : g(⋅)
                    # qfunc((xgrid[iap] - bp[nx*(iss-1)+ia] - T) / (ω*H*sgrid[isp]): (a'-(1+r)exp(-γ)(a-c_t(a,s) - T))/ w*s'*H
                    # f[iss, isp]/ sgrid[isp]                                      : p(s'|s)/s'
                    # β*R*exp(-γ))/ω*H  (line 429)                                 β*(1+r)*exp(-z)/ω_{t+1}*H{t+1}
                    # When you're trying to match code to paper, be creative in the ordering of stuff as the variables in the paper are commutatively multipled in different orders
                    sumn += (xswts[nx*(isp-1)+iap]/c[nx*(isp-1)+iap]) *
                        qfunc((xgrid[iap] - bp[nx*(iss-1)+ia] - T) /
                              (ω*H*sgrid[isp])) * f[iss,isp] ./ sgrid[isp]
                end
            end
            l_out[nx*(iss-1)+ia] = (β*R*(exp(-γ))/ω*H)*sumn
        end
    end
    return l_out
end


# This function is the exprression for m_{t+1} in the middle of page 7 of the paper. See the key in lines 417-424 of this file to help you in matching stuff up
@inline function kolmogorov_fwd_hetdsgegovdebt(nx::Int, ns::Int, ω::S, H::S, T::S, R::S,
                                               γ::S, qfunc::Function, xgrid::Vector{S},
                                               sgrid::Vector{S}, bp::Vector{S},
                                               f::Matrix{S}) where {S<:AbstractFloat}
    tr = zeros(nx*ns,nx*ns)
    for iss=1:ns
        for ia=1:nx
            for isp=1:ns
                for iap=1:nx
                    tr[nx*(isp-1)+iap, nx*(iss-1)+ia] = qfunc((xgrid[iap] -
                                                               bp[nx*(iss-1)+ia] - T) /
                                                              (ω*H*sgrid[isp])) * f[iss,isp] ./
                                                              (ω * H * sgrid[isp])
                end
            end
        end
    end
    return tr
end

# also known as qfunction of g(⋅)
@inline function mollifier_hetdsgegovdebt(z::S, ehi::S, elo::S) where {S<:AbstractFloat}
    In = 0.443993816237631
    if z<ehi && z>elo
        temp = -1.0 + 2.0 * (z - elo) / (ehi - elo)
        return (2.0 / (ehi - elo)) * exp(-1.0 / (1.0 - temp^2)) / In
    end
    return 0.0
end
# also known as qfunction of g(⋅) but the derivative of it
@inline function dmollifier_hetdsgegovdebt(x::S, ehi::S, elo::S) where {S<:AbstractFloat}
    In = 0.443993816237631
    if x<ehi && x>elo
        temp = (-1.0 + 2.0*(x-elo)/(ehi-elo))
        out  = -(2*temp./((1 - temp.^2).^2)).*(2/(ehi-elo)) .*
            mollifier_hetdsgegovdebt(x, ehi, elo)
    else
        out = 0.0
    end
    return out
end

function calibrate_pLH_pHL(m::HetDSGEGovDebt)
    target_mpc, target_pc0 = get_setting(m, :targets)
    σt_mpc, σt_pc0         = get_setting(m, :target_σt)

    function minimize_penalty(m, pLHpHL::Vector{Float64})
        m[:pLH] = pLHpHL[1]
        m[:pHL] = pLHpHL[2]
        steadystate!(m)
        if get_setting(m, :auto_reject)
            m <= Setting(:auto_reject, false)
            @show pLHpHL, Inf
            return Inf
        end
        @show pLHpHL, m[:mpc].value, m[:pc0].value
        println(0.5*((log(m[:mpc])-log(target_mpc))^2/σt_mpc +
                    (log(m[:pc0].value)-log(target_pc0))^2/σt_pc0^2))
        return 0.5*((log(m[:mpc])-log(target_mpc))^2/σt_mpc +
                    (log(m[:pc0].value)-log(target_pc0))^2/σt_pc0^2)
    end
    minimize_me(x) = minimize_penalty(m, x)
    res = optimize(minimize_me, [0.005, 0.005], [0.095, 0.095], [0.01125, 0.03],
                   Fminbox(Optim.NelderMead()), Optim.Options(f_calls_limit=300))
    println(res)
    return res
end
