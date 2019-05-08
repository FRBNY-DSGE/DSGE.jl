function persistent_skill_process(sH_over_sL::AbstractFloat, pLH::AbstractFloat,
                                  pHL::AbstractFloat, ns::Int)
    f1 = [[1-pLH pLH];[pHL 1-pHL]] # f1[i,j] is prob of going from i to j
    ss_skill_distr = [pHL/(pLH+pHL); pLH/(pLH+pHL)]
    slo    = 1.0 / (ss_skill_distr'*[1;sH_over_sL])
    sgrid  = slo*[1;sH_over_sL]
    sscale = sgrid[2] - sgrid[1]
    swts   = (sscale/ns)*ones(ns) # Quadrature weights
    f      = f1 ./ repeat(swts',ns,1)
    return f, sgrid, swts, sscale
end

function cash_grid(sgrid::AbstractArray, ω::AbstractFloat, H::AbstractFloat,
                   r::AbstractFloat, η::AbstractFloat, γ::AbstractFloat,
                   T::AbstractFloat, zlo::AbstractFloat, na::Int)
    smin = minimum(sgrid)*zlo                                   # lowest possible skill
    xlo_ss = ω*smin*H - (1+r)*η*exp(-γ) + T + sgrid[1]*ω*H*0.05 # lowest SS possible cash on hand

    xlo = xlo_ss                        # lower bound on cash on hand - could be < xlo_ss
    xhi = max(xlo*2, xlo + 12.0) # TODO # upper bound on cash on hand
    xscale = (xhi-xlo)                  # size of w grids

    # Make grids
    xgrid = collect(range(xlo,stop = xhi, length = na)) # Evenly spaced grid
    xwts  = (xscale/na)*ones(na)          # Quadrature weights
    return xgrid, xwts, xlo, xhi, xscale
end
