# constructs f which is the transition matrix between different skill levels
# also contstructs the sgrid (agin this uses uniform weights where the weights sum to the difference between shi and slo, not 1)
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

# constructs the cash on hand grid (known as both agrid and xgrid and kinda used interchangebly throughout the code
function cash_grid(sgrid::AbstractArray, ω::AbstractFloat, H::AbstractFloat,
                   r::AbstractFloat, η::AbstractFloat, γ::AbstractFloat,
                   T::AbstractFloat, zlo::AbstractFloat, na::Int)
    smin = minimum(sgrid)*zlo                                   # lowest possible skill
    xlo_ss = ω*smin*H - (1+r)*η*exp(-γ) + T + sgrid[1]*ω*H*0.05 # lowest SS possible cash on hand

    xlo = xlo_ss                        # lower bound on cash on hand - could be < xlo_ss
    xhi = max(xlo*2, xlo + 12.0)         # upper bound on cash on hand
    xscale = (xhi-xlo)                  # size of w grids

    # Make grids
    xgrid = collect(range(xlo,stop = xhi, length = na)) # Evenly spaced grid
    xwts  = (xscale/na)*ones(na)          # Quadrature weights
    return xgrid, xwts, xlo, xhi, xscale
end

"""
```
generate_us_and_zs(ni, nz)
```

There is no need to recall this function, unless one wants to undo the seeding
in all past saved output.
"""
function generate_us_and_zs(ni, nz)
    us = rand(ni, 8)
    uz = rand(ni, 8)

    zgrid  = collect(range(0., stop = 2., length = nz))
    zprob  = [2*mollifier_hetdsgegovdebt(zgrid[i], 2., 0.) / nz for i=1:nz]
    zprob /= sum(zprob)

    zcdf = cumsum(zprob)
    zave = 0.5 * zgrid[1:nz-1] + 0.5 * zgrid[2:nz]
    zs   = zsample(uz, zgrid, zcdf, ni, nz)

    return us, zs
    #=
    JLD2.jldopen("$HETDSGEGOVDEBT/reference/us_zs.jld2", true, true, true, IOStream) do file
        file["us"] = us
        file["zs"] = zs
    end
    =#
end
