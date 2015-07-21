module FinancialFrictionsFunctions

using Distributions

export ζ_spb_fn, ζ_bω_fn, ζ_zω_fn, nk_fn, μ_fn, ω_fn, G_fn, Γ_fn, dG_dω_fn, d2G_dω2_fn,
       dΓ_dω_fn, d2Γ_dω2_fn, dG_dσ_fn, d2G_dωdσ_fn, dΓ_dσ_fn, d2Γ_dωdσ_fn

# This file contains functions that are used to compute financial frictions steady-state
# values from parameter values. Specifically, it is called by models/m990/parameters.jl

function ζ_spb_fn(z, σ, sprd)
    zetaratio = ζ_bω_fn(z, σ, sprd)/ζ_zω_fn(z, σ, sprd)
    nk = nk_fn(z, σ, sprd)
    return -zetaratio/(1-zetaratio)*nk/(1-nk)
end

function ζ_bω_fn(z, σ, sprd)
    nk          = nk_fn(z, σ, sprd)
    μstar       = μ_fn(z, σ, sprd)
    ωstar       = ω_fn(z, σ)
    Γstar       = Γ_fn(z, σ)
    Gstar       = G_fn(z, σ)
    dΓ_dωstar   = dΓ_dω_fn(z)
    dG_dωstar   = dG_dω_fn(z, σ)
    d2Γ_dω2star = d2Γ_dω2_fn(z, σ)
    d2G_dω2star = d2G_dω2_fn(z, σ)
    return ωstar*μstar*nk*(d2Γ_dω2star*dG_dωstar - d2G_dω2star*dΓ_dωstar)/
        (dΓ_dωstar - μstar*dG_dωstar)^2/sprd/(1 - Γstar + dΓ_dωstar*(Γstar - μstar*Gstar)/
            (dΓ_dωstar - μstar*dG_dωstar))
end

function ζ_zω_fn(z, σ, sprd)
    μstar = μ_fn(z, σ, sprd)
    return ω_fn(z, σ)*(dΓ_dω_fn(z) - μstar*dG_dω_fn(z, σ))/
        (Γ_fn(z, σ) - μstar*G_fn(z, σ))
end

nk_fn(z, σ, sprd) = 1 - (Γ_fn(z, σ) - μ_fn(z, σ, sprd)*G_fn(z, σ))*sprd
μ_fn(z, σ, sprd)  = (1 - 1/sprd)/(dG_dω_fn(z, σ)/dΓ_dω_fn(z)*(1 - Γ_fn(z, σ)) + G_fn(z, σ))
ω_fn(z, σ)        = exp(σ*z - σ^2/2)
G_fn(z, σ)        = cdf(Normal(), z-σ)
Γ_fn(z, σ)        = ω_fn(z, σ)*(1 - cdf(Normal(), z)) + cdf(Normal(), z-σ)
dG_dω_fn(z, σ)    = pdf(Normal(), z)/σ
d2G_dω2_fn(z, σ)  = -z*pdf(Normal(), z)/ω_fn(z, σ)/σ^2
dΓ_dω_fn(z)       = 1 - cdf(Normal(), z)
d2Γ_dω2_fn(z, σ)  = -pdf(Normal(), z)/ω_fn(z, σ)/σ
dG_dσ_fn(z, σ)    = -z*pdf(Normal(), z-σ)/σ
d2G_dωdσ_fn(z, σ) = -pdf(Normal(), z)*(1 - z*(z-σ))/σ^2
dΓ_dσ_fn(z, σ)    = -pdf(Normal(), z-σ)
d2Γ_dωdσ_fn(z, σ) = (z/σ-1)*pdf(Normal(), z)
end # module
