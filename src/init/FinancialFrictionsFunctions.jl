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
    nk = nk_fn(z, σ, sprd)
    μstar = μ_fn(z, σ, sprd)
    ωstar = ω_fn(z, σ)
    Γstar = Γ_fn(z, σ)
    Gstar = G_fn(z, σ)
    dΓ_dωstar = dΓ_dω_fn(z)
    dG_dωstar = dG_dω_fn(z, σ)
    d2Γ_dω2star = d2Γ_dω2_fn(z, σ)
    d2G_dω2star = d2G_dω2_fn(z, σ)
    return ωstar*μstar*nk*(d2Γ_dω2star*dG_dωstar - d2G_dω2star*dΓ_dωstar)/(dΓ_dωstar - μstar*dG_dωstar)^2/sprd/(1 - Γstar + dΓ_dωstar*(Γstar - μstar*Gstar)/(dΓ_dωstar - μstar*dG_dωstar))
end

function ζ_zω_fn(z, σ, sprd)
    μstar = μ_fn(z, σ, sprd)
    return ω_fn(z, σ)*(dΓ_dω_fn(z) - μstar*dG_dω_fn(z, σ))/(Γ_fn(z, σ) - μstar*G_fn(z, σ))
end

function nk_fn(z, σ, sprd)
    return 1 - (Γ_fn(z, σ) - μ_fn(z, σ, sprd)*G_fn(z, σ))*sprd
end

function μ_fn(z, σ, sprd)
    return (1 - 1/sprd)/(dG_dω_fn(z, σ)/dΓ_dω_fn(z)*(1 - Γ_fn(z, σ)) + G_fn(z, σ))
end

function ω_fn(z, σ)
    return exp(σ*z - σ^2/2)
end

function G_fn(z, σ)
    return cdf(Normal(), z-σ)
end

function Γ_fn(z, σ)
    return ω_fn(z, σ)*(1 - cdf(Normal(), z)) + cdf(Normal(), z-σ)
end

function dG_dω_fn(z, σ)
    return pdf(Normal(), z)/σ
end

function d2G_dω2_fn(z, σ)
    return -z*pdf(Normal(), z)/ω_fn(z, σ)/σ^2
end

function dΓ_dω_fn(z)
    return 1 - cdf(Normal(), z)
end

function d2Γ_dω2_fn(z, σ)
    return -pdf(Normal(), z)/ω_fn(z, σ)/σ
end

function dG_dσ_fn(z, σ)
    return -z*pdf(Normal(), z-σ)/σ
end

function d2G_dωdσ_fn(z, σ)
    return -pdf(Normal(), z)*(1 - z*(z-σ))/σ^2
end

function dΓ_dσ_fn(z, σ)
    return -pdf(Normal(), z-σ)
end

function d2Γ_dωdσ_fn(z, σ)
    return (z/σ-1)*pdf(Normal(), z)
end

end # module
