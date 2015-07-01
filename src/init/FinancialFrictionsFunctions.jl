module FinancialFrictionsFunctions

using Distributions

export ζ_spb_fn, ζ_bω_fn, ζ_zω_fn, nk_fn, μ_fn, ω_fn, G_fn, Γ_fn, dG_dω_fn, d2G_dω2_fn, dΓ_dω_fn, d2Γ_dω2_fn, dG_dσ_fn, d2G_dωdσ_fn, dΓ_dσ_fn, d2Γ_dωdσ_fn


# This file contains functions that are used to compute financial frictions steady-state values from parameter values.
# Specifically, it is called by models/m990/parameters.jl

function ζ_spb_fn(z, σ, sprd)
    zetaratio = ζ_bω_fn(z, σ, sprd)/ζ_zω_fn(z, σ, sprd)
    nk = nk_fn(z, σ, sprd)
    return -zetaratio/(1-zetaratio)*nk/(1-nk)
end

function ζ_bω_fn(z, σ, sprd)
    nk = nk_fn(z, σ, sprd)
    mustar = μ_fn(z, σ, sprd)
    omegastar = ω_fn(z, σ)
    Gammastar = Γ_fn(z, σ)
    Gstar = G_fn(z, σ)
    dGammadomegastar = dΓ_dω_fn(z)
    dGdomegastar = dG_dω_fn(z, σ)
    d2Gammadomega2star = d2Γ_dω2_fn(z, σ)
    d2Gdomega2star = d2G_dω2_fn(z, σ)
    return omegastar*mustar*nk*(d2Gammadomega2star*dGdomegastar-d2Gdomega2star*dGammadomegastar)/(dGammadomegastar-mustar*dGdomegastar)^2/sprd/(1-Gammastar+dGammadomegastar*(Gammastar-mustar*Gstar)/(dGammadomegastar-mustar*dGdomegastar))
end

function ζ_zω_fn(z, σ, sprd)
    mustar = μ_fn(z, σ, sprd)
    return ω_fn(z, σ)*(dΓ_dω_fn(z) - mustar*dG_dω_fn(z, σ))/(Γ_fn(z, σ) - mustar*G_fn(z, σ))
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
