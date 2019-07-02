# functions that are used to compute financial frictions
# steady-state values from parameter values
@inline function ζ_spb_fn(z, σ, spr)
    zetaratio = ζ_bω_fn(z, σ, spr)/ζ_zω_fn(z, σ, spr)
    nk = nk_fn(z, σ, spr)
    return -zetaratio/(1-zetaratio)*nk/(1-nk)
end

@inline function ζ_bω_fn(z, σ, spr)
    nk          = nk_fn(z, σ, spr)
    μstar       = μ_fn(z, σ, spr)
    ω_star       = ω_fn(z, σ)
    Γstar       = Γ_fn(z, σ)
    Gstar       = G_fn(z, σ)
    dΓ_dω_star   = dΓ_dω_fn(z)
    dG_dω_star   = dG_dω_fn(z, σ)
    d2Γ_dω2star = d2Γ_dω2_fn(z, σ)
    d2G_dω2star = d2G_dω2_fn(z, σ)
    return ω_star*μstar*nk*(d2Γ_dω2star*dG_dω_star - d2G_dω2star*dΓ_dω_star)/
        (dΓ_dω_star - μstar*dG_dω_star)^2/spr/(1 - Γstar + dΓ_dω_star*(Γstar - μstar*Gstar)/
            (dΓ_dω_star - μstar*dG_dω_star))
end

@inline function ζ_zω_fn(z, σ, spr)
    μstar = μ_fn(z, σ, spr)
    return ω_fn(z, σ)*(dΓ_dω_fn(z) - μstar*dG_dω_fn(z, σ))/
        (Γ_fn(z, σ) - μstar*G_fn(z, σ))
end

nk_fn(z, σ, spr) = 1 - (Γ_fn(z, σ) - μ_fn(z, σ, spr)*G_fn(z, σ))*spr
μ_fn(z, σ, spr)  =
    (1 - 1/spr)/(dG_dω_fn(z, σ)/dΓ_dω_fn(z)*(1 - Γ_fn(z, σ)) + G_fn(z, σ))
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
