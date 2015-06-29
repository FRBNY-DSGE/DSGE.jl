module FinancialFrictionsFunctions

using Distributions

export zetaspbfcn, zetabomegafcn, zetazomegafcn, nkfcn, mufcn, omegafcn, Gfcn, Gammafcn, dGdomegafcn, d2Gdomega2fcn, dGammadomegafcn, d2Gammadomega2fcn, dGdsigmafcn, d2Gdomegadsigmafcn, dGammadsigmafcn, d2Gammadomegadsigmafcn


# This file contains functions that are used to compute financial frictions steady-state values from parameter values.
# Specifically, it is called by models/m990/parameters.jl

function zetaspbfcn(z, σ, sprd)
    zetaratio = zetabomegafcn(z, σ, sprd)/zetazomegafcn(z, σ, sprd)
    nk = nkfcn(z, σ, sprd)
    return -zetaratio/(1-zetaratio)*nk/(1-nk)
end

function zetabomegafcn(z, σ, sprd)
    nk = nkfcn(z, σ, sprd)
    mustar = mufcn(z, σ, sprd)
    omegastar = omegafcn(z, σ)
    Gammastar = Gammafcn(z, σ)
    Gstar = Gfcn(z, σ)
    dGammadomegastar = dGammadomegafcn(z)
    dGdomegastar = dGdomegafcn(z, σ)
    d2Gammadomega2star = d2Gammadomega2fcn(z, σ)
    d2Gdomega2star = d2Gdomega2fcn(z, σ)
    return omegastar*mustar*nk*(d2Gammadomega2star*dGdomegastar-d2Gdomega2star*dGammadomegastar)/(dGammadomegastar-mustar*dGdomegastar)^2/sprd/(1-Gammastar+dGammadomegastar*(Gammastar-mustar*Gstar)/(dGammadomegastar-mustar*dGdomegastar))
end

function zetazomegafcn(z, σ, sprd)
    mustar = mufcn(z, σ, sprd)
    return omegafcn(z, σ)*(dGammadomegafcn(z) - mustar*dGdomegafcn(z, σ))/(Gammafcn(z, σ) - mustar*Gfcn(z, σ))
end

function nkfcn(z, σ, sprd)
    return 1 - (Gammafcn(z, σ) - mufcn(z, σ, sprd)*Gfcn(z, σ))*sprd
end

function mufcn(z, σ, sprd)
    return (1 - 1/sprd)/(dGdomegafcn(z, σ)/dGammadomegafcn(z)*(1 - Gammafcn(z, σ)) + Gfcn(z, σ))
end

function omegafcn(z, σ)
    return exp(σ*z - σ^2/2)
end

function Gfcn(z, σ)
    return cdf(Normal(), z-σ)
end

function Gammafcn(z, σ)
    return omegafcn(z, σ)*(1 - cdf(Normal(), z)) + cdf(Normal(), z-σ)
end

function dGdomegafcn(z, σ)
    return pdf(Normal(), z)/σ
end

function d2Gdomega2fcn(z, σ)
    return -z*pdf(Normal(), z)/omegafcn(z, σ)/σ^2
end

function dGammadomegafcn(z)
    return 1 - cdf(Normal(), z)
end

function d2Gammadomega2fcn(z, σ)
    return -pdf(Normal(), z)/omegafcn(z, σ)/σ
end

function dGdsigmafcn(z, σ)
    return -z*pdf(Normal(), z-σ)/σ
end

function d2Gdomegadsigmafcn(z, σ)
    return -pdf(Normal(), z)*(1 - z*(z-σ))/σ^2
end

function dGammadsigmafcn(z, σ)
    return -pdf(Normal(), z-σ)
end

function d2Gammadomegadsigmafcn(z, σ)
    return (z/σ-1)*pdf(Normal(), z)
end

end # module
