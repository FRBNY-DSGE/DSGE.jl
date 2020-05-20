# This script holds measurement error, which figures
# out what EE and MM matrices to return based on the subspec of the model
function measurement_error(m::DSGEVECM{T}) where {T <: Real}
    ss = subspec(m)
    if ss == "no subspecs with measurement error yet"
    else
        EE = zeros(T, n_observables(m), n_observables(m))
        MM = zeros(T, n_observables(m), n_shocks(m))
    end

    return EE, MM
end
