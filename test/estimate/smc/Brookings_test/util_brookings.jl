"""
```
standard_spec!(m, vint, directory_path; cdvt = vint,
    dsid = data_id(m), cdid = cond_id(m), fcast-date = Dates.lastdayofquarter(Dates.today()),
    fcast_end = DSGE.quartertodate("2030-Q1"))
```
adds the standard settings to the model object `m` for
the Brookings Phillips Curve project.
"""
function standard_spec!(m::AbstractDSGEModel, vint::String,
                        directory_path::String = dirname(@__FILE__); cdvt::String = vint,
                        dsid::Int = data_id(m),
                        cdid::Int = cond_id(m),
                        fcast_date::Dates.Date = Dates.lastdayofquarter(Dates.today()),
                        fcast_end::Date = DSGE.quartertodate("2024-Q4"),
                        four_folders_down::Bool = false,
                        save_orig::Bool = false)
    m <= Setting(:data_vintage, vint)
    if four_folders_down
        m <= Setting(:saveroot, joinpath(directory_path, "../../../../save/"))
        m <= Setting(:dataroot, joinpath(directory_path, "../../../../save/input_data"))
    else
        m <= Setting(:saveroot, joinpath(directory_path, "../../../save/"))
        m <= Setting(:dataroot, joinpath(directory_path, "../../../save/input_data"))
    end
    if save_orig
        if four_folders_down
            m <= Setting(:saveroot, joinpath(directory_path, "../../../../save_orig/"))
            m <= Setting(:dataroot, joinpath(directory_path, "../../../../save/input_data"))
        else
            m <= Setting(:saveroot, joinpath(directory_path, "../../../save_orig/"))
            m <= Setting(:dataroot, joinpath(directory_path, "../../../save/input_data"))
        end
    end
    m <= Setting(:data_id, dsid)
    m <= Setting(:cond_id, cdid)
    m <= Setting(:hours_first_observable, true)

    m <= Setting(:cond_vintage, cdvt)
    m <= Setting(:date_forecast_start,  fcast_date)
    m <= Setting(:forecast_horizons, subtract_quarters(fcast_end, fcast_date) + 1)
    m <= Setting(:impulse_response_horizons, 20)
    m <= Setting(:date_conditional_end, fcast_date)
    m <= Setting(:use_population_forecast, true)
    m <= Setting(:forecast_uncertainty_override, Nullable{Bool}(false))

    m <= Setting(:sampling_method, :SMC)
    m <= Setting(:forecast_block_size, 500)
    m <= Setting(:forecast_jstep, 1)

    if cdid == 1
        m <= Setting(:cond_semi_names, [])
        m <= Setting(:cond_full_names, [:obs_hours])
    end
end

function standard_spec200129!(m::AbstractDSGEModel, vint::String,
                        directory_path::String = dirname(@__FILE__); cdvt::String = vint,
                        dsid::Int = data_id(m),
                        cdid::Int = cond_id(m),
                        fcast_date::Dates.Date = Dates.lastdayofquarter(Dates.today()),
                        fcast_end::Date = DSGE.quartertodate("2024-Q4"),
                        four_folders_down::Bool = false)
    m <= Setting(:data_vintage, vint)
    if four_folders_down
        m <= Setting(:saveroot, joinpath(directory_path, "../../../../save/200129"))
        m <= Setting(:dataroot, joinpath(directory_path, "../../../../save/input_data"))
    else
        m <= Setting(:saveroot, joinpath(directory_path, "../../../save/200129"))
        m <= Setting(:dataroot, joinpath(directory_path, "../../../save/input_data"))
    end
    m <= Setting(:data_id, dsid)
    m <= Setting(:cond_id, cdid)
    m <= Setting(:hours_first_observable, true)

    m <= Setting(:cond_vintage, cdvt)
    m <= Setting(:date_forecast_start,  fcast_date)
    m <= Setting(:forecast_horizons, subtract_quarters(fcast_end, fcast_date) + 1)
    m <= Setting(:impulse_response_horizons, 20)
    m <= Setting(:date_conditional_end, fcast_date)
    m <= Setting(:use_population_forecast, true)
    m <= Setting(:forecast_uncertainty_override, Nullable{Bool}(false))

    m <= Setting(:sampling_method, :SMC)
    m <= Setting(:forecast_block_size, 500)
    m <= Setting(:forecast_jstep, 1)

    m <= Setting(:date_regime2_start_text, "900331", true, "reg2start",
                 "The text version to be saved of when regime 2 starts")
    m <= Setting(:preZLB, "false", true, "preZLB", "")

    if cdid == 1
        m <= Setting(:cond_semi_names, [])
        m <= Setting(:cond_full_names, [:obs_hours])
    end
end

function standard_spec200204_dsgevar!(m::AbstractDSGEModel, vint::String,
                        directory_path::String = dirname(@__FILE__); cdvt::String = vint,
                        dsid::Int = data_id(m),
                        cdid::Int = cond_id(m),
                        fcast_date::Dates.Date = Dates.lastdayofquarter(Dates.today()),
                        fcast_end::Date = DSGE.quartertodate("2024-Q4"),
                        four_folders_down::Bool = false)
    m <= Setting(:data_vintage, vint)
    if four_folders_down
        m <= Setting(:saveroot, joinpath(directory_path, "../../../../save/"))
        m <= Setting(:dataroot, joinpath(directory_path, "../../../../save/input_data"))
    else
        m <= Setting(:saveroot, joinpath(directory_path, "../../../save/"))
        m <= Setting(:dataroot, joinpath(directory_path, "../../../save/input_data"))
    end
    m <= Setting(:data_id, dsid)
    m <= Setting(:cond_id, cdid)
    m <= Setting(:hours_first_observable, false)

    m <= Setting(:cond_vintage, cdvt)
    m <= Setting(:date_forecast_start,  fcast_date)
    m <= Setting(:forecast_horizons, subtract_quarters(fcast_end, fcast_date) + 1)
    m <= Setting(:impulse_response_horizons, 20)
    m <= Setting(:date_conditional_end, fcast_date)
    m <= Setting(:use_population_forecast, false)
    m <= Setting(:forecast_uncertainty_override, Nullable{Bool}(false))

    m <= Setting(:sampling_method, :SMC)
    m <= Setting(:forecast_block_size, 500)
    m <= Setting(:forecast_jstep, 1)

    m <= Setting(:date_regime2_start_text, "900331", true, "reg2start",
                 "The text version to be saved of when regime 2 starts")
    m <= Setting(:preZLB, "false", true, "preZLB", "")
    m <= Setting(:dsgevar, "true", true, "dsgevar", "indicates output from DSGEVAR estimation")
end

"""
```
get_κp(m::AbstractDSGEModel, θ::Matrix{Float64})
```
computes the slope of the price Phillips curve, given a
`n_draws × n_parameters` matrix.
"""
function get_κp(m::AbstractDSGEModel, θ::Matrix{Float64})
    κ = Vector{Float64}(undef, size(θ, 1))
    for i = 1:size(θ, 1)
        κ[i] = get_κp(m, θ[i, :])
    end
    return κ
end

function get_κp(m::AbstractDSGEModel, θ::Vector{Float64})
    DSGE.update!(m, θ)
    return ((1 - m[:ζ_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))*
            (1 - m[:ζ_p]))/(m[:ζ_p]*((m[:Φ]- 1)*m[:ϵ_p] + 1))/(1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
end

"""
```
get_κw(m::AbstractDSGEModel, θ::Matrix{Float64})
```
computes the slope of the wage Phillips curve, given a
`n_draws × n_parameters` matrix.
"""
function get_κw(m::AbstractDSGEModel, θ::Matrix{Float64})
    κ = Vector{Float64}(undef, size(θ, 1))
    for i = 1:size(θ, 1)
        κ[i] = get_κw(m, θ[i, :])
    end
    return κ
end

function get_κw(m::AbstractDSGEModel, θ::Vector{Float64})
    DSGE.update!(m, θ)
    return (1 - m[:ζ_w]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))*
    (1 - m[:ζ_w])/(m[:ζ_w]*((m[:λ_w] - 1)*m[:ϵ_w] + 1))/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
end

"""
```
draw_κ(m::AbstractDSGEModel; N::Int = 100000,
       kernel_density::Bool = true) where {S <: Real}
```
draws from the implied prior for κp and κw, given priors
for the parameters which specify these slopes (e.g. ζ_p).
"""
function draw_κ(m::AbstractDSGEModel; N::Int = 100000,
                kernel_density::Bool = true) where {S <: Real}
    paras = zeros(N, length(m.parameters))
    for i = 1:N
        success = false
        while !success
            try
                paras[i, :] = rand(m.parameters, 1)
            catch e
                continue
            end
            success = true
        end
    end
    κp = get_κp(m, paras)
    κw = get_κw(m, paras)

    if kernel_density
        kdep = kde(κp)
        kdew = kde(κw)

        return kdep, kdew
    else

        return κp, κw
    end
end
