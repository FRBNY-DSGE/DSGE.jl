"""
meansbands_matrix_all{S<:AbstractString}(m::AbstractModel, input_type::Symbol.
    output_vars = Vector{Symbol}; subset_string = "", verbose::Symbol = :low)

Reformat `MeansBands` object into matrices, and save to individual files.

### Inputs

- `m::AbstractModel`: model object
- `input_type`: same as input into `forecast_all`
- `output_vars`: same as input into `forecast_all`
- `cond_type`: same as input into `forecast_all`

### Keyword Arguments

- `subset_string`: if `input_type == :subset`, this is a subset identified
  string. Otherwise, it should be left empty.
- `verbose::Symbol`: desired frequency of function progress messages printed to
  standard out. One of `:none`, `:low`, or `:high`
"""
function meansbands_matrix_all(m::AbstractModel, input_type::Symbol,
                               output_vars::Vector{Symbol}, cond_type::Symbol;
                               subset_string = "", verbose::Symbol = :low)

    ## Step 0: Determine full set of output_vars necessary for plotting desired results
    #          Specifically, if output_vars contains shockdecs but not
    #          trend or deterministic trends, add those

    output_vars = add_requisite_output_vars(output_vars)
    output_vars = [symbol("mb", x) for x in output_vars]
    output_dir = workpath(m, "forecast", "")
    outfiles = DSGE.get_output_files(m, "forecast", input_type, output_vars, cond_type,
                                     pathfcn = workpath, subset_string = subset_string)
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println()
        info("Converting means and bands to matrices for input_type = $input_type, cond_type = $cond_type...")
        println("Start time: $(now())")
        println("Means and bands matrices will be saved in $output_dir")
    end

    mbs = Dict{Symbol,MeansBands}()
    for input in keys(outfiles)
        fn = outfiles[input]
        mbs[input] = read_mb(fn)
    end

    meansbands_matrix_all(m, mbs; verbose = verbose)

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("\nConversion of means and bands complete: $(now())")
    end
end

function meansbands_matrix_all(m::AbstractModel, mbs::Dict{Symbol,MeansBands};
                               verbose::Symbol = :low)

    for output_var in keys(mbs)

        mb = mbs[output_var]
        prod     = product(mb)                    # history? forecast? shockdec?
        cl       = class(mb)                      # pseudo? obs? state?
        condtype = cond_type(mb)

        ## Get name of file to write
        outfile = DSGE.get_output_files(m, "forecast", para(mb), [output_var], cond_type(mb),
                    pathfcn = workpath, subset_string = mb.metadata[:subset_string],
                    fileformat = :h5)[output_var]

        base = basename(outfile)[3:end]
        outfile = joinpath(dirname(outfile), "mb_matrix_"*base)

        ## Convert MeansBands objects to matrices
        print("* Extracting means and bands matrices for $prod...")
        means, bands = meansbands_matrix(mb)

        ## Save to file
        h5open(outfile, "w") do f
            f["means"] = means
            f["bands"] = bands
        end
        println("wrote matrix-form means and bands for ($prod$cl, $condtype) to $outfile\n")
    end
end

"""
```
meansbands_matrix{S<:AbstractString}(mb::MeansBands)
```

Convert a `MeansBands` object to matrix form.
"""
function meansbands_matrix(mb::MeansBands)

    # extract useful info
    vars     = get_vars_means(mb)             # get names of variables
    nvars    = length(vars)                   # get number of variables

    tmp      = setdiff(names(mb.means), [:date])[1]
    T        = eltype(mb.means[tmp])          # type of elements of mb structure

    nperiods = n_periods_means(mb)            # get number of periods
    prod     = product(mb)                    # history? forecast? shockdec?
    inds     = mb.metadata[:indices]          # mapping from names of vars to indices
    bands_list = get_density_bands(mb)        # which bands are stored
    nbands   = length(bands_list)             # how many bands are stored

    # extract  matrices from MeansBands structure
    means, bands = if prod in [:hist, :forecast, :dettrend]

        # construct means and bands arrays
        means = Array{T,2}(nvars, nperiods)
        bands = Array{T}(nbands, nvars, nperiods)

        # extract each series and place in arrays
        for series in setdiff(names(mb.means), [:date])
            ind = inds[series]
            means[ind,:] = convert(Array{T},mb.means[series])

            for (i,band) in enumerate(bands_list)  # these are ordered properly already
                bands[i,ind,:] = convert(Array{T},mb.bands[series][symbol(band)])
            end
        end

        means, bands

    elseif prod in [:shockdec]

        shock_inds = mb.metadata[:shock_indices]
        nshocks = length(shock_inds)

        # construct means and bands arrays
        means = Array{T}(nvars, nperiods, nshocks)
        bands = Array{T}(nbands, nvars, nperiods, nshocks)

        for series in setdiff(names(mb.means), [:date])

            var, shock = DSGE.parse_mb_colname(series)
            ind = inds[var]

            for shock in keys(shock_inds)
                shock_ind = shock_inds[shock]
                means[ind,:,shock_ind] = convert(Array{T}, mb.means[series])

                for (band_ind, band) in enumerate(bands_list)
                    bands[band_ind, ind, :, shock_ind] =
                        convert(Array{T}, mb.bands[series][symbol(band)])
                end
            end
        end

        means, bands

    elseif prod in [:trend]

        means = convert(Matrix{Float64}, mb.means)

        bands = Array{T}(nbands, nvars)
        for var in keys(mb.bands)
            ind = inds[var]
            bands[:, ind] = convert(Matrix{T}, mb.bands[var])
        end

        means, bands
    end

    # return matrix
    means, bands
end

