"""
meansbands_matrix_all{S<:AbstractString}(m::AbstractModel, input_type::Symbol.
                     output_vars = Vector{Symbol};
                               subset_string = "")

Reformat MeansBands object in Julia into a matrix.

### Inputs
- `m::AbstractModel`: model object
- `output_vars`: same as input into `forecast_all`
- `subset_string`: if `input_type==:subset`, this is a subset
  identified string. Otherwise, it should be left empty.
"""
function meansbands_matrix_all(m::AbstractModel, input_var::Symbol,
                               output_vars::Vector{Symbol}, cond_type::Symbol;
                               subset_string = "")

    output_vars = [symbol("mb", x) for x in output_vars]
    outfiles = DSGE.get_output_files(m, "forecast", input_var, output_vars, cond_type,
                                     pathfcn = workpath, subset_string = subset_string)

    mbs = Dict{Symbol,MeansBands}()
    for input in keys(outfiles)

        fn = outfiles[input]
        mbs[input] = jldopen(fn, "r") do f
            read(f, "mb")
        end
    end


    meansbands_matrix_all(m, mbs)
end

function meansbands_matrix_all(m::AbstractModel, mbs::Dict{Symbol,MeansBands})

    for mb in values(mbs)
        meansbands_matrix(m, mb)
    end

end

function meansbands_matrix(m::AbstractModel, mb::MeansBands)

    output_var = symbol("mb_matrix_", product(mb), class(mb))
    outfile = DSGE.get_output_files(m, "forecast", para(mb), [output_var], cond_type(mb),
                                    pathfcn = workpath, subset_string = mb.metadata[:subset_string],
                                    fileformat = :h5)

    meansbands_matrix(mb, outfile[output_var])
end

function meansbands_matrix{S<:AbstractString}(mb::MeansBands, outfile::S)

    # extract useful info
    vars     = get_vars_means(mb)             # get names of variables
    nvars    = length(vars)                   # get number of variables
    T        = eltype(mb.means[vars[1]])      # type of elements of mb structure
    nperiods = n_periods_means(mb)            # get number of periods
    prod     = product(mb)                    # history? forecast? shockdec?
    cl       = class(mb)                      # pseudo? obs? state?
    inds     = mb.metadata[:indices]          # mapping from names of vars to indices
    bands_list = get_density_bands(mb)        # which bands are stored
    nbands   = length(bands_list)             # how many bands are stored
    condtype = cond_type(mb)

    print("* Extracting means and bands matrices for $prod...")

    # extract  matrices from MeansBands structure
    means, bands = if prod in [:hist, :forecast]

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
        # TODO
        shock_inds = mb.metadata[:shock_inds]
        nshocks = length(shock_inds)

        # construct means and bands arrays
        means = Array{T,3}(nvars, nperiods, nshocks)
        bands = Array{T,4}(nbands, nvars, nperiods, nshocks)

        for series in setdiff(names(mb.means), [:date])
            ind = inds[series]

            for shock in keys(shock_inds)
                means[ind,:,shock] = mb.means[series]
            end
        end

        means, bands
    end

    # write to file

    h5open(outfile, "w") do f
        f["means"] = means
        f["bands"] = bands
    end

    println("wrote matrix-form means and bands for ($prod$cl, $condtype) to $outfile\n")

    # return matrix
    means
end

