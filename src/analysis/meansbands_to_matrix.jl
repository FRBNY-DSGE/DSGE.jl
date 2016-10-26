# Reformat MeansBands object in Julia into a matrix

"""
meansbands_matrix_all{S<:AbstractString}(m::AbstractModel, input_type::Symbol.
                     output_vars = Vector{Symbol};
                               subset_string = "")

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

class(mb::MeansBands) = mb.metadata[:class]
product(mb::MeansBands) = mb.metadata[:product]
cond_type(mb::MeansBands) = mb.metadata[:cond_type]
para(mb::MeansBands) = mb.metadata[:para]

function Base.show(io::IO, mb::MeansBands)
    @printf io "MeansBands\n"
    @printf io "  class: %s\n" class(mb)
    @printf io "  product: %s\n" product(mb)
    @printf io "  cond: %s\n" cond_type(mb)
    @printf io "  para: %s\n" para(mb)
    @printf io "  dates: %s - %s\n" startdate_means(mb) enddate_means(mb)
    @printf io "  # of variables: %s\n" n_vars_means(mb)
    @printf io "  bands: %s\n" get_density_bands(mb, uniqueify=true)
end


###################################
## MEANS
###################################

"""
```
n_vars_means(mb::MeansBands)
````

Get number of variables (`:y_t`, `:OutputGap`, etc) in `mb.means`
"""
n_vars_means(mb::MeansBands) = length(get_vars_means(mb))

"""
```
get_vars_means(mb::MeansBands)
````

Get variables (`:y_t`, `:OutputGap`, etc) in `mb.means`
"""
get_vars_means(mb::MeansBands) = setdiff(names(mb.means), [:date])


"""
```
n_periods_means(mb::MeansBands)
```

Get number of periods in `mb.means`
"""
n_periods_means(mb::MeansBands) = size(mb.means,1)

"""
```
startdate_means(mb::MeansBands)
```

Get first period in`mb.means`. Assumes `mb.means[product]` is already sorted by date.
"""
startdate_means(mb::MeansBands) = mb.means[:date][1]

"""
```
enddate_means(mb::MeansBands)
```

Get last period for which `mb` stores means. Assumes `mb.means[product]` is already sorted by date.
"""
enddate_means(mb::MeansBands) = mb.means[:date][end]


###################################
## BANDS
###################################

"""
```
n_vars_bands(mb::MeansBands)
```

Get number of variables (`:y_t`, `:OutputGap`, etc) for which `mb`
stores bands for the specified `product` (`hist`, `forecast`, `shockdec`, etc).
"""
n_vars_bands(mb::MeansBands) = length(mb.bands)


"""
```
n_periods_bands(mb::MeansBands)
```

Get number of periods for which `mb` stores bands for the specified
`product` (`hist`, `forecast`, `shockdec`, etc).
"""
function n_periods_bands(mb::MeansBands)
    size(mb.bands[collect(keys(mb.bands))[1]],1)
end

"""
```
startdate_bands(mb::MeansBands)
```

Get first period for which `mb` stores bands. Assumes `mb.bands` is already sorted by date.
"""
startdate_bands(mb::MeansBands) = mb.bands[collect(keys(mb.bands))][:date][1]

"""
```
enddate_bands(mb::MeansBands)
```

Get last period in `mb.bands`. Assumes `mb.bands` is already sorted by date.
"""
enddate_bands(mb::MeansBands) = mb.bands[:date]

"""
```
get_density_bands(mb)
```

Return a list of the bands stored in mb.bands.
"""
function get_density_bands(mb::MeansBands; uniqueify=false, ordered=true)

    # extract one of the keys in mb.bands
    var  = collect(keys(mb.bands))[1]

    # get all the columns in the corresponding dataframe that aren't dates
    strs = map(string,names(mb.bands[var]))
    strs = setdiff(strs, ["date"])

    lowers = strs[map(ismatch, repmat([r"LB"], length(strs)), strs)]
    uppers = strs[map(ismatch, repmat([r"UB"], length(strs)), strs)]

    # sort
    if ordered
        sort!(lowers, rev=true)
        sort!(uppers)
    end

    # return both upper and lower bands, or just percents, as desired
    strs = if uniqueify
        unique([split(x, " ")[1] for x in [lowers; uppers]])
    else
        [lowers; uppers]
    end

    return strs
end