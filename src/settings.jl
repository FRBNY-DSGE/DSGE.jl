"""
Setting{T<:Any}
- `key::Symbol`: Name of setting
- `value::T`: Value of setting 
- `savestring::Bool`: Indicates whether to append this setting's code
and value to output file names. If true, output file names will
include a suffix of the form _code1=val1_code2=val2_etc. where codes
are listed in alphabetical order.
- `code::AbstractString`: string of <=4 characters to print to output
file suffixes when `savestring` is true.
- `description::AbstractString`: Short description of what the setting
is used for.
"""
immutable Setting{T} 
    key::Symbol                  # name of setting
    value::T                     # whatever the setting is
    savestring::Bool             # whether or not to add this setting to the savestring
    code::AbstractString         # what gets printed to the savestring
    description::AbstractString  # description of what the setting is for
end

# for printing codes to filename string 
Base.convert{T<:Number, U<:Number}(::Type{T}, s::Setting{U}) = convert(T, s.value)
Base.convert{T<:AbstractString, U<:AbstractString}(::Type{T}, s::Setting{U}) = convert(T, s.value)

Base.promote_rule{T<:Number,U<:Number}(::Type{Setting{T}}, ::Type{U}) = promote_rule(T,U)
Base.promote_rule{T<:AbstractString,U<:AbstractString}(::Type{Setting{T}}, ::Type{U}) = promote_rule(T,U)
Base.promote_rule(::Type{Setting{Bool}}, ::Type{Bool}) = promote_rule(Bool, Bool)

Base.string(s::Setting{AbstractString}) = string(s.value)

to_filename(s::Setting) = "$(s.code)=$(s.value)"


# key, value constructor
function Setting(key, value)
    Setting(key, value, false, "", "")
end

function Setting(key, value, description)
    Setting(key, value, false, "", description)
end


"""
`(<=){T}(m::AbstractDSGEModel{T}, s::Setting)`

Syntax for adding a setting to a model/overwriting a setting: m <= setting
"""
function (<=){T}(m::AbstractDSGEModel{T}, s::Setting)
    if s.savestring 
        # Add to a sorted dictionary of things to print
        insert!(m._filestrings, s.key, to_filename(s))
    end

    m.settings[s.key] = s
end


"""
`get_setting(m::AbstractDSGEModel, setting::Symbol)`

Returns the value of the setting
"""
function get_setting(m::AbstractDSGEModel, s::Symbol)

    if m.testing && in(s, keys(m.test_settings))
        return m.test_settings[s].value
    end

    return m.settings[s].value
end


"""
`default_settings(m::AbstractDSGEModel)`

The following Settings are constructed, initialized and added to
`m.settings`:

### I/O 

- `datapathroot::Setting{ASCIIString}`: The root directory for
  model input data.
- `savepathroot::Setting{ASCIIString}`: The root directory for model output.
- `data_vintage`::Setting{ASCIIString}`: Data vintage identifier,
  formatted YYMMDD (e.g. data from October 30, 2015 is identified by
  the string "151030".)

### Anticipated Shocks
- `num_anticipated_shocks::Setting{Int}`: Number of anticipated policy shocks.
- `num_anticipated_shocks_padding::Setting{Int}`: Padding for
  `num_anticipated_shocks`.
- `num_anticipated_lags::Setting{Int}`: Number of periods back to
incorporate zero bound expectations.
- `num_presample_periods::Setting{Int}`: Number of periods in the
  presample

### Estimation 

- `reoptimize::Setting{Bool}`: Whether to reoptimize the posterior
mode. If `false` (the default), `estimate()` reads in a previously
found mode.
- `recalculate_hessian::Setting{Bool}`: Whether to reecalculate the
hessian at the mode. If `false` (the default), `estimate()` reads in
a previously computed Hessian.

#### Metropolis-Hastings 

- `num_mh_simulations::Setting{Int}`: Number of draws from the
posterior distribution per block.
- `num_mh_blocks::Setting{Int}`: Number of blocks to run
Metropolis-Hastings.
- `num_mh_burn::Setting{Int}`: Number of blocks to discard as burn-in
for Metropolis-Hastings
- `mh_thinning_step::Setting{Int}`: Save every `mh_thinning_step`-th
draw in Metropolis-Hastings.
"""
function default_settings(m::AbstractDSGEModel)

    # I/O File locations
    modelpath = normpath(joinpath(dirname(@__FILE__), "..","save","output_data",spec(m),subspec(m)))
    datapath = normpath(joinpath(dirname(@__FILE__), "..","save","input_data"))

    
    m <= Setting(:modelpathroot, modelpath, "Root of data directory structure")
    m <= Setting(:datapathroot, datapath, "Input data directory path")

    # Anticipated shocks
    m <= Setting(:num_anticipated_shocks,         6, "Number of anticipated policy shocks")
    m <= Setting(:num_anticipated_shocks_padding, 20, "Padding for anticipated policy shocks")
    m <= Setting(:num_anticipated_lags,  24, "Number of periods back to incorporate zero bound expectations")

    # TODO: should be set when data are read in
    m <= Setting(:num_presample_periods, 2, "Number of periods in the presample")

    # General computation
    m <= Setting(:use_parallel_workers, true, "Use available parallel workers in computations")

    # Estimation
    m <= Setting(:reoptimize,          false, "Reoptimize the posterior mode")
    m <= Setting(:recalculate_hessian, false, "Recalculate the hessian at the mode")
    m <= Setting(:max_hessian_free_params, typemax(Int), "Max number of free params for which to calculate Hessian")
    m <= Setting(:num_mh_simulations,  10000, "Number of draws per block in Metropolis-Hastings")
    m <= Setting(:num_mh_blocks,       22   , "Number of blocks for Metropolis-Hastings")
    m <= Setting(:num_mh_burn,         2    , "Number of blocks to use as burn-in in Metropolis-Hastings")
    m <= Setting(:mh_thinning_step,    5    , "How often to write draw to file in Metropolis-Hastings")

    # Data vintage. Default behavior: choose the most recent data file
    input_files = readdir(inpath(m, "data", "")) 
    vint = 0
    for file in input_files
        if ismatch(r"^\s*data*", file)
            regmatch = match(r"^\s*data_(?P<vint>\d{6})", file)
            vint = max(vint, parse(Int,regmatch[:vint]))
        end
    end
    vint = "$vint"
    m <= Setting(:data_vintage, vint, true, "vint", "Date of data")
    
end


"""
`default_test_settings(m::AbstractDSGEModel)`

The following Settings are constructed, initialized and added to
`m.test_settings`. Their purposes are identical to those in
`m.settings`, but these values are used to test DSGE.jl.


### I/O Locations and identifiers
- `modelpathroot::Setting{ASCIIString}`: A temporary directory in /tmp/
- `datapathroot::Setting{ASCIIString}`: DSGEroot/test/reference/
- `data_vintage::Setting{ASCIIString}`: "_REF"

### Metropolis-Hastings 
- `num_mh_simulations::Setting{Int}`: 100 
- `num_mh_blocks::Setting{Int}`: 1
- `num_mh_burn::Setting{Int}`: 0 
- `mh_thinning_step::Setting{Int}`: 1
"""
function default_test_settings(m::AbstractDSGEModel)
    
    test = m.test_settings

    # I/O
    datapathroot = normpath(joinpath(dirname(@__FILE__), "..","test","reference"))
    modelpathroot = mktempdir()
    
    test[:modelpathroot] = Setting(:modelpathroot, modelpathroot,
                                       "Where to write files when in test mode")

    test[:datapathroot] = Setting(:datapathroot, datapathroot,
                                       "Location of input files when in test mode" )

    test[:data_vintage] = Setting(:data_vintage, "REF", true, "vint", "Reference data identifier")

    test[:use_parallel_workers] = Setting(:use_parallel_workers, false, false, "parw", 
                                            "Use available parallel workers in computations")

    test[:max_hessian_free_params] = Setting(:max_hessian_free_params, 3, false, "mhfp",
                                            "Max number of free params for which to calculate Hessian")
    
    # Metropolis-Hastings
    test[:num_mh_simulations] = Setting(:num_mh_simulations, 100, false, "nsim",
                                        "Number of parameter draws per block for testing Metropolis-Hastings") 
    
    test[:num_mh_blocks]      = Setting(:num_mh_blocks, 1, false, "nblc",
                                        "Number of blocks to draw parameters for testing Metropolis-Hastings")
    
    test[:num_mh_burn]        = Setting(:num_mh_burn,   0, false, "nbrn",
                                        "Number of burn-in blocks for testing Metropolis-Hastings")
    
    test[:mh_thinning_step]   = Setting(:mh_thinning_step, 1, false, "thin",
                                        "Thinning step for testing Metropolis-Hastings")
end
