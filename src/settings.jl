"""
```
Setting{T}
```

The `Setting` type is an interface for computational settings that affect how the code runs
without affecting the mathematical definition of the model. It also provides support for
non-conflicting file names for output of 2 models that differ only in the values of their
computational settings.

### Fields
- `key::Symbol`: Name of setting
- `value::T`: Value of setting
- `print::Bool`: Indicates whether to append this setting's code and value to output file
  names. If true, output file names will include a suffix of the form _code1=val1_code2=val2
  etc. where codes are listed in alphabetical order.
- `code::AbstractString`: string of <=4 characters to print to output file suffixes when
  `print=true`.
- `description::AbstractString`: Short description of what the setting is used for.
"""
type Setting{T}
    key::Symbol                  # name of setting
    value::T                     # whatever the setting is
    print::Bool                  # whether or not to add this setting to the print
    code::AbstractString         # what gets printed to the print
    description::AbstractString  # description of what the setting is for
end

# for printing codes to filename string
Base.convert{T<:Number, U<:Number}(::Type{T}, s::Setting{U}) = convert(T, s.value)
Base.convert{T<:AbstractString, U<:AbstractString}(::Type{T}, s::Setting{U}) = convert(T, s.value)

Base.promote_rule{T<:Number,U<:Number}(::Type{Setting{T}}, ::Type{U}) = promote_rule(T,U)
Base.promote_rule{T<:AbstractString,U<:AbstractString}(::Type{Setting{T}}, ::Type{U}) = promote_rule(T,U)
Base.promote_rule(::Type{Setting{Bool}}, ::Type{Bool}) = promote_rule(Bool, Bool)

Base.string(s::Setting{AbstractString}) = string(s.value)

to_filestring(s::Setting) = "$(s.code)=$(s.value)"

# key, value constructor
Setting(key, value) = Setting(key, value, false, "", "")
Setting(key, value, description) = Setting(key, value, false, "", description)

function Base.show(io::IO, s::Setting)
    @printf io "%s: %s" s.key s.value
end

"""
```
(<=)(m::AbstractModel, s::Setting)
```
Syntax for adding a setting to a model/overwriting a setting via `m <= Setting(...)`
"""
function (<=)(m::AbstractModel, s::Setting)
    if !m.testing
        setting_field_name = :settings
    else
        setting_field_name = :test_settings
    end

    if !haskey(getfield(m, setting_field_name), s.key)
        getfield(m, setting_field_name)[s.key] = s
    else
        update!(getfield(m, setting_field_name)[s.key], s)
    end
end

function update!(a::Setting, b::Setting)
    # Make sure Setting a can appropriately be updated by b.
    if a.key ≠ b.key
        return
    end

    a.value = b.value

    # b overrides the print boolean of a if:
    # - a.print is false and b.print is true.
    # - a.print is true, b.print is false, and b.code is non-empty.
    # We must be able to tell if b was created via a constructor like `Setting(:key,
    # value)`, in which case the print, code, and description values are set to defaults. We
    # do not overwrite if we can't determine whether or not those fields are just the
    # defaults.
    if (a.print == false && b.print == true) ||
       (a.print == true && b.print == false && !isempty(b.code))
       a.print = b.print
   end
    if !isempty(b.code) && b.code ≠ a.code
        a.code = b.code
    end
    if !isempty(b.description) && b.description ≠ a.description
        a.description = b.description
    end
end

"""
```
get_setting(m::AbstractModel, setting::Symbol)
```

Returns the value of the setting
"""
function get_setting(m::AbstractModel, s::Symbol)

    if m.testing && in(s, keys(m.test_settings))
        return m.test_settings[s].value
    end

    return m.settings[s].value
end

"""
```
default_test_settings!(m::AbstractModel)
```

The following Settings are constructed, initialized and added to
`m.test_settings`. Their purposes are identical to those in
`m.settings`, but these values are used to test DSGE.jl.

### I/O Locations and identifiers
- `saveroot::Setting{String}`: A temporary directory in /tmp/
- `dataroot::Setting{String}`: dsgeroot/test/reference/
- `data_vintage::Setting{String}`: "_REF"

### Metropolis-Hastings
- `n_mh_simulations::Setting{Int}`: 100
- `n_mh_blocks::Setting{Int}`: 1
- `n_mh_burn::Setting{Int}`: 0
- `mh_thin::Setting{Int}`: 1
"""
function default_test_settings!(m::AbstractModel)

    test = m.test_settings

    # I/O
    dataroot = normpath(joinpath(dirname(@__FILE__), "..","test","reference"))
    saveroot = mktempdir()

    #General
    test[:saveroot] = Setting(:saveroot, saveroot,
        "Where to write files when in test mode")
    test[:dataroot] = Setting(:dataroot, dataroot,
        "Location of input files when in test mode" )
    test[:data_vintage] = Setting(:data_vintage, "REF", true, "vint",
        "Reference data identifier")
    test[:date_mainsample_end] = Setting(:date_mainsample_end, quartertodate("2015-Q3"),
        "End date of main sample")
    test[:use_parallel_workers] = Setting(:use_parallel_workers, false, false, "parw",
        "Use available parallel workers in computations")
    test[:n_hessian_test_params] = Setting(:n_hessian_test_params, 3, false, "mhfp",
        "Max number of free params for which to calculate Hessian")

    # Metropolis-Hastings
    test[:n_mh_simulations] = Setting(:n_mh_simulations, 100, false, "nsim",
        "Number of parameter draws per block for testing Metropolis-Hastings")
    test[:n_mh_blocks] = Setting(:n_mh_blocks, 1, false, "nblc",
        "Number of blocks to draw parameters for testing Metropolis-Hastings")
    test[:n_mh_burn] = Setting(:n_mh_burn, 0, false, "nbrn",
        "Number of burn-in blocks for testing Metropolis-Hastings")
    test[:mh_thin] = Setting(:mh_thin, 1, false, "thin",
        "Thinning step for testing Metropolis-Hastings")

    return test
end
