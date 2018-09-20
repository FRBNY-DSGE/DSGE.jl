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
  names. If true, output file names will include a suffix of the form `_code1=val1_code2=val2`
  etc. where codes are listed in alphabetical order.
- `code::String`: string of <=4 characters to print to output file suffixes when
  `print=true`.
- `description::String`: Short description of what the setting is used for.
"""
mutable struct Setting{T}
    key::Symbol                  # name of setting
    value::T                     # whatever the setting is
    print::Bool                  # whether or not to add this setting to the print
    code::String         # what gets printed to the print
    description::String  # description of what the setting is for
end

# for printing codes to filename string
Base.convert(::Type{T}, s::Setting{U}) where {T<:Number, U<:Number} = convert(T, s.value)
Base.convert(::Type{String}, s::Setting{String}) = convert(String, s.value)

Base.promote_rule(::Type{Setting{T}}, ::Type{U}) where {T<:Number,U<:Number} = promote_rule(T,U)
Base.promote_rule(::Type{Setting{String}}, ::Type{String}) = String # promote_rule(String, String)
Base.promote_rule(::Type{Setting{Bool}}, ::Type{Bool})  = Bool #promote_rule(Bool, Bool)

Base.string(s::Setting{String}) = string(s.value)

to_filestring(s::Setting) = string(s.code) * "=" * string(s.value)

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

    if !haskey(getfield(m, setting_field_name), s.key) ||
        isa(s.value, Function)
        getfield(m, setting_field_name)[s.key] = s
    else
        update!(getfield(m, setting_field_name)[s.key], s)
    end
end

"""
```
update!(a::Setting, b::Setting)
```

Update `a` with the fields of `b` if:

- The `key` field is updated if `a.key == b.key`
- The `print` boolean and `code` string are overwritten if `a.print` is false and
  `b.print` is true, or `a.print` is true, `b.print` is false, and
  b.code is non-empty.
- The `description` field is updated if `b.description` is nonempty
"""
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
    else
        return m.settings[s].value
    end
end
