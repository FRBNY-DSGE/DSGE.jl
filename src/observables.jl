"""
```
mutable struct Observable
```

### Fields

- `key::Symbol`
- `input_series::Vector{Symbol}`: vector of mnemonics, each in the form
  `:MNEMONIC__SOURCE` (e.g. `:GDP__FRED`). This vector is parsed to determine
  source (e.g. per-capita consumption gets population and consumption).
- `fwd_transform::Function`: Extracts appropriate `input_series` from a
  DataFrame of levels, and transforms data to model units (for example, computes
  per-capita growth rates from levels).
- `rev_transform::Function`: Transforms a series from model units into
  observable units. May take kwargs.
- `name::String`: e.g. \"Real GDP growth\"
- `longname::String`: e.g. \"Real GDP growth per capita\"
"""
mutable struct Observable
    key
    input_series::Vector{Symbol} # (vector of Mnemonics, e.g. GDPC1@FRED)
                                 # parse this to determine source
                                 # (eg consumption per cap gets population and consumption)

    fwd_transform::Function      # corresponds to data_transforms
                                 # field. Operates on dataframe of levels, that has every data series
                                 # requested by any observable

    rev_transform::Function      # corresponds to getyp, or getyp4q. May take kwargs (untransformed,
                                 # 1q transforms, 4q transforms, levels/growth rates, etc)

    name::String
    longname::String
end


"""
```
mutable struct PseudoObservable
```

### Fields

- `key::Symbol`
- `name::String`: e.g. \"Flexible Output Growth\"
- `longname::String`: e.g. \"Output that would prevail in a flexible-price
  economy\"
- `rev_transform::Function`: Transforms a series from model units into
  observable units. May take kwargs.
"""
mutable struct PseudoObservable
    key::Symbol
    name::String
    longname::String
    rev_transform::Function
end

function PseudoObservable(k::Symbol)
    PseudoObservable(k, "", "", identity)
end

function check_mnemonics(levels::DataFrame, mnemonics::Symbol)
    for mnemonic in mnemonics
        @assert in(mnemonic, names(levels)) "Dataframe is missing $(mnemonic)"
    end
end
