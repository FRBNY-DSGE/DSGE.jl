"""
```
type Observable
```

### Fields
- `key::Symbol`
- `input_series::Vector{Symbol}`: (vector of Mnemonics, e.g. GDPC1@FRED). This vector is parsed to determine source (eg consumption per capita gets population and consumption).
- `fwd_transform::Function`: Extracts appropriate `input_series` from a DataFrame of levels, and transforms data to model inputs (for example, computes per capita growth rates from levels). 
- `rev_transform::Function`: Transforms a series from model units into observable units. May take kwargs. 
- `name::UTF8String`: eg. "Real GDP growth"
- `longname::UTF8String`: eg. "Real GDP growth per capita"
- `units::Symbol`: For untransformed variables, units would be "thousands" or "percent" or "levels". Transformed, could be "1q" or "4q" or something. This needs some thought.
"""
type Observable
    key
    input_series::Vector{Symbol} # (vector of Mnemonics, e.g. GDPC1@FRED)
                                 # parse this to determine source
                                 # (eg consumption per cap gets population and consumption)

    fwd_transform::Function      # corresponds to data_transforms
                                 # field. Operates on dataframe of levels, that has every data series
                                 # requested by any observable

    rev_transform::Function      # corresponds to getyp, or getyp4q. May take kwargs (untransformed,
                                 # 1q transforms, 4q transforms, levels/growth rates, etc)

    name::UTF8String
    longname::UTF8String
    # units::Symbol
end


"""
```
type PseudoObservable
```

### Fields
`key::Symbol` 
`name::UTF8String`
`longname::UTF8String`
`rev_transform::Function`
"""
type PseudoObservable
    key::Symbol 
    name::UTF8String
    longname::UTF8String
    rev_transform::Function
end

function PseudoObservable(k::Symbol)
    PseudoObservable(k,"","", identity)
end

type PseudoObservableMapping{T}
    inds::Dict{Symbol,Int}
    ZZ::Matrix{T}
    DD::Array{T}
end