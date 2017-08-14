"""
```
type ShockGroup
```

The `ShockGroup` type is used only for `plot_shock_decompositions`, in which
each group of shocks gets its own bar.

### Fields

- `name::String`
- `shocks::Vector{Symbol}`
- `color::Colorant`
"""
type ShockGroup
    name::String
    shocks::Vector{Symbol}
    color::Colorant
end

function ShockGroup(name::String, shocks::Vector{Symbol}, color_name::Symbol)
    color = parse(Colorant, color_name)
    return ShockGroup(name, shocks, color)
end