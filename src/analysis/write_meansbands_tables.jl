"""
```
write_meansbands_table(filename::AbstractString, mb::MeansBands, var:Symbol)
```

```
write_meansbands_table(filename::AbstractString, mb::MeansBands, [vars:Vector{Symbol} = Vector{Symbol}()])
```

Write means and bands (ordered) to a csv file in `tablespath(m, "forecast")`

### Inputs
- `m`: Model object
- `mb`: a meansbands object

## Optional inputs
- `var`: which variable in `mb` you want to print means and bands
  for. If left out or empty, all variables will be printed.
"""
function write_meansbands_table(m::AbstractModel, mb::MeansBands, var::Symbol)

    cond  = mb.metadata[:cond_type]
    para  = mb.metadata[:para]

    filename = get_forecast_filename(m, para, cond, var,
                                     pathfcn = tablespath, fileformat = :csv)

    means = mb.means[[:date, var]]
    bands = mb.bands[var][[:date; map(symbol, get_density_bands(mb))]]

    df = join(bands, means, on = :date)
    rename!(df, var, symbol("mean$var"))

    writetable(filename, df)
    println(" * Wrote means and bands for $var to $filename")

end

function write_meansbands_table(m::AbstractModel, mb::MeansBands, vars::Vector{Symbol} = Vector{Symbol}())
    if isempty(vars)
        vars = setdiff(names(mb.means), [:date])
    end
    for var in vars
        write_meansbands_table(m, mb, var)
    end
end
