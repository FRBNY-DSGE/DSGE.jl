"""
`transform_data(m::AbstractModel, levels::DataFrame)`

Transform data loaded in levels and order columns appropriately for the DSGE model.
"""
function transform_data(m::AbstractModel, levels::DataFrame)

    transformed = DataFrame()
    transformed[:date] = levels[:date]
    
    for (i, series) in enumerate(keys(m.data_transforms))
        transformed[series] = call(m.data_transforms[series], levels)
    end
    
    sort!(transformed, cols = :date)
end

