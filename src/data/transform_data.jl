
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
    
    # Order the matrix columns appropriately
    transformed
end

function order_data!(m::AbstractModel, transformed::DataFrame)

    
end
