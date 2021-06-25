"""
`init_subspec!(m::AnSchorfheide)`

Initializes a model subspecification by overwriting parameters from
the original model object with new parameter objects. This function is
called from within the model constructor.
"""
function init_subspec!(m::AnSchorfheide) # unless you want to add subspecifications, do not edit this function definition.
    if subspec(m) == "ss0"               # For examples of subspecs, see `smets_wouters` or `smets_wouters_orig`
        return
    else
        error("This subspec should be a 0")
    end
end
