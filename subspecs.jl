"""
`init_subspec!(m::Model1010)`

Initializes a model subspecification by overwriting parameters from
the original model object with new parameter objects. This function is
called from within the model constructor.
"""
function init_subspec!(m::Model1010)
    if subspec(m) == "ss1"
        return
    else
        error("This subspec is not defined.")
    end
end

