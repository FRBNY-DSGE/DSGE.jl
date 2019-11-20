"""
`init_subspec!(m::HetDSGELag)`

Initializes a model subspecification by overwriting parameters from
the original model object with new parameter objects. This function is
called from within the model constructor.
"""
function init_subspec!(m::HetDSGELag)
    if subspec(m) == "ss0"
        return
    elseif subspec(m) == "ss1"
        return ss1!(m)
    else
        error("This subspec should be a 0")
    end
end

"""
```
ss1!(m::HetDSGELag)
```

Initializes subspec 1 of `HetDSGELag`.
This shuts down all shocks except the technology process (z shock).
"""
function ss1!(m::HetDSGELag)

end
