"""
`init_subspec!(m::RepDSGEGovDebt)`

Initializes a model subspecification by overwriting parameters from
the original model object with new parameter objects. This function is
called from within the model constructor.
"""
function init_subspec!(m::RepDSGEGovDebt, het::HetDSGEGovDebt)
    m.subspec = subspec(het)
    return
end
