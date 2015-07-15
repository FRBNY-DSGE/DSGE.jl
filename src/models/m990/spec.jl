# Model-specific specifications
function spec_vars()
    dict = Dict{String, Any}()
    
    # Number of anticipated policy shocks
    dict["nant"] = 6

    # Padding for nant
    dict["nantpad"] = 20

    # Number of periods back we should start incorporating zero bound expectations
    # ZLB expectations should begin in 2008 Q4
    dict["antlags"] = 24

    return dict
end
