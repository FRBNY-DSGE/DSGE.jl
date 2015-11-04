#using Match

function initialize_subspec(m::Model990)

    if subspec(m) == "ss0"
        return
    elseif subspec(m) == "ss5"
        return ss5(model)
    else
        error("This subspec is not defined.")
    end
    
    ## This is neat pattern matching we can do, but has deprication warnings...
    ## @match subspec(m) begin
    ##     "ss0"          => return
    ##     "ss5"          => ss5(m)
    ##     _              => error("This subspec is not defined.")
    ## end
end

function ss5(m::Model990)
    m <= parameter(m[:ι_p], 0.000, valuebounds=(0.0,0.0), prior=PointMass(0.0), fixed=true)
    m <= parameter(m[:ι_w], 0.000, valuebounds=(0.0,0.0), prior=PointMass(0.0), fixed=true)
end


