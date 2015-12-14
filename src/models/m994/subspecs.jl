#using Match

function initialize_subspec(m::Model994)

    if subspec(m) == "ss8"
        return
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

