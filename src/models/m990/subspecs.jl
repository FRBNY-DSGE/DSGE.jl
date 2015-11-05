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

    m <= parameter(:ι_p, 0.0, fixed=true,
                   description= "ι_p: The persistence of last period's inflation in
                   the equation that describes the intertemporal
                   change in prices for intermediate goods producers
                   who cannot adjust prices. The change in prices is a
                   geometric average of steady-state inflation
                   (π_star, with weight (1-ι_p)) and last period's
                   inflation (π_{t-1})).",
                   texLabel="\\iota_p")


    m <= parameter(:ι_w,   0.0, fixed=true,
                   description="ι_w: This is the something something.",
                   texLabel="\\iota_w")
    
end


