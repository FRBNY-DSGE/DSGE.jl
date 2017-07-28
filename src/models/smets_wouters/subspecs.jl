"""
`init_subspec!(m::SmetsWouters)`

Initializes a model subspecification by overwriting parameters from
the original model object with new parameter objects. This function is
called from within the model constructor.
"""
function init_subspec!(m::SmetsWouters)

    if subspec(m) == "ss0"
        return
    elseif subspec(m) == "ss1"
        return ss1!(m,5)
    elseif subspec(m) == "ss2"
        return ss1!(m,10)
    elseif subspec(m) == "ss3"
        return ss1!(m,2)
    else
        error("This subspec is not defined.")
    end
end

function ss1!(m::SmetsWouters,divisor::Int)
    scale = (1/divisor)
    #The following values come from 0.1*stdev(_obs)
    m <= parameter(:e_y, scale*0.868241996, fixed=true, 
                   description = "e_y: Measurement error on GDP", tex_label="e_y")
    m <= parameter(:e_L, scale*0.283703727, fixed=true, 
                   description = "e_L: Measurement error on hours worked", tex_label="e_L")
    m <= parameter(:e_w, scale*0.600015046, fixed=true, 
                   description = "e_w: Measurement error on wages", tex_label = "e_w")
    m <= parameter(:e_π, scale*0.600386299, fixed=true, 
                   description = "e_π: Measurement error on GDP deflator", tex_label= "e_π")
    m <= parameter(:e_R, scale*0.846722459, fixed=true,
                   description = "e_R: Measurement error on nonominal rate of interest", tex_label="e_R")
    m <= parameter(:e_c, scale*0.710297826, fixed=true, 
                   description = "e_c: Measurement error on consumption", tex_label="e_c")
    m <= parameter(:e_i, scale*2.405946712, fixed=true, 
                   description = "e_i: Measurement error on investment", tex_label="e_i")

end
