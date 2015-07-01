# A type that bundles together the model indices dictionaries from
#   models/m$(spec)/modelinds.jl
type ModelInds
    endostates::Dict{String, Int64}
    exoshocks::Dict{String, Int64}
    expshocks::Dict{String, Int64}
    equations::Dict{String, Int64}
end



# Given an array of names, return a dictionary mapping names to indices
function makedict{T<:String}(names::Array{T, 1})
    return [names[i] => i for i = 1:length(names)]
end
