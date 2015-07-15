# A type that bundles together the model indices dictionaries from
#   models/m$(spec)/modelinds.jl
type ModelInds
    endostates::Dict{String, Int64}
    exoshocks::Dict{String, Int64}
    expshocks::Dict{String, Int64}
    eqconds::Dict{String, Int64}
    endostates_postgensys::Dict{String, Int64}
    observables::Dict{String, Int64}
end



# Given an array of names, return a dictionary mapping names to indices
# When optional field `start` provided, first index is start+1
function makedict{T<:String}(names::Array{T, 1}; start::Int = 0)
    return [names[i] => start+i for i = 1:length(names)]
end
