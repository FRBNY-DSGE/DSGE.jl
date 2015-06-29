# A type that bundles together the model indices functions from models/m$(spec)/modelinds.jl, as well as the number of states, etc.
type ModelInds
    endostates::Dict{ASCIIString, Int64}
    exoshocks::Dict{ASCIIString, Int64}
    expshocks::Dict{ASCIIString, Int64}
    equations::Dict{ASCIIString, Int64}
end

