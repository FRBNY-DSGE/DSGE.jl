using DataStructures: SortedDict, insert!


T = gensym()
immutable Setting{T} 
    key::Symbol                  # name of setting
    value::T                     # whatever the setting is
    savestring::Bool             # whether or not to add this setting to the savestring
    code::AbstractString         # what gets printed to the savestring
    description::AbstractString  # description of what the setting is for
end

# for printing codes to filename string 
Base.convert{T<:Number, U<:Number}(::Type{T}, s::Setting{U}) = convert(T, s.value)
#Base.convert(::Type{Bool}, s::Setting{Bool}) = s.value
Base.convert{T<:AbstractString, U<:AbstractString}(::Type{T}, s::Setting{U}) = convert(T, s.value)
#Base.convert{T<:AbstractString}(::Type{T}, s::Setting{ASCIIString}) = s.value

Base.promote_rule{T<:Number,U<:Number}(::Type{Setting{T}}, ::Type{U}) = promote_rule(T,U)
Base.promote_rule{T<:AbstractString,U<:AbstractString}(::Type{Setting{T}}, ::Type{U}) = promote_rule(T,U)
Base.promote_rule(::Type{Setting{Bool}}, ::Type{Bool}) = promote_rule(Bool, Bool)

Base.string(s::Setting{AbstractString}) = string(s.value)

filename_string(s::Setting) = "$(s.code)=$(s.value)"




# key, value constructor
function Setting(key, value)
    Setting(key, value, false, "", "")
end

function Setting(key, value, description)
    Setting(key, value, false, "", description)
end


"""
(<=){T}(m::AbstractDSGEModel{T}, s::Setting)

Syntax for adding a setting to a model/overwriting a setting: m <= setting
"""
function (<=){T}(m::AbstractDSGEModel{T}, s::Setting)
    if s.savestring 
        # Add to a sorted dictionary of things to print
        insert!(m._filestrings, s.key, filename_string(s))
    end

    m.settings[s.key] = s
end


"""
get_setting(m::AbstractDSGEModel, setting::Symbol)

Returns the value of the setting
"""
function get_setting(m::AbstractDSGEModel, setting::Symbol)
    s_test = symbol(setting,"_test")
    
    if m.testing && in(s_test, keys(m.test_settings))
        return m.test_settings[s_test].value
    end
    
    m.settings[setting].value
end


"""
default_settings(m::AbstractDSGEModel)

Add default settings to the model's settings dictionary
"""
function default_settings(m::AbstractDSGEModel)

    # Subspec number
    m <= Setting(:subspec, "ss0", "Model sub-specification number")

    # I/O File locations
    modelpath = normpath(joinpath(dirname(@__FILE__), "..","save",spec(m)))
    datapath = normpath(joinpath(dirname(@__FILE__), "..","save",spec(m),"input_data"))

    
    m <= Setting(:modelpathroot, modelpath, "Root of data directory structure")
    m <= Setting(:datapathroot, datapath, "Input data directory path")

    # Anticipated shocks
    m <= Setting(:num_anticipated_shocks,         6, "Number of anticipated policy shocks")
    m <= Setting(:num_anticipated_shocks_padding, 20, "Padding for anticipated policy shocks")
    m <= Setting(:num_anticipated_lags,  24, "Number of periods back to incorporate zero bound expectations")

    # TODO: should be set when data are read in
    m <= Setting(:num_presample_periods, 2, "Number of periods in the presample")

    # Estimation
    m <= Setting(:reoptimize,          false, "Reoptimize the posterior mode")
    m <= Setting(:recalculate_hessian, false, "Recalculate the hessian at the mode")
    m <= Setting(:num_mh_simulations,  10000, "Number of draws per block in Metropolis-Hastings")
    m <= Setting(:num_mh_blocks,       22   , "Number of blocks for Metropolis-Hastings")
    m <= Setting(:num_mh_burn,         2    , "Number of blocks to use as burn-in in Metropolis-Hastings")
    m <= Setting(:mh_thinning_step,    5    , "How often to write draw to file in Metropolis_Hastings")

    # Data vintage
    m <= Setting(:data_vintage,        "REF", "Date of data")
    
    # Test settings
    default_test_settings(m)
    
end


"""
default_test_settings(m::AbstractDSGEModel)

Add default testing settings to the model's settings dictionary
"""
function default_test_settings(m::AbstractDSGEModel)

    test = Dict{Symbol,Setting}()

    test[:num_mh_simulations_test] = Setting(:num_mh_simulations_test, 100, false, "nsim",
                                            "Number of parameter draws per block in Metropolis-Hastings") 
    
    test[:num_mh_blocks_test]      = Setting(:num_mh_blocks_test, 1, false, "nblc",
                                             "Number of blocks to draw parameters in Metropolis-Hastings")
    
    test[:num_mh_burn_test]        = Setting(:num_mh_burn_test,   0, false, "nbrn",
                                             "Number of burn-in blocks in Metropolis-Hastings")
    
    test[:mh_thinning_step_test]   = Setting(:mh_thinning_step_test, 1, false, "thin",
                                             "Thinning step in Metropolis-Hastings")


    # I/O
    datapathroot_test = normpath(joinpath(dirname(@__FILE__), "..","test","reference"))
    modelpathroot_test = ""
    
    test[:modelpathroot_test] = Setting(:modelpathroot_test, modelpathroot_test,
                                       "Where to write files when in test mode")

    test[:datapathroot_test] = Setting(:datapathroot_test, datapathroot_test,
                                       "Location of input files when in test mode" )

    m.test_settings = test
end
    
