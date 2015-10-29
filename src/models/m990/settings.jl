

type Setting{T<:Any}
    key::Symbol                  # name of setting
    value::T                     # whatever the setting is
    savestring::Bool             # whether or not to add this setting to the savestring
    code::AbstractString         # what gets printed to the savestring
    description::AbstractString  # description of what the setting is for
end

typealias SettingDict{T} Dict{Setting{T}}

# for printing codes to filename string - maybe?
Base.convert{T<:AbstractString, U<:Integer}(::Type{T}, s::Setting{U}) = @sprintf "%s=%d" s.code s.value
Base.convert{T<:AbstractString, U::Bool}(::Type{T}, s::Setting{U}) = @sprintf "%s=%d" s.code s.value

# key, value constructor
function Setting(key::Symbol, value)
    Setting(key, value, false, "", "")
end

function Setting(key::Symbol, value, description::AbstractString)
    Setting(key, value, false, "", description)
end


"""
(<=){T}(m::AbstractDSGEModel{T}, s::Setting{T})

Syntax for adding a setting to a model/overwriting a setting: m <= setting
"""
function (<=){T}(m::AbstractDSGEModel{T}, s::Setting)
    if s.savestring 
        # Add to a sorted dictionary of things to print
    end
    setindex!(m.settings, s, s.key)
end


"""
default_settings(m::AbstractDSGEModel)

Add default settings to the model's settings dictionary
"""
function default_settings(m::AbstractDSGEModel)

    # spec and subspec numbers
    spec = split(basename(@__FILE__),'.')[1]   
    m <= Setting(:spec, spec, "Model specification number")
    m <= Setting(:subspec, 0, "Model sub-specification number")

    # I/O File locations
    savepath = normpath(joinpath(dirname(@__FILE__), *("../../../save/",spec)))
    m <= Setting(:savepath, savepath,                                 "Root of data directory structure")
    m <= Setting(:inpath,   joinpath(dirname(savepath),"input_data"), "Input data directory path")

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

    
    test[:savepath_test]           = Setting(:savepath_test, mktempdir(m.settings[:savepath])) 

    m <= Setting(:test_settings, test)
end
    
