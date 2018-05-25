function get_decomp_filename(m_new::AbstractModel, m_old::AbstractModel,
                             input_type::Symbol, cond_new::Symbol, cond_old::Symbol,
                             product::Symbol, class::Symbol;
                             pathfcn::Function = rawpath,
                             fileformat::Symbol = :jld)
    output_var = Symbol(product, class)

    fn_new = get_forecast_filename(m_new, input_type, cond_new, output_var,
                                   pathfcn = pathfcn, fileformat = fileformat)
    fn_old = get_forecast_filename(m_old, input_type, cond_old, output_var)

    dir = dirname(fn_new)
    base_new, ext = splitext(basename(fn_new))
    base_old, _   = splitext(basename(fn_old))
    base_old = replace(base_old, string(output_var), "")
    base_old = spec(m_old) * "_" * subspec(m_old) * base_old

    return joinpath(dir, base_new * "__" * base_old * ext)
end

function get_decomp_output_files(m_new::AbstractModel, m_old::AbstractModel,
                                 input_type::Symbol, cond_new::Symbol, cond_old::Symbol,
                                 classes::Vector{Symbol})
    output_files = Dict{Symbol, String}()
    for comp in [:state, :shock, :data, :param, :total]
        for class in classes
            product = Symbol(:decomp, comp)
            output_var = Symbol(product, class)
            output_files[output_var] =
                get_decomp_filename(m_new, m_old, input_type, cond_new, cond_old, product, class,
                                    pathfcn = rawpath, fileformat = :h5)
        end
    end
    return output_files
end

function get_decomp_mean_file(m_new::AbstractModel, m_old::AbstractModel,
                              input_type::Symbol, cond_new::Symbol, cond_old::Symbol,
                              class::Symbol)
    get_decomp_filename(m_new, m_old, input_type, cond_new, cond_old, :decomp, class,
                        pathfcn = workpath, fileformat = :jld)
end

function write_forecast_decomposition(m_new::AbstractModel, m_old::AbstractModel,
                                      input_type::Symbol, classes::Vector{Symbol},
                                      hs::Union{Int, UnitRange{Int}},
                                      decomp_output_files::Dict{Symbol, String},
                                      decomps::Dict{Symbol, Matrix{Float64}};
                                      block_number::Nullable{Int} = Nullable{Int}(),
                                      block_inds::Range{Int} = 1:0,
                                      verbose::Symbol = :low)
    for comp in [:state, :shock, :data, :param, :total]
        for class in classes
            prod = Symbol(:decomp, comp)
            var = Symbol(prod, class)
            filepath = decomp_output_files[var]

            if isnull(block_number) || get(block_number) == 1
                jldopen(filepath, "w") do file
                    # Write metadata
                    # Pass in m_old because its date_mainsample_end is used to calculate dates
                    DSGE.write_forecast_metadata(m_old, file, prod, class, hs = hs)

                    # Pre-allocate HDF5 dataset which will contain all draws
                    if !isnull(block_number) && get(block_number) == 1
                        # Determine forecast output size
                        dims = DSGE.get_forecast_output_dims(m_new, input_type, Symbol(:forecast, class))
                        dims = (dims[1], dims[2], length(hs)) # dims is ndraws x nvars x nperiods
                        block_size = forecast_block_size(m_new)
                        chunk_dims = collect(dims)
                        chunk_dims[1] = block_size

                        # Initialize dataset
                        pfile = file.plain
                        HDF5.d_create(pfile, "arr", datatype(Float64), dataspace(dims...), "chunk", chunk_dims)
                    end
                end
            end

            jldopen(filepath, "r+") do file
                if isnull(block_number)
                    write(file, "arr", decomps[var])
                else
                    DSGE.write_forecast_block(file, decomps[var], block_inds)
                end
            end

            if DSGE.VERBOSITY[verbose] >= DSGE.VERBOSITY[:high]
                println(" * Wrote $(basename(filepath))")
            end
        end # of loop over classes
    end # of loop over comps
end