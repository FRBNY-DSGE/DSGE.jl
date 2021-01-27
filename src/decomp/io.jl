function get_decomp_filename(m_new::M, m_old::M, input_type::Symbol,
                             cond_new::Symbol, cond_old::Symbol,
                             product::Symbol, class::Symbol;
                             pathfcn::Function = rawpath,
                             fileformat::Symbol = :jld2,
                             forecast_string_new = "", forecast_string_old = "") where M<:AbstractDSGEModel
    output_var = Symbol(product, class)

    fn_new = get_forecast_filename(m_new, input_type, cond_new, output_var,
                                   pathfcn = pathfcn, fileformat = fileformat, forecast_string = forecast_string_new)
    fn_old = get_forecast_filename(m_old, input_type, cond_old, output_var, forecast_string = forecast_string_old)

    dir = dirname(fn_new)
    base_new, ext = splitext(basename(fn_new))
    base_old, _   = splitext(basename(fn_old))
    base_old = replace(base_old, string(output_var) => "")
    base_old = spec(m_old) * "_" * subspec(m_old) * base_old

    return joinpath(dir, base_new * "__" * base_old * ext)
end

function get_decomp_output_files(m_new::M, m_old::M, input_type::Symbol,
                                 cond_new::Symbol, cond_old::Symbol,
                                 classes::Vector{Symbol}; forecast_string_new = "", forecast_string_old = "") where M<:AbstractDSGEModel
    output_files = Dict{Symbol, String}()
    for comp in [:data, :news, :shockdec, :dettrend, :para, :total]
        for class in classes
            product = Symbol(:decomp, comp)
            output_var = Symbol(product, class)
            output_files[output_var] =
                get_decomp_filename(m_new, m_old, input_type, cond_new, cond_old,
                                    product, class, pathfcn = rawpath, fileformat = :jld2,
                                    forecast_string_new = forecast_string_new, forecast_string_old = forecast_string_old)
        end
    end
    return output_files
end

function get_decomp_mean_file(m_new::M, m_old::M, input_type::Symbol,
                              cond_new::Symbol, cond_old::Symbol,
                              class::Symbol; forecast_string_new = "", forecast_string_old = "") where M<:AbstractDSGEModel
    get_decomp_filename(m_new, m_old, input_type, cond_new, cond_old, :decomp, class,
                        pathfcn = workpath, fileformat = :jld2,
                        forecast_string_new = forecast_string_new, forecast_string_old = forecast_string_old)
end

function write_forecast_decomposition(m_new::M, m_old::M, input_type::Symbol,
                                      classes::Vector{Symbol},
                                      decomp_output_files::Dict{Symbol, String},
                                      decomps::Dict{Symbol, Array{Float64}};
                                      block_number::Nullable{Int} = Nullable{Int}(),
                                      block_inds::AbstractRange{Int} = 1:0,
                                      verbose::Symbol = :low,
                                      forecast_string_new = "", forecast_string_old = "") where M<:AbstractDSGEModel
    for comp in [:data, :news, :shockdec, :dettrend, :para, :total]
        for class in classes
            prod = Symbol(:decomp, comp)
            var = Symbol(prod, class)
            filepath = decomp_output_files[var]

            if isnull(block_number) || get(block_number) == 1
                # Write forecast metadata to a jld2 and the raw forecast output
                # to an h5. The data in the HDF5 will be transferred to the jld2
                # and the h5 file will be deleted when combine_raw_forecast_output_and_metadata
                # is executed.
                JLD2.jldopen(filepath, true, true, true, IOStream) do file
                    # Write metadata
                    # Pass in m_old because its historical and forecast dates are used
                    write_forecast_metadata(m_old, file, var)
                end

                h5open(replace(filepath, "jld2" => "h5"), "w") do file
                    # Pre-allocate HDF5 dataset which will contain all draws
                    if !isnull(block_number) && get(block_number) == 1
                        # Determine forecast output size
                        jstep = get_jstep(m_new, n_forecast_draws(m_new, input_type))
                        dims = collect(size(decomps[var]))
                        dims[1] = convert(Int, n_forecast_draws(m_new, input_type) / jstep)

                        # Determine chunk (block) size
                        block_size = convert(Int, forecast_block_size(m_new) / jstep)
                        chunk_dims = copy(dims)
                        chunk_dims[1] = block_size

                        # Initialize dataset
                        # pfile = file.plain
                        if isdefined(HDF5, :create_dataset)
                            HDF5.create_dataset(file, "arr", datatype(Float64), dataspace(dims...); chunk = chunk_dims)
                        else
                            HDF5.d_create(file, "arr", datatype(Float64), dataspace(dims...), "chunk", chunk_dims)
                        end
                    end
                end
            end

            h5open(replace(filepath, "jld2" => "h5"), "r+") do file
                if isnull(block_number)
                    write(file, "arr", decomps[var])
                else
                    write_forecast_block(file, decomps[var], block_inds)
                end
            end

            println(verbose, :high, " * Wrote $(basename(filepath))")
        end # of loop over classes
    end # of loop over comps
end
