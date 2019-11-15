# block is for if there are multiple blocks of parameters being produce in succession
# then it is meant to indicate the ordering of blocks for visual comparison.
# I.e. if param_range 1:10 corresponds to the first block, and 11:20 corresponds to the
# second block then the files will be aptly named ...block=1.pdf, ...block=2.pdf
function plot_posterior_intervals(m::AbstractDSGEModel; cloud::ParticleCloud = ParticleCloud(m, 0),
                                  df::DataFrame = DataFrame(),
                                  plotroot::String = figurespath(m, "estimate"),
                                  verbose::Symbol = :low, label::String = "", title::String
                                  = isempty(label) ? "" : "$label Posterior Intervals",
                                  block::Int64 = 0, include_fixed::Bool = false,
                                  param_range::UnitRange = include_fixed ? UnitRange(1, n_parameters(m)) : UnitRange(1, n_parameters_free(m)),
                                  excl_list::Vector{Symbol} = Vector{Symbol}(undef, 0))

    # If the moments to be plotted are not pre-provided
    if isempty(df)
        df = load_posterior_moments(m; cloud = cloud, load_bands = true, include_fixed =
                                    include_fixed, excl_list = excl_list)
    end

    # Indexing out a block of the parameters
    df = df[param_range, :]
    df[!,:param_inds] = param_range
    parameter_labels = df[!,:param]

    # Plot non-outliers
    p = @df df plot(:param_inds, :post_ub, label = "", linealpha = 0, marker =
                    :hline, markerstrokecolor = :black, title = title, series_annotations =
                    parameter_labels)
    @df df plot!(p, :param_inds, :post_lb, label = "", linealpha = 0, marker =
                 :hline, markerstrokecolor = :black)
    @df df plot!(:param_inds, :post_mean, label = "", linealpha = 0, marker =
                 :hline, markerstrokecolor = :black)

    rescale_annotations_fontsize!(p, 6)
    relocate_annotations!(p, :right, :top)

    if !isempty(plotroot)
        vint = data_vintage(m)
        output_file = joinpath(plotroot, "posterior_intervals_" * "vint=$vint")
        if !iszero(block)
            output_file = output_file * "_block=$block"
        end
        output_file = output_file * ".pdf"
        save_plot(p, output_file, verbose = verbose)
    else
        return p
    end
end

# Calculate the posterior intervals averaged over the estimations in the
# clouds vector
function plot_posterior_intervals(m::AbstractDSGEModel, clouds::Vector{ParticleCloud};
                                  plotroot::String = figurespath(m, "estimate"),
                                  verbose::Symbol = :low, label::String = "",
                                  title::String = isempty(label) ? "" : "$label Posterior Intervals",
                                  block::Int64 = 0, include_fixed::Bool = false,
                                  param_range::UnitRange = include_fixed ? UnitRange(1, n_parameters(m)) : UnitRange(1, n_parameters_free(m)),
                                  excl_list::Vector{Symbol} = Vector{Symbol}(undef, 0))

    df = load_posterior_moments(m, clouds, load_bands = true,
                                include_fixed = include_fixed, excl_list = excl_list)

    plot_posterior_intervals(m; df = df, plotroot = plotroot, verbose = verbose,
                             label = label, title = title, param_range = param_range,
                             block = block, include_fixed = include_fixed, excl_list = excl_list)
end

# block is for if there are multiple blocks of parameters being produce in succession
# then it is meant to indicate the ordering of blocks for visual comparison.
# I.e. if param_range 1:35 corresponds to the first block, and 36:68 corresponds to the
# second block then the files will be aptly named ...b1.pdf, ...b2.pdf
function plot_posterior_interval_comparison(m_baseline::AbstractDSGEModel,
                                            m_comparison::AbstractDSGEModel;
                                            cloud_baseline::ParticleCloud = ParticleCloud(m_baseline, 0),
                                            cloud_comparisons::Vector{ParticleCloud} = [ParticleCloud(m_comparison, 0)],
                                            df_baseline::DataFrame = DataFrame(),
                                            df_comparisons::Vector{DataFrame} = [DataFrame()],
                                            plotroot::String = figurespath(m_baseline, "estimate"),
                                            verbose::Symbol = :low,
                                            in_deviations::Bool = false,
                                            scale_by_std::Bool = false,
                                            base_label::String = "",
                                            comp_labels::Vector{String} = [""],
                                            title::String = in_deviations ? "$comp_label deviations from $base_label Interval Comparisons" : "$base_label and $comp_label Interval Comparisons",
                                            block::Int64 = 0, include_fixed::Bool = false,
                                            param_range::UnitRange = include_fixed ? UnitRange(1, n_parameters(m_baseline)) : UnitRange(1, n_parameters_free(m_baseline)),
                                            excl_list::Vector{Symbol} = Vector{Symbol}(undef, 0),
                                            filename_tag::String = "")

    if isempty(df_baseline)
        df_baseline = load_posterior_moments(m_baseline; cloud = cloud_baseline, load_bands = true,
                                             include_fixed = include_fixed, excl_list = excl_list)
    end

    if isempty(df_comparisons[1])
        for i=1:length(cloud_comparisons)
            df_comparisons[i] = load_posterior_moments(m_comparison; cloud = cloud_comparisons[i],
                                                      load_bands = true, include_fixed = include_fixed,
                                                      excl_list = excl_list)
        end
    end
    for i=1:length(df_comparisons)
        df_comparisons[i] = df_comparisons[i][param_range, :]
        df_comparisons[i][:param_inds] = param_range
    end

    # Indexing out a block of the parameters
    df_baseline = df_baseline[param_range, :]
  #  df_comparison = df_comparison[param_range, :]
    df_baseline[:param_inds] = param_range
  #  df_comparison[:param_inds] = param_range

    parameter_labels = df_baseline[:param]

    if in_deviations
        df_deviations = DataFrame()
        df_deviations[:param] = df_baseline[:param]
        df_deviations[:param_inds] = df_baseline[:param_inds]

        # Baseline model mean normalized to 0
        df_deviations[:base_mean] = zeros(size(df_baseline, 1))

        # Comparison model mean and bands normalized to deviations from baseline mean and
        # bands (hence the mean deviation may be higher than either of the bands' deviations
        # in the resulting plot)
        df_deviations[:comp_mean] = df_comparison[:post_mean] - df_baseline[:post_mean]
        df_deviations[:comp_ub] = df_comparison[:post_ub] - df_baseline[:post_ub]
        df_deviations[:comp_lb] = df_comparison[:post_lb] - df_baseline[:post_lb]

        # Scale the comparison mean and lower and upper bounds deviations by the std. of the
        # parameter under the baseline
        if scale_by_std
            df_deviations[:comp_mean] = df_deviations[:comp_mean]./df_baseline[:post_std]
            df_deviations[:comp_ub] = df_deviations[:comp_ub]./df_baseline[:post_std]
            df_deviations[:comp_lb] = df_deviations[:comp_lb]./df_baseline[:post_std]
        end

        p = @df df_deviations plot(:param_inds, :comp_ub, label = "", linealpha = 0, marker =
                                :hline, markerstrokecolor = :blue, title = "",
                                series_annotations = parameter_labels, legend=:bottomright)
        @df df_deviations plot!(p, :param_inds, :comp_lb, label = "", linealpha = 0, marker =
                                :hline, markerstrokecolor = :blue, legend=:bottomright)
        @df df_deviations plot!(p, :param_inds, :base_mean, label = "", marker = :hline,
                                markerstrokecolor = :black, linecolor = :black, legend=:bottomright)
        @df df_deviations plot!(p, :param_inds, :comp_mean, label = "", linealpha = 0, marker =
                                :hline, markerstrokecolor = :red, legend=:bottomright)

        rescale_annotations_fontsize!(p, 6)
        relocate_annotations!(p, :right, :top)
    else
        p = @df df_baseline plot(:param_inds, :post_ub, label = base_label, linealpha = 0, marker =
                              :hline, markerstrokecolor = :black, title = "", series_annotations =
                              parameter_labels, legend=:bottomright)
        @df df_baseline plot!(p, :param_inds, :post_lb, label = "", linealpha = 0, marker =
                              :hline, markerstrokecolor = :black, legend=:bottomright)
        @df df_baseline plot!(p, :param_inds, :post_mean, label = "", linealpha = 0, marker =
                              :hline, markerstrokecolor = :black, legend=:bottomright)
        colors = [:red, :blue, :green, :orange, :yellow]
        for i = 1:length(df_comparisons)
            @df (df_comparisons[i]) plot!(p, :param_inds, :post_ub, label = comp_labels[i], linealpha = 0, marker =
                                    :hline, markerstrokecolor = colors[i], legend=:bottomright)
            @df (df_comparisons[i]) plot!(p, :param_inds, :post_lb, label = "", linealpha = 0, marker =
                                    :hline, markerstrokecolor = colors[i], legend=:bottomright)
            @df (df_comparisons[i]) plot!(p, :param_inds, :post_mean, label = "", linealpha = 0, marker =
                                    :hline, markerstrokecolor = colors[i], legend=:bottomright)
        end
        rescale_annotations_fontsize!(p, 6)
        relocate_annotations!(p, :right, :top)
    end

    if !isempty(plotroot)
        vint = data_vintage(m_baseline)
        output_file = joinpath(plotroot, "posterior_comparison_intervals_" *
                               "deviat=$(in_deviations)_" * "vint=$vint" * "_" *filename_tag)
        if !iszero(block)
            output_file = output_file * "_block=$block"
        end
        output_file = output_file * ".pdf"
        save_plot(p, output_file, verbose = verbose)
    else
        return p
    end
end

# Calculate the posterior intervals averaged over the estimations in the
# clouds vectors
function plot_posterior_interval_comparison(m_baseline::AbstractDSGEModel,
                                            m_comparison::AbstractDSGEModel,
                                            clouds_baseline::Vector{ParticleCloud},
                                            clouds_comparisons::Vector{Vector{ParticleCloud}};
                                            plotroot::String = figurespath(m_baseline, "estimate"),
                                            verbose::Symbol = :low,
                                            in_deviations::Bool = false,
                                            scale_by_std::Bool = false,
                                            base_label::String = "",
                                            comp_labels::Vector{String} = [""],
                                            title::String = in_deviations ? "$(comp_labels[1]) deviations from $base_label Interval Comparisons" : "$base_label and $(comp_labels[1]) Interval Comparisons",
                                            block::Int64 = 0, include_fixed::Bool = false,
                                            param_range::UnitRange = include_fixed ? UnitRange(1, n_parameters(m_baseline)) : UnitRange(1, n_parameters_free(m_baseline)),
                                            excl_list::Vector{Symbol} = Vector{Symbol}(undef, 0),
                                            filename_tag::String = "")

    df_baseline = load_posterior_moments(m_baseline, clouds_baseline, load_bands = true,
                                         include_fixed = include_fixed, excl_list = excl_list)
    df_comparisons = Vector{DataFrame}(undef, length(clouds_comparisons))
    for i=1:length(clouds_comparisons)
        df_comparisons[i] = load_posterior_moments(m_comparison, clouds_comparisons[i], load_bands = true,
                                         include_fixed = include_fixed, excl_list = excl_list)
    end

    plot_posterior_interval_comparison(m_baseline, m_comparison, df_baseline = df_baseline, df_comparisons = df_comparisons,
                                       plotroot = plotroot, verbose = verbose, in_deviations = in_deviations, scale_by_std = scale_by_std,
                                       base_label = base_label, comp_labels = comp_labels, title = "", param_range = param_range, block =
                                       block, excl_list = excl_list, filename_tag = filename_tag)
end
