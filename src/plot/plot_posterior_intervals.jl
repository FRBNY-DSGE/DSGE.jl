# block is for if there are multiple blocks of parameters being produce in succession
# then it is meant to indicate the ordering of blocks for visual comparison.
# I.e. if param_range 1:10 corresponds to the first block, and 11:20 corresponds to the
# second block then the files will be aptly named ...block=1.pdf, ...block=2.pdf
function plot_posterior_intervals(m::AbstractModel; plotroot::String = figurespath(m, "estimate"),
                                  verbose::Symbol = :low, label::String = "", title::String
                                  = isempty(label) ? "" : "$label Posterior Intervals",
                                  param_range::UnitRange = UnitRange(1, n_parameters_free(m)),
                                  block::Int64 = 0,
                                  excl_list::Vector{Symbol} = Vector{Symbol}(0))

    df = load_posterior_moments(m; load_bands = true, include_fixed = false, excl_list =
                                excl_list)

    # Indexing out a block of the parameters
    df = df[param_range, :]
    df[:param_inds] = param_range
    parameter_labels = df[:param]

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
    end
end

# block is for if there are multiple blocks of parameters being produce in succession
# then it is meant to indicate the ordering of blocks for visual comparison.
# I.e. if param_range 1:35 corresponds to the first block, and 36:68 corresponds to the
# second block then the files will be aptly named ...b1.pdf, ...b2.pdf
function plot_posterior_interval_comparison(m_baseline::AbstractModel,
                                            m_comparison::AbstractModel;
                                            plotroot::String = figurespath(m, "estimate"),
                                            verbose::Symbol = :low,
                                            in_deviations::Bool = false,
                                            scale_by_std::Bool = false,
                                            base_label::String = "",
                                            comp_label::String = "",
                                            title::String = in_deviations ? "$comp_label deviations from $base_label Interval Comparisons" : "$base_label and $comp_label Interval Comparisons",
                                            param_range::UnitRange =
                                            UnitRange(1, n_parameters_free(m_baseline)),
                                            block::Int64 = 0,
                                            excl_list::Vector{Symbol} = Vector{Symbol}(0))

    df_baseline = load_posterior_moments(m_baseline; load_bands = true,
                                         include_fixed = false, excl_list = excl_list)
    df_comparison = load_posterior_moments(m_comparison; load_bands = true,
                                           include_fixed = false, excl_list = excl_list)

    # Indexing out a block of the parameters
    df_baseline = df_baseline[param_range, :]
    df_comparison = df_comparison[param_range, :]
    df_baseline[:param_inds] = param_range
    df_comparison[:param_inds] = param_range

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
                                :hline, markerstrokecolor = :blue, title = title,
                                series_annotations = parameter_labels)
        @df df_deviations plot!(p, :param_inds, :comp_lb, label = "", linealpha = 0, marker =
                                :hline, markerstrokecolor = :blue)
        @df df_deviations plot!(p, :param_inds, :base_mean, label = "", marker = :hline,
                                markerstrokecolor = :black, linecolor = :black)
        @df df_deviations plot!(p, :param_inds, :comp_mean, label = "", linealpha = 0, marker =
                                :hline, markerstrokecolor = :red)

        rescale_annotations_fontsize!(p, 6)
        relocate_annotations!(p, :right, :top)
    else
        p = @df df_baseline plot(:param_inds, :post_ub, label = base_label, linealpha = 0, marker =
                              :hline, markerstrokecolor = :black, title = title, series_annotations =
                              parameter_labels)
        @df df_baseline plot!(p, :param_inds, :post_lb, label = "", linealpha = 0, marker =
                              :hline, markerstrokecolor = :black)
        @df df_comparison plot!(p, :param_inds, :post_ub, label = comp_label, linealpha = 0, marker =
                                :hline, markerstrokecolor = :red)
        @df df_comparison plot!(p, :param_inds, :post_lb, label = "", linealpha = 0, marker =
                                :hline, markerstrokecolor = :red)
        @df df_baseline plot(:param_inds, :post_mean, label = "", linealpha = 0, marker = :hline,
                                 markerstrokecolor = :black)
        @df df_comparison plot!(p, :param_inds, :post_mean, label = "", linealpha = 0, marker =
                                :hline, markerstrokecolor = :red)

        rescale_annotations_fontsize!(p, 6)
        relocate_annotations!(p, :right, :top)
    end

    if !isempty(plotroot)
        vint = data_vintage(m_baseline)
        output_file = joinpath(plotroot, "posterior_comparison_intervals_" *
                               "deviat=$(in_deviations)_" * "vint=$vint")
        if !iszero(block)
            output_file = output_file * "_block=$block"
        end
        output_file = output_file * ".pdf"
        save_plot(p, output_file, verbose = verbose)
    end
end
