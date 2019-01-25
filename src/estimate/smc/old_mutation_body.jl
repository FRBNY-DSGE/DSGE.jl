            like_new = -Inf
            post_new = -Inf
            like_old_data = -Inf

            try
                update!(m, para_new)
                like_new = likelihood(m, data; sampler = true,
                                      use_chand_recursion = use_chand_recursion,
                                      verbose = verbose)
                if like_new == -Inf
                    post_new = like_old_data = -Inf
                end
                post_new = ϕ_n*like_new + prior(m) - para_new_density
                like_old_data = isempty(old_data) ? 0. : likelihood(m, old_data; sampler = true,
                                                         use_chand_recursion = use_chand_recursion,
                                                         verbose = verbose)
            catch err
                if isa(err, ParamBoundsError)
                    println("ParamBoundsError!")
                    post_new = like_new = like_old_data = -Inf
                else
                    throw(err)
                end
            end

            # Accept Reject
            post_old = post - para_old_density
            η        = exp(post_new - post_old)

            if η > 1.
                η = 1.0
            end
            #@show  η, alp[1]
            #@show step_pr, step_lik
            #@show post_new, post_old, post, para_old_density

            if step_prob < η # accept
                para      = para_new
                like      = like_new
                post      = post_new + para_new_density # Have to add it back so as to not accumulate
                like_prev = like_old_data               # para_new_density throughout the iterations
                accept   += length(block_a) #true
            end

            @show η, alp[1], my_alp
            @assert abs(η - alp[1]) < 1e-4

            @show like_new, lik1
            @show post_new, para_new_density, pr1

            if !(like_new == -Inf && lik1 == -Inf) @assert abs(like_new - lik1) < 1e-4 end
            if !(post_new == -Inf && pr1 == -Inf) @assert abs(post_new + para_new_density - pr1) < 10 end
