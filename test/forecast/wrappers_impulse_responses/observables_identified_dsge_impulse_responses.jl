path = dirname(@__FILE__)

# Set up arguments
global m = AnSchorfheide(testing = true)
m <= Setting(:impulse_response_horizons, 10)
m[:σ_z].value = sqrt(.0201332) # need to change parameters or you get a PositiveDefiniteException
m[:σ_g].value = sqrt(1.09725)  # in the long run identification
m[:σ_R].value = sqrt(.0590387)
m[:ρ_z].value = 0.99
m[:ρ_g].value = 0.99
m[:ρ_R].value = 0.99
θ = map(x -> x.value, m.parameters)
DSGE.update!(m, θ)
system = compute_system(m)
θmat = convert(Matrix{Float64}, hcat(θ, θ)')
obs_shock = zeros(n_observables(m))
obs_shock[1] = 1.

@testset "Wrapper for Cholesky-identified DSGE impulse responses" begin
    for do_parallel in [false, true]
        for do_flip in [false, true]
            for method in [:cholesky, :choleskyLR, :cholesky_long_run]
                states_chol, obs_chol, pseudo_chol =
                    impulse_responses(system, impulse_response_horizons(m),
                                      Matrix{Float64}(I, n_observables(m), n_observables(m)),
                                      obs_shock, flip_shocks = do_flip, restriction = method)

                states1, obs1, pseudo1 = impulse_responses(m, θ, :mode, method, 1;
                                                           flip_shocks = do_flip,
                                                           parallel = do_parallel)
                states2, obs2, pseudo2 = impulse_responses(m, θmat, :full, method, 1,
                                                           flip_shocks = do_flip,
                                                           parallel = do_parallel)
                @test @test_matrix_approx_eq states_chol states1[:, :, 1]
                @test @test_matrix_approx_eq obs_chol obs1[:, :, 1]
                @test @test_matrix_approx_eq pseudo_chol pseudo1[:, :, 1]
                @test @test_matrix_approx_eq states_chol states2[:, :, 1]
                @test @test_matrix_approx_eq obs_chol obs2[:, :, 1]
                @test @test_matrix_approx_eq pseudo_chol pseudo2[:, :, 1]
                @test @test_matrix_approx_eq states_chol states2[:, :, 2]
                @test @test_matrix_approx_eq obs_chol obs2[:, :, 2]
                @test @test_matrix_approx_eq pseudo_chol pseudo2[:, :, 2]
            end
        end
    end

    # No tests but check that there is no error if we create a MeansBands
    mb = impulse_responses(m, θmat, :full, :cholesky, 1,
                           create_meansbands = true,
                           test_meansbands = true)
    mb = impulse_responses(m, θmat, :full, :choleskyLR, 1,
                           create_meansbands = true,
                           test_meansbands = true)
    mb = impulse_responses(m, θmat, :full, :maxBC, 1,
                           create_meansbands = true,
                           test_meansbands = true)
end

@testset "Wrapper for `maxBC`-identified DSGE impulse responses" begin
    for do_parallel in [false, true]
        for do_flip in [false, true]
            for method in [:maxBC, :maximum_business_cycle_variance]
                states_maxBC, obs_maxBC, pseudo_maxBC =
                    impulse_responses(system, impulse_response_horizons(m),
                                      (2 * π / 32, 2 * π / 6), 1;
                                      flip_shocks = do_flip)

                states1, obs1, pseudo1 = impulse_responses(m, θ, :mode, method, 1,
                                                           flip_shocks = do_flip,
                                                           parallel = do_parallel)
                states2, obs2, pseudo2 = impulse_responses(m, θmat, :full, method, 1,
                                                           flip_shocks = do_flip,
                                                           parallel = do_parallel)
                @test @test_matrix_approx_eq states_maxBC states1[:, :, 1]
                @test @test_matrix_approx_eq obs_maxBC obs1[:, :, 1]
                @test @test_matrix_approx_eq pseudo_maxBC pseudo1[:, :, 1]
                @test @test_matrix_approx_eq states_maxBC states2[:, :, 1]
                @test @test_matrix_approx_eq obs_maxBC obs2[:, :, 1]
                @test @test_matrix_approx_eq pseudo_maxBC pseudo2[:, :, 1]
                @test @test_matrix_approx_eq states_maxBC states2[:, :, 2]
                @test @test_matrix_approx_eq obs_maxBC obs2[:, :, 2]
                @test @test_matrix_approx_eq pseudo_maxBC pseudo2[:, :, 2]
            end
        end
    end

    # No tests but check that there is no error if we create a MeansBands
    mb = impulse_responses(m, θmat, :full, :cholesky, 1,
                           create_meansbands = true,
                           test_meansbands = true)
    mb = impulse_responses(m, θmat, :full, :choleskyLR, 1,
                           create_meansbands = true,
                           test_meansbands = true)
    mb = impulse_responses(m, θmat, :full, :maxBC, 1,
                           create_meansbands = true,
                           test_meansbands = true)
end

nothing
