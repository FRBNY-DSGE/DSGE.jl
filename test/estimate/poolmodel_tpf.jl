using JLD2, StateSpaceRoutines, BenchmarkTools

# Verifies that tpf works properly with a PoolModel object
tpf_main_input = load(joinpath(dirname(@__FILE__), "../reference/tpf_poolmodel.jld2"))
tuning = tpf_main_input["tuning"]
lpd_805 = tpf_main_input["lpd_805"]
lpd_904 = tpf_main_input["lpd_904"]
data = exp.([reshape(lpd_904, 1, length(lpd_904)); reshape(lpd_805, 1, length(lpd_805))])
m = DSGE.PoolModel("ss0")

# Define state space
Φpm, Ψpm, F_ϵpm, F_upm, ~ = compute_system(m)

# Load in test inputs and outputs. The saved TPF outputs are RNG-dependent, so they live in a
# version-tagged reference file (Julia 1.7+ switched the default RNG to Xoshiro256++):
#   "111" → <1.5,  "150" → 1.5–1.6,  "1126" → ≥1.7.  Deterministic inputs stay unversioned.
ver = VERSION < v"1.5" ? "111" : (VERSION < v"1.7" ? "150" : "1126")
test_file_inputs = load(joinpath(dirname(@__FILE__), "../reference/tpf_aux_inputs.jld2"))
outputs_path = joinpath(dirname(@__FILE__),
                        "../reference/tpf_aux_outputs_poolmodel_version=$(ver).jld2")

# Flip to true to regenerate the saved reference outputs (RNG-dependent), then flip back.
writing_output = false

# Guard the load: on a first regen run the version-tagged file may not exist yet (the
# writing_output blocks below create it). When writing, comparisons are self-consistent.
test_file_outputs_pm = isfile(outputs_path) ? load(outputs_path) : Dict{String, Any}()

φ_old = test_file_inputs["phi_old"]
norm_weights = test_file_inputs["norm_weights"]
coeff_terms = test_file_inputs["coeff_terms"]
log_e_1_terms = test_file_inputs["log_e_1_terms"]
log_e_2_terms = test_file_inputs["log_e_2_terms"]
inc_weights = test_file_inputs["inc_weights"]
HH = ones(1,1)

s_t_nontemp = [.5; .5] * ones(1,1000) # this might be desirable since all particles same then
Random.seed!(1793)
Ψ47_pm(x) = Ψpm(x,data[:,47])
weight_kernel!(coeff_terms, log_e_1_terms, log_e_2_terms, φ_old, Ψ47_pm, data[:, 47],
               s_t_nontemp, det(HH), inv(HH);
               initialize = false,
               poolmodel = true)
φ_new = next_φ(φ_old, coeff_terms, log_e_1_terms, log_e_2_terms, length(data[:,47]), tuning[:r_star], 2)
correction!(inc_weights, norm_weights, φ_new, coeff_terms, log_e_1_terms, log_e_2_terms, length(data[:,47]))
true_inc_wt = φ_new^(1/2) * coeff_terms[1] * exp(log_e_1_terms[1]) * exp(φ_new * log_e_2_terms[1])

@testset "Corection and Auxiliary Tests" begin
    @test coeff_terms[1] ≈ (φ_old)^(-1)
    @test log_e_1_terms[1] ≈ log(Ψpm(s_t_nontemp[:,1], data[:,47])) * -φ_old
    @test log_e_2_terms[1] ≈ log(Ψpm(s_t_nontemp[:,1], data[:,47]))
    @test φ_new ≈ 1.0
    @test inc_weights[1] ≈ true_inc_wt
end

φ_old = test_file_inputs["phi_old"]
norm_weights = test_file_inputs["norm_weights"]
coeff_terms = test_file_inputs["coeff_terms"]
log_e_1_terms = test_file_inputs["log_e_1_terms"]
log_e_2_terms = test_file_inputs["log_e_2_terms"]
inc_weights = test_file_inputs["inc_weights"]
s_t_nontemp0 = [.5; .5] * ones(1, 900)
s_t_nontemp = hcat(s_t_nontemp0, repeat([.49; .51], 1, 100))
weight_kernel!(coeff_terms, log_e_1_terms, log_e_2_terms, φ_old, Ψ47_pm, data[:, 47],
               s_t_nontemp, det(HH), inv(HH);
               initialize = false,
               poolmodel = true)
φ_new = next_φ(φ_old, coeff_terms, log_e_1_terms, log_e_2_terms, length(data[:,47]), tuning[:r_star], 2)
correction!(inc_weights, norm_weights, φ_new, coeff_terms, log_e_1_terms, log_e_2_terms, length(data[:,47]))
true_inc_wt = φ_new^(1/2) * (coeff_terms[1] * exp(log_e_1_terms[1])
                             * exp(φ_new * log_e_2_terms[1]) * 900/1000
                             + coeff_terms[end] * exp(log_e_1_terms[end])
                             * exp(φ_new * log_e_2_terms[end]) * 100/1000)

@testset "Ensure different states lead to different measurement errors" begin
    @test norm_weights[1] != norm_weights[end]
    @test log_e_1_terms[1] != log_e_1_terms[end]
    @test log_e_2_terms[1] != log_e_2_terms[end]
    @test inc_weights[1] != inc_weights[end]
    @test true_inc_wt ≈ mean(inc_weights)  # analytic mixture weight == empirical mean (by construction)
end

## Mutation Tests
s_t1_temp = [.49; .51] * ones(1,1000)
ϵ_t = reshape(test_file_inputs["eps_t"][1,:], 1, 1000)
QQ = F_ϵpm.σ * ones(1,1)
accept_rate = test_file_inputs["accept_rate"]
c = test_file_inputs["c"]

c = update_c(c, accept_rate, tuning[:target_accept_rate])
Random.seed!(47)
StateSpaceRoutines.mutation!(Φpm, Ψ47_pm, QQ, det(HH), inv(HH), φ_new, data[:,47],
                             s_t_nontemp, s_t1_temp, ϵ_t, c, tuning[:n_mh_steps];
                             poolmodel = true)
if writing_output
    test_file_outputs_pm["s_t_nontemp"] = copy(s_t_nontemp)
    test_file_outputs_pm["eps_t"]       = copy(ϵ_t)
end
@testset "Mutation Tests" begin
    @test s_t_nontemp[1] ≈ test_file_outputs_pm["s_t_nontemp"][1]
    @test ϵ_t[1] ≈ test_file_outputs_pm["eps_t"][1]
end

## Whole TPF Tests
Random.seed!(47)
s_init = reshape(rand(Uniform(.4, .6), 1000), 1, 1000)
s_init = [s_init; 1 .- s_init]
out_no_parallel = tempered_particle_filter(data, Φpm, Ψpm, F_ϵpm, F_upm, s_init;
                                           tuning..., verbose = :none, fixed_sched = [1.],
                                           parallel = false, poolmodel = true)
Random.seed!(47)
out_parallel_one_worker = tempered_particle_filter(data, Φpm, Ψpm, F_ϵpm, F_upm, s_init;
                                                   tuning..., verbose = :none, fixed_sched = [1.],
                                                   parallel = true, poolmodel = true)
if writing_output
    test_file_outputs_pm["out_no_parallel"]         = out_no_parallel
    test_file_outputs_pm["out_parallel_one_worker"] = out_parallel_one_worker
end
@testset "TPF tests" begin
    @test out_no_parallel[1] ≈ test_file_outputs_pm["out_no_parallel"][1]
    # Parallel path (DistributedArrays); RNG-stream-sensitive, so compare within tolerance to
    # the version-tagged saved value (see the `ver` tier above) rather than a hand-set literal.
    @test abs(out_parallel_one_worker[1] - test_file_outputs_pm["out_parallel_one_worker"][1]) < 0.05
end

if writing_output
    JLD2.jldopen(outputs_path, "w") do file
        for (k, v) in test_file_outputs_pm
            file[k] = v
        end
    end
end

################
# Benchmarking #
################
run_benchmarks = false
if run_benchmarks
    results = Tuple{String, Any}[]

    # correction! / weight_kernel! pipeline (rebuild inputs each sample; both mutate in place)
    b_corr = @benchmark begin
        correction!(iw, nw, $φ_new, ct, le1, le2, length($(data[:,47])))
    end setup = (iw  = copy($inc_weights); nw  = copy($norm_weights);
                 ct  = copy($coeff_terms);  le1 = copy($log_e_1_terms); le2 = copy($log_e_2_terms))
    push!(results, ("correction!", b_corr))

    # whole tempered particle filter (serial)
    b_tpf = @benchmark tempered_particle_filter($data, $Φpm, $Ψpm, $F_ϵpm, $F_upm, $s_init;
                                                $tuning..., verbose = :none, fixed_sched = [1.],
                                                parallel = false, poolmodel = true)
    push!(results, ("tempered_particle_filter (serial)", b_tpf))

    println("\n===== estimate/poolmodel_tpf benchmark results =====")
    for (name, b) in results
        println(rpad(name, 34), " time: ", BenchmarkTools.prettytime(median(b).time),
                "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end

nothing
