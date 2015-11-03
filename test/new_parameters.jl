import DSGE: Exponential

#TODO: add tests and specify behavior for fixed params
# transformations
for T in subtypes(Transform)
    u = parameter(:σ_pist, 2.5230, (1e-8, 5.), (1e-8, 5.), T(), fixed=false)
    @test ( toreal(u) |> x -> tomodel(u,x) ) == u.value
    
    if !isa(T,Type{Untransformed})
        # check toreal and tomodel to different things if T is not Untransformed
        @test toreal(u,u.value) != tomodel(u,u.value)
    end
end

# probability
N = 10^2
u = parameter(:bloop, 2.5230, (1e-8, 5.), (1e-8, 5.), SquareRoot(); fixed = true)
v = parameter(:cat, 2.5230, (1e-8, 5.), (1e-8, 5.), Exponential(), Gamma(2.00, 0.1))

pvec = @compat ParameterVector{Float64}(N)
for i in 1:length(pvec)
	pvec[i] = (i%2 == 0) ? u : v
end
@test_approx_eq logpdf(pvec) 50*logpdf(v)
@test_approx_eq pdf(pvec) exp(50*logpdf(v))

updated = update(pvec, ones(length(pvec)))
update!(pvec, ones(length(pvec)))

@test all(updated .== pvec)
@test logpdf(pvec) == logpdf(updated)

# test we only update unfixed parameters
for p in pvec
	if p.fixed
		@test p.value == 2.5230
	elseif isa(p, Parameter)
		@test p.value == one(Float64)
	end
end

# vector of new values must be the same length
@test_throws AssertionError update!(pvec, ones(length(pvec)-1))

for w in [parameter(:moop, 3.0, fixed=false), parameter(:moop, 3.0; scaling = log, fixed=false)]
	# new values must be of the same type
	@test_throws MethodError parameter(w, one(Int))

	# new value is out of bounds
	@test_throws ParamBoundsError parameter(w, -1.)
end

