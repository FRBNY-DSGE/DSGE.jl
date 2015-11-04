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

# subspecs
function sstest(m::Model990)
    # Change all the fields of an unfixed parameter
    m <= parameter(m[:ι_w], 0.000, valuebounds=(0.0, .9999), transform_parameterization=(0.0,0.9999), transform=Untransformed(), prior=Normal(0.0,1.0))

    # Change the value to something outside the bounds
    @test_throws ParamBoundsError m <= parameter(m[:ζ_p], 0.0)

    # Change an unfixed parameter to be fixed
    m <= parameter(m[:ι_p], 0.000, valuebounds=(0.0,0.0), prior=PointMass(0.0), fixed=true)

    # Change a fixed parameter
    m <= parameter(m[:δ],   0.02,  valuebounds=(0.02,0.02), prior=PointMass(0.02), fixed=true)

    # incomplete change for fixed parameter - shouldnt go through
    m <= parameter(m[:ϵ_p], 11.0, valuebounds=(11.0,11.0), fixed=true)

    # incomplete change for unfixed parameter to fixed - shouldnt go through
    m <= parameter(m[:ψ1], 0.0, fixed=true)
    
    steadystate!(m)
end

m = Model990()
sstest(m)

@test m[:ι_w].value == 0.0
@test m[:ι_w].valuebounds == (0.0, .9999)
@test m[:ι_w].transform == Untransformed()
@test m[:ι_w].transform_parameterization == (0.0,0.9999)
@test isa(m[:ι_w].prior.value, Normal)

@test m[:ι_p].value == 0.0
@test m[:ι_p].valuebounds == (0.0, 0.0)
@test isa(m[:ι_p].prior.value, PointMass)
@test m[:ι_p].fixed == true

@test m[:δ].value == 0.02
@test m[:δ].valuebounds == (0.02, 0.02)
@test isa(m[:δ].prior.value, PointMass)
@test m[:δ].prior.value.μ == 0.02
@test m[:δ].fixed == true

@test m[:ϵ_p].value == 10.0
@test m[:ϵ_p].valuebounds == (10.0,10.0)
@test m[:ϵ_p].fixed == true

@test m[:ψ1].value == 1.3679
@test m[:ψ1].transform == DSGE.Exponential()
@test isa(m[:ψ1].prior.value, Normal)
@test m[:ψ1].fixed==false
