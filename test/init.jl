include("init.jl")
m990 = init("990")

@test typeof(m990) == Model

# Test parameter promotion
alp = m990.Θ.alp
@test typeof(promote(alp,int8(3)) == (Float64, Float64)
@test typeof(promote(alp,float16(3)) == (Float64, Float64)

# Futher tests
