include("init.jl")
m990 = init("990")

@test typeof(m990) == Model

# Test parameter promotion
alp = m990.Θ.alp
@test promote_type(Param, Float16) = Float64
@test promote_type(Param, Int8) = Float64

# Futher tests
