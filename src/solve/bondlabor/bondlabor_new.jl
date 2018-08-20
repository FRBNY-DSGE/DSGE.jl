using DSGE
using BenchmarkTools
using Distributions
using JLD
using Roots

include("klein_solve.jl")

test_output = true

m = BondLabor()

# Steady-state computation
steadystate!(m)
@btime steadystate!(m)

test_output && include("test/steady_state.jl")

# Jacobian computation
JJ = DSGE.jacobian(m)
@btime JJ = DSGE.jacobian(m)

test_output && include("test/jacobian.jl")

throw("stop!")

# Create PPP matrix
P1 = kron(eye(ns),ones(nx,1))
Ptemp = eye(nx)
Ptemp = Ptemp[:,2:end]
P2 = kron(eye(ns),Ptemp)
P = [P1 P2]

(Q,R)=qr(P)

S         = Q[:,ns+1:end]'; #
Qleft     = cat([1 2],S,[1],eye(nx*ns),[1])
Qx        = cat([1 2],S,[1])
Qy        = cat([1 2],eye(nx*ns),[1])

# Solve
tic()
gx2, hx2, gx, hx =  klein_solve(JJ, Qleft, Qx, Qy)
toc()

test_output && include("test/solve.jl")
