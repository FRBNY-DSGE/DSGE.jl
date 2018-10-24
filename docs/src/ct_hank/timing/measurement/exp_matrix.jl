using EHANK

### Scalar error
# Given dt, what size growth rate Q is needed such that
# a kth-order approximation of power series is no longer okay
dt = 1. / 3 # monthly for quarterly
tol = 1e-5
Q_low = 0 # run bisection search
Q_high = 1
kvec = [1; 2; 3]
Q_final = zeros(length(kvec))

for ind = 1:length(kvec)
    Q = 0.
    for i = 1:100000
        Q += .01
        truth = exp(Q * dt)/Q - 1 / Q
        power = kvec[ind]
        approx = dt
        for j = 1:power
            approx += Q^j * dt^(j+1) / factorial(j+1)
        end
        if abs(truth - approx) >= tol
            Q_final[ind] = Q
            break
        end
    end
end

### Matrix error w/KrusellSmith
# Given T, for what size dt is
# a kth-order approximation of power series no longer okay
tol = 1e-8
dt_vec = collect(linspace(1/90, .5, 10))
power_vec = zeros(length(dt_vec))

# Solve model
m = KrusellSmith()
T, ~, ~, ~, ~ = solve(m)

# Bisection search
for ind = 1:length(dt_vec)
    dt = dt_vec[ind]
    power = 0
    while power <= 10
        power += 1
        approx = dt
        for j = 1:power
            approx += T^j * dt^(j+1) / factorial(j+1)
        end
        if maximum(abs.(T^power * dt^(power + 1) / factorial(power + 1))) <= tol
            power_vec[ind] = power-1
            break
        end
        if power == 10
            power_vec[ind] = NaN
        end
    end
end
