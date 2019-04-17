using JLD
using DifferentialEquations
using Plots
using Printf
include("../../../../src/estimate/kalman.jl")
include("../../../../src/estimate/ct_filters/ct_kalman_filter_mdn.jl")



# start with AR(1) model

# set up transition/meas eqs matrixes
n_vars = 1
n_states = 1
n_shocks = 1
T = Array{Float64}(undef,n_states,n_states)
R = Array{Float64}(undef,n_states,n_shocks)
Z = Array{Float64}(undef,n_vars,n_states)
T = -10.0
R = 1.0
Z = 1.0


# Set up SDE
f(u,p,t) = T*u
g(u,p,t) = R
dt = 1.0/90.0 # daily frequency when one period is one quarter
tspan = (0.0,200.) # 50 years of data
u0 = zeros(Float64, n_states) #initial values
W = WienerProcess(0.0, 0., 0.)
prob = SDEProblem(f, g, u0, tspan, W)

# generate data (solve using Euler-maruyama)
     sol = solve(prob, EM(), dt = dt)

   gr()
    default(show = true)
p = plot(sol)

# created 50 years worth of daily data -> transform data sampled at quarterly frequency for discrete observations
Delta_t = 1.0
del_t = dt./Delta_t
data_d = sol.u
data = Array{Float64}(undef,round.(length(data_d)*dt./Delta_t;digits = 0))
for i = 1:length(data)
    data[i] =data_d[1+ (i-1)*90]
end

E = zeros(n_vars, n_vars) * Delta_t # Shock is Browian motion -> variance is Δ_t b/c this is quarterly observation data
Q = Matrix{Float64}(I,n_shocks,n_shocks) * dt    # Shock is Brownian motion -> variance is dt b/c this is state data



prob_out = ct_kalman_simple(T, Z, R*Q*R', E, mean_0, var_0, data_y, dt)

#out = ct_kalman_filter(data, T, R, C, Q, Z, D, E, Delta_t)
#@show true_lik = sum(out[1])




 #   p = plot(ind_yy,yy[ind_yy], color = :red)#, label = "High vol",xlims=(glimpdf[1],glimpdf[2]))


#=
SeHyoun's function

% Test Kalman Filer <kalman_ct>
n_simul = 1000;
n_step = 50;
make_plots = false;
F = [-1, 0; 0, -1];
F_sim = F;

Q = 0.02*speye(2);
% Q = [0.1, 0.02; 0.02, 0.1];
R = 0.01*speye(2);
% R = [0.1, 0.02; 0.02, 0.1];
H = speye(2);

mean_init = [2; -2];
var_init = 0.01*speye(2);

mean_now = mean_init;
data_x = mean_now + sqrt(0.1)*randn(2,1);
var_now = Var_init;

dt = 0.25;

F_val = 0.1:0.1:1.9;
data_y_iter = zeros(2, n_simul);
for iter_time = 1:n_simul
  for i = 1:n_step
    data_x = data_x + F*data_x*dt/n_step + sqrt(dt/n_step)*sqrtm(Q)*randn(2, 1);
  end
  data_y = H*data_x + sqrtm(R)*randn(2, 1);
  data_y_iter(:, iter_time) = data_y;
end
figure(101);
% plot(0:dt:dt*(n_simul-1), data_y_iter);
plot(data_y_iter(1,:), data_y_iter(2, :), '.-');

time_step = dt*ones(n_simul, 1);
for iter_outer = 1:length(F_val)

  F_sim = F_val(iter_outer)*F;

  [mean_now, var_now, log_like] = kalman_ct(F_sim, H, Q, R, mean_init, var_init, data_y_iter, time_step);

  log_like_iter(iter_outer) = log_like;

  figure(4);
  clf;
  plot(F_val(1:iter_outer), log_like_iter(1:iter_outer));
  xlim([0.1, 1.9]);
  drawnow
end

=#
