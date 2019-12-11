using JLD
using DifferentialEquations
using Plots
using Printf
include("../../../../src/estimate/ct_filters/ct_kalman_simple.jl")
#include("../../../../src/estimate/kalman.jl")
#include("../../../../src/estimate/ct_filters/ct_kalman_filter_mdn.jl")


EX =2
# set up transition/meas eqs matrixes
if EX == 1
    #  AR(1) model
    n_vars = 1
    n_states = 1
    n_shocks = 1
    T = Array{Float64}(undef,n_states,n_states)
    T .= -10.0
    R = Array{Float64}(undef,n_states,n_shocks)
    R = 10.0
    Z = Array{Float64}(undef,n_vars,n_states)
    Z .= 1.0
    
elseif EX == 2
    #  independent AR(1) models -- load only on the first one
    n_vars = 1
    n_states = 2
    n_shocks = 1
    T = Array{Float64}(undef,n_states,n_states)
    T[1,1] = -10.0
    T[2,2] = -10.0    
    R = Array{Float64}(undef,n_states,n_shocks)
    R .= 1.0
    @show R    
    Z = Array{Float64}(undef,n_vars,n_states)
    Z[1,1] = 1.0
    Z[1,2] = 0.0
end

# generate data (solve using Euler-maruyama)

n_quarters = 200.
global s_t = zeros(Float64, n_states) #initial values

# #use Differential Equations
# # Set up SDE
# f(u,p,t) = T*u
# g(u,p,t) = R
# dt = 1.0/90.0 # daily frequency when one period is one quarter
# tspan = (0.0,n_quarters) # 50 years of data
# W = WienerProcess(0.0, 0., 0.)
# prob = SDEProblem(f, g, s_t, tspan, W)
#      sol = solve(prob, EM(), dt = dt)
# s_1T = vcat(sol.u...)
# data_d = s_1T*Z'

#    gr()
#     default(show = true)
# p = plot(sol.t, data_d)
# sleep(2)

data_d = Array{Float64}(undef,Int(round(n_quarters/dt)),n_vars)
time_d = Array{Float64}(undef,Int(round(n_quarters/dt)))
n_step = 10.0
for i_t = 1:Int(round(n_quarters/dt))
    for i = 1:Int(n_step)
        s_t = s_t + T*s_t*dt/n_step + sqrt(dt/n_step)*R*randn(n_shocks, 1);
    end
    data_d[i_t,:] = Z*s_t
    time_d[i_t] = i_t*dt
end

gr()
default(show = true)
p = plot(time_d, data_d, color = :red)

# created 50 years worth of daily data -> transform data sampled at quarterly frequency for discrete observations
Delta_t = 1.0#1.0/90.0##1.0
del_t = dt./Delta_t



@show n_y = Int(round.(length(data_d)*dt/Delta_t;digits = 0))
data = Array{Float64}(undef,n_y,n_vars)
for i = 1:n_y
    data[i,:] .= data_d[1+(i-1)*Int(floor(Delta_t/dt)),:]
end

E = zeros(n_vars, n_vars) # meas error
Q = Matrix{Float64}(I,n_shocks,n_shocks) # variance of shocks
mean_0 = zeros(n_states,1) # intial state mean
var_0 = 0.1*Matrix{Float64}(I,n_states,n_states) # initial state variance
for sig = 0.1:0.1:2.0
    prob_out = ct_kalman_simple(T, Z, R*(sig*Q)*R', E, mean_0, var_0, data, [Delta_t])
    @show [sig prob_out]
end

# not sure the Delta_t abd dt are correct
#E = zeros(n_vars, n_vars) * Delta_t # Shock is Browian motion -> variance is Δ_t b/c this is quarterly observation data
#Q = Matrix{Float64}(I,n_shocks,n_shocks) * dt    # Shock is Brownian motion -> variance is dt b/c this is state data
#out = ct_kalman_filter(data, T, R, C, Q, Z, D, E, Delta_t)
#@show true_lik = sum(out[1])




 #   #, label = "High vol",xlims=(glimpdf[1],glimpdf[2]))


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
