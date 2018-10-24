using StateSpaceRoutines
using StateSpaceRoutines.kalman_filter
using Distributions
using DataFrames
using Plots

#=

# Written by William Chen, April 2018
# Implements a Lucas economy driven by
# productivity shocks in continuous time
# and tests the ability of the temperered
# particle filter to recover parameters
# from simulated observations of asset prices

function generate_ar1_cont_time_data(year_ct,sigma,gA,rho,sd_measure_error,obs_freq)

# Written by William Chen, April 2018
# Generates data from a Lucas economy driven by
# productivity shocks in continuous time
# and tests the ability of the temperered
# particle filter to recover parameters
# from simulated observations of asset prices

# parameters
#sigma = 0.03 # volatility
#gA = 0.4 # rate of mean reversion
#rho = 0.05 # discount factor
years = year_ct # number of years
#sd_measure_error = 0.001 # standard deviation of measurement error
#obs_freq = "quarter" # frequency of observable data fed to particle filter

# dictionary mapping observation frequency to number of elapsed hours
observation_frequency = Dict("day" => 24, "week" => 8*24, "everyotherweek" => 16*24, "month" => 30*24, "quarter" => 90*24)

# generate data
#T_data = years * 360 * 24 # hours in each year
#dt = 1/360/24 # size of "dt" increment
T_data = years*360
dt = 1/360
Brownian = Distributions.Normal(0,sqrt(dt)) # Brownian motion driving productivity
shocks = rand(Brownian, T_data-1)
logA_data = zeros(T_data) # zero vector to be filled with logA observations
logA_data[1] = 0
for t = 2:T_data
    logA_data[t] = logA_data[t-1] - gA*logA_data[t-1]*dt+ sigma*shocks[t-1]
end
logq_data = logA_data - log(rho)
F_obs = Distributions.Normal(0,sd_measure_error)
noisy_logq_data = logq_data + rand(F_obs,length(logq_data)) # scramble with some noise
obs_freq_hrs = observation_frequency[obs_freq] # frequency of observable data in hours
obs_freq_days = convert(Int,observation_frequency[obs_freq]/24) # frequency of observable data in days
obs_data = logq_data[collect(1:obs_freq_days:T_data)] # actually seen observables given user chosen frequency
#noisy_obs_data = noisy_logq_data[collect(1:obs_freq_hrs:T_data)]
noisy_obs_data = obs_data + rand(F_obs,length(obs_data))
writedlm("ar1_ct_data.txt", noisy_obs_data)
end

=#

###

function get_oneassethank_statespace()

NTreps = 50
num_years = 40
max_sigma_lik = zeros(NTreps)
max_gA_lik = zeros(NTreps)

#for NT = 1:NTreps
# define true parameters and test vectors
n_reps = 15 # how many reps for each parameter
sigma_true = 1 # volatility 0.02,.1
sigma_vec = linspace(0.8,1.2,n_reps)
gA_true = 1-.8^4 # rate of mean reversion 0.3,0.5
gA_vec = linspace(1-.95^4,1-.4^4,n_reps)
rho_true = 0.05 # discount factor 0.03, 0.1
rho_vec = linspace(0.045,0.055,10)
sd_measure_error = sqrt(sigma_true^2/10) # standard deviation of measurement error
sim_freq = "day" # frequency of simulated data b/n observable data points

# dictionary mapping simulation frequency to number of simulated states
# Note that this simulated dictionary assumes quarterly data. Need to modify accordingly if data is not quarterly.
simulated_frequencies = Dict("hours"=> 90*24,"quarter_day" => 90*6,"day" => 90, "week" => 12, "everyotherweek" => 6, "month" => 3)
n_sim_states = simulated_frequencies[sim_freq]

# dictionary mapping simulation frequency to number of sub-divisions per period (assumed to be a year)
observation_frequency = Dict("hours" => 360*24,"quarter_day"=> 360*6, "day" => 360, "week" => 48, "everyotherweek" => 24, "month" => 12)
sim_dt = 1/observation_frequency[sim_freq] # dt of simulated states

# get data; only need to comment out if you need to create the data
#generate_ar1_cont_time_data(num_years,sigma_true,gA_true,rho_true,sd_measure_error,"quarter")

df = readtable("ar1_ct_data.txt", header = false, separator = ' ')
data = convert(Matrix{Float64}, df)'

# likelihood vectors
lik_sigma = zeros(length(sigma_vec),1)
lik_gA = zeros(length(gA_vec),1)
lik_rho = zeros(length(rho_vec),1)

# tuning of particle filter
tuning = Dict(:verbose => :none, :r_star => 2., :c => 0.3, :accept_rate => 0.4, :target => 0.4,
              :xtol => 0., :resampling_method => :systematic, :N_MH => 1,
              :n_particles => 1000, :n_presample_periods => 0,
              :adaptive => true, :allout => false, :parallel => false)


sigma = sigma_true
sigmaA = sigma
rho = rho_true
gA = gA_true

# Define transition and measurement functions
F_u = Distributions.MvNormal(sd_measure_error^2*eye(1,1))
F_ϵ = Distributions.MvNormal(sigma^2*sim_dt*eye(n_sim_states,n_sim_states))

# create state & measurement matrices
TTT = zeros(n_sim_states,n_sim_states) # state transition matrix
RRR = zeros(n_sim_states,n_sim_states) # matrix on state transition errors
CCC = zeros(n_sim_states) # constant term in state transiiton
DD = zeros(1) # constant term in measurement transition
ZZ = zeros(1,n_sim_states) # observable given states
QQ = zeros(n_sim_states,n_sim_states) # Covariance matrix for state transition
EE = zeros(1,1) # observable covariance matrix

for i = 1:n_sim_states
    if i == 1
        TTT[i,n_sim_states] = (1-gA*sim_dt)
        RRR[i,i] = 1
    else
        TTT[i,n_sim_states] = (1-gA*sim_dt)^i
        RRR[i,i] = 1
        for j = 1:(i-1)
            RRR[i,i-j] = (1-gA*sim_dt)^j
        end
    end
    QQ[i,i] = sigmaA^2*sim_dt
end
DD[1] = -log(rho)
ZZ[1,n_sim_states] = 1
EE[1,1] = sd_measure_error^2
# Generation of the initial state draws
#lik_true = kalman_filter(data, TTT, RRR, CCC, QQ, ZZ, DD, EE)
#lik_true = sum(lik_true[1])



return  data, TTT, RRR, CCC, QQ, ZZ, DD, EE, n_sim_states
end
#=
for j = 1:length(lik_sigma)
    sigma = sigma_vec[j]
    sigmaA = sigma
    rho = rho_true
    gA = gA_true

    # Define transition and measurement functions
    F_u = Distributions.MvNormal(sd_measure_error^2*eye(1,1))
    F_ϵ = Distributions.MvNormal(sigma^2*sim_dt*eye(n_sim_states,n_sim_states))

    # create state & measurement matrices
    TTT = zeros(n_sim_states,n_sim_states) # state transition matrix
    RRR = zeros(n_sim_states,n_sim_states) # matrix on state transition errors
    CCC = zeros(n_sim_states) # constant term in state transiiton
    DD = zeros(1) # constant term in measurement transition
    ZZ = zeros(1,n_sim_states) # observable given states
    QQ = zeros(n_sim_states,n_sim_states) # Covariance matrix for state transition
    EE = zeros(1,1) # observable covariance matrix

    for i = 1:n_sim_states
        if i == 1
            TTT[i,n_sim_states] = (1-gA*sim_dt)
            RRR[i,i] = 1
        else
            TTT[i,n_sim_states] = (1-gA*sim_dt)^i
            RRR[i,i] = 1
            for k = 1:(i-1)
                RRR[i,i-k] = (1-gA*sim_dt)^k
            end
        end
        QQ[i,i] = sigma^2*sim_dt
    end
    DD[1] = -log(rho)
    ZZ[1,n_sim_states] = 1
    EE[1,1] = sd_measure_error^2
    # Generation of the initial state draws
    out = kalman_filter(data, TTT, RRR, CCC, QQ, ZZ, DD, EE)
    lik_sigma[j] = sum(out[1])
end
for j = 1:length(lik_gA)
    sigma = sigma_true
    sigmaA = sigma
    rho = rho_true
    gA = gA_vec[j]

    # Define transition and measurement functions
    F_u = Distributions.MvNormal(sd_measure_error^2*eye(1,1))
    F_ϵ = Distributions.MvNormal(sigma^2*sim_dt*eye(n_sim_states,n_sim_states))

    # create state & measurement matrices
    TTT = zeros(n_sim_states,n_sim_states) # state transition matrix
    RRR = zeros(n_sim_states,n_sim_states) # matrix on state transition errors
    CCC = zeros(n_sim_states) # constant term in state transiiton
    DD = zeros(1) # constant term in measurement transition
    ZZ = zeros(1,n_sim_states) # observable given states
    QQ = zeros(n_sim_states,n_sim_states) # Covariance matrix for state transition
    EE = zeros(1,1) # observable covariance matrix

    for i = 1:n_sim_states
        if i == 1
            TTT[i,n_sim_states] = (1-gA*sim_dt)
            RRR[i,i] = 1
        else
            TTT[i,n_sim_states] = (1-gA*sim_dt)^i
            RRR[i,i] = 1
            for k = 1:(i-1)
                RRR[i,i-k] = (1-gA*sim_dt)^k
            end
        end
        QQ[i,i] = sigma^2*sim_dt
    end
    DD[1] = -log(rho)
    ZZ[1,n_sim_states] = 1
    EE[1,1] = sd_measure_error^2
    # Generation of the initial state draws
    out = kalman_filter(data, TTT, RRR, CCC, QQ, ZZ, DD, EE)
    lik_gA[j] = sum(out[1])
end
#max_gA_lik[NT] = gA_vec[indmax(lik_gA)]
#max_sigma_lik[NT] = sigma_vec[indmax(lik_sigma)]
# for j = 1:length(lik_rho)
#     println(j)
#     sigma = sigma_true
#     rho = rho_vec[j]
#     gA = gA_true
#
#     # Define transition and measurement functions
#     F_obs = Distributions.MvNormal(sd_measure_error^2*sim_dt*eye(1,1))
#     F_state = Distributions.MvNormal(sigma^2*sim_dt*eye(n_sim_states,n_sim_states))
#
#     function state_fnct(s_t::Vector{Float64}, eps_t::Vector{Float64})
#         s_tNew= zeros(length(s_t)) # out vector
#         s_tNew[1] = -gA*s_t[length(s_t)]*sim_dt + eps_t[1]
#         for i = 2:length(s_t)
#             s_tNew[i] = s_tNew[i-1] - gA*s_tNew[i-1]*sim_dt + eps_t[i]
#         end
#         return s_tNew
#     end
#
#     obs_fnct(s_t::Vector{Float64},u_t::Vector{Float64}) = s_t[1] - log(rho) + u_t # gives log q given current period state
#
#     # Generation of the initial state draws
#     s0 = rand(F_state,1)[:,1]
#     s_init = initialize_state_draws(s0, F_state,state_fnct,tuning[:n_particles])
#     lik_rho[j,1] = tempered_particle_filter(data, state_fnct, obs_fnct, F_state, F_obs, s_init; tuning...)
# end
#
# # make Plots
plot(sigma_vec, lik_sigma, xlab = "sigma", ylab = "log lik", title = "log lik plot varying sigma")
plot(gA_vec, lik_gA, xlab = "gA", ylab = "log lik", title = "log lik plot varying gA")
#plot(rho_vec, lik_rho, xlab = "rho", ylab = "log lik", title = "log lik plot varying rho")
#end
=#