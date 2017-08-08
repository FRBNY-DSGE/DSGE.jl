"""
```
function mutation(c::Float64, N_MH::Int64, deterministic::Bool, system::System{Float64}, y_t::Array{Float64,1}, s_init::Array{Float64,1}, ε_init::Array{Float64,1}, nonmissing::Array{Bool,1})
```
Runs random-walk Metropolis Hastings for ONE particle. The caller should loop through all particles, calling this method on each. 
m: The model. Used to get the setting tpf_c (which is VERY important it turns out...) 
y_t: A 1 x num_measurement vector at time t for all observed y (GDP, PCE, FFR, etc.)
s_init: The starting state before mutation.
ε_init: The starting epsilon (state error) before mutation.
"""
function mutation_testing(system::System{Float64}, y_t::Array{Float64,1}, s_init::Array{Float64,1}, ε_init::Array{Float64,1}, c::Float64, N_MH::Int64, deterministic::Bool, nonmissing::Array{Bool,1}, rand_mat::Array{Float64})
    #------------------------------------------------------------------------
    # Setup
    #------------------------------------------------------------------------
    DD     = system.measurement.DD[nonmissing]
    ZZ     = system.measurement.ZZ[nonmissing,:]
    EE     = system.measurement.EE[nonmissing,nonmissing]
    RRR    = system.transition.RRR[:,nonmissing]
    TTT    = system.transition.TTT
    QQ     = system.measurement.QQ[nonmissing,nonmissing]
    sqrtS2 = RRR*Matrix(chol(nearestSPD(QQ)))'

    # Initialize ind_s and ind_ε
    ind_s = s_init
    ind_ε = ε_init

    # Store length of y_t
    n_obs = length(y_t)

    # Initialize acceptance counter to zero
    accept = 0
    accept_rate=1
    #------------------------------------------------------------------------
    # Metropolis-Hastings Steps
    #------------------------------------------------------------------------
    for i=1:N_MH
        
        if !deterministic rand_mat = randn(size(QQ,1)) end
        #else rand_mat = get_setting(m,:tpf_rand_mat) end
        
        # Generate new draw of ε from a N(ε_init, c²I) distribution, c tuning parameter, I identity
        ε_new = ε_init + c*rand_mat
        
        # Use the state equation to calculate the corresponding state from that ε 
        s_new_fore = TTT*s_init + sqrtS2*ε_new

        # Use the state equation to calculate the state corresponding to ε_init
        s_init_fore = TTT*s_init + sqrtS2*ε_init

        # Calculate difference between data and expected y from measurement equation and calculated states from above for both the new draw of ε_new and the old ε_init (we do this to calculate probabilities below. Since the error is still drawn from a Normal and every_thing is still linear, we know that y will also be normal. See equation of multivariate normal for how error_new and error_init enter into the pdf).
        error_new  = y_t - ZZ*s_new_fore - DD
        error_init = y_t - ZZ*s_init_fore - DD
        # Calculate the top and bottom probabilities for the α ratio.
        post_new = log(pdf(MvNormal(zeros(n_obs),EE),error_new)[1]*pdf(MvNormal(zeros(length(ε_new)),eye(length(ε_new),length(ε_new))),ε_new)[1])
        post_init = log(pdf(MvNormal(zeros(n_obs),EE),error_init)[1]*pdf(MvNormal(zeros(length(ε_init)),eye(length(ε_init),length(ε_init))),ε_init)[1])

        # α is the probability of accepting the new particle 
        α = exp(post_new - post_init)

        # Accept the particle with probability α
        num = (deterministic) ? 0.5 : rand()
        
        if num < α 
            # Accept
            ind_s = s_new_fore
            ind_ε = ε_new
            accept += 1
        else 
            # Reject and keep the old particle unchanged
            ind_s = s_init_fore
            ind_ε = ε_init
            
            # THIS IS DIFFERENT THAN ORIGINAL MATLAB BUT WE THINK MATLAB IS WRONG 
            #accept = 0
        end
        ε_init = ind_ε
        accept_rate = accept/N_MH
        c = c*(0.95+0.1*exp(20*(accept_rate-0.4))/(1+exp(20*(accept_rate-0.4))))
    end

    return ind_s, ind_ε, accept_rate 
end

