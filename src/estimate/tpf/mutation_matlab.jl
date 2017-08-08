"""
```
function mutation{S<:AbstractFloat}(system::System{S}, y_t::Array{S,1}, s_init::Array{S,1}, 
        ε_init::Array{S,1}, c::S, N_MH::Int, deterministic::Bool, nonmissing::Array{Bool,1})
```
Runs random-walk Metropolis Hastings for single particle. The caller should loop through 
all particles, calling this method on each. 

### Inputs

- `system::System{S}`: state-space system matrices
- `y_t`: 1 x num_measurement vector at time t for all observed y (GDP, PCE, FFR, etc.)
- `s_init`: starting state before mutation
- `ε_init`: starting state error before mutation

### Outputs

- `s_hat`:
- `ε_hat`: 
- `accept_rate`:  

"""
function mutation_matlab{S<:AbstractFloat}(system::System{S}, y_t::Array{S,1}, s_init::Array{S,1}, 
    ε_init::Array{S,1}, c::S, N_MH::Int, deterministic::Bool, nonmissing::Array{Bool,1}, 
    rand_mat::Array{S}, mutation_rand_mat1, mutation_rand_mat2)
    #------------------------------------------------------------------------
    # Setup
    #------------------------------------------------------------------------
    DD     = system[:DD][nonmissing]
    ZZ     = system[:ZZ][nonmissing,:]
    EE     = system[:EE][nonmissing,nonmissing]
    RRR    = system[:RRR][:,nonmissing]
    TTT    = system[:TTT]
    QQ     = system[:QQ][nonmissing,nonmissing]
    sqrtS2 = RRR*Matrix(chol(nearestSPD(QQ)))'

    # Initialize s_hat and ε_hat
    s_hat = s_init
    ε_hat = ε_init

    # Store length of y_t
    n_obs = length(y_t)

    # Initialize acceptance counter to zero
    accept = 0

    #------------------------------------------------------------------------
    # Metropolis-Hastings Steps
    #------------------------------------------------------------------------
    for i=1:N_MH
        
        if !deterministic rand_mat = randn(size(QQ,1)) end
        #else rand_mat = get_setting(m,:tpf_rand_mat) end
        
        rand_mat = (i == 1) ? mutation_rand_mat1 : mutation_rand_mat2

        # Generate new draw of ε from a N(ε_init, c²I) distribution, c tuning parameter, I identity
        ε_new = ε_init + c*rand_mat
        
        # Use the state equation to calculate the corresponding state from that ε 
        s_new_e = TTT*s_init + sqrtS2*ε_new

        # Use the state equation to calculate the state corresponding to ε_init
        s_init_e = TTT*s_init + sqrtS2*ε_init

        # Calculate difference between data and expected y from measurement equation and calculated states from above for both the new draw of ε_new and the old ε_init (we do this to calculate probabilities below. Since the error is still drawn from a Normal and every_thing is still linear, we know that y will also be normal. See equation of multivariate normal for how error_new and error_init enter into the pdf).
        error_new  = y_t - ZZ*s_new_e - DD
        error_init = y_t - ZZ*s_init_e - DD

        # Calculate the top and bottom probabilities for the α ratio.
        post_new = log(pdf(MvNormal(zeros(n_obs),EE),error_new)[1]*pdf(MvNormal(zeros(length(ε_new)),eye(length(ε_new),length(ε_new))),ε_new)[1])
        post_init = log(pdf(MvNormal(zeros(n_obs),EE),error_init)[1]*pdf(MvNormal(zeros(length(ε_init)),eye(length(ε_init),length(ε_init))),ε_init)[1])

        # α is the probability of accepting the new particle 
        α = exp(post_new - post_init)

        # Accept the particle with probability α
        num = (deterministic) ? 0.5 : rand()
        
        if num < α 
            # Accept
            s_hat = s_new_e
            ε_hat = ε_new
            accept += 1
        else 
            # Reject and keep particle unchanged
            s_hat = s_init_e
            ε_hat = ε_init
            
            # THIS IS DIFFERENT THAN ORIGINAL MATLAB BUT WE THINK MATLAB IS WRONG 
            #accept = 0
        end
        ε_init = ε_hat
    end
    accept_rate = accept/N_MH
    return s_hat, ε_hat, accept_rate 
end

