"""
```
mutation{S<:AbstractFloat}(system::System{S}, y_t::Array{S,1}, s_init::Array{S,1}, 
        ε_init::Array{S,1}, c::S, N_MH::Int, nonmissing::Array{Bool,1})
```
Runs random-walk Metropolis Hastings for single particle. The caller should loop through 
all particles, calling this method on each. 

### Inputs

- `system`: state-space system matrices
- `y_t`: vector of observables at time t
- `s_init`: vector of starting state before mutation (ŝ in paper) 
- `ε_init`: vector of starting state error before mutation
- `c`: scaling factor used to achieve a desired acceptance rate, adjusted via:
    
    cₙ = cₙ₋₁f(1-R̂ₙ₋₁(cₙ₋₁))

    Where c₁ = c_star and R̂ₙ₋₁(cₙ₋₁) is the emprical rejection rate based on mutation
    phase in iteration n-1. Average is computed across previous N_MH RWMH steps.

- `N_MH`: number of Metropolis Hastings steps
- `nonmissing`: vector of booleans used to remove NaN values from matrices in system object

### Outputs

- `s_out`: mutated state vector
- `ε_out`: output ε shock corresponding to state vector 
- `accept_rate`: acceptance rate across N_MH steps

"""
function mutation{S<:AbstractFloat}(system::System{S}, y_t::Array{S,1}, s_init::Array{S,1}, 
    ε_init::Array{S,1}, c::S, N_MH::Int, nonmissing::Array{Bool,1})
    
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

    # Initialize s_out and ε_out
    s_out = s_init
    ε_out = ε_init

    # Store length of y_t, ε
    n_obs    = length(y_t)
    n_states = length(ε_init)

    # Initialize acceptance counter to zero
    accept = 0

    #------------------------------------------------------------------------
    # Metropolis-Hastings Steps
    #------------------------------------------------------------------------
    for i=1:N_MH
        
        # Generate new draw of ε from a N(ε_init, c²I) distribution, c tuning parameter, I identity
        ε_new = ε_init + c*randn(size(QQ, 1))
        
        # Use the state equation to calculate the corresponding state from that ε 
        s_new_e = TTT*s_init + sqrtS2*ε_new

        # Use the state equation to calculate the state corresponding to ε_init
        s_init_e = TTT*s_init + sqrtS2*ε_init

        # Calculate difference between data and expected y from measurement equation
        error_new  = y_t - ZZ*s_new_e - DD
        error_init = y_t - ZZ*s_init_e - DD

        # Calculate posteriors
        post_new = log(pdf(MvNormal(zeros(n_obs), EE), error_new)[1] * 
                       pdf(MvNormal(zeros(n_states), eye(n_states, n_states)), ε_new)[1])
        post_init = log(pdf(MvNormal(zeros(n_obs), EE), error_init)[1] * 
                        pdf(MvNormal(zeros(n_states), eye(n_states, n_states)), ε_init)[1])

        # Calculate α, probability of accepting the new particle 
        α = exp(post_new - post_init)

        # Accept the particle with probability α
        if rand() < α 
            # Accept and update particle
            s_out = s_new_e
            ε_out = ε_new
            accept += 1
        else 
            # Reject and keep particle unchanged
            s_out = s_init_e
            ε_out = ε_init
        end
        ε_init = ε_out
    end

    # Calculate acceptance rate
    accept_rate = accept/N_MH

    return s_out, ε_out, accept_rate 
end