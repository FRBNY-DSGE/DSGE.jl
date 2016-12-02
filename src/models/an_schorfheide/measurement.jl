"""
```
measurement{T<:AbstractFloat}(m::AnSchorfheide{T}, TTT::Matrix{T}, RRR::Matrix{T},
                              CCC::Matrix{T}; shocks::Bool = true)
```

Assign measurement equation
```
X_t = ZZ*S_t + DD + u_t
```
where
```
u_t = eta_t + MM*eps_t
var(eta_t) = EE
var(u_t) = HH = EE + MM*QQ*MM'
cov(eps_t,u_t) = VV = QQ*MM'
```
"""
function measurement{T<:AbstractFloat}(m::AnSchorfheide{T},
                                       TTT::Matrix{T},
                                       RRR::Matrix{T},
                                       CCC::Matrix{T};
                                       shocks::Bool = true)
endo = m.endogenous_states
exo  = m.exogenous_shocks
obs  = m.observables

# If shocks = true, then return measurement equation matrices with rows and columns for
# anticipated policy shocks
if shocks
    _n_observables = n_observables(m)
    _n_states = n_states_augmented(m)
    _n_shocks_exogenous = n_shocks_exogenous(m)
    endo_new = m.endogenous_states_augmented
else
    _n_observables = n_observables(m) - n_anticipated_shocks(m)
    _n_states = n_states_augmented(m) - n_anticipated_shocks(m)
    _n_shocks_exogenous = n_shocks_exogenous(m) - n_anticipated_shocks(m)
    endo_new = Dict(
        [(key,m.endogenous_states_augmented[key] - n_anticipated_shocks(m)) for key in keys(m.endogenous_states_augmented)])
end

    ZZ = zeros(_n_observables, _n_states)
    DD = zeros(_n_observables, 1)
    MM = zeros(_n_observables, _n_shocks_exogenous)
    EE = zeros(_n_observables, _n_observables)
    QQ = zeros(_n_shocks_exogenous, _n_shocks_exogenous)

    ## Output growth - Quarterly!

    tau = m.parameters[1] 
    kappa = m.parameters[2] 
    psi1 = m.parameters[3] 
    psi2 = m.parameters[4] 
    rA = m.parameters[5] 
    piA = m.parameters[6] 
    gammaQ = m.parameters[7] 
    rho_R = m.parameters[8] 
    rho_g = m.parameters[9] 
    rho_z = m.parameters[10] 
    sigma_R = m.parameters[11] 
    sigma_g = m.parameters[12] 
    sigma_z = m.parameters[13]

    eq_y = 1
    eq_pi = 2
    eq_ffr = 3

    #  number of observation variables 

    ny = 3

    #  model variable indices 

    y_t   = 1
    pi_t   = 2
    R_t   = 3
    y1_t  = 4
    g_t   = 5
    z_t   = 6
    Ey_t1   = 7
    Epi_t1  = 8

    #  shock indices 

    z_sh = 1
    g_sh = 2
    R_sh = 3

    DD[eq_y,1] = gammaQ
    DD[eq_pi,1] = piA
    DD[eq_ffr,1] = piA+rA+4*gammaQ

    ZZ[eq_y,y_t] =  1
    ZZ[eq_y,y1_t] =  -1 
    ZZ[eq_y, z_t] = 1

    ZZ[eq_pi,pi_t] =  4

    ZZ[eq_ffr,R_t] = 4

    # with measurement errors (from dsge1_me.yaml)
    HH = zeros(ny,ny)  
    HH[eq_y, y_t] = (0.20*0.579923)^2
    HH[eq_pi, pi_t] = (0.20*1.470832)^2
    HH[eq_ffr, R_t] = (0.20*2.237937)^2

    #Don't square h
    QQ[z_sh,z_sh] = (sigma_z)^2
    QQ[g_sh,g_sh] = (sigma_g)^2
    QQ[R_sh,R_sh] = (sigma_R)^2
    
    VV    = QQ*MM'
    VVall = [[RRR*QQ*RRR' RRR*VV];
             [VV'*RRR'    HH]]

    #VVall *= 0

    return Measurement(ZZ, DD, QQ, EE, MM, VVall)
end
