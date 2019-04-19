"""
```
ct_kalman_filter(y, T, R, C, Q, Z, D, E, tspan, s_0 = Vector(), P_0 = Matrix();
    outputs = [:loglh, :pred, :filt], Nt0 = 0)
tspan - Number of time units betwen each observation
"""
function ct_kalman_filter(y::Matrix{S},
    T::Matrix{S}, R::Matrix{S}, C::Vector{S},
    Q::Matrix{S}, Z::Matrix{S}, D::Vector{S}, E::Matrix{S}, tspan::Float64,
    s_0::Vector{S} = Vector{S}(undef, 0), P_0::Matrix{S} = Matrix{S}(undef, 0, 0);
    outputs::Vector{Symbol} = [:loglh, :pred, :filt],
    Nt0::Int = 0, method = Tsit5(),
    reltol::Float64 = 1e-8, abstol::Float64 = 1e-8) where {S<:AbstractFloat}

    # Determine outputs
    return_loglh = :loglh in outputs
    return_pred  = :pred  in outputs
    return_filt  = :filt  in outputs

    # Dimensions
    Ns = size(T, 1) # number of states
    Nt = size(y, 2) # number of periods of data

    # Initialize inputs and outputs
    k = StateSpaceRoutines.KalmanFilter(T, R, C, Q, Z, D, E, s_0, P_0)

    mynan = convert(S, NaN)
    loglh  = return_loglh ? fill(mynan, Nt)         : Vector{S}(0)
    s_pred = return_pred  ? fill(mynan, Ns, Nt)     : Matrix{S}(0, 0)
    P_pred = return_pred  ? fill(mynan, Ns, Ns, Nt) : Array{S, 3}(0, 0, 0)
    s_filt = return_filt  ? fill(mynan, Ns, Nt)     : Matrix{S}(0, 0)
    P_filt = return_filt  ? fill(mynan, Ns, Ns, Nt) : Array{S, 3}(0, 0, 0)

    # Populate initial states
    s_0 = k.s_t
    P_0 = k.P_t

    # Loop through periods t
    for t = 1:Nt
        # Forecast
        forecast!(k, tspan; method = method, reltol = reltol, abstol = abstol)
        if return_pred
            s_pred[:,    t] = k.s_t
            P_pred[:, :, t] = k.P_t
        end

        # Update and compute log-likelihood
        StateSpaceRoutines.update!(k, y[:, t]; return_loglh = return_loglh)
        if return_filt
            s_filt[:,    t] = k.s_t
            P_filt[:, :, t] = k.P_t
        end
        if return_loglh
            loglh[t]        = k.loglh_t
        end

        # Update s_0 and P_0 if Nt0 > 0
        if t == Nt0
            s_0 = k.s_t
            P_0 = k.P_t
        end
    end

    # Populate final states
    s_T = k.s_t
    P_T = k.P_t

    # Remove presample periods
    loglh, s_pred, P_pred, s_filt, P_filt =
        StateSpaceRoutines.remove_presample!(Nt0, loglh, s_pred, P_pred, s_filt, P_filt; outputs = outputs)

    return loglh, s_pred, P_pred, s_filt, P_filt, s_0, P_0, s_T, P_T
end


"""
```
ct_kalman_filter(y, T, R, C, Q, Z, D, E, tspan, s_0 = Vector(), P_0 = Matrix();
    outputs = [:loglh, :pred, :filt], Nt0 = 0)
tspan - Number of time units betwen each observation
n_subinterval - Number of subintervals to use ODE integration with
"""
function ct_kalman_filter(y::Matrix{S},
                          T::Matrix{S}, R::Matrix{S}, C::Vector{S}, Q::Matrix{S},
                          Z::Matrix{S}, D::Vector{S}, E::Matrix{S}, tspan::Float64, n_subinterval::Int64,
                          s_0::Vector{S} = Vector{S}(0), P_0::Matrix{S} = Matrix{S}(0, 0);
                          outputs::Vector{Symbol} = [:loglh, :pred, :filt],
                          Nt0::Int = 0, method = Tsit5(),
                          reltol::Float64 = 1e-8, abstol::Float64 = 1e-8) where {S<:AbstractFloat}

    # Determine outputs
    return_loglh = :loglh in outputs
    return_pred  = :pred  in outputs
    return_filt  = :filt  in outputs

    # Dimensions
    Ns = size(T, 1) # number of states
    Nt = size(y, 2) # number of periods of data

    # Initialize inputs and outputs
    k = StateSpaceRoutines.KalmanFilter(T, R, C, Q, Z, D, E, s_0, P_0)

    mynan = convert(S, NaN)
    loglh  = return_loglh ? fill(mynan, Nt)         : Vector{S}(0)
    s_pred = return_pred  ? fill(mynan, Ns, Nt)     : Matrix{S}(0, 0)
    P_pred = return_pred  ? fill(mynan, Ns, Ns, Nt) : Array{S, 3}(0, 0, 0)
    s_filt = return_filt  ? fill(mynan, Ns, Nt)     : Matrix{S}(0, 0)
    P_filt = return_filt  ? fill(mynan, Ns, Ns, Nt) : Array{S, 3}(0, 0, 0)

    # Populate initial states
    s_0 = k.s_t
    P_0 = k.P_t

    # Loop through periods t
    for t = 1:Nt
        # Forecast
        forecast!(k, tspan, n_subinterval; method = method, reltol = reltol, abstol = abstol)
        if return_pred
            s_pred[:,    t] = k.s_t
            P_pred[:, :, t] = k.P_t
        end

        # Update and compute log-likelihood
        StateSpaceRoutines.update!(k, y[:, t]; return_loglh = return_loglh)
        if return_filt
            s_filt[:,    t] = k.s_t
            P_filt[:, :, t] = k.P_t
        end
        if return_loglh
            loglh[t]        = k.loglh_t
        end

        # Update s_0 and P_0 if Nt0 > 0
        if t == Nt0
            s_0 = k.s_t
            P_0 = k.P_t
        end
    end

    # Populate final states
    s_T = k.s_t
    P_T = k.P_t

    # Remove presample periods
    loglh, s_pred, P_pred, s_filt, P_filt =
        remove_presample!(Nt0, loglh, s_pred, P_pred, s_filt, P_filt; outputs = outputs)

    return loglh, s_pred, P_pred, s_filt, P_filt, s_0, P_0, s_T, P_T
end

"""
Forecast functions
"""
function forecast!(k::StateSpaceRoutines.KalmanFilter, tspan::Float64; method = Tsit5(),
                   reltol::Float64 = 1e-8, abstol::Float64 = 1e-8) where {S<:AbstractFloat}
    T, R, Q = k.T, k.R, k.Q
    s0, P0 = k.s_t, k.P_t
    f(u,p,t) = T*u
    g(u,p,t) = T*u + u*T' + R*Q*R'
    tspan_int = (0., tspan)
    s_prob = ODEProblem(f, s0, tspan_int)
    P_prob = ODEProblem(g, P0, tspan_int)
    s_sol = DifferentialEquations.solve(s_prob, method, reltol = reltol, abstol = abstol)
    P_sol = DifferentialEquations.solve(P_prob, method, reltol = reltol, abstol = abstol)

    k.s_t = s_sol.u[end]
    k.P_t = P_sol.u[end]
end

# INCOMPLETE
# Not sure how to write down the ODE system correctly since I don't have the same
# kind of recursion that I can do otherwise. Moreover, I effectively will get the
# same predicted states regardless, so I'm not sure if this will improve the
# measurement equation or accuracy.
function forecast!(k::StateSpaceRoutines.KalmanFilter, tspan::Float64, n_subinterval::Int64; method = Tsit5(),
                   reltol::Float64 = 1e-8, abstol::Float64 = 1e-8) where {S<:AbstractFloat}
    T, R, Q = k.T, k.R, k.Q
    s_filt, P0 = k.s_t, k.P_t
    indiv_dim = size(T, 1)
    s0 = s_filt[1:indiv_dim]
    f(u,p,t) = T*u
    g(u,p,t) = T*u + u*T' + R*Q*R'
    tspan_int = (0., span)
    s_prob = ODEProblem(f, s0, tspan_int)
    P_prob = ODEProblem(g, P0, tspan_int)
    s_sol = solvect(s_prob, method, reltol = reltol, abstol = abstol)
    P_sol = solvect(P_prob, method, reltol = reltol, abstol = abstol)
    k.s_t = s_sol.u
    k.P_t = P_sol.u
end
