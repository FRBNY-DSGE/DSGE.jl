function ct_kalman_simple(T, Z, Q, E, mean_0, var_0, data_y, dt)
    # Compute Kalman filter of a continuous time process sampled at discrete time
    # ported from SeHyoun Ahn
    #
    #   & dx_t = Tx_t dt + R dW_t \\
    #   & y_{t_k} = Zx_{t_k} + v_{t_k} \\
    #   & x_0 \sim N(mean_0, var_0) \\
    #   & k\in\{1, 2, \ldots, n\} \\
    #   & W_t = \text{Brownian Motion} \\
    #   & v_k \sim N(0, E) \\
    #   & dt_i = t_{i} - t_{i-1}
    #
    # Args:
    #   T (matrix): see equation above
    #   Z (matrix): see equation above
    #   Q (matrix): :math:`Q = RR'` i.e., variance
    #   E (matrix): see equation above
    #   mean_0 (vector): see equation above
    #   var_0 (matrix): see equation above
    #   data_y (vector): :math:`y_k` above
    #   dt (vector): time period between data measurement (see above)
    #
    # Note:
    #   - The data does not have to come with uniform step size.
    #   - There is a `kalman <https://www.mathworks.com/help/control/ref/kalman.html>`_ function in MATLAB as well, but it is under control system toolbox.
    #
    # Returns:
    #   : [mean_out, var_out, prob_out]
    #
    #   - **mean_out** (vector): :math:`E(x_{t_n}| y_{t_1},\ldots, y_{t_n})`
    #   - **var_out** (matrix): :math:`var(x_{t_n}| y_{t_1},\ldots, y_{t_n})`
    #   - **prob_out** (double): log-likelihood of the data
    #
    # Reference:
    #
    #   - `Wikipedia article on Kalman Filtering <https://en.wikipedia.org/wiki/Kalman_filter#Hybrid_Kalman_filter>`_
    #   - :cite:`grewal2011kalman`
    #
    # .. Syntax:
    #
    # ..   [mean_out, var_out, prob_out] = kalman_ct(T, Z, Q, E, mean_0, var_0, data_y, dt)
    #

    # Get data size
    n_data = size(data_y, 1)
    if (length(dt) < n_data) 
        dt = repeat(dt, n_data)
    end
    n_dim = length(mean_0)
    prob_out = 0
    mean_out = mean_0
    var_out = var_0

    # Step through data points

    for iter_data = 1:n_data

        #       tspan = (0.0,dt[iter_data][1]) 
        # s0 = mean_0 #initial values
        # P0 = var_out
        #        f(u,p,t) = T*u
        #     g(u,p,t) = T*u + u*T' + Q
        #           s_prob = ODEProblem(f, mean_out, tspan)
        #     P_prob = ODEProblem(g, P0, tspan)
        # sol = solve(s_prob,Tsit5(),reltol=1e-4,abstol=1e-4)

        ## forecasting state (this is different for ct)
        mean_out = exp(dt[iter_data][1]*T)*mean_out
        AB = exp(dt[iter_data][1]*[hcat(T, Q); hcat(zeros(n_dim,n_dim), -T')])*[var_out; Matrix{Float64}(I,n_dim,n_dim)]
        var_inter = AB[1:n_dim, :]/AB[n_dim+1:end, :]
        ##########

        ## the rest is the same
        var_obs = Z*var_inter*Z' .+ E

        v_t = data_y[iter_data,:] .- Z*mean_out
        prob_out = prob_out - .5*n_dim*log(2*pi) -.5*log(abs(det(var_obs))) - 0.5* v_t'*(var_obs\v_t)
        mean_out = mean_out .+ var_inter*Z'*(var_obs\v_t)
        var_out = var_inter .- var_inter*Z'*(var_obs\Z)*var_inter

    end

    return prob_out
end



#=

function [mean_out, var_out, prob_out] = kalman_ct(F, H, Q, R, mean_0, var_0, data_y, dt)
  % Compute Kalman filter of a continuous time process sampled at discrete time
  %
  % .. math::
  %
  %   & dx_t = Fx_t dt + J dW_t \\
  %   & y_{t_k} = Hx_{t_k} + v_{t_k} \\
  %   & x_0 \sim N(mean_0, var_0) \\
  %   & k\in\{1, 2, \ldots, n\} \\
  %   & W_t = \text{Brownian Motion} \\
  %   & v_k \sim N(0, R) \\
  %   & dt_i = t_{i} - t_{i-1}
  %
  % Args:
  %   F (matrix): see equation above
  %   H (matrix): see equation above
  %   Q (matrix): :math:`Q = JJ'` i.e., variance
  %   R (matrix): see equation above
  %   mean_0 (vector): see equation above
  %   var_0 (matrix): see equation above
  %   data_y (vector): :math:`y_k` above
  %   dt (vector): time period between data measurement (see above)
  %
  % Note:
  %   - The data does not have to come with uniform step size.
  %   - There is a `kalman <https://www.mathworks.com/help/control/ref/kalman.html>`_ function in MATLAB as well, but it is under control system toolbox.
  %
  % Returns:
  %   : [mean_out, var_out, prob_out]
  %
  %   - **mean_out** (vector): :math:`E(x_{t_n}| y_{t_1},\ldots, y_{t_n})`
  %   - **var_out** (matrix): :math:`var(x_{t_n}| y_{t_1},\ldots, y_{t_n})`
  %   - **prob_out** (double): log-likelihood of the data
  %
  % Reference:
  %
  %   - `Wikipedia article on Kalman Filtering <https://en.wikipedia.org/wiki/Kalman_filter#Hybrid_Kalman_filter>`_
  %   - :cite:`grewal2011kalman`
  %
  % .. Syntax:
  %
  % ..   [mean_out, var_out, prob_out] = kalman_ct(F, H, Q, R, mean_0, var_0, data_y, dt)
  %

  % Get data size
  n_data = size(data_y, 2);
  if (length(dt) < n_data) & (length(n_data) == 1)
    dt = repmat(dt, n_data, 1);
  end
  n_dim = length(mean_0);

  prob_out = 0;
  mean_out = mean_0;
  var_out = var_0;

  % Step through data points
  for iter_data = 1:n_data
    mean_out = expm(dt(iter_data)*F)*mean_out;
    AB = expm(dt(iter_data)*[F, Q; H'*(R\H), -F])*[var_out; speye(n_dim)];
    var_inter = AB(1:n_dim, :)*inv(AB((n_dim+1):end, :));

    var_obs = H*var_inter*H' + R;
    v_t = data_y(:, iter_data) - H*mean_out;

    prob_out = prob_out - n_dim*log(2*pi)/2 ...
               -(log(abs(det(var_obs))) + v_t'*(var_obs\v_t))/2;

    mean_out = mean_out + var_inter*H'*(var_obs\v_t);
    var_out = var_inter - (var_inter*H')*(var_obs\H)*var_inter;
  end
end




=#
