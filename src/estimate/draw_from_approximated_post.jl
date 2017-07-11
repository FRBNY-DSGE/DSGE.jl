function draw_from_approximated_post(m::AbstractModel, data::Array{Float64}, dist::Distribution)
    success=false
    n_params = n_parameters(m)
    T = typeof(m.parameters[1].value)
    draw = Array{T}(n_params)
    logpost = 0
    loglh = 0
    while !success
        for i in 1:n_params
            draw[i] = if !m.parameters[i].fixed
                prior = rand(dist, m)
                while !(m.parameters[i].valuebounds[1] <prior <m.parameters[i].valuebounds[2])
                    prior = rand(dist,m)
                end
                prior
            else
                m.parameters[i].value
            end
         end
         try 
             logpost = posterior!(m,convert(Array{Float64,1},draw), data)
             loglh = logpost - prior(m)
             success = true
         end
    end
    return (draw, logpost, loglh)
end
