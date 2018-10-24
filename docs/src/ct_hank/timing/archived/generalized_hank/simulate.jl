"""
simulate(T, R, periods, steps, shocks; method = :explicit)

"""
function simulate(T::Matrix, R::Array, periods::Int64, steps::Int64, shocks::Matrix;
                  method::Symbol = :explicit, transformation::Matrix = eye(size(T,1)),
                  to_extract::Vector{Int64} = collect(1:size(transformation, 1)))

    @assert steps >= periods "number of steps must exceed number of periods"
    dt = periods / (steps - 1)

    n_vars = size(T, 1)
    values = zeros(n_vars, steps+1)

    if method == :explicit
        T_scaled = sparse(eye(T)) + T*dt
        for n in 1:steps
            values[:,n+1] = T_scaled * (values[:,n] + sqrt(dt) * R * shocks[:,n])
        end
        values = values[:,2:end]
    elseif method == :implicit
        T_inv = inv(sparse(eye(T)) - T*dt)
        for n in 1:steps
            values[:,n+1] = T_inv * (values[:,n] + sqrt(dt) * R * shocks[:,n])
        end
        values = values[:,2:end]
    else
        error("method must be :explicit or :implicit")
    end
    init_mat = transformation*values
    if length(to_extract) == size(init_mat, 1)
        out_mat = init_mat
    else
        out_mat = zeros(length(to_extract), size(init_mat, 2)) # preallocate memory
        for i in 1:length(to_extract)
            out_mat[i, :] = init_mat[to_extract[i], :]
        end
    end
    return out_mat
end