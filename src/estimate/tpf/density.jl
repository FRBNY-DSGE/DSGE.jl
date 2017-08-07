"""
```
density{S<:Float64}(φ_new::S, φ_old::S, y_t::Array{S,1}, p_error::Array{S,1}, 
    EE::Array{S,2}; initialize::Bool=false)
```
### Input



### Output

Returns the appropriate probability (right now from a Normal distribution) evaluated at p_error.
"""
function density{S<:Float64}(φ_new::S, φ_old::S, y_t::Array{S,1}, p_error::Array{S,1}, 
                                   EE::Array{S,2}; initialize::Bool=false)

    # Initialization step (using 2π instead of φ_old)
    if initialize
        return (φ_new/(2*pi))^(length(y_t)/2) * (det(EE)^(-1/2)) * exp((-1/2)*p_error'*φ_new*inv(EE)*p_error)[1]
    # Non-initialization step (tempering and final iteration)
    else
        return (φ_new/φ_old)^(length(y_t)/2) * exp(-1/2*p_error'*(φ_new - φ_old)*inv(EE)*p_error)[1]
    end
 end