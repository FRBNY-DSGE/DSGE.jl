# Dispatching on max and min to allow comparison of Complex numbers to other numbers

import Base: max, min, >, <, =>, <=

function max(a::Complex, b::Complex)
    if iszero(a.im) && iszero(b.im)
        re_a = a.re
        re_b = b.re
        return re_a > re_b ? a : b
    else
        abs_a = abs(a)
        abs_b = abs(b)
        if abs_a == abs_b
            angle_a = angle(a)
            angle_b = angle(b)
            if max(angle_a, angle_b) == angle_a
                return a
            else
                return b
            end
        else
            if max(abs_a, abs_b) == abs_a
                return a
            else
                return b
            end
        end
    end
end

function min(a::Complex, b::Complex)
    if iszero(a.im) && iszero(b.im)
        re_a = a.re
        re_b = b.re
        return re_a < re_b ? a : b
    else
        abs_a = abs(a)
        abs_b = abs(b)
        if abs_a == abs_b
            angle_a = angle(a)
            angle_b = angle(b)
            if min(angle_a, angle_b) == angle_a
                return a
            else
                return b
            end
        else
            if min(abs_a, abs_b) == abs_a
                return a
            else
                return b
            end
        end
    end
end

function max{T<:Real}(a::Complex, b::T)
    abs_a = abs(a)
    if max(abs_a, b) == abs_a
        return a
    else
        return convert(Complex{Float64}, b)
    end
end

function min{T<:Real}(a::Complex, b::T)
    abs_a = abs(a)
    if min(abs_a, b) == abs_a
        return a
    else
        return convert(Complex{Float64}, b)
    end
end

function max{T<:Real}(a::T, b::Complex)
    abs_b = abs(b)
    if max(a, abs_b) == a
        return convert(Complex{Float64}, a)
    else
        return b
    end
end

function min{T<:Real}(a::T, b::Complex)
    abs_b = abs(b)
    if min(a, abs_b) == a
        return convert(Complex{Float64}, a)
    else
        return b
    end
end

# Greater than
function >(a::Complex, b::Complex)
    return a.re > b.re
end

function >{T<:Real}(a::T, b::Complex)
    return a > b.re
end

function >{T<:Real}(a::Complex, b::T)
    return a.re > b
end

# Less than
function <(a::Complex, b::Complex)
    return a.re < b.re
end

function <{T<:Real}(a::T, b::Complex)
    return a < b.re
end

function <{T<:Real}(a::Complex, b::T)
    return a.re < b
end

# Greater than or equal to
function =>(a::Complex, b::Complex)
    return a.re => b.re
end

function =>{T<:Real}(a::T, b::Complex)
    return a => b.re
end

function =>{T<:Real}(a::Complex, b::T)
    return a.re => b
end

# Less than or equal to
function <=(a::Complex, b::Complex)
    return a.re <= b.re
end

function <={T<:Real}(a::T, b::Complex)
    return a <= b.re
end

function <={T<:Real}(a::Complex, b::T)
    return a.re <= b
end
