using Base.Test, Compat, HDF5, DSGE
#using MATLAB      #Comment out MATLAB to suppress warnings, still need package to convert files 

# minusnan(x, y) evaluates x-y in a way that treats NaN like Inf and sets Inf - Inf = 0
minusnan{S<:@compat(AbstractFloat), T<:@compat(AbstractFloat)}(x::S, y::T) =  minusnan(complex(x), complex(y))

function minusnan{S<:AbstractFloat, T<:AbstractFloat}(x::Complex{S}, y::Complex{T})
    x = isnan(x) || isinf(x) ? Inf : x
    y = isnan(y) || isinf(y) ? Inf : y
    return isinf(x) && isinf(y) ? 0 : x - y
end

# percenterr(x,y) returns the absolute value of the percentage error of y wrt x
function percenterr{S<:AbstractFloat, T<:AbstractFloat}(x::Complex{S}, y::Complex{T})
    return abs((x-y)/x)
end

# checksigns(x,y) checks to see if x and y have the same sign. returns 1 if they do, 0 if not, and -1 if one is zero but the other isnt
function checksigns(x::Number,y::Number)
    sx = sign(x)
    sy = sign(y)

    if sx == sy 
        return 1
    elseif isnan(x) && isnan(y)
        
        return 1
        
    elseif sx == 0 || sy == 0 # one is zero
        return -1
        
    else
        return 0
    end
end

# TODO: decide what a sensible default ε value is for our situation
# Compares matrices, reports absolute differences, returns true if all entries close enough
function test_matrix_eq{R<:AbstractFloat, S<:AbstractFloat, T<:AbstractFloat}(expected::Array{R}, actual::Array{S}; ε::T = 1e-4, noisy::Bool = false)

    n_entries = length(expected)
    same_sign = map(checksigns, expected, actual)
    diff_sign_zeros = count(x-> x < 0, same_sign) # # that are different because 1 sign is 0
    n_diff_sign = count(x -> x != 1, same_sign)

    if noisy
        println("$n_diff_sign of $n_entries entries have opposite signs")
        println("$diff_sign_zeros of $n_entries entries have different signs, but one of the entries is 0")
    end

    return test_matrix_eq(complex(expected), complex(actual); ε=ε, noisy=noisy)
end

# Complex-valued input matrices

function test_matrix_eq{R<:AbstractFloat, S<:AbstractFloat, T<:AbstractFloat}(expected::Array{Complex{R}}, actual::Array{Complex{S}}; ε::T = 1e-4, noisy::Bool = false)

    # Matrices of different sizes return false
    if size(expected) != size(actual)
        if noisy
            println("Size expected $(size(expected)), actual $(size(actual))\n")
        end
        return false
    end

    n_dims = ndims(expected)
    n_entries = length(expected)

    # Count differences and find max
    abs_diff = abs(map(minusnan, expected, actual))
    n_neq = countnz(abs_diff)
    n_not_approx_eq = count(x -> x > ε, abs_diff)
    max_abs_diff = maximum(abs_diff)
    max_inds = ind2sub(size(abs_diff), indmax(abs_diff))

    # Count percent differences and find max
    pct_err = map(percenterr, expected, actual) #percenterr returns an absolute value
    max_pct_err = maximum(pct_err)
    max_pct_inds = ind2sub(size(pct_err), indmax(pct_err))
    
    # Print output
    if noisy
        println("$n_neq of $n_entries entries with abs diff > 0")
        println("$n_not_approx_eq of $n_entries entries with abs diff > $ε")
        # println("$n_diff_sign of $n_entries entries have opposite signs")
        if n_neq != 0
            println("Max abs diff of $max_abs_diff at entry $max_inds")
            if isa(expected, Matrix)
                println("The entries at $max_inds are $(expected[max_inds[1],max_inds[2]]) (expected) and $(actual[max_inds[1], max_inds[2]]) (actual)")
            end
            println("Max percent error of $max_pct_err at entry $max_pct_inds")
            if isa(expected, Matrix)
                println("The entries at $max_pct_inds are $(expected[max_pct_inds[1],max_pct_inds[2]]) (expected) and $(actual[max_pct_inds[1], max_pct_inds[2]]) (actual)\n")
            end
        else
            println("Max abs diff of 0\n")
        end
    end

    # Return true if all entries within ε
    return n_not_approx_eq == 0
end

function test_matrix_abs_eq{R<:AbstractFloat, S<:AbstractFloat, T<:AbstractFloat}(expected::
             Array{R}, actual::Array{S}; ε::T = 1e-4, noisy::Bool = false)
    # Matrices of different sizes return false
    if size(expected) != size(actual)
        if noisy
            println("Size expected $(size(expected)), actual $(size(actual))\n")
        end
        return false
    end
    
    n_dims = ndims(expected)
    n_entries = length(expected)

    expected = abs(expected)
    actual = abs(actual)
    
    # Count differences and find max
    abs_diff = abs(map(minusnan, expected, actual))
    n_neq = countnz(abs_diff)
    n_not_approx_eq = count(x -> x > ε, abs_diff)
    max_abs_diff = maximum(abs_diff)
    max_inds = ind2sub(size(abs_diff), indmax(abs_diff))
    
    # Print output
    if noisy
        println("$n_neq of $n_entries entries with abs diff > 0")
        println("$n_not_approx_eq of $n_entries entries with abs diff > $ε")
        if n_neq != 0
            println("Max abs diff of $max_abs_diff at entry $max_inds\n")
        else
            println("Max abs diff of 0\n")
        end
    end

    # Return true if all entries within ε
    return n_not_approx_eq == 0
end

#function test_eigs_eq{R<:AbstractFloat, S<:AbstractFloat, T<:AbstractFloat, U<:AbstractFloat, V<:AbstractFloat}(eigvals_expected::Array{R}, eigvals_actual::Array{S}, eigvecs_expected::Array{U},  eigvecs_actual::Array{V}; ε::T = 1e-12, noisy::Bool = true)

function test_eigs_eq(eigvals_expected::Array, eigvals_actual::Array, eigvecs_expected::Array,  eigvecs_actual::Array; ε = 1e-12, noisy::Bool = true)
   
    # check to see if sizes correspond
    if size(eigvecs_expected) != size(eigvecs_actual)
        if noisy
            println("Eigenvectors matrix size expected $(size(eigvecs_expected)), actual $(size(eigvecs_actual))\n")
        end
        return false
    end

    if size(eigvals_expected) != size(eigvals_actual)
        if noisy
            println("Eigenvalues vector size expected $(size(eigvals_expected)), actual $(size(eigvals_actual))\n")
        end
        return false
    end
    
    # Take absolute values of everything
    abs_eigvals_expected = abs(eigvals_expected);
    abs_eigvals_actual   = abs(eigvals_actual);

    abs_eigvecs_expected = abs(eigvecs_expected);
    abs_eigvecs_actual   = abs(eigvecs_actual);

    ## @test abs_eigvecs_actual .> 0
    ## @test abs_eigvecs_expected .> 0
    
    
    # Stick the eigenvalues on top of the eigenvectors, sort the columns of the matrix, and then separate
    # the eigenvals and vecs back out again
    
    expected = [abs_eigvals_expected'; abs_eigvecs_expected]
    actual   = [abs_eigvals_actual'; abs_eigvecs_actual]

    expected_sorted = sortcols(expected)
    actual_sorted   = sortcols(actual)

    eigvals_expected_sorted = expected_sorted[1,:] 
    eigvals_actual_sorted   = actual_sorted[1,:]
    
    eigvecs_expected_sorted = expected_sorted[2:end,:]
    eigvecs_actual_sorted = actual_sorted[2:end,:]
    
    # See if expected eigenvalues/vectors == actual eigenvalues/vectors
    vals_eq = test_matrix_eq(eigvals_expected_sorted, eigvals_actual_sorted,ε = ε, noisy=true);
    vecs_eq = test_matrix_eq(eigvecs_expected_sorted, eigvecs_actual_sorted, ε = ε, noisy=true)

    # return vals_eq && vecs_eq
    return eigvals_expected_sorted, eigvals_actual_sorted, eigvecs_expected_sorted, eigvecs_actual_sorted
    
end

function test_test_eigs_eq()
    A_eigvals = [0.4, -6.2, 3.4, 4.7]
    B_eigvals = [7.3, 2.4, -0.3, 10.2 ]

    A_eigvecs = [.2   -1.0   .4 -100.0;
                 -.5   0.0  6.4    8.5;
                 0.0   2.0  -1.0   1.4;
                 5.0   5.0  7.0    8.0]
    
    B_eigvecs = [-.5   3.7   .4 -100.0;
                 0.0   2.0  -1.0   1.4;
                 5.0  -8.0  7.0   -4.0;
                 -.5   0.0  6.4    8.5]

    @test !test_eigs_eq(A_eigvals, B_eigvals, A_eigvecs, B_eigvecs)

    E_eigvals = [.1, .4, .3, .2]
    F_eigvals = [.2, .4, .1, .3]
 

    E_eigvecs = [0.2   0.4   -0.9    -0.7;
                -0.3   0.5    0.8     0.0;
                 0.1   0.1    0.1     0.1;
                 3.7   4.2    8.8    -0.3]

    F_eigvecs = [-0.7  0.4    0.2    -0.9;
                  0.0  0.5   -0.3     0.8;
                  0.1  0.1    0.1     0.1;
                 -0.3  4.2    3.7     8.8 ]

    @test test_eigs_eq(E_eigvals, F_eigvals, E_eigvecs, F_eigvecs)

    C_eigvals = [.10, .20]
    D_eigvals = [.20, .10]

    C_eigvecs = [.10 -.30;
                 2.3 0.5]

    D_eigvecs = [-.30    .10 ;
                  0.5    2.3]

    
    @test test_eigs_eq(C_eigvals, D_eigvals, C_eigvecs, D_eigvecs)

    println("test_eigs_eq passed")
end


function find_matrix_diffs{T<:Number}(m0::Matrix{T}, m1::Matrix{T}; ε=1e-4)

    if size(m0) != size(m1)
        println("The matrices are not the same size.")
        return
    end

    n_entries = length(m0)
    
    abs_diff = abs(map(minusnan, m0, m1))
    rows, cols = size(abs_diff)


    # Get row, column indices for absolute differences
    v_indices = find(x-> x>ε, abs_diff)
    m_indices = ind2sub(size(abs_diff), v_indices)

    abs_diff_entries = m_indices
    #abs_diff_entries = collect(zip(m_indices[1],m_indices[2]))
    
    # Figure out which entries have diff signs
    same_sign = map(checksigns, m0, m1)
    
    # Get row, column indices for entries with opposite signs
    v_indices = find(x-> x != 1, same_sign)
    m_indices = ind2sub(size(abs_diff), v_indices)

    opp_sign_entries = m_indices
    #opp_sign_entries = collect(zip(m_indices[1], m_indices[2]))
    
    ## for (i,k) in enumerate(zipped)
    ##     println(k)
    ## end
    
    return abs_diff_entries, opp_sign_entries
    
end


function get_diagonal{T<:Number}(matrix::Matrix{T})
    rows,cols = size(matrix)

    if rows != cols
        println("Not a square matrix")
        return
    end

    diag = zeros(rows)
    for i in 1:rows
        diag[i] = matrix[i,i]
    end

    return diag
end

# Complex numbers are parsed weirdly from readcsv, so we build a complex array using regex
function readcsv_complex{T<:@compat(AbstractString)}(file::T)
    matrix_str = readcsv(file)
    rows, cols = size(matrix_str)
    matrix = complex(zeros(size(matrix_str)))

    for j = 1:cols
        for i = 1:rows
            value = matrix_str[i, j]
            if isa(value, AbstractString)
                m = match(r"(\d+\.?\d*e?-?\d*) ?([+-]) ?(\d+\.?\d*+e?-?\d*)im?", value)
                if m == nothing
                    error("Regex didn't match anything at ($i, $j) entry")
                end
                a = parse(m.captures[1])
                b = parse(m.captures[2] * m.captures[3])
                matrix[i, j] = complex(a, b)
            elseif isa(value, Number)
                matrix[i, j] = complex(value)
            else
                error("($i, $j) entry not a string or number")
            end
        end
    end
    return matrix
end



function test_util()
    test_test_matrix_eq()
    test_readcsv_complex()
    test_test_eigs_eq()
    println("### All tests in util.jl tests passed\n")
end


function test_test_matrix_eq()
    # Matrices of different sizes returns false
    m0 = zeros(2, 3)
    m1 = zeros(2, 2)
    @test !test_matrix_eq(m0, m1)

    # Returns true
    m2 = [0.0001 0.0; 0.0 -0.0001]
    @test test_matrix_eq(m1, m2)

    # Returns false
    m3 = [0.0001 0.0; 0.0 -0.0002]
    @test !test_matrix_eq(m1, m3)
    @test test_matrix_eq(m1, m3; ε=2e-4) # but true with larger ε

    # Arguments of different float type returns true
    m2_float16 = convert(Matrix{Float16}, m2)
    ε = convert(Float32, 0.0001)
    @test test_matrix_eq(m1, m2)

    # Complex-valued matrices
    @test test_matrix_eq(complex(m1), complex(m2))
    @test !test_matrix_eq(complex(m1), complex(m3))

    # 3D arrays
    m4 = zeros(2, 2, 2)
    m5 = ones(2, 2, 2)
    @test !test_matrix_eq(m4, m5)

    println("test_matrix_eq tests passed\n")
end


function test_readcsv_complex()
    readcsv_complex("readcsv_complex.csv")
    println("readcsv_complex tests passed\n")
end

# Convert a mat file to an hdf5 file. No compression is used.
# Caller must be using HDF5 and MATLAB packages
function mat_to_hdf5(matfileName, h5fileName)
    mf = MatFile(matfileName)
    vars = MATLAB.variable_names(mf)
    
    # Write each variable to h5 file
    h5 = HDF5.h5open(h5fileName,"w") do h5
    
        for name in vars
            var = MATLAB.get_variable(mf, name)
            HDF5.write(h5, name, var)
        end

    end

    MATLAB.close(mf)
end

function matvar_to_hdf5var(matfileName,matVarName,h5fileName,h5VarName)

    # Open mat file and get the desired variable
    mf = MatFile(matfileName)
    var = MATLAB.get_variable(mf, matVarName)
    
    # Open the HDF5 file and write the variable
    h5 = HDF5.h5open(h5fileName,"w")
    HDF5.write(h5, h5VarName, var)

    # Close the files
    HDF5.close(h5)
    MATLAB.close(mf)
end

function compare_mat_hdf5(matfile,h5file)
    mf = MatFile(matfile)
    mvars = MATLAB.variable_names(mf)
    
    # Compare variables
    h5 = HDF5.h5open(h5file,"w") do h5
        h5vars = HDF5.names(h5)
      
        for name in mvars
            if(in(h5vars,name))
                mvar  = MATLAB.get_variable(mf, name)
                h5var = HDF5.read(h5f, name)

                test_matrix_eq(mvar, h5var,noisy=true)
                
            else
                println("Missing: h5 file doesn't contain variable $name")
            end
            
            
        end

    end

    MATLAB.close(mf)
   
end

function compare_matvar_hdf5var(matfile, mname, h5file, h5name; ε=1e-4, noisy=true)
    #Get matlab variable
    mf = MatFile(matfile);
    mvar = MATLAB.get_variable(mf, mname);
    MATLAB.close(mf);    

    # get hdf5 variable
    h5 = HDF5.h5open(h5file,"r"); 
    h5var = HDF5.read(h5, h5name);
    HDF5.close(h5);
    

    if isa(h5var,Vector) || isa(h5var,Matrix)
        return mvar, h5var, test_matrix_eq(mvar, h5var, ε=ε, noisy=noisy);
    elseif isa(h5var,Float64) || isa(h5var,Float32)
        return mvar, h5var, mvar == h5var;
    end
    
end

function get_diff_symbols{S<:AbstractFloat, T<:AbstractDSGEModel}(m::T, expected::Matrix{S}, actual::Matrix{S}; ε=1e-4, noisy=true)

    diff_entries, opp_sign_entries = find_matrix_diffs(expected, actual; ε=ε)

    diff_states = Set(diff_entries[2])
    diff_eqcond = Set(diff_entries[1])

    diff_state_inds  = Dict{Int64,Symbol}()
    diff_eqcond_inds = Dict{Int64,Symbol}()

    println("States:")
    for (i,k) in enumerate(m.endogenous_states)
        j= k[2]
        if in(j, diff_states)
            get!(diff_state_inds, j, k[1])
            println(k)
        end
    end

    println("\n")
    println("Equilibrium conditions:")
    for (i,k) in enumerate(m.equilibrium_conditions)
        j= k[2]

        if in(j, diff_eqcond)
            get!(diff_eqcond_inds, j, k[1])
            println(k)
        end
    end
    
    println("\nEntries:")
    zipped = collect(zip(diff_entries[1], diff_entries[2]))
    for entry in zipped
        diff = abs(expected[entry[1],entry[2]] - actual[entry[1],entry[2]])
        @printf "%s, %s, %f, %f, %f\n" diff_eqcond_inds[entry[1]] diff_state_inds[entry[2]] diff expected[entry[1],entry[2]] actual[entry[1],entry[2]]
    end
    
end
