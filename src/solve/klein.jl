function klein(m::AbstractModel)

    nw::Int64 = get_setting(m, :nw)
    normalize_distr::Bool = get_setting(m, :normalize_distr_variables)

    #################
    # Linearization:
    #################
    JJ = jacobian(m)

    #################
    # Normalization:
    #################
    Jac1 = normalize_matrices(nw, JJ, normalize_distr)

    ##########################
    # Klein Solution Method---apply generalized Schur decomposition a la Klein (2000)
    ##########################

    # NK is number of predetermined variables
    NK = 2*nw + 1
    # n is number of variables (predet + non-predet)
	n = size(Jac1,1)

    # A and B matrices
    A = Array{Float64, 2}(size(Jac1,2), n)
	A::Array{Float64, 2} = Jac1[:,1:n]
	B = -Jac1[:,n+1:2*n]

    # Apply generlaized Schur decomposition
    # A ≈ QZ[:Q]*QZ[:S]*QZ[:Z]'
    # B ≈ QZ[:Q]*QZ[:T]*QZ[:Z]'
	QZ = schurfact(A,B)

    # Reorder so that stable comes first
    alpha::Vector{complex(promote_type(eltype(A), eltype(B)))} = QZ.alpha
    beta::Vector{complex(promote_type(eltype(A), eltype(B)))} = QZ.beta
	eigs = real(QZ.beta ./ QZ.alpha)
	eigselect::AbstractArray{Bool} = abs.(eigs) .< 1 # returns false for NaN gen. eigenvalue which is correct here bc they are > 1
	ordschur!(QZ, eigselect)

    # Check that number of stable eigenvalues equals the number of predetermined state variables
	nk = sum(eigselect)
	if nk>NK
	    warn("Equilibrium is locally indeterminate")
	elseif nk<NK
	    warn("No local equilibrium exists")
	end


	U::Array{Float64, 2}= QZ[:Z]'
	T::Array{Float64, 2} = QZ[:T]
	S::Array{Float64, 2} = QZ[:S]

    U11 = Array{Float64, 2}(NK, NK)
    U12 = Array{Float64, 2}(NK, NK-2)
    U21 = Array{Float64, 2}(NK-2, NK)
    U22 = Array{Float64, 2}(NK-2, NK-2)
    S11 = Array{Float64, 2}(NK, NK)
    T11 = Array{Float64, 2}(NK, NK)

    U11 = U[1:NK,1:NK]
	U12 = U[1:NK,NK+1:end]
	U21 = U[NK+1:end,1:NK]
	U22 = U[NK+1:end,NK+1:end]

	S11 = S[1:NK,1:NK]
	T11 = T[1:NK,1:NK]

    # Find minimum norm solution to U₂₁ + U₂₂*g_x = 0 (more numerically stable than -U₂₂⁻¹*U₂₁)
    gx_coef = Array{Float64, 2}(n-NK, NK)
	gx_coef = -U22'*pinv(U22*U22')*U21

    #  Solve for h_x (in a more numerically stable way)
	S11invT11 = S11\T11;
	Ustuff = (U11 + U12*gx_coef);
	invterm = pinv(eye(NK)+gx_coef'*gx_coef);
    #hx_coef = Array{Float64, 2}(NK, NK)
	hx_coef = invterm*Ustuff'*S11invT11*Ustuff;

	# Ensure that hx and S11invT11 should have same eigenvalues
	#(eigst,valst) = eig(S11invT11);
	#(eighx,valhx) = eig(hx_coef);
	eigst = eigvals(S11invT11)
    eighx = eigvals(hx_coef)
    if abs.(maximum(abs.(eighx))-maximum(abs.(eigst)))>1e-4
		warn("max abs eigenvalue of S11invT11 and hx are different!")
	end

	# next, want to represent policy functions in terms of meaningful things
	# gx_fval = Qy'*gx_coef*Qx
	# hx_fval = Qx'*hx_coef*Qx

	return gx_coef, hx_coef
end

function normalize_matrices(nw::Int64, JJ::Array{Float64, 2}, normalize_distr::Bool)
if normalize_distr
    P = eye(nw)
    P[:,1] = ones(nw)

    Q = Array{Float64, 2}(nw, nw)
    Q::Array{Float64, 2}, _ = qr(P)

    S = Array{Float64, 2}(nw, nw-1)
    S::Array{Float64,2} = Q[:, 2:end]

    Qf = zeros(4*nw, 4*nw+2)
    Qf[1:nw, 1:nw] = eye(nw)
    Qf[nw+1:2*nw-1, nw+1:2*nw] = S'
    Qf[2*nw:3*nw-2, 2*nw+1:3*nw] = S'
    Qf[3*nw-1:4*nw, 3*nw+1:4*nw+2] = eye(nw+2)
    #Qf = cat([1 2],eye(nw),S,S,eye(nw),[1],[1])

    # Qx and Qy are for normalizing any variables that represent distributions
    # the S component is for normalizing a distribution, the other identity portions
    # are for non-distributional variables
    # The distribution in the states, x, is the lagged μ
    # The distribution in the jumps, y, is the current μ
    # Qx = cat([1 2], S, eye(nw), [1], [1]). 161 x 162
    # Qy = cat([1 2], S, eye(nw)). 159 x 160

    # whenever outputs and/or inputs are densities,
    # pre/postmultiply the jacobians from step 2 by an appropriate
    # matrix so things integrate to 1

    Qright = zeros(8*nw+4, 8*nw)
    Qright[1:nw, 1:nw-1] = S
    Qright[nw+1:2*nw+2, nw:2*nw+1] = eye(nw+2)
    Qright[2*nw+3:3*nw+2, 2*nw+2:3*nw] = S
    Qright[3*nw+3:4*nw+2, 3*nw+1:4*nw] = eye(nw)
    Qright[4*nw+3:5*nw+2, 4*nw+1:5*nw-1] = S
    Qright[5*nw+3:6*nw+4, 5*nw:6*nw+1] = eye(nw+2)
    Qright[6*nw+5:7*nw+4, 6*nw+2:7*nw] = S
    Qright[7*nw+5:8*nw+4, 7*nw+1:8*nw] = eye(nw)
    #Qright = cat([1,2],Qx',Qy',Qx',Qy')
    Jac1 = Qf*JJ*Qright
else
    Jac1 = JJ
end
    return Jac1
end

# Need an additional transition_equation function to properly stack the
# individual state and jump transition matrices/shock mapping matrices to
# a single state space for all of the model_states
function klein_transition_matrices(m::AbstractModel,
                                   TTT_state::Matrix{Float64}, TTT_jump::Matrix{Float64})
    TTT = zeros(n_model_states(m), n_model_states(m))

    # Loading mapping time t states to time t+1 states
    TTT[1:n_backward_looking_states(m), 1:n_backward_looking_states(m)] = TTT_state

    # Loading mapping time t jumps to time t+1 states
    TTT[1:n_backward_looking_states(m), n_backward_looking_states(m)+1:end] = 0.

    # Loading mapping time t states to time t+1 jumps
    TTT[n_backward_looking_states(m)+1:end, 1:n_backward_looking_states(m)] = TTT_jump*TTT_state

    # Loading mapping time t jumps to time t+1 jumps
    TTT[n_backward_looking_states(m)+1:end, n_backward_looking_states(m)+1:end] = 0.

    RRR = shock_loading(m, TTT_jump)

    return TTT, RRR
end
