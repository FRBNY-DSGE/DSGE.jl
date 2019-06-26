
function solve(Jac0, Qleft, Qx, Qy)
	# step756: whenever outputs and/or inputs are densities,
	# pre/postmultiply the jacobians from step 2 by an appropriate
	# matrix so things integrate to 1
	Qright = cat([1,2],Qx',Qy',Qx',Qy')
	Jac1 = Qleft*Jac0*Qright
	#println("finished step 3")

	# apply generalized Schur decomposition a la Klein (2000)
	n = size(Jac1,1)
	A =  Jac1[:,1:n]
	B = -Jac1[:,n+1:2*n]
	QZ        = schur(A,B)
	eigs      = QZ.β ./ QZ.α
	eigselect = abs.(eigs) .< 1 # returns false for NaN gen. eigenvalues
	                           # which is correct here bc they are > 1
	QZ = ordschur(QZ,eigselect)
	# note: because we did (A,B), beta gives eigs corresponding to B

	NK = size(Qx,1) # number of predetermined vars. NOTE: this is after
					# removing elements corresponding to densities
	nk = sum(convert(Array{Int64,1},eigselect)) # number of stable eigs
	println("NK = ",NK)
	println("nk = ",nk)
	if nk>NK
	    warn("equilibrium is locally indeterminate")
	elseif nk<NK
	    warn("no local equilibrium exists")
	end

	U = QZ[:Z]'
	T = QZ[:T]
	S = QZ[:S]
	U11 = U[1:NK,1:NK]
	U12 = U[1:NK,NK+1:end]
	U21 = U[NK+1:end,1:NK]
	U22 = U[NK+1:end,NK+1:end]
	S11 = S[1:NK,1:NK]
	T11 = T[1:NK,1:NK]

	gx_coef = - U22'*pinv(U22*U22')*U21
	# this is the minimum norm solution to U21 + U22*gx = 0
	# in principle one could calculate gx as -pinv(U22)*U21
	# however, the minimum norm solution is more numerically stable
	S11invT11 = S11\T11;
	Ustuff = (U11 + U12*gx_coef);
	invUstuff = Ustuff\eye(NK);
	invterm = pinv(eye(NK)+gx_coef'*gx_coef);
	hx_coef = invterm*Ustuff'*S11invT11*Ustuff;
	# again, in principle there are other ways to calculate hx
	# but this one is the most numerically stable
	# check: hx and S11invT11 should have same eigenvalues
	(eigst,valst) = eig(S11invT11);
	(eighx,valhx) = eig(hx_coef);
	if abs.(maximum(abs.(eighx))-maximum(abs.(eigst)))>1e-4
		warn("max abs eigenvalue of S11invT11 and hx are different!")
	end

	# next, want to represent policy functions in terms of meaningful things

	gx_fval = Qy'*gx_coef*Qx
	hx_fval = Qx'*hx_coef*Qx

	return gx_fval, hx_fval, gx_coef, hx_coef
end
