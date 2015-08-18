function hessdiag = hessizero_saveall(fcn,x,Verbose,varargin)
% computes hessian of function fcn (string) evaluated at x (vector)
% varargin are the other inputs of fcn
% if Verbose, display error messages , results, etc.
% 11/12/01 translated by Marco DelNegro in matlab from Frank Schorfheide's program in gauss

%% index of free parameters
para_free = 1-x(:,2);
fpara_free = find(para_free);
nfree = length(fpara_free);

%% actual max
x = x(:,1);

npara = length(x);
ndx = 6;
dx =  exp(-(6:2:(6+(ndx-1)*2))');
hessdiag = zeros(npara, npara, ndx);



% Compute Diagonal elements first
for seli = fpara_free'

	if Verbose; fprintf(1,'\n\n Hessian Element: (%2.2g %2.2g)',[seli,seli]); end;
	for i=1:ndx;
		paradx = x;
		parady = x;
		paradx(seli) = paradx(seli) + dx(i);
		parady(seli) = parady(seli) - dx(i);
        
        fx  = eval([fcn '(x,varargin{:})']);
		fdx = eval([fcn '(paradx,varargin{:})']);
		fdy = eval([fcn '(parady,varargin{:})']);
		hessdiag(seli, seli, i) = -( 2*fx - fdx - fdy)/(dx(i))^2; 
    end
  	if Verbose == 2; fprintf(1,'\n Values: %2.6f',hessdiag(seli, seli, :));     end;
end

% Now compute off-diagonal elements
for II = 1:(nfree-1);
   seli = fpara_free(II);	
   for JJ = II+1:nfree;
	selj = fpara_free(JJ);
    	if Verbose; fprintf(1,'\n\n Hessian Element: (%2.2g %2.2g)',[seli,selj]); end;
	    for i=1:ndx;
			paradx = x;
			parady = x;
			paradx(seli) = paradx(seli) + dx(i);
			parady(selj) = parady(selj) - dx(i);
			paradxdy = paradx;
			paradxdy(selj) = paradxdy(selj) - dx(i);
    	
            fx  = eval([fcn '(x,varargin{:})']);
    		fdx = eval([fcn '(paradx,varargin{:})']);
    		fdy = eval([fcn '(parady,varargin{:})']);
            fdxdy = eval([fcn '(paradxdy,varargin{:})']);
			hessdiag(seli, selj, i) = -( fx - fdx - fdy + fdxdy )/(dx(i)*dx(i));
            hessdiag(selj, seli, i) = hessdiag(seli, selj, i);
        end
    	if Verbose == 2; fprintf(1,'\n Values: %2.6f',hessdiag(seli, selj, :)); end;
    end
end

