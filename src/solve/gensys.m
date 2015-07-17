function [G1,C,impact,fmat,fwt,ywt,gev,eu]=gensys(g0,g1,c,psi,pi,div)
%function [G1,C,impact,fmat,fwt,ywt,gev,eu]=gensys(g0,g1,c,psi,pi,div)
%System given as
%        g0*y(t)=g1*y(t-1)+c+psi*z(t)+pi*eta(t),
%with z an exogenous variable process and eta being endogenously determined
%one-step-ahead expectational errors.  Returned system is
%       y(t)=G1*y(t-1)+C+impact*z(t)+ywt*inv(I-fmat*inv(L))*fwt*z(t+1) .
% If z(t) is inlags.i.d., the last term drops out.
% If div is omitted from argument list, a div>1 is calculated.
% eu(1)=1 for existence, eu(2)=1 for uniqueness.  eu(1)=-1 for
% existence only with not-s.c. z; eu=[-2,-2] for coincident zeros.
% By Christopher A. Sims
% Updated by Marco Del Negro, Vasco Curdia, Daria Finocchiaro and others.

  qzerror_flag = 0;
  invwarn_flag = 0;
  
  eu=[0;0];
  realsmall=1e-6;
  fixdiv=(nargin==6);
  n=size(g0,1);

  if any(any(isinf(g0))) | any(any(isnan(g0))) | ...
     any(any(isinf(g1))) | any(any(isnan(g1)));

    %We have inf or nan in g0, g1
    [inang0,jnang0] = find(~isfinite(g0))    
    [inang1,jnang1] = find(~isfinite(g1))
    eu=[-4;-4];
    G1=[];
    C=[];
    impact=[];
    fmat=[];
    fwt=[];
    ywt=[];
    gev=[];
    return
  end  

  try;
    [a b q z v]=qz(g0,g1);
  catch exception;
    disp(exception.identifier); 
    disp(exception.message); 
    disp('QZ did not work ');
    qzerror_flag = 1;
  end;

  if qzerror_flag;
    eu=[-4;-4];
    G1=[];
    C=[];
    impact=[];
    fmat=[];
    fwt=[];
    ywt=[];
    gev=[];
    return
  end  

  if ~fixdiv
    div=1.01; 
  end

  nunstab=0;
  zxz=0;

  for i=1:n
    %Div calc
    if ~fixdiv
      if abs(a(i,i)) > 0
        divhat=abs(b(i,i))/abs(a(i,i));
        if 1+realsmall<divhat && divhat<div
          div=.5*(1+divhat);
        end
      end
    end
    %End of div calc

    nunstab=nunstab+(abs(b(i,i))>div*abs(a(i,i)));
    if abs(a(i,i))<realsmall && abs(b(i,i))<realsmall
      zxz=1;
    end
  end

  div;
  nunstab;
  if ~zxz
    [a b q z]=qzdiv(div,a,b,q,z);
  end

  gev=[diag(a) diag(b)];

  if zxz
    disp('Coincident zeros.  Indeterminacy and/or nonexistence.')
    eu=[-2;-2];
    return
  end

  q1=q(1:n-nunstab,:);
  q2=q(n-nunstab+1:n,:);
  z1=z(:,1:n-nunstab)';
  z2=z(:,n-nunstab+1:n)';
  a2=a(n-nunstab+1:n,n-nunstab+1:n);
  b2=b(n-nunstab+1:n,n-nunstab+1:n);
  etawt=q2*pi;
  zwt=q2*psi;
  [ueta,deta,veta]=svd(etawt);
  md=min(size(deta));
  bigev=find(diag(deta(1:md,1:md))>realsmall);
  ueta=ueta(:,bigev);
  veta=veta(:,bigev);
  deta=deta(bigev,bigev);
  [uz,dz,vz]=svd(zwt);
  md=min(size(dz));
  bigev=find(diag(dz(1:md,1:md))>realsmall);
  uz=uz(:,bigev);
  vz=vz(:,bigev);
  dz=dz(bigev,bigev);
  if isempty(bigev)
    exist=1;
  else
    exist=norm(uz-ueta*ueta'*uz) < realsmall*n;
  end
  if ~isempty(bigev)
    zwtx0=b2\zwt;

    zwtx=zwtx0;
    M=b2\a2;
    for i=2:nunstab
      zwtx=[M*zwtx zwtx0];
    end
    zwtx=b2*zwtx;
    
    if any(any(isinf(zwtx))) | any(any(isnan(zwtx)));
      disp('inf or nan in zwtx')
      eu=[-2;-2];
      G1=[];
      C=[];
      impact=[];
      fmat=[];
      fwt=[];
      ywt=[];
      gev=[];
      return          
    end  

    [ux,dx,vx]=svd(zwtx);
    md=min(size(dx));
    bigev=find(diag(dx(1:md,1:md))>realsmall);
    ux=ux(:,bigev);
    vx=vx(:,bigev);
    dx=dx(bigev,bigev);
    existx=norm(ux-ueta*ueta'*ux) < realsmall*n;
  else
    existx=1;
  end
  %----------------------------------------------------
  % Note that existence and uniqueness are not just matters of comparing
  % numbers of roots and numbers of endogenous errors.  These counts are
  % reported below because usually they point to the source of the problem.
  %------------------------------------------------------


  [ueta1,deta1,veta1]=svd(q1*pi);
  md=min(size(deta1));
  bigev=find(diag(deta1(1:md,1:md))>realsmall);
  ueta1=ueta1(:,bigev);
  veta1=veta1(:,bigev);
  deta1=deta1(bigev,bigev);

  if existx | nunstab==0
    %Solution exists
    eu(1)=1;
  else
    %Does solution exist for unforecastable z only? Else, there is
    %no solution due to unstable roots and/or endogenous errors.
    if exist
      eu(1)=-1;
    end
  end

  if isempty(veta1)
    unique=1;
  else
    unique=norm(veta1-veta*veta'*veta1)<realsmall*n;
  end

  %Unique solution?
  if unique
    eu(2)=1;
  else
    fprintf(1,'Indeterminacy.  %d loose endog errors.\n',size(veta1,2)-size(veta,2));
  end

  tmat = [eye(n-nunstab) -(ueta*(deta\veta')*veta1*deta1*ueta1')'];
  G0= [tmat*a; zeros(nunstab,n-nunstab) eye(nunstab)];
  G1= [tmat*b; zeros(nunstab,n)];
  %----------------------
  % G0 is always non-singular because by construction there are no zeros on
  % the diagonal of a(1:n-nunstab,1:n-nunstab), which forms G0's ul corner.
  %-----------------------
  G0I=inv(G0);
  G1=G0I*G1;
  usix=n-nunstab+1:n;
  C=G0I*[tmat*q*c;(a(usix,usix)-b(usix,usix))\q2*c];
  impact=G0I*[tmat*q*psi;zeros(nunstab,size(psi,2))];
  fmat=b(usix,usix)\a(usix,usix);
  fwt=-b(usix,usix)\q2*psi;
  ywt=G0I(:,usix);
  %-------------------- above are output for system in terms of z'y -------
  G1=real(z*G1*z');
  C=real(z*C);
  impact=real(z*impact);
  ywt=z*ywt;
