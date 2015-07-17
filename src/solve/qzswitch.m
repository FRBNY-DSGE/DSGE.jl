function [A,B,Q,Z] = qzswitch(i,A,B,Q,Z)
%function [A,B,Q,Z] = qzswitch(i,A,B,Q,Z)
%
% Takes U.T. matrices A, B, orthonormal matrices Q,Z, interchanges
% diagonal elements i and i+1 of both A and B, while maintaining
% Q'AZ' and Q'BZ' unchanged.  If diagonal elements of A and B
% are zero at matching positions, the returned A will have zeros at both
% positions on the diagonal.  This is natural behavior if this routine is used
% to drive all zeros on the diagonal of A to the lower right, but in this case
% the qz transformation is not unique and it is not possible simply to switch
% the positions of the diagonal elements of both A and B.
 realsmall=sqrt(eps)*10;
%realsmall=1e-3;
a = A(i,i); d = B(i,i); b = A(i,i+1); e = B(i,i+1);
c = A(i+1,i+1); f = B(i+1,i+1);
		% A(i:i+1,i:i+1)=[a b; 0 c];
		% B(i:i+1,i:i+1)=[d e; 0 f];
if (abs(c)<realsmall & abs(f)<realsmall)
	if abs(a)<realsmall
		% l.r. coincident 0's with u.l. of A=0; do nothing
		return
	else
		% l.r. coincident zeros; put 0 in u.l. of a
		wz=[b; -a];
		wz=wz/sqrt(wz'*wz);
		wz=[wz [wz(2)';-wz(1)'] ];
		xy=eye(2);
	end
elseif (abs(a)<realsmall & abs(d)<realsmall)
	if abs(c)<realsmall
		% u.l. coincident zeros with l.r. of A=0; do nothing
		return
	else
		% u.l. coincident zeros; put 0 in l.r. of A
		wz=eye(2);
		xy=[c -b];
		xy=xy/sqrt(xy*xy');
		xy=[[xy(2)' -xy(1)'];xy];
	end
else
	% usual case
	wz = [c*e-f*b, (c*d-f*a)'];
	xy = [(b*d-e*a)', (c*d-f*a)'];
	n = sqrt(wz*wz');
	m = sqrt(xy*xy');
	if m<eps*100
		% all elements of A and B proportional
		return
	end
   wz = n\wz;
   xy = m\xy;
   wz = [wz; -wz(2)', wz(1)'];
   xy = [xy;-xy(2)', xy(1)'];
end
A(i:i+1,:) = xy*A(i:i+1,:);
B(i:i+1,:) = xy*B(i:i+1,:);
A(:,i:i+1) = A(:,i:i+1)*wz;
B(:,i:i+1) = B(:,i:i+1)*wz;
Z(:,i:i+1) = Z(:,i:i+1)*wz;
Q(i:i+1,:) = xy*Q(i:i+1,:);