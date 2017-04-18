function [A,qs] = fcnAddDirichlet(A,qs,meshObj)

nx = meshObj.nc(1) ;
nz = meshObj.nc(2) ;

% Which nodes are on outside
sx = zeros(nx+1,1) ;
sz = zeros(nz+1,1) ;
sx([1 end]) = 1 ;
sz([1 end]) = 1 ;
[sx,sz] = ndgrid(sx,sz) ;
sx = sx(:) ; sz = sz(:) ;
s = bsxfun(@max,sx,sz) ;

% Replace rows in ficticious source and matrix
qs(s==1) = 0 ;
S = speye(length(s)) ;
A(s==1,:) = S(s==1,:) ;









