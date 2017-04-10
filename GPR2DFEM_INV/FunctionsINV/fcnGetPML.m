function [PML,diagPML] = fcnGetPML(meshObj,f)

% INPUTS
% meshObj: mesh object
%       f: frequency
%
% OUTPUTS:
% [ex,ez]: corresponding values at cell centers for envoking PML

nx = meshObj.nc(1) ;
nz = meshObj.nc(2) ;
a0 = 1.76 ;
lambda = 3e8/f ;          % maximum physicall possible wavelength

% PML Thickness
Xn = meshObj.r0(1) + cumsum([0; meshObj.hx]) ;          % Node vector 1D
kx = find(lambda > Xn-meshObj.r0(1),1,'last') ;         % #PML
Lx = Xn(kx+1) - meshObj.r0(1) ;                         % Length PML
xLB = Xn(1) + Lx ;
xUB = Xn(end) - Lx ;

Zn = meshObj.r0(2) + cumsum([0; meshObj.hz]) ;
kz =  find(lambda > Zn-meshObj.r0(2),1,'last') ; 
Lz = Zn(kz+1) - meshObj.r0(2) ;
zLB = Zn(1) + Lz ;
zUB = Zn(end) - Lz ;

% PML ID
xID = zeros(nx,1) ; xID([1:kx end+1-kx:end]) = 1 ;
zID = zeros(nz,1) ; zID([1:kz end+1-kz:end]) = 1 ;
[xID, zID] = ndgrid(xID, zID) ;

Xc = meshObj.r0(1) + cumsum([meshObj.hx]) - meshObj.hx/2 ;
Zc = meshObj.r0(2) + cumsum([meshObj.hz]) - meshObj.hz/2 ;
[Xc, Zc]   = ndgrid(Xc, Zc) ;
Xc = Xc(:) ; Zc = Zc(:) ;

% PML VALUES
ex = ones(nx*nz,1) ;
lx = bsxfun(@min,abs(Xc(xID==1)-xLB),abs(Xc(xID==1)-xUB)) ;
ex(xID==1) = 1 - 1i*a0*(lx/Lx).^2 ;

ez = ones(nx*nz,1) ;
lz = bsxfun(@min,abs(Zc(zID==1)-zLB),abs(Zc(zID==1)-zUB)) ;
ez(zID==1) = 1 - 1i*a0*(lz/Lx).^2 ;

% FINALIZE
PML = [ex,ez] ;
diagPML = spdiags(ex.*ez,0,length(ex),length(ez)) ;




