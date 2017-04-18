function opObj = fcnGetDiffOperators(meshObj)

nx = meshObj.nc(1) ;
nz = meshObj.nc(2) ;
hx = meshObj.hx ;
hz = meshObj.hz ;

% NODE TO EDGE GRADIENT OPERATOR
Gx = spdiags(1./hx,0,nx,nx)*spdiags(ones(nx+1,1)*[-1 1],[0,1],nx,nx+1) ;
Gx = kron(speye(nz+1),Gx) ;
Gz = spdiags(1./hz,0,nz,nz)*spdiags(ones(nz+1,1)*[-1 1],[0,1],nz,nz+1) ;
Gz = kron(Gz,speye(nx+1)) ;

G = [Gx; Gz] ;

% DIVERGENCE OPERATOR
D = -G' ;       % It needed the minus for some reason. Likely since right-handed coordinates.

% NODE TO CENTER AVERAGE
Anc  = kron(spdiags(ones(nz+1,1)*[1/2 1/2],[0,1],nz,nz+1),spdiags(ones(nx+1,1)*[1/2 1/2],[0,1],nx,nx+1)) ;
Aecx = kron(spdiags(ones(nz+1,1)*[1/2 1/2],[0,1],nz,nz+1),speye(nx)) ;
Aecz = kron(speye(nz),spdiags(ones(nx+1,1)*[1/2 1/2],[0,1],nx,nx+1)) ;
Aec  = [Aecx Aecz] ;
 
% COLLECT OPERATORS
opObj.Grad = G ;
opObj.Div  = D ;
opObj.Anc  = Anc ;
opObj.Aecx = Aecx ;
opObj.Aecz = Aecz ;
opObj.Aec  = Aec ;


% TESTING OPERATORS
% xn = meshObj.r0(1) + cumsum([0; meshObj.hx]) ;
% zn = meshObj.r0(2) + cumsum([0; meshObj.hz]) ;
% [Xn, Zn] = ndgrid(xn, zn) ;
% Xn = Xn(:) ; Zn = Zn(:) ;
% 
% xc = meshObj.r0(1) + cumsum(meshObj.hx) - meshObj.hx/2 ;
% zc = meshObj.r0(2) + cumsum(meshObj.hz) - meshObj.hz/2 ;
% [Xex, Zex] = ndgrid(xc, zn) ;
% Xex = Xex(:) ; Zex = Zex(:) ;
% [Xez, Zez] = ndgrid(xn, zc) ;
% Xez = Xez(:) ; Zez = Zez(:) ;
% 
% F   = exp(-(Xn/3).^2 - (Zn/5).^2) ;
% dFx = -(2*Xex/9).*exp(-(Xex/3).^2 - (Zex/5).^2) ;
% dFz = -(2*Zez/25).*exp(-(Xez/3).^2 - (Zez/5).^2) ;
% ddF = -(2/9 + 2/25)*exp(-(Xn/3).^2 - (Zn/5).^2) + ((2*Xn/9).^2).*exp(-(Xn/3).^2 - (Zn/5).^2) + ((2*Zn/25).^2).*exp(-(Xn/3).^2 - (Zn/5).^2) ;
% 
% dQ  = (opObj.Grad)*F ;
% ddQ = ((opObj.Grad)'*opObj.Grad)*F ;
% 
% subplot(1,3,1)
% imagesc(reshape(dFx,602,603)) ; colorbar ;
% subplot(1,3,2)
% imagesc(reshape(dQ(1:end/2),602,603)) ; colorbar ;
% subplot(1,3,3)
% imagesc(reshape(dFx-dQ(1:end/2),602,603)) ; colorbar ;
% 
% subplot(1,3,1)
% imagesc(reshape(dFz,602,603)) ; colorbar ;
% subplot(1,3,2)
% imagesc(reshape(dQ(end/2+1:end),602,603)) ; colorbar ;
% subplot(1,3,3)
% imagesc(reshape(dFz-dQ(end/2+1:end),602,603)) ; colorbar ;
% 
% subplot(1,3,1)
% imagesc(reshape(ddF,603,603)) ; colorbar ;
% subplot(1,3,2)
% imagesc(reshape(ddQ,603,603)) ; colorbar ;
% subplot(1,3,3)
% imagesc(reshape(ddF-ddQ,603,603)) ; colorbar ;





