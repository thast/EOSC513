function A = fcnComputeA(meshObj,opObj,m,f)

% NODE LOCATIONS
nx = meshObj.nc(1) ;
nz = meshObj.nc(2) ;

% PHYSICAL PROPERTIES
if isempty(m) == 1 ;
    S = (1/3e8)^2*ones(nx*nz,1) ;
else
    S  = (m.value).^-2 ;
end

if isempty(m) == 0 ;
    % GET PML PARAMETERS
    a0 = 1.76 ;
    lambda = 3e8/f ;          % 2 wavelength of PML

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
    % pmlID = bsxfun(@max,xID,zID) ; clear xID zID ;          % ID for in or out of PML
    % pmlID = pmlID(:) ;

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

    % FORM OPERATOR
    Grad = opObj.Grad ;
    Div  = opObj.Div ;
    Anc  = opObj.Anc ;
    Aecx = opObj.Aecx ;
    Aecz = opObj.Aecz ;
    
    w = 2*pi*f ;
    B = spdiags(Anc'*(S.*ex.*ez),0,(nx+1)*(nz+1),(nx+1)*(nz+1)) ;
    M = spdiags([Aecx'*(ez./ex); Aecz'*(ex./ez)],0,nx*(nz+1)+nz*(nx+1),nx*(nz+1)+nz*(nx+1)) ;
    
    A  = Div*(M*Grad) + w^2*B ;

else
    
    Grad = opObj.Grad ;
    Div  = opObj.Div ;
    Anc  = opObj.Anc ;
    
    B = spdiags(Anc'*S,0,(nx+1)*(nz+1),(nx+1)*(nz+1)) ;
    w = 2*pi*f ;
    A = (Div*Grad) + w^2*B ;
    
end










