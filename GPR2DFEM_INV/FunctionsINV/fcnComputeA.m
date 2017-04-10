function [A,L] = fcnComputeA(meshObj,opObj,m,f,PML)

ex = PML(:,1) ;
ez = PML(:,2) ;
nx = meshObj.nc(1) ;
nz = meshObj.nc(2) ;

% PHYSICAL PROPERTIES
if isempty(m) == 1 ;
    S = (1/3e8)^2*ones(nx*nz,1) ;
else
    S  = m.^-2 ;
end

% FORM OPERATORS
Grad = opObj.Grad ;
Div  = opObj.Div ;
Anc  = opObj.Anc ;
Aecx = opObj.Aecx ;
Aecz = opObj.Aecz ;

w = 2*pi*f ;
B = spdiags(Anc'*(S.*ex.*ez),0,(nx+1)*(nz+1),(nx+1)*(nz+1)) ;
M = spdiags([Aecx'*(ez./ex); Aecz'*(ex./ez)],0,nx*(nz+1)+nz*(nx+1),nx*(nz+1)+nz*(nx+1)) ;
L = Div*(M*Grad) ;

A  = L + w^2*B ;












