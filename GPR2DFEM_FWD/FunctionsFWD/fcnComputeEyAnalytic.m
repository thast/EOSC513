function Ey = fcnComputeEyAnalytic(surveyObj,meshObj,TxID,f)

% NODE LOCATIONS
X0 = meshObj.r0(1) ;
Z0 = meshObj.r0(2) ;
Xn = X0 + cumsum([0; meshObj.hx]) ;
Zn = Z0 + cumsum([0; meshObj.hz]) ;
[Xn,Zn] = ndgrid(Xn,Zn) ;
Xn = Xn(:) ; Zn = Zn(:) ;

% ANALYTIC E FIELD FROM SOURCE
x0 = surveyObj.Tx{TxID}(1) ;
z0 = surveyObj.Tx{TxID}(2) ;
 m = surveyObj.Tx{TxID}(3) ;
 w = 2*pi*f ;
 R = sqrt((Xn-x0).^2 + (Zn-z0).^2) ;

mu0 = 4*pi*1e-7 ;
ep0 = 8.8542e-12 ;
% tao = 0.05*min([meshObj.hx; meshObj.hz]) ;

Ey = (m./(1i*w*ep0*4*pi))*((R).^-3).*exp(-1i*w*sqrt(mu0*ep0)*R).*(w^2*mu0*ep0*R.^2 - 1i*w*sqrt(mu0*ep0)*R - 1) ;
% Ey = (isnan(Ey)==1) == 1 ;







