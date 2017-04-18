function qs = fcnComputeSourceNearest(meshObj,surveyObj,ii)

% Transmitter Properties
x0 = surveyObj.Tx{ii}(1) ;
z0 = surveyObj.Tx{ii}(2) ;
 m = surveyObj.Tx{ii}(3) ;

% Mesh Properties
Xn = meshObj.r0(1) + cumsum([0; meshObj.hx]) ;
Zn = meshObj.r0(2) + cumsum([0; meshObj.hz]) ;
[Xn, Zn] = ndgrid(Xn, Zn) ;
Xn = Xn(:) ; Zn = Zn(:) ;

% Move Source
qs = zeros(length(Xn),1) ;
 R = (Xn-x0).^2 + (Zn-z0).^2 ;
 k = find(R == min(R),1,'first') ;
qs(k) = m ;