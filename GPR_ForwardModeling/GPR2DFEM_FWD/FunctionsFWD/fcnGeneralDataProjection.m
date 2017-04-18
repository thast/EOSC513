function P = fcnGeneralDataProjection(surveyObj,meshObj,ii)

rec_loc = unique(surveyObj.data(surveyObj.data(:,1)==ii,2:3),'rows','stable') ;
x0  = rec_loc(:,1) ;
z0  = rec_loc(:,2) ;
ncx = meshObj.nc(1) ;
ncz = meshObj.nc(2) ;

P = sparse(size(rec_loc,1),(ncx+1)*(ncz+1)) ;

% NODES

Xn = meshObj.r0(1) + [0; cumsum(meshObj.hx)] ;
Zn = meshObj.r0(2) + [0; cumsum(meshObj.hz)] ;

for ii = 1:size(rec_loc,1) ;
    
    j = find(Xn<rec_loc(ii,1)+1e-7,1,'last') ;
    k = find(Zn<rec_loc(ii,2)+1e-7,1,'last') ;
    
    % Data Vector Locations
    p1 = (k-1)*(ncx+1) + j ;
    p2 = (k-1)*(ncx+1) + j + 1 ;
    p3 = k*(ncx+1) + j ;
    p4 = k*(ncx+1) + j + 1 ;
    
    A = 1/(meshObj.hx(j)*meshObj.hz(k)) ;
    
    P(ii,[p1 p2 p3 p4]) = A*[(Xn(j+1)-x0(ii))*(Zn(k+1)-z0(ii)) (x0(ii)-Xn(j))*(Zn(k+1)-z0(ii)) (Xn(j+1)-x0(ii))*(z0(ii)-Zn(k)) (x0(ii)-Xn(j))*(z0(ii)-Zn(k))] ;
    
end

