function [gradF,H] = fcnComputeGradHessian(surveyObj,meshObj,opObj,m,P,qs,lambda) 

% INITIALIZE SELECTED PROPERTIES
nx = meshObj.nc(1) ;
nz = meshObj.nc(2) ;
Anc = opObj.Anc ;

% INITIALIZE GRADIENT AND HESSIAN
gradF = zeros(numel(m),1) ;
    H = sparse(numel(m),numel(m)) ;

f = unique(surveyObj.data(:,4)) ;    

fprintf('%s','UPDATING GRADIENT AND HESSIAN')

% LOOP OVER FREQUENCIES
for pp = 1:numel(f) ;
    
    fprintf('\n%s %i %s','AT FREQUENCY',pp,'FOR TRANSMITTER')
    
    w = 2*pi*f(pp) ;

    % Get PML
    [PML,diagPML] = fcnGetPML(meshObj,f(pp)) ;
    % Compute A Matrix and PML Laplacian
    [A,L] = fcnComputeA(meshObj,opObj,m,f(pp),PML) ;
    
    for qq = 1:numel(surveyObj.Tx) ;
        
        fprintf('% i',qq) ;
        
        % Get data
        d = surveyObj.data(surveyObj.data(:,1)==qq & surveyObj.data(:,4)==f(pp),[5,6])*[1 1]' ; 
        % Solve for "adjusted field" (u_bar in Eq. 4)
        u = (P{qq}'*P{qq} + lambda.^2*(A'*A))\(P{qq}'*d + lambda^2*A'*qs{qq}) ;
        diagU = spdiags(u,0,(nx+1)*(nz+1)) ;
        % Add to Gradient and Hessian
        gradF = gradF + lambda^2*real(w^2*(diagU*(Anc'*diagPML))'*(L*u + w^2*(diagU*(Anc'*diagPML))*m - qs{qq})) ;
        H = H + lambda^2*w^4*real((diagU*(Anc'*diagPML))'*(diagU*(Anc'*diagPML))) ;
        
    end
    
    
end

fprintf('\n%s\n','GRADIENT AND HESSIAN UPDATED') 










