function x1 = fcnSolveWaveEqPCG(A,qs)

% Compute Preconditionner and Apply to Both Sides (Jacobi Preconditionner)
N = size(A,1) ;
Pc = spdiags(1./diag(A),0,N,N) ;

 A = Pc*A ;
qs = Pc*qs ;

x1 = A\qs ;





% % Initiate for Conjugate Gradient
% x0 = zeros(N,1) ;
% r0 = qs - A*x0 ;
% p0 = r0 ;
% count  = 0 ;
% 
% while count < length(x0) ;
%     
%     count = count + 1 ;
%     alpha = (r0'*r0)/(p0'*(A*p0)) ;
%        x1 = x0 + alpha*p0 ;
%        r1 = r0 - alpha*(A*p0) ;
%        
%        if norm(r1) < 1e-12 ;
%            return
%        end
%        
%        B  = (r1'*r1)/(r0'*r0) ;
%        p1 = r1 + B*p0 ;
%        
%        % Update
%        x0 = x1 ;
%        r0 = r1 ;
%        p0 = p1 ;
%        
% end



