function x1 = fcnSolveAdjU(A,P,qs,d,lambda)

% DEFINE SYSTEM AND RHS
B = P'*P + A'*A ;
RHS = P'*d + lambda^2*A'*qs ;

% ADD JACOBI PRECONDITIONNER
N = size(B,1) ;
Pc = spdiags(1./diag(B),0,N,N) ;
B = Pc*B ;
RHS = Pc*RHS ;

% INITIALIZE PCG
x0 = zeros(N,1) ;
r0 = RHS - B*x0 ;
p0 = r0 ;
count  = 0 ;

while count < length(x0) ;
    
    count = count + 1 ;
    alpha = (r0'*r0)/(p0'*(B*p0)) ;
       x1 = x0 + alpha*p0 ;
       r1 = r0 - alpha*(B*p0) ;
       
       if norm(r1) < 1e-12 ;
           return
       end
       
       M  = (r1'*r1)/(r0'*r0) ;
       p1 = r1 + M*p0 ;
       
       % Update
       x0 = x1 ;
       r0 = r1 ;
       p0 = p1 ;
       
end











