function u = fcnSolveE(surveyObj,meshObj,opObj,modelObj,qs,PML)

for pp = 1:numel(surveyObj.Tx) ;
    
    f = unique(surveyObj.data(surveyObj.data(:,1)==pp,4)) ;
    
    for qq = 1:numel(f) ;
        
        A = fcnComputeA(meshObj,opObj,modelObj,f(qq),PML) ;
        u{pp,qq} = A\qs{pp} ;
        
    end

end










