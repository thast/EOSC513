function writeVRMpred(surveyObj,workDir)

Filename = strcat(workDir,'\VRMdata_pred.txt') ;

FID = fopen(Filename,'wt') ;
fclose(FID) ;

IgnFlag = surveyObj.flag ;
TxIDvec = surveyObj.data(:,1) ;

% ORDER DOES NOT NEED TO BE REVERSED SINCE I THINK KRIS FIXED THIS
for jj = 1:numel(surveyObj.Tx) ;
    
    % Obtain Number of Receivers and Times
        t = unique(surveyObj.data(TxIDvec==jj,5)) ;
     nRec = numel(surveyObj.data(TxIDvec==jj,1))/numel(t) ;
    nTime = numel(t) ;

    % Write File    
    FID = fopen(Filename,'a') ;
    fprintf(FID,'\n') ;
    fclose(FID) ;
    
    Datajj = [surveyObj.data(TxIDvec==jj,2:5) IgnFlag*ones(nRec*nTime,3) surveyObj.data(TxIDvec==jj,6:end)] ;
    
    dlmwrite(Filename,Datajj,'delimiter','\t','precision','%.8e','-append') ;
    
end




