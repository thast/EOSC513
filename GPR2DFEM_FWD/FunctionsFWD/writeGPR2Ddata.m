function writeGPR2Ddata(surveyObj,dataFilename)

% Open File
fid = fopen(dataFilename,'w') ;

fprintf(fid,'%s %s\n%s %i\n','IGNORE',surveyObj.flag,'N_TRX',numel(surveyObj.Tx)) ;

for pp = 1:numel(surveyObj.Tx) ;
    
    fprintf(fid,'\n\n%s\n','TRX_EY_DIPOLE') ;
    fprintf(fid,' %g %g %g\n',surveyObj.Tx{pp}) ;
    
    dataTemp = surveyObj.data(surveyObj.data(:,1)==pp,2:end) ;
    COUNT = 0 ;

    fprintf(fid,'\n%s %i\n%s %i','N_RECV',surveyObj.n_rec(pp),'N_FREQ',surveyObj.n_freq(pp)) ;
    
    for rr = 1:surveyObj.n_rec(pp) ;
        for ff = 1:surveyObj.n_freq(pp) ;
            COUNT = COUNT + 1 ;
            fprintf(fid,'\n%.8e\t%.8e\t%.8e\t%.8e\t%.8e',dataTemp(COUNT,:)) ;
        end
    end

end
    
fclose(fid) ;


















