function fcnCreateTxRxFileGPR2D(workDir)

% This function take the parameters defined in parametersTxRxLines.m and
% create the TxRxFile for transmitters of type TxRxLines
% 
% INPUTS
% workDir : Working Directory as a string
%
% PARAMETERS REQUIRED FOR FILE parametersTxRxVRM
%       Tx : Cell array where each Tx{i} is an 1 X 3 array with [x z m]
%       Rx : Cell array where each Rx{i} is an N X 2 array with unique
%            [xj zj] observation locations.
%        f : Frequency vector for all observed frequencies (should be in Hz)
%     Flag : String for IGNORE data (I like NaN)
% 
% OUTPUTS
% TxRxFileGPR2D.txt

%% WRITE THE FILE

% Obtain parameters from file
% parametersTxRxLines

paramFile = strcat(workDir,'\parametersTxRxGPR2D') ;

paramFileAssert = strcat(workDir,'\parametersTxRxGPR2D.m') ;
assert(exist(paramFileAssert)==2,'Parameters file does not exist or has incorrect name') ;
run(paramFile)

Filename = strcat(workDir,'\TxRxFileGPR2D.txt') ;

% Open File
fid = fopen(Filename,'w') ;

fprintf(fid,'%s %s\n%s %i\n','IGNORE',Flag,'N_TRX',numel(Tx)) ;

for pp = 1:numel(Tx) ;
    
    fprintf(fid,'\n\n%s\n','TRX_EY_DIPOLE') ;
    fprintf(fid,' %g %g %g\n',Tx{pp}) ;

    fprintf(fid,'\n%s %i\n%s %i','N_RECV',size(Rx{pp},1),'N_FREQ',numel(f)) ;

    for rr = 1:size(Rx{pp},1) ;

        for ff = 1:numel(f) ;
            fprintf(fid,'\n%.8e\t%.8e\t%.8e',Rx{pp}(rr,1),Rx{pp}(rr,2),f(ff)) ;
        end

    end

end
    
fclose(fid) ;



















