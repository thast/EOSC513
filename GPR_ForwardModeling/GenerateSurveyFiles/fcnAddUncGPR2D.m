function fcnAddUncGPR2D(workDir,dpreFile,Pct,Floor)

% LOAD DATA

dpreFile = strcat(workDir,'\',dpreFile) ;

FID = fopen(dpreFile,'r') ;

A = textscan(FID,'%s','delimiter','\n') ;

fclose(FID) ;

A = A{1} ;
A = A(cellfun(@isempty,A)==0) ;

% Ignore flag and number of transmitters
temp1 = strsplit(A{1},' ') ;
temp2 = strsplit(A{2},' ') ;
assert(strcmp(temp1{1},'IGNORE')==1,'Ignore flag not found') ;
assert(strcmp(temp2{1},'N_TRX')==1,'Number of transmitters missing') ;
numFlag = str2num(temp1{2}) ;
    nTx = str2num(temp2{2}) ;
    

% FIND LOCATION OF TRANSMITTER TYPES
TxID = cellfun(@(x) regexp(x,'TRX_EY_DIPOLE'),A,'UniformOutput',false) ;
TxID = cellfun(@(x) isempty(x),TxID) ;
TxID = find(TxID == 0) ;
assert(isempty(TxID)==0,'No transmitters found') ;
TxID = [TxID; size(A,1)+1] ;

% Pre-define data array
data = [] ;


for ii = 1:numel(TxID)-1 ;
    
    type{ii} = 'TRX_EY_DIPOLE' ;
    k = TxID(ii) ;

    p = sscanf(A{TxID(ii)+1},'%g %g %g') ;
    Tx{ii} = p' ;

    % Organize receivers
    temp1 = strsplit(A{k+2},' ') ;
    temp2 = strsplit(A{k+3},' ') ;
    assert(strcmp(temp1{1},'N_RECV')==1,'Number of receivers missing or incorrectly formatted') ;
    assert(strcmp(temp2{1},'N_FREQ')==1,'Number of tme channels missing or incorrectly formatted') ;
    n_rec(ii) = str2num(temp1{2}) ;
    n_freq(ii) = str2num(temp2{2}) ;

    temprx = A(k+4:k+3+n_rec(ii)*n_freq(ii)) ;                          % Finds the right rows of A
    temprx = cellfun(@(x) strsplit(x,'\t'),temprx,'UniformOutput',false) ;  % Splits and converts to cell array
    temprx = cellfun(@(x) str2double(x),temprx,'UniformOutput',false) ;     % Turns strings into double
    temprx = [ii*ones(n_rec(ii)*n_freq(ii),1) cell2mat(temprx)] ;             % Turns cell array into matrix and assigns transmitter ID

    data = [data; temprx] ;

end



% ADD UNCERTAINTIES

data = [data Pct*abs(data(:,[5 6]))+Floor] ;

% PRINT OBS FILE

dobsFile = strcat(workDir,'/dobsGPR2D.txt') ;

% Open File
fid = fopen(dobsFile,'w') ;

fprintf(fid,'%s %s\n%s %i\n','IGNORE',numFlag,'N_TRX',nTx) ;

for pp = 1:numel(Tx) ;
    
    fprintf(fid,'\n\n%s\n','TRX_EY_DIPOLE') ;
    fprintf(fid,' %g %g %g\n',Tx{pp}) ;
    
    dataTemp = data(data(:,1)==pp,2:end) ;
    COUNT = 0 ;

    fprintf(fid,'\n%s %i\n%s %i','N_RECV',n_rec(pp),'N_FREQ',n_freq(pp)) ;
    
    for rr = 1:n_rec(pp) ;
        for ff = 1:n_freq(pp) ;
            COUNT = COUNT + 1 ;
            fprintf(fid,'\n%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e',dataTemp(COUNT,:)) ;
        end
    end

end
    
fclose(fid) ;

