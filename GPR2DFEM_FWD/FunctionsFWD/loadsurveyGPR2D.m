function surveyObj = loadsurveyGPR2D(TxRxFile)

FID = fopen(TxRxFile,'r') ;

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


surveyObj.flag = numFlag ;
surveyObj.Tx = Tx ;
surveyObj.type = type ;
surveyObj.data = data ;
surveyObj.n_rec = n_rec ;
surveyObj.n_freq = n_freq ;











