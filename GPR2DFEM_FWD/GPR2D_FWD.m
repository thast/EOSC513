function GPR2D_FWD(workDir)

% PARAMETERS FOR FORWARD MODEL
paramFile = strcat(workDir,'\parametersGPR2D_FWD') ;

paramFileAssert = strcat(workDir,'\parametersGPR2D_FWD.m') ;
assert(exist(paramFileAssert)==2,'Parameters file does not exist or has incorrect name') ;
run(paramFile);

% Add full path to files
codeDir   = strcat(codeDir,'\FunctionsFWD') ;
TxRxFile  = strcat(workDir,'\',TxRxFile) ;
meshFile  = strcat(workDir,'\',meshFile) ;
modelFile = strcat(workDir,'\',modelFile) ;

addpath(codeDir)

%% LOAD FILES

fprintf('\n%s\n','LOADING FILES')

% Load Mesh File
meshObj = loadmeshGPR2D(meshFile) ;

% Load Model File
modelObj = loadmodelGPR2D(modelFile) ;

% Load Survey File
surveyObj = loadsurveyGPR2D(TxRxFile) ;

fprintf('\n%s\n','FILES LOADED')

%% PREDICT DATA

surveyObj.data = [surveyObj.data zeros(size(surveyObj.data,1),2)] ;

% Get operators
opObj = fcnGetDiffOperators(meshObj) ;

% Compute Each Source

% FOR EACH TRANSMITTER
for ii = 1:numel(surveyObj.Tx) ;
    
    TxID = ii ;
    
    % Generate Data Projection Matrix
    P = fcnGeneralDataProjection(surveyObj,meshObj,ii) ;
    
    % SOURCE (NEAREST NODE)
    qs = fcnComputeSourceNearest(meshObj,surveyObj,ii) ;
    
    % FOR EACH FREQUENCY
    for jj = 1:numel(unique(surveyObj.data(:,4))) ;
        
        % Transmitter ID and frequency
        f = unique(surveyObj.data(:,4)) ;
        f = f(jj) ;
        
        % Compute A Matrix for Model
         A = fcnComputeA(meshObj,opObj,modelObj,f) ;
        % Add Dirichlet Boundary Condition
%          [A,qs] = fcnAddDirichlet(A,qs,meshObj) ;
        % Solve System with Preconditionned Conjugate Gradient
        fprintf('\n%s\n','SOLVE SYSTEM')
        fprintf('%s %i\n','TRANSMITTER:',ii)
        fprintf('%s %.3e %s\n','FREQUENCY:',f,'HZ')
         E = fcnSolveWaveEqPCG(A,qs) ;

        % Predict Data
        surveyObj.data(surveyObj.data(:,1)==ii & surveyObj.data(:,4)==f,[5 6]) = P*[real(E) imag(E)] ;
        
        % Plot Data
        fcnPlotFieldGPR2D(E,meshObj) ;   
        
    end
end

% WRITE OUTPUT FILE

dataFilename = strcat(workDir,'\dpreGPR2DFWD.txt') ;
writeGPR2Ddata(surveyObj,dataFilename) ;

fprintf('\n%s\n','DATA PREDICTED')

rmpath(codeDir)








