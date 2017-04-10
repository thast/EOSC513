function GPR2D_INV(workDir)

% PARAMETERS FOR FORWARD MODEL
paramFile = strcat(workDir,'\parametersGPR2D_INV') ;

paramFileAssert = strcat(workDir,'\parametersGPR2D_INV.m') ;
assert(exist(paramFileAssert)==2,'Parameters file does not exist or has incorrect name') ;
run(paramFile);

% Add full path to files
codeDir       = strcat(codeDir,'\FunctionsINV') ;
dobsFile      = strcat(workDir,'\',dobsFile) ;
meshFile      = strcat(workDir,'\',meshFile) ;
modelFile     = strcat(workDir,'\',modelFile) ;

addpath(codeDir)

%% LOAD FILES

fprintf('\n%s\n','LOADING FILES')

% Load Mesh File
meshObj = loadmeshGPR2D(meshFile) ;

% Load Model File (Starting model)
modelObj = loadmodelGPR2D(modelFile) ;

% Load Survey File
surveyObj = loadsurveyGPR2D(dobsFile) ;

fprintf('\n%s\n','FILES LOADED')

%% OPERATORS, PROJECTION MATRICIES AND SOURCES

fprintf('\n%s\n','COMPUTE DIFF OPERATORS, PROJECTION MATRICIES AND SOURCES')

% Get operators
opObj = fcnGetDiffOperators(meshObj) ;

% Projection Matricies and Sources
for ii = 1:numel(surveyObj.Tx) ;
    
    % Generate Data Projection Matrix for Tx{ii}
    P{ii} = fcnGeneralDataProjection(surveyObj,meshObj,ii) ;
    
    % SOURCE (NEAREST NODE)
    qs{ii} = fcnComputeSourceNearest(meshObj,surveyObj,ii) ;
    
end

%% OPTIMIZATION

% INITIALIZE
m = modelObj.value.^(-2) ;
lambda = 1 ;



% while

[gradF,H] = fcnComputeGradHessian(surveyObj,meshObj,opObj,m,P,qs,lambda) ;


% PML = fcnGetPML(meshObj,f) ;
% A = fcnComputeA(meshObj,opObj,modelObj,f,PML) ;
% u = fcnSolveE(surveyObj,meshObj,opObj,modelObj,qs,PML) ;
% gradF = fcnComputeGradient(surveyObj,meshObj,opObj,PML,u,qs,m,lambda) ;











%% >>>>>>>>>>>>0LD<<<<<<<<<<<<<<<

surveyObj.data = [surveyObj.data zeros(size(surveyObj.data,1),2)] ;

% Get operators
opObj = fcnGetDiffOperators(meshObj) ;



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
        
%         % FICTITIOUS SOURCES
%         % Compute Freespace Field
%         Ey = fcnComputeEyAnalytic(surveyObj,meshObj,TxID,f) ;
%         % Compute A Matrix in Freespace
%          A = fcnComputeA(meshObj,opObj,[],f) ;
%         % Compute Fictitious Source
%         qs = A*Ey ;
        
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
%         fcnPlotFieldGPR2D(E,meshObj) ;   
        
    end
end

% WRITE OUTPUT FILE

dataFilename = strcat(workDir,'\dpreGPR2DFWD.txt') ;
writeGPR2Ddata(surveyObj,dataFilename) ;

fprintf('\n%s\n','DATA PREDICTED')

rmpath(codeDir)








