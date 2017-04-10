function fcnViewModelGRPR2D(workDir,meshFile,modelFile)

meshFile  = strcat(workDir,'\',meshFile) ;
modelFile = strcat(workDir,'\',modelFile) ;

FID = fopen(meshFile,'r') ;

A = textscan(FID,'%s','delimiter','\n') ;
A = A{1} ;

fclose(FID) ;

for ii = 1:length(A) ;
    temp = strsplit(A{ii},' !') ;
    temp = strsplit(temp{1},' ') ;
    A{ii} = temp ;
end

% GET MESH PROPERTIES

% Number of cells
nc = cellfun(@str2num,A{1}) ;

% Origin
r0 = cellfun(@str2num,A{2}) ;

% x cells
k = cellfun(@isempty,strfind(A{3},'*')) ;
pad = cellfun(@str2num,A{3}(k==1)) ;
core = A{3}(k==0) ;
core = str2double(strsplit(core{1},'*')) ;
core = core(1)*ones(1,core(2)) ;
np = length(pad)/2 ;
hx = [pad(1:np) core pad(end-np+1:end)]' ;

% z cells
k = cellfun(@isempty,strfind(A{4},'*')) ;
pad = cellfun(@str2num,A{4}(k==1)) ;
core = A{4}(k==0) ;
core = str2double(strsplit(core{1},'*')) ;
core = core(1)*ones(1,core(2)) ;
np = length(pad)/2 ;
hz = [pad(1:np) core pad(end-np+1:end)]' ;

meshObj = struct('r0',r0,'nc',nc,'hx',hx,'hz',hz) ;

% LOAD MODEL

m = dlmread(modelFile) ;
m = reshape(m,meshObj.nc(1),meshObj.nc(2)) ;
m = fliplr(m) ;

% PLOT MODEL

imagesc(m') ; colorbar ;
set(gca,'fontsize',16) ;
xlabel('X [m]') ;
ylabel('Z [m]') ;

Xtick = linspace(1,meshObj.nc(1),5) ;
Ztick = linspace(1,meshObj.nc(2),5) ;
Xticklabel = linspace(meshObj.r0(1), meshObj.r0(1) + sum(meshObj.hx),5) ;
Zticklabel = linspace(meshObj.r0(2) + sum(meshObj.hz),meshObj.r0(2), 5) ;

set(gca,'xtick',Xtick,'xticklabel',sprintf('%.2f\n',Xticklabel)) ;
set(gca,'ytick',Ztick,'yticklabel',sprintf('%.2f\n',Zticklabel)) ;




