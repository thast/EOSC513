function meshObj = loadmeshGPR2D(meshFile)

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



% ASSIGN MESH OBJECT PROPERTIES
meshObj = struct('r0',r0,'nc',nc,'hx',hx,'hz',hz) ;


