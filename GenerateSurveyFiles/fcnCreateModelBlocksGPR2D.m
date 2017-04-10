function fcnCreateModelBlocksGPR2D(meshFile,geoInt,geoProps)

% INPUTS
% meshFile : Mesh file including full path
%   geoInt : N X 4 array with rows [xmin xmax zmin zmax] to give block ends
% geoProps : N X 1 array with propagation velocity 
%
% OUTPURS
% modelGPR2D.txt : Model output in same directory as mesh file


% LOAD FILE

[workDir,fileName,ext] = fileparts(meshFile) ;

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

% Center Locations
xc = r0(1) + cumsum(hx)-hx/2 ;
zc = r0(2) + cumsum(hz)-hz/2 ;

[Xc,Zc] = ndgrid(xc,zc) ;
Xc = Xc(:) ; Zc = Zc(:) ;

% CREATE MODEL

m = zeros(nc(1)*nc(2),1) ;

for pp = 1:size(geoInt,1) ;
   
    xmin = geoInt(pp,1) ;
    xmax = geoInt(pp,2) ;
    zmin = geoInt(pp,3) ;
    zmax = geoInt(pp,4) ;
    
    m(Xc>xmin & Xc<xmax & Zc>zmin & Zc<zmax) = geoProps(pp) ;
    
end

% PRINT MODEL

fileName = strcat(workDir,'\modelGPR2D.txt.') ;

FID = fopen(fileName,'w') ;
fprintf(FID,'%.8e\n',m) ;
fclose(FID) ;




















% Fill in blocks



