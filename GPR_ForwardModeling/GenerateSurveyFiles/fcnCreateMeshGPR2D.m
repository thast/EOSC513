function fcnCreateMeshGPR2D(workDir)

% This function take the parameters defined in parametersMeshVRM.m and
% create the mesh for the VRM code (different from UBC GIF)
% 
% INPUTS
% workDir : Working Directory as a string
%
% PARAMETERS REQUIRED FOR FILE parametersTxRxVRM
%       X0 : 1X2 array containing bottom west corner   
%       dh : 1X2 array containing minimum cell size in x,z
%       nc : 1X2 array containing number of cells in core mesh
%       np : 1X2 array containing number of padding cells
%      fac : 1X2 array containing padding factor
% 
% OUTPUTS
% MeshVRM.txt : Mesh file specific for VRM code (right-handed coordinates)

%% WRITE THE FILE

% Obtain parameters from the parameter file

paramFile = strcat(workDir,'\parametersMeshGPR2D') ;

paramFileAssert = strcat(workDir,'\parametersMeshGPR2D.m') ;
assert(exist(paramFileAssert)==2,'Parameters file does not exist or has incorrect name') ;
run(paramFile)

% Calculate outputs
hxc = strcat(num2str(dh(1),'%g'),'*',num2str(nc(1),'%g')) ;
hzc = strcat(num2str(dh(2),'%g'),'*',num2str(nc(2),'%g')) ;

Filename = strcat(workDir,'\meshGPR2D.txt') ;

% Open File
fid = fopen(Filename,'w') ;

fprintf(fid,'%i %i %s\n',nc+2*np,'! Number of Cells [x z]') ;
fprintf(fid,'%g %g %s\n',X0,'! Bottom Southwest Corner [x z]') ;

if sum(np) == 0 ;   % No padding
    fprintf(fid,'%s %s\n',hxc,'! # Cells in x') ;
    fprintf(fid,'%s %s',hzc,'! # Cells in z') ;
elseif sum(np) ~= 0 ;   % Padding
    fprintf(fid,'%g ',dh(1)*fac.^linspace(np(1),1,np(1))) ;
    fprintf(fid,'%s',hxc) ;
    fprintf(fid,' %g',dh(1)*fac.^linspace(1,np(1),np(1))) ;
    fprintf(fid,' %s\n','! Cell widths x-direction') ;
    fprintf(fid,'%g ',dh(2)*fac.^linspace(np(2),1,np(2))) ;
    fprintf(fid,'%s',hzc) ;
    fprintf(fid,' %g',dh(2)*fac.^linspace(1,np(2),np(2))) ;
    fprintf(fid,' %s','! Cell widths z-direction') ;
end

fclose(fid) ;



















