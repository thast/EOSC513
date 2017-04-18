function fcnPlotFieldGPR2D(E,meshObj)

% REORDER

E = reshape(E,meshObj.nc(1)+1,meshObj.nc(2)+1) ;
E = fliplr(E) ;

subplot(1,2,1) ;
imagesc(real(E')) ; colorbar ;
set(gca,'fontsize',16) ;
xlabel('X [m]') ;
ylabel('Z [m]') ;

Xtick = linspace(1,meshObj.nc(1)+1,5) ;
Ztick = linspace(1,meshObj.nc(2)+1,5) ;
Xticklabel = linspace(meshObj.r0(1), meshObj.r0(1) + sum(meshObj.hx),5) ;
Zticklabel = linspace(meshObj.r0(2) + sum(meshObj.hz),meshObj.r0(2), 5) ;

set(gca,'xtick',Xtick,'xticklabel',sprintf('%.2f\n',Xticklabel)) ;
set(gca,'ytick',Ztick,'yticklabel',sprintf('%.2f\n',Zticklabel)) ;

subplot(1,2,2) ;
imagesc(imag(E')) ; colorbar ;
set(gca,'fontsize',16) ;
xlabel('X [m]') ;
ylabel('Z [m]') ;

Xtick = linspace(1,meshObj.nc(1)+1,5) ;
Ztick = linspace(1,meshObj.nc(2)+1,5) ;
Xticklabel = linspace(meshObj.r0(1), meshObj.r0(1) + sum(meshObj.hx),5) ;
Zticklabel = linspace(meshObj.r0(2) + sum(meshObj.hz),meshObj.r0(2), 5) ;

set(gca,'xtick',Xtick,'xticklabel',sprintf('%.2f\n',Xticklabel)) ;
set(gca,'ytick',Ztick,'yticklabel',sprintf('%.2f\n',Zticklabel)) ;











