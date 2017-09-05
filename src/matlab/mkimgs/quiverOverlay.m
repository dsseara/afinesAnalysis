% Script for generating a quiver overlay of the Afines data
% Makes one folder with an overlay over images of the actin and one folder with an overlay over a heatmap of the velocity divergence

clear, close all

fname = 'interpedData.mat';
myVars = {'grid', 'binParams', 'divV'};
load(fname, myVars{:});
fname = 'simdata.mat';
myVars = {'adata', 'params'};
load(fname, myVars{:})

if isunix
    mkdir('Overlay/heatmap');
    mkdir('Overlay/filaments');
elseif ispc
    mkdir('Overlay\heatmap');
    mkdir('Overlay\filaments');
end

xmin = params.xRange(1);
ymin = params.xRange(1);
xmax = params.yRange(2);
ymax = params.yRange(2);

timestep = 10; % Plot every timestep-th frame
count = 1;
frames = 1:timestep:size(grid.vx,3);

for ii=frames
    disp(count)
    adat=adata(:,:,ii);
    
    figure;
    set(gcf,'Visible', 'off');
    for k=1:1:params.npoly
        idx=find(adat(:,4)==k-1);
        X1=adat(idx(:),1);Y1=adat(idx(:),2);
        ind=rangesearch([X1 Y1],[X1(1) Y1(1)],10);
        plot(X1(ind{1}(:)),Y1(ind{1}(:)),'-r','MarkerSize',4,'MarkerFaceColor','r','LineWidth',2);
        hold on
    end

    quiver(grid.x, grid.y, grid.vx(:,:,ii), grid.vy(:,:,ii), 'k');

    xlim([xmin xmax])
    ylim([ymin ymax])
    pbaspect([1 1 1])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    axis tight 

    if isunix
        saveas(gcf, [pwd, '/Overlay/filaments/' ,num2str(count), '.tif']);
    elseif ispc
        saveas(gcf, [pwd, '\Overlay\filaments\' ,num2str(count), '.tif']);
    end

    clf;

    pcolor(grid.x, grid.y, divV(:,:,ii))
    shading interp
    colormap jet
    caxis([-2, 2])
    colorbar;
    hold on;
    quiver(grid.x, grid.y, grid.vx(:,:,ii), grid.vy(:,:,ii),'k')
    
    xlim([params.xRange(1) params.xRange(2)])
    ylim([params.yRange(1) params.yRange(2)])
    pbaspect([1 1 1])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    axis tight 

    if isunix
        saveas(gcf, [pwd, '/Overlay/heatmap/' ,num2str(count), '.tif'])
    elseif ispc
        saveas(gcf, [pwd, '\Overlay\heatmap\' ,num2str(count), '.tif'])
    end

    count=count+1;  
    close all
end
