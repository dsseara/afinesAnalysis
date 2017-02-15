% Script for generating a quiver overlay of the Afines data

clear, close all

fname = 'interpedData.mat';
myVars = {'grid', 'binParams'};
load(fname, myVars{:});
myVars = {'adata', 'params'};
load('simdata.mat', myVars{:})

if ~exist('Overlay','dir')
    mkdir('Overlay')
end

xmin = params.xRange(1);
ymin = params.xRange(1);
xmax = params.yRange(2);
ymax = params.yRange(2);

timestep = 1;
count = 1;
actinFrames = 1:binParams.numTimeSteps:size(adata,3);

for ii=1:timestep:size(grid.vx,3)
    disp(ii)
    adat=adata(:,:,actinFrames(ii));
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
        saveas(gcf, [pwd, '/Overlay/' ,num2str(count), '.tif']);
    elseif ispc
        saveas(gcf, [pwd, '\Overlay\' ,num2str(count), '.tif']);
    end

    count=count+1;  
    close all
end
