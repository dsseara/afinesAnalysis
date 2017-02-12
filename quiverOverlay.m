% Script for generating a quiver overlay of the Afines data

clear, close all

fname = 'interpedData3.mat';
myVars = {'grid', 'binParams'};
load(fname, myVars{:});
myVars = {'adata', 'params'};
load('simdata.mat', myVars{:})


xmin = params.xRange(1);
ymin = params.xRange(1);
xmax = params.yRange(2);
ymax = params.yRange(2);



timestep = 1;
count = 1;
actinFrames = 1:binParams.numTimeSteps:size(adata,3);

for i=1:timestep:size(grid.vx,3)
    disp(i)
    adat=adata(:,:,actinFrames(i));
    figure;
    set(gcf,'Visible', 'off');
    for k=1:1:params.npoly
        idx=find(adat(:,4)==k-1);
        X1=adat(idx(:),1);Y1=adat(idx(:),2);
        ind=rangesearch([X1 Y1],[X1(1) Y1(1)],10);
        plot(X1(ind{1}(:)),Y1(ind{1}(:)),'-r','MarkerSize',4,'MarkerFaceColor','r','LineWidth',2);
        hold on
    end

    quiver(grid.x, grid.y, grid.vx(:,:,i), grid.vy(:,:,i), 'k');

    count=count+1;  
    xlim([xmin xmax])
    ylim([ymin ymax])
    pbaspect([1 1 1])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    axis tight 
    saveas(gcf, [pwd, '/Overlay/interp_timestep10/' ,num2str(count), '.tif'])
    close all
end
