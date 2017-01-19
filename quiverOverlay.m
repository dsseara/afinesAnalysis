% Script for generating a quiver overlay of the Afines data

clear, close all

fname = 'interpedData.mat';
myVars = {'adata', 'adataVelInterped', 'xGrid', 'yGrid'};
load(fname, myVars{:});

npoly = 500;
xmin = -25;
ymin = -25;
xmax = 25;
ymax = 25;
totalTime = size(adata,3);

timestep = 25;
count = 1;

for i=1:timestep:totalTime
    disp(i)
    adat=adata(:,:,i);
    figure;
    set(gcf,'Visible', 'off');
    for k=1:1:npoly
        idx=find(adat(:,4)==k-1);
        X1=adat(idx(:),1);Y1=adat(idx(:),2);
        ind=rangesearch([X1 Y1],[X1(1) Y1(1)],10);
        plot(X1(ind{1}(:)),Y1(ind{1}(:)),'-r','MarkerSize',4,'MarkerFaceColor','r','LineWidth',2);
        hold on
    end

    quiver(xGrid, yGrid, adataVelInterped(:,:,1,i), adataVelInterped(:,:,2,i), 'k');

    count=count+1;  
    xlim([xmin xmax])
    ylim([ymin ymax])
    pbaspect([1 1 1])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    axis tight 
    saveas(gcf, [pwd, '/Overlay/highThreshold/' ,num2str(count), '.tif'])
    close all
end
