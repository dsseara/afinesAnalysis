load([pwd,'\simdata.mat']);

if exist([pwd,'\actin_tiffs'])==0, mkdir([pwd,'\actin_tiffs']); end

%specify size of the plot (1 inch = 100 pixels)
units = 'inches';
imx = 5;
imy = 5;

%enter number of polymers
npoly=500;

nmotors=mnum(1);
pmotors=pnum(1);
nt=size(adata,3);
h_old=[];

%enter the x and y ranges below
xmin=-25;
xmax=25;
ymin=-25;
ymax=25;

timestep=20;
count=1;  %added by mm

for tm=1:timestep:nt

    adat=adata(:,:,tm);
    
    hfig = fig('units',units,'width',imx,'height',imy,'border','off');
    
    for k=1:1:npoly
        idx=find(adat(:,4)==k-1);
        X1=adat(idx(:),1);Y1=adat(idx(:),2);
        ind=rangesearch([X1 Y1],[X1(1) Y1(1)],10);
        h=plot(X1(ind{1}(:)),Y1(ind{1}(:)),'-r','MarkerSize',4,'MarkerFaceColor','r','LineWidth',2);
        hold on
    end

    count=count+1;   %added by mm
    xlim([xmin xmax])
    ylim([ymin ymax])
    pbaspect([1 1 1])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    set(gca,'box','off')
    axis square
    axis off
    set(gcf,'PaperUnits',units,'PaperPosition',[0 0 imx imy])
    print('-dtiff',[pwd,'\actin_tiffs\',num2str([count]),'.tif'], '-r100')

    close all;
end

