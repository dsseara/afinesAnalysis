clear; close all;

if ispc
    load([pwd,'\simdata.mat']);
    if exist([pwd,'\actin_tiffs'],'dir')==0, mkdir([pwd,'\actin_tiffs']); end
elseif isunix
    load([pwd,'/simdata.mat']);
    if exist([pwd,'/actin_tiffs'],'dir')==0, mkdir([pwd,'/actin_tiffs']); end
end

%specify size of the plot (1 inch = 100 pixels)
units = 'inches';
imx = 5;
imy = 5;

%nmotors=params.mnum(1); % Not plotting motors or xlinkers here, so don't need
%pmotors=params.pnum(1);
nt=size(adata,3);
h_old=[];

%enter the x and y ranges below
xmin = params.xRange(1);
xmax = params.xRange(2);
ymin = params.yRange(1);
ymax = params.yRange(2);

deltaT=20;
count=0;  %added by mm

for tm=1:deltaT:nt

    adat=adata(:,:,tm);

    %Find number of polymers in this frame
    np = params.npoly(tm);

    hfig = fig('units',units,'width',imx,'height',imy,'border','off');
    
    for k=1:1:np
        idx=find(adat(:,4)==k-1);
        X1=adat(idx(:),1);Y1=adat(idx(:),2);
        ind=rangesearch([X1 Y1],[X1(1) Y1(1)],10);
        h=plot(X1(ind{1}(:)),Y1(ind{1}(:)),'-r','MarkerSize',4,'MarkerFaceColor','r','LineWidth',2);
        hold on
    end

    count=count+1;   %added by mm
    xlim([xmin xmax])
    ylim([ymin ymax])
    str = sprintf('t = %0.2f', params.timestep(tm));
    text(xmax - params.L/7, ymin+params.L/20, str, 'FontWeight', 'bold', 'FontSize', 12);
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

